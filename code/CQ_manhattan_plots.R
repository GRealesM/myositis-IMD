#########################################
##                                     ##
##     Preparing Manhattan plots       ##
##                                     ##
#########################################

# Author: Guillermo Reales 
# Date last updated: 2024/05/02

# Background: We want to see the Manhattan plots of some good colocalisation examples.


## Load packages and required data sets
library(data.table)
library(magrittr)
library(ggplot2)
library(cowplot)

setDTthreads(20)
setwd("/home/gr440/rds/rds-cew54-basis/Projects/myositis-IMD/code")

coloc <- fread("../data/coloc_results-v3.tsv")
# Import mapped genes 
mg <- fread("../data/mapped.genes-v3.tsv") %>% unique
mg <- mg[,.(SNPID, pid, nearestGene)]
# Add some info to the coloc table, and extract rsids to map
coloc <- merge(coloc, mg, by="pid", all.x = TRUE) # First on driver SNPs
setnames(coloc, c("SNPID", "nearestGene"), c("driver.rsid", "driver.nearestgene"))
coloc <- merge(coloc, mg, by.x="bestsnp", by.y = "pid", all.x = TRUE) # Then on candidate snps. Bear in mind that we only mapped candidate SNPs with H4 > 0.5
setnames(coloc, c("SNPID", "nearestGene"), c("bestsnp.rsid", "bestsnp.nearestgene"))
coloc[, c("CHR38", "BP38"):= tstrsplit(pid, split = ":")][, CHR38:= as.numeric(CHR38)][, BP38:= as.numeric(BP38)]

# List of diseases with their file names
st4 <- fread("../tables/ST4_coloc_diseases.tsv")
st4 <- st4[, .(Trait, coloc_Label)]

# We'll list here our selected examples
# - rs2476601 - PTPN22, good signal. We can use the IIM (R) example and its colocalising diseases.
# - rs5754217 - UBE2L3, strong coloc signal, only SSc. We can use IIM (R) to try.

snp <- "rs2476601" # Choose example
iim <- "IIM (R)" # Choose IIM dataset to include
cls <- coloc[ H4 > 0.4 & driver.rsid == snp & trait.myos == iim, ]
dpid <- cls$pid[1]
chr <- cls$CHR38[1]
bp <- cls$BP38[1]
w = 1e6 # window
at <- c(cls$trait.myos, cls$trait.other)
tf <- st4[coloc_Label %in% at, Trait]

# Prepare data
# Since we need dense data, and we don't have the full files for FinnGen and PanUKBBR2, we need to download them and process them.
# I've been doing this for different files, so most of them will be in ../data/fg_sumstats/, but I added the piece of code below to
# download and process files automatically, if some is missing.
data <- lapply(tf, function(x){
        message("Working on ", x, ".")
        if(grepl("FinnGenR7", x)){
			fn <- gsub("_FinnGen_FinnGenR7_1", "", x)
			if(!file.exists(paste0("../data/fg_sumstats/finngen_R7_", fn,".gz"))) system(paste0("wget ***REMOVED***finngen_R7_", fn, ".gz -O ../data/fg_sumstats/finngen_R7_",fn,".gz"))
			y = fread(paste0("../data/fg_sumstats/finngen_R7_", fn,".gz"), tmpdir = "tmp")
			y[ , pid38:=paste(`#chrom`, pos, sep = ":")]
			y <- y[, .(`#chrom`, pos, ref, alt, beta, sebeta, pval)]
			setnames(y, c("#chrom", "pos", "ref", "alt", "beta", "sebeta", "pval"), c("CHR38", "BP38", "REF", "ALT", "BETA", "SE", "P"))
        } else if(grepl("PanUKBBR2", x)){
			if(!file.exists(paste0("../data/fg_sumstats/",x,"-hg38.tsv.gz"))){
				# PanUKBB requires a bit more processing
				# Transform into .gz so we can import it in fread
				mf <- fread("../data/Metadata_20230906-v1.tsv")
				url <- mf[Trait == x, URL]
				system(paste0("wget  ", url, " -O ../data/fg_sumstats/", x, ".bgz"))
				system(paste0("zcat ../data/fg_sumstats/",x,".bgz > ../data/fg_sumstats/",x,".tsv ; gzip ../data/fg_sumstats/", x, ".tsv"))
				f <- fread(paste0("../data/fg_sumstats/",x, ".tsv.gz"), tmpdir = "tmp")
				f <- f[, .(chr, pos, ref, alt, beta_meta, se_meta, neglog10_pval_meta)]
				names(f) <- c("CHR", "BP", "REF", "ALT", "BETA", "SE", "P")
				f <- f[complete.cases(f)]
				f[, P:=10^-P] # P-values are stored as -log10(P), so we'll revert them -- even though we'll likely don't need the P-values

				# Save it
				fwrite(f, paste0("../data/fg_sumstats/",x,"-hr.tsv.gz"), sep = "\t")

				# Pipeline it 
				system(paste0("cd ../data/fg_sumstats/; ~/rds/rds-cew54-basis/GWAS_tools/01-Pipeline/pipeline_v5.3.2_beta.sh -f ", x,"-hr.tsv.gz"))
				unlink(paste0("../data/fg_sumstats/", x, ".bgz"))
				unlink(paste0("../data/fg_sumstats/", x, ".tsv.gz"))
				unlink(paste0("../data/fg_sumstats/", x, "-hr.tsv.gz"))
			}
				y <- fread(paste0("../data/fg_sumstats/",x,"-hg38.tsv.gz"), tmpdir = "tmp")
		} else{
            y  <- fread(paste0("~/rds/rds-cew54-basis/02-Processed/", x, "-hg38.tsv.gz"), tmpdir = "tmp")

        }
    
        y <- y[ CHR38 == chr & BP38 >= bp - w & BP38 <= bp + w & !is.na(SE) & !is.na(BETA)]
        y[ , pid:=paste(CHR38, BP38, sep =":")]
        y <- y[ !duplicated(pid), .(pid, CHR38, BP38, REF, ALT, BETA, SE, P)]
        y[, study:=st4[ Trait == x, coloc_Label]]

})  %>% rbindlist
commonpid <- data[, .(count = .N), by = pid][ count == max(count), pid] # Keep SNPs common to all datasets (Note: bear in mind that we're not checking alleles here).
data <- data[pid %in% commonpid]
data

 ggplot(data, aes(x=BP38, y=-log10(P))) + 
            geom_point() + 
            geom_point(data = data[pid == dpid], col = "red", size = 3)+
            facet_grid(study~.,scales="free_y") + 
            # geom_vline(xintercept=bp,col="red")+
            geom_hline(yintercept = -log10(5e-8), col = "red")+
            ggtitle(paste0(dpid, " / ", snp))+
            theme_cowplot() + 
            theme(axis.title.x = element_blank(),
                  strip.background = element_rect(colour="black", fill = "white"),
                  )



sessionInfo()
# R version 4.3.3 (2024-02-29)
# Platform: x86_64-redhat-linux-gnu (64-bit)
# Running under: Rocky Linux 8.9 (Green Obsidian)

# Matrix products: default
# BLAS/LAPACK: /usr/lib64/libopenblaso-r0.3.15.so;  LAPACK version 3.9.0

# locale:
#  [1] LC_CTYPE=en_GB.UTF-8       LC_NUMERIC=C               LC_TIME=en_GB.UTF-8        LC_COLLATE=en_GB.UTF-8     LC_MONETARY=en_GB.UTF-8    LC_MESSAGES=en_GB.UTF-8   
#  [7] LC_PAPER=en_GB.UTF-8       LC_NAME=C                  LC_ADDRESS=C               LC_TELEPHONE=C             LC_MEASUREMENT=en_GB.UTF-8 LC_IDENTIFICATION=C       

# time zone: GB
# tzcode source: system (glibc)

# attached base packages:
# [1] stats     graphics  grDevices utils     datasets  methods   base     

# other attached packages:
# [1] cowplot_1.1.3     ggplot2_3.5.0     magrittr_2.0.3    data.table_1.15.4

# loaded via a namespace (and not attached):
#  [1] vctrs_0.6.5       cli_3.6.2         rlang_1.1.3       generics_0.1.3    jsonlite_1.8.8    labeling_0.4.3    glue_1.7.0        colorspace_2.1-0  scales_1.3.0     
# [10] fansi_1.0.6       grid_4.3.3        munsell_0.5.1     tibble_3.2.1      lifecycle_1.0.4   compiler_4.3.3    dplyr_1.1.4       pkgconfig_2.0.3   httpgd_2.0.1     
# [19] R.oo_1.26.0       farver_2.1.1      systemfonts_1.0.6 R.utils_2.12.3    R6_2.5.1          unigd_0.1.1       tidyselect_1.2.1  utf8_1.2.4        pillar_1.9.0     
# [28] R.methodsS3_1.8.2 tools_4.3.3       withr_3.0.0       gtable_0.3.4 