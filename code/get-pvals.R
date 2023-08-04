### Getting p-values of top candidate SNPs in their respective myositis datasets.
### We'll also take the chance to made some ladder plots of SNPs of interest

# Date: 2023/08/01
# Author: Guillermo Reales

setwd("/home/gr440/rds/rds-cew54-basis/Projects/myositis-IMD/code")

# Load pkgs

library(data.table)
setDTthreads(15)
library(magrittr)

# Set variables and load summary statistics


files = list(   dmy.m = "DMY_Miller_26291516_1-hg38.tsv.gz",
                dmy.r = "DMY_Rothwell_up_1-hg38.tsv.gz",
                jdm.m = "JDM_Miller_26291516_1-hg38.tsv.gz",
                jdm.r = "JDM_Rothwell_up_1-hg38.tsv.gz",
                jo1m.r = "JO1M_Rothwell_up_1-hg38.tsv.gz",
                myo.m = "MYO_Miller_26291516_1-hg38.tsv.gz",
                myo.r = "IIM_Rothwell_up_1-hg38.tsv.gz",
                pm.m =   "PM_Miller_26291516_1-hg38.tsv.gz",
                pm.r = "PM_Rothwell_up_1-hg38.tsv.gz") %>%
    lapply(.,  function(f) (file.path("~/rds/rds-cew54-basis/02-Processed",f)))


man=fread("~/rds/rds-cew54-basis/03-Bases/IMD_basis/SNP.manifest.38.tsv")
snps=fread("~/rds/rds-cew54-basis/03-Bases/IMD_basis/Manifest_build_translator.tsv")
snps[, pid38:=paste(CHR38, BP38, sep=":")][, pid19:=paste(CHR19, BP19, sep=":")]
snps <- merge(snps, man, by="pid38")

data=lapply(files, fread)


# Import mapped genes 
mg <- fread("../data/mapped.genes.tsv") %>% unique
mg <- mg[,.(SNPID, pid, nearestGene)]

# Import coloc table and add some info to the coloc table, and extract rsids to map
coloc <- fread("../tables/coloc_results_dfilt.tsv")
coloc <- merge(coloc, mg, by="pid") # First on driver SNPs
setnames(coloc, c("SNPID", "nearestGene"), c("driver.rsid", "driver.nearestgene"))
coloc <- merge(coloc, mg, by.x="bestsnp", by.y = "pid", all.x = TRUE) # Then on candidate snps. Bear in mind that we only mapped candidate SNPs with H4 > 0.5
setnames(coloc, c("SNPID", "nearestGene"), c("bestsnp.rsid", "bestsnp.nearestgene"))


# List top SNPs
#top <- c("rs11692867","rs2736340","rs819991","rs13236009","rs3821236","rs2476601","rs597808","rs597808")
top <- c("rs2736340", "rs13236009","rs2476601") 

topids <- coloc[ bestsnp.rsid %in% top, .(bestsnp, bestsnp.rsid)]  %>% unique
studies <- names(data)

# Extract required SNPs
d2 <- lapply(1:9, function(i){
        x <- data[[i]]
        x[, pid := paste(CHR38, BP38, sep=":")]
        x  <- x[ pid %in% topids$bestsnp , .(pid, SNPID, REF, ALT, BETA, SE, P)]
        x[, study := studies[i]]
        x

}) %>% rbindlist

d2 <- merge(d2, topids, by.x="pid", by.y="bestsnp")

d2[ bestsnp.rsid %in% top[1:2], .(bestsnp.rsid, pid, REF,ALT, P, study)]
#    bestsnp.rsid         pid REF ALT           P  study
#  1:   rs13236009 7:129023119   T   G 2.01914e-03  dmy.m
#  2:   rs13236009 7:129023119   T   G 1.00165e-02  dmy.r
#  3:   rs13236009 7:129023119   T   G 7.47704e-01  jdm.m
#  4:   rs13236009 7:129023119   T   G 2.05076e-01  jdm.r
#  5:   rs13236009 7:129023119   T   G 1.16643e-03 jo1m.r
#  6:   rs13236009 7:129023119   T   G 1.34854e-03  myo.m
#  7:   rs13236009 7:129023119   T   G 6.41047e-05  myo.r
#  8:   rs13236009 7:129023119   T   G 2.29058e-02   pm.m
#  9:   rs13236009 7:129023119   T   G 5.17526e-02   pm.r
# 10:    rs2736340  8:11486464   C   T 2.69300e-03  dmy.m
# 11:    rs2736340  8:11486464   C   T 1.43846e-02  dmy.r
# 12:    rs2736340  8:11486464   C   T 4.18800e-03  jdm.m
# 13:    rs2736340  8:11486464   C   T 1.85708e-03  jdm.r
# 14:    rs2736340  8:11486464   C   T 2.28270e-01 jo1m.r
# 15:    rs2736340  8:11486464   C   T 1.52900e-04  myo.m
# 16:    rs2736340  8:11486464   C   T 1.92214e-04  myo.r
# 17:    rs2736340  8:11486464   C   T 2.96900e-02   pm.m
# 18:    rs2736340  8:11486464   C   T 1.03374e-02   pm.r

## Extract info on rs2476601 to create a ladder plot
d3  <- d2[bestsnp.rsid == "rs2476601"]

ml <- data.table(mlabel = c("DM (M)", "DM (R)", "JDM (M)","JDM (R)", "Anti-Jo1+ (R)", "IIM (M)", "IIM (R)", "PM (M)", "PM (R)"),
                 trait.myos = c("dmy.m", "dmy.r", "jdm.m", "jdm.r", "jo1m.r", "myo.m", "myo.r", "pm.m", "pm.r"))
d3  <- merge(d3, ml, by.x = "study", by.y = "trait.myos")
myoc <- c(`PM (R)` = "#CF000F", `PM (M)` = "#CF000F", `PM (FG)` = "#CF000F", `DM (R)` = "#2E8856", `DM (M)` = "#2E8856", `IIM (R)` = "#1460AA", `IIM (M)` = "#1460AA", `JDM (M)` = "#B8860B", `JDM (R)` = "#B8860B", `IBM (R)` = "#E65722", `Anti-Jo1+ (R)` ="#1C2833", `DPM (FG)` = "#053061")

library(ggplot2)
library(cowplot)
dpm <- ggplot(d3, aes(x = BETA, y = mlabel, xmin=BETA-SE, xmax=BETA+SE, colour = mlabel))+
  geom_pointrange()+
  geom_vline(xintercept = 0, col="red", lty=2)+
  scale_colour_manual(values = myoc)+
  xlab("Beta")+
  theme_cowplot(12)+
  theme(legend.position = "none", axis.title.y = element_blank())
dpm
ggsave("../figures/ladder_plot_rs2476601.png", dpm, height = 3, width = 5, bg="white")


sessionInfo()
# R version 4.1.3 (2022-03-10)
# Platform: x86_64-pc-linux-gnu (64-bit)
# Running under: Rocky Linux 8.7 (Green Obsidian)

# Matrix products: default
# BLAS/LAPACK: /usr/local/software/spack/spack-rhel8-20210927/opt/spack/linux-centos8-icelake/gcc-11.2.0/intel-oneapi-mkl-2021.4.0-s2cksi33smowj5zlqvmew37cufvztdkc/mkl/2021.4.0/lib/intel64/libmkl_gf_lp64.so.1

# locale:
#  [1] LC_CTYPE=en_GB.UTF-8       LC_NUMERIC=C              
#  [3] LC_TIME=en_GB.UTF-8        LC_COLLATE=en_GB.UTF-8    
#  [5] LC_MONETARY=en_GB.UTF-8    LC_MESSAGES=en_GB.UTF-8   
#  [7] LC_PAPER=en_GB.UTF-8       LC_NAME=C                 
#  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
# [11] LC_MEASUREMENT=en_GB.UTF-8 LC_IDENTIFICATION=C       

# attached base packages:
# [1] stats     graphics  grDevices utils     datasets  methods   base     

# other attached packages:
# [1] cowplot_1.1.1     ggplot2_3.4.1     magrittr_2.0.3    data.table_1.14.8

# loaded via a namespace (and not attached):
#  [1] Rcpp_1.0.10       tidyselect_1.2.0  munsell_0.5.0     colorspace_2.1-0 
#  [5] R6_2.5.1          ragg_1.2.5        rlang_1.0.6       fansi_1.0.4      
#  [9] dplyr_1.1.0       tools_4.1.3       grid_4.1.3        gtable_0.3.1     
# [13] R.oo_1.25.0       utf8_1.2.3        cli_3.6.0         withr_2.5.0      
# [17] systemfonts_1.0.4 tibble_3.1.8      lifecycle_1.0.3   textshaping_0.3.6
# [21] farver_2.1.1      later_1.3.0       vctrs_0.5.2       R.utils_2.12.2   
# [25] glue_1.6.2        labeling_0.4.2    compiler_4.1.3    pillar_1.8.1     
# [29] generics_0.1.3    scales_1.2.1      R.methodsS3_1.8.2 httpgd_1.3.1     
# [33] jsonlite_1.8.4    pkgconfig_2.0.3  