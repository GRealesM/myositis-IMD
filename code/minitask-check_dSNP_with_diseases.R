################################################
### Check SNPs for key genes and diseases    ###
################################################

# Author: Guillermo Reales
# Date: 2022-01-07

# Background: Upon doing pathway analysis and extracting the exclusive pathways for PC8, 
# two pathways stood out (ie. GO:0042267 and GO:0038110, associated with NK toxicity and IL-2, respectively).
# These pathways are associated with 8 genes (described below). We want to see if the driver SNPs 
# associated with these genes are associated with diseases at one and the other side of the component.


## Load stuff
library(data.table)

# Genes
pw1 <- c("GZMB", "RASGRP1", "NCR1","IL12A", "IL12B")
pw2 <- c("IL2RA", "IL2RB", "SYK")


# SNP manifest
manifest.translator <- fread("../data/Manifest_build_translator.tsv")
manifest.translator[,pid19:=paste(CHR19, BP19, sep = ":")][, pid38:=paste(CHR38, BP38, sep = ":")]
manifest.translator <- manifest.translator[,c(1,8:9)]
SNP.manifest <- merge(manifest.translator, copy(cupcake::SNP.manifest), by.x = "pid19", by.y = "pid")

# Annotations. Same treatment as in the report
otg <- readRDS("../data/imdbasis-driver-snps-otg-annotations.RDS")
otg <- merge(SNP.manifest[,.(pid38, SNPID, ref_a1, ref_a2)], otg, by.x="pid38", by.y="pid")
# Check concordance in alleles
all(otg$ref_a1 == otg$ref_allele)
## [1] TRUE
all(otg$ref_a2 == otg$alt_allele) # FALSE, let's see
## [1] FALSE
nrow(otg[otg$ref_a2 != otg$alt_allele]) # 152 instances of alt allele discordance. We'll weed them out
## [1] 152
otg <- otg[ref_a2 == alt_allele]

# Remove redundant columns
otg[, c("source_list", "source_score_list", "wh.vep") := NULL]

# Include rotation values in otg
otg <- merge(otg, manifest.translator[, c("pid19", "pid38")], by = "pid38")

rotmat <- as.data.table(cupcake::rot.pca, keep.rownames = TRUE)
names(rotmat)[1] <- "pid19"
names(rotmat) <- gsub("PC", "rot", names(rotmat))

otg <- merge(otg, rotmat, by="pid19")

# Add eQTL data
eqtl <- fread("../data/Causal_eQTLs_full_20210915.tsv") 
# Noticed a bug in the procesing of snp_study and study_tissue. Let's briefly fix it
eqtl[, snp_study:=gsub(pattern = "(_[A-Za-z0-9\\+ -]+)$", "", snp_study, perl = TRUE)][, study_tissue:=gsub(pattern = "(_[A-Za-z0-9\\+ -]+)$", "", study_tissue, perl = TRUE)][,study:=gsub(pattern = "(_[A-Za-z0-9\\+ -]+)$", "", study, perl = TRUE)][,quant_method:=gsub("([A-Za-z0-9_]+)_([a-z]+)", "\\2", study, perl = TRUE)][,study:=gsub(pattern = "(_[a-z]+)$", "", study, perl = TRUE)]


# Take a look at the driver SNPs associated with the genes
otg[PC8==TRUE & tss_score == 1 & hgnc_symbol %in% pw1, .(pid19, pid38, SNPID, ref_allele, alt_allele, source, hgnc_symbol, rot8 )]
#          pid19       pid38      SNPID ref_allele alt_allele
# 1: 14:25102160 14:24632954  rs8192917          C          T
# 2: 15:38838264 15:38546063  rs8035957          T          C
# 3: 19:55380214 19:54868759 rs17771967          A          G
# 4: 3:159728878 3:160011091  rs6441286          T          G
# 5: 5:158759900 5:159332892  rs2546890          A          G
# 6: 5:158787385 5:159360377  rs6556412          G          A
#                         source hgnc_symbol        rot8
# 1: eqtl,pqtl,canonical_tss,vep        GZMB  0.05962021
# 2:      vep,eqtl,canonical_tss     RASGRP1  0.04079473
# 3:               canonical_tss        NCR1 -0.08594753
# 4:          eqtl,canonical_tss       IL12A -0.05335533
# 5:               canonical_tss       IL12B -0.04110094
# 6:               canonical_tss       IL12B -0.04523441

eqtl[rot8 !=0 & gene_name %in% pw1, .(pid, SNPID.basis, ref, alt, gene_name, snp_study)]
#            pid SNPID.basis ref alt gene_name
# 1: 15:38546063   rs8035957   T   C   RASGRP1
# 2: 15:38546063   rs8035957   T   C   RASGRP1
#                                       snp_study
# 1:       rs8035957_T_C_Schmiedel_2018_ge_B cell
# 2: rs8035957_T_C_Fairfax_2012_microarray_B cell

otg[PC8==TRUE & tss_score == 1 & hgnc_symbol %in% pw2, .(pid38, SNPID, ref_allele, alt_allele, hgnc_symbol)]
#          pid38     SNPID ref_allele alt_allele hgnc_symbol
# 1:  10:6068912 rs7090530          C          A       IL2RA
# 2: 22:37185445  rs229527          C          A       IL2RB
# 3:  9:90801254  rs290986          A          G         SYK

eqtl[rot8 !=0 & gene_name %in% pw2, .(pid, SNPID.basis, ref, alt, gene_name, snp_study)]
# Empty data.table (0 rows and 6 cols): pid,SNPID.basis,ref,alt,gene_name,snp_study


# Join all
o1 <- otg[PC8==TRUE & tss_score == 1 & hgnc_symbol %in% pw1, .(pid38, SNPID, ref_allele, alt_allele, hgnc_symbol)]
names(o1) <- c("pid38", "SNPID", "REF","ALT","gene_name")
e1 <- unique(eqtl[rot8 !=0 & gene_name %in% pw1, .(pid, SNPID.basis, ref, alt, gene_name)])
names(e1) <- c("pid38", "SNPID", "REF","ALT", "gene_name")
dsnp1  <- rbind(o1,e1)
dsnp1[, pathway:="pw1"]

o2 <- otg[PC8==TRUE & tss_score == 1 & hgnc_symbol %in% pw2, .(pid38, SNPID, ref_allele, alt_allele, hgnc_symbol)]
names(o2) <- c("pid38", "SNPID", "REF","ALT","gene_name")
o2[, pathway:="pw2"]

smr <- rbind(dsnp1, o2)
smr
#           pid38      SNPID REF ALT gene_name pathway
#  1: 14:24632954  rs8192917   C   T      GZMB     pw1
#  2: 15:38546063  rs8035957   T   C   RASGRP1     pw1
#  3: 19:54868759 rs17771967   A   G      NCR1     pw1
#  4: 3:160011091  rs6441286   T   G     IL12A     pw1
#  5: 5:159332892  rs2546890   A   G     IL12B     pw1
#  6: 5:159360377  rs6556412   G   A     IL12B     pw1
#  7: 15:38546063  rs8035957   T   C   RASGRP1     pw1
#  8:  10:6068912  rs7090530   C   A     IL2RA     pw2
#  9: 22:37185445   rs229527   C   A     IL2RB     pw2
# 10:  9:90801254   rs290986   A   G       SYK     pw2


# Sounds about right. Let's now check the diseases
# We have Asthma and MS on the positive side, and PSC, T1D, and UC on the other side.
# We'll look for those representatives of the diseases that were significant at the component
redir <- "../../../03-Bases/IMD_basis/reduced_datasets"

dsets <- grep("MS_IMSGC_31604244|AST_Moffatt_20860503|T1D_Chiou_34012112|UC_Liu_26192919_1" , dir(redir), perl=TRUE, value= TRUE)

diss <- c("MS", "AST", "T1D", "UC")

dis.driver <- data.table()

for (i in 1:length(dsets)){

	x <- fread(paste0(redir, "/", dsets[i]))
	x[, pid38:=paste(CHR38, BP38, sep=":")]
	x  <- x[pid38 %in% smr$pid38]
	x[, trait:=diss[i]]
	dis.driver <- rbindlist(list(dis.driver, x), fill=TRUE)

}
dis.driver  <- merge(dis.driver, smr[, .(pid38, gene_name, pathway)], by="pid38")

# Let's take a look at the pathways
dis.driver[pathway == "pw1" & P < 10e-4]
#          pid38      SNPID CHR38      BP38 REF ALT       BETA         SE         P trait gene_name pathway
# 1: 19:54868759 rs17771967    19  54868759   A   G  0.0832947 0.02220000 1.696e-04    UC      NCR1     pw1
# 2: 5:159332892  rs2546890     5 159332892   A   G -0.1169830 0.01641812 1.039e-12   AST     IL12B     pw1
# 3: 5:159360377  rs6556412     5 159360377   G   A -0.0775164 0.01734212 7.828e-06   AST     IL12B     pw1
# 4: 5:159360377  rs6556412     5 159360377   G   A  0.0935996 0.02240000 2.992e-05    UC     IL12B     pw1
dis.driver[pathway == "pw2" & P < 10e-4]
#          pid38     SNPID CHR38     BP38 REF ALT      BETA         SE         P trait gene_name pathway
# 1:  10:6068912 rs7090530    10  6068912   C   A  0.192842 0.01450100 2.350e-40   T1D     IL2RA     pw2
# 2: 22:37185445  rs229527    22 37185445   C   A  0.104077 0.01413800 1.820e-13   T1D     IL2RB     pw2
# 3:  9:90801254  rs290986     9 90801254   A   G -0.100207 0.02078889 1.434e-06   AST       SYK     pw2





