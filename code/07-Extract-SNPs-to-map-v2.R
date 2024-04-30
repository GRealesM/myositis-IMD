#########################################
##                                     ##
##  EXTRACTING DRIVER AND CANDIDATES   ##
##                                     ##
#########################################

# Author: Guillermo Reales 
# Date last updated: 2024/04/29

# Background: After finishing coloc, we'll need to extract the chromosomes, rsids, alleles, etc. from the relevant driver and bestsnps.
# Then we'll use 08-Mapping-genes.py on the generated file to extract the names of the nearest genes.

# This script is meant to be run at the HPC
setwd("/home/gr440/rds/rds-cew54-basis/Projects/myositis-IMD/code")


#########################################


# Load packages
library(data.table)
setDTthreads(15)
library(magrittr)

# Load data
cc <- fread("../data/coloc_results-v3.tsv")

snps <- c(cc[ H4 > 0.5, unique(pid)],cc[ H4 > 0.5, unique(bestsnp)]) %>% unique

# We'll extract all info from the DM (M) dataset
mm <- fread("~/rds/rds-cew54-basis/02-Processed/MYO_Miller_26291516_1-hg38.tsv.gz")
mr <- fread("~/rds/rds-cew54-basis/02-Processed/IIM_Rothwell_up_1-hg38.tsv.gz")
mm[, pid:=paste(CHR38, BP38, sep=":")]
mr[, pid:=paste(CHR38, BP38, sep=":")]
mm <- mm[pid %in% snps, .(pid, CHR38, BP38, SNPID, REF, ALT)]
mr <- mr[pid %in% snps, .(pid, CHR38, BP38, SNPID, REF, ALT)]
all(snps %in% mm$pid) # Not all SNPs present
setdiff(snps, mm$pid) # Locate SNP
mr[pid %in% setdiff(snps, mm$pid)] # is it in the alternative dataset? It is
mm <- rbindlist(list(mm, mr[pid %in% setdiff(snps, mm$pid)])) # add missing snp to dataset
mm[pid == "8:11484361", SNPID:="rs2736336"] # add missing rsid

names(mm)[2:4] <- c("CHR", "BP", "SNP")


# Save
fwrite(mm, "../data/snp.to.map-v3.tsv", sep="\t")

#########################################

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
# [1] magrittr_2.0.3    data.table_1.15.4

# loaded via a namespace (and not attached):
# [1] compiler_4.3.3    cli_3.6.2         tools_4.3.3       R.methodsS3_1.8.2 jsonlite_1.8.8    R.utils_2.12.3    rlang_1.1.3       R.oo_1.26.0    