## Extracting driver and best SNPs

# Author: Guillermo Reales
# Date: 2023-08-07

# After finishing coloc, we'll need to extract the chromosomes, rsids, alleles, etc. from the relevant driver and bestsnps.
# Then we'll use mapping-genes.py on the generated file to extract the names of the nearest genes.

setwd("/home/gr440/rds/rds-cew54-basis/Projects/myositis-IMD/code")

library(data.table)
setDTthreads(15)
library(magrittr)

cc <- fread("../data/coloc_results_dfilt-v3.tsv")

snps <- c(cc[ H4 > 0.5, unique(pid)],cc[ H4 > 0.5, unique(bestsnp)]) %>% unique

# We'll extract all info from the DM (M) dataset
dm <- fread("~/rds/rds-cew54-basis/02-Processed/DMY_Miller_26291516_1-hg38.tsv.gz")
dm[, pid:=paste(CHR38, BP38, sep=":")]
dm <- dm[pid %in% snps, .(pid, CHR38, BP38, SNPID, REF, ALT)]
names(dm)[2:4] <- c("CHR", "BP", "SNP")


fwrite(dm, "../data/snp.to.map.tsv", sep="\t")

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
# [1] magrittr_2.0.3    data.table_1.14.8

# loaded via a namespace (and not attached):
# [1] compiler_4.1.3    cli_3.6.0         R.methodsS3_1.8.2 jsonlite_1.8.4   
# [5] R.utils_2.12.2    rlang_1.0.6       R.oo_1.25.0  
