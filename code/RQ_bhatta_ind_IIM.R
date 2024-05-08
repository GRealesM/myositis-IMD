#########################################
##                                     ##
##     Computing Bhatta one by one     ##
##                                     ##
#########################################

# Author: Guillermo Reales 
# Date last updated: 2024/05/07

# Background: Responding to a reviewer query, we wanted to run the Bhattacharyya distance process but using one myositis dataset at a time.


### Set paths

# This script is meant to be run at the HPC
setwd("/home/gr440/rds/rds-cew54-basis/Projects/myositis-IMD/code")
bpath <- "/home/gr440/rds/rds-cew54-basis/03-Bases/IMD_basis/"

### Load required packages

library(data.table)
library(magrittr)
library(stringr)
library(pheatmap)
library(reshape2)

# Load key data about datasets
bhd <- fread("../data/bh_dist.tsv", sep="\t")

bhd[, logbhat:=log(1 + bhat.dist)]
bhd[, T1:=str_trunc(T1, 50, ellipsis = " [...]")][, T2:=str_trunc(T2, 50, ellipsis = " [...]")]

# Call clusters in Bhattacharyya. First for the full dataset, so we extract the clustered diseases
bh.cl <- acast(bhd, T1~T2, value.var = "bhat.dist") %>% as.dist %>% hclust(., method = "complete") # Cluster by raw
bhatta.d <-  acast(bhd, T1~T2, value.var = "logbhat") # Show log

bh.ph <- pheatmap(bhatta.d,
                  cluster_cols = bh.cl, 
                  cluster_rows = bh.cl,
                  annotation_names_row = FALSE,
                  show_colnames = FALSE,
                  annotation_legend = FALSE,
                  cutree_rows = 9,
                  treeheight_col=0,
                  fontsize_row=8)

tbh <- bh.ph$tree_row
tbhc <- cutree(tbh, k=9) # manually adjusted to capture visually-selected IMD
bhatta.sel <- names(tbhc[tbhc ==5]) # Myositis cluster is 9
myosets <- bhatta.sel[grepl("Miller|Rothwell", bhatta.sel)] # Myositis datasets
bhatta.sel <- bhatta.sel[!grepl("Miller|Rothwell", bhatta.sel)] # keep close IMD according to Bhatta only.

bhind <- lapply(myosets, function(x){
    myoexcl <- myosets[!myosets %in% x]
    b1 <- bhd[ !T1 %in% myoexcl & !T2 %in% myoexcl] # Remove excluded datasets, keep the focus only
    b1.cl <- acast(b1, T1~T2, value.var = "bhat.dist") %>% as.dist %>% hclust(., method = "complete") # Cluster by raw
    b1.c <- cutree(b1.cl, k=9) # manually adjusted to capture visually-selected IMD
    cl.n <- b1.c[x] # Extract cluster number
    b1.sel <- names(b1.c[b1.c == cl.n]) # Extract IMD in the cluster 
    b1.sel <- b1.sel[!b1.sel %in% x]
    r1 <- data.table(Myositis = x,
    N.match = length(b1.sel[b1.sel %in% bhatta.sel]), # Number of IMD in cluster also present in the all myositis Bhatta result
    N.nomatch = length(setdiff(b1.sel, bhatta.sel))) # Number of IMD in cluster not present in the all myo Bhatta result
    r1

})  %>% rbindlist
bhind

qf <- fread("../data/qf.tsv") %>% .[grepl("Miller|Rothwell", Label), .(Label, N0, N1, N)]
bhind <- merge(bhind, qf, by.x = "Myositis", by.y = "Label")
bhind <- bhind[order(N, decreasing = TRUE)]

fwrite(bhind, "../tables/RQ_ind_bh_clustering.tsv", sep = "\t")

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
# [1] reshape2_1.4.4    pheatmap_1.0.12   stringr_1.5.1     magrittr_2.0.3    data.table_1.15.4

# loaded via a namespace (and not attached):
#  [1] RColorBrewer_1.1-3 R6_2.5.1           httpgd_2.0.1       gtable_0.3.4       glue_1.7.0         lifecycle_1.0.4    cli_3.6.2          unigd_0.1.1       
#  [9] scales_1.3.0       grid_4.3.3         vctrs_0.6.5        systemfonts_1.0.6  compiler_4.3.3     plyr_1.8.9         tools_4.3.3        munsell_0.5.1     
# [17] Rcpp_1.0.12        colorspace_2.1-0   rlang_1.1.3        jsonlite_1.8.8     stringi_1.8.3 