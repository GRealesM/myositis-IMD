#########################################
##                                     ##
##          RUNNING DPMUnc             ##
##                                     ##
#########################################

# Author: Guillermo Reales 
# Date last updated: 2022/08/24


# This script will
# * Run DPMUnc on files supplied by command line.

# NOTE: This script is meant to be run at the HPC via sbatch (slurm_02_RunDPMUnc_myo)

##########################################

# Load packages
library(data.table)
#devtools::install_github("nichollskc/DPMUnc") # install latest version
library(DPMUnc)

# Helper functions
read_matrix <- function(x){
                  ds <- fread(x)
                  ds[, Label:=NULL] # Remove Trait and label columns, keep just the numbers!
                  mat <- as.matrix(ds)
                  colnames(mat) <- NULL
                  return(mat)
}

# Load datasets
args <- commandArgs(trailingOnly = TRUE) # First argument should be the "experiment name", second the seed
exp <- args[1]
seed <- as.numeric(args[2])
delta <- paste0("../data/", exp, "_Delta.tsv")
var <- paste0("../data/", exp, "_Var.tsv")

message("Running DPMUnc on ", exp, " using seed = ", seed, ".")

obsData <- read_matrix(delta)
obsVars <- read_matrix(var)

# Standard run
DPMUnc(obsData, obsVars, saveFileDir = paste0("../data/DPMUnc_results/", exp, "_", seed), seed = seed, kappa0 = 0.01, alpha0 = 2, beta0 = 0.1, nIts = 5000000, scaleData=TRUE)
message("Done!")


########


sessionInfo()
# R version 4.1.3 (2022-03-10)
# Platform: x86_64-pc-linux-gnu (64-bit)
# Running under: Rocky Linux 8.8 (Green Obsidian)

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
# [1] DPMUnc_0.0.0.9000 data.table_1.14.8

# loaded via a namespace (and not attached):
# [1] compiler_4.1.3 cli_3.6.0      Rcpp_1.0.10    jsonlite_1.8.4 rlang_1.0.6   