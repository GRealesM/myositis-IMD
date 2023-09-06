#########################################
##                                     ##
##   COMPUTING BHATTACHARYYA DISTANCE  ##
##                                     ##
#########################################

# Author: Guillermo Reales 
# Date last updated: 2023/09/06


# This script will
# * Import projections and metadata files
# * Access the list of reduced datasets (at the HPC)
# * Compute Bhattacharyya distance across each of the selected 66 datasets and 7 key myositis PCs


# NOTE: This script is meant to be run at the HPC

##########################################


### Set paths

# This script is meant to be run at the HPC
setwd("/home/gr440/rds/rds-cew54-basis/Projects/myositis-IMD/code")
bpath <- "/home/gr440/rds/rds-cew54-basis/03-Bases/IMD_basis/"

### Load required packages

library(data.table)
library(magrittr)
#install.packages("fpc")
library(fpc)
library(cupcake)


# Load key data about datasets
ps <- fread("../data/ps2.tsv")
qs <- fread("../data/qs2.tsv")

# Key features
pcs <- paste0("PC", c(1,2,3,8,9,12,13))

# Since files are in hg38 and the basis is in hg19, we use a dirty fix to liftover back to hg19 (easier to do at this stage than in the previous one, due to the huge size of files).
build_dict <- fread(paste0(bpath, "Manifest_build_translator.tsv"))
build_dict[, `:=`(pid38=paste(CHR38,BP38,sep=":"), pid=paste(CHR19,BP19,sep=":"))]
build_dict <- build_dict[,c("pid38", "pid")] 

fs <- paste0(qs$Trait, "-ft.tsv")

covar.m <- lapply(fs, function(x){

        fss <- fread(paste0(bpath, "reduced_datasets/", x))
        # Some checks
        fss <- unique(fss)
        fss[fss == ""] <- NA # Some missing data might pass as empty string. This will fix that	
        fss <- na.omit(fss, cols = c("pid38", "BETA", "SE", "P"))
        dups <- fss$pid38[duplicated(fss$pid38)]

        if(length(dups) > 0){
            dupmessage= "This file has duplicated pids. I removed them prior to projection. You might want to check it."
            message(dupmessage)
            fss <- fss[!pid38 %in% dups] # Remove all duplicated instances, to be safe
        }
        fss <- merge(fss, build_dict)

        pids <- fss$pid
        seb <- fss$SE

        v <- seb * shrinkage[pids] * rot.pca[pids,]
        var.proj  <- t(v) %*% LD[pids,pids] %*% v

        var.proj <- var.proj[pcs, pcs]
        var.proj
})


bh.list <- lapply(seq_along(qs$Label), function(i){

                dd <- data.table(T1 = qs$Label[i], T2 = qs$Label) # Create a data.table to store combinations

                # Data for focus trait (ie. T1)
                rr1  <- ps[ Label == qs$Label[i] & PC %in% pcs, Delta]
                names(rr1) <- ps[ Label == qs$Label[i] & PC %in% pcs, PC]

                cv1 <- as.matrix(covar.m[[i]])
                stopifnot(row.names(cv1) == names(rr1))

                # Then, take each trait in turn, extract necessary data, and compute the Bhattacharyya distance
                bhd <- sapply(seq_along(qs$Label), function(x){
                        rr2 <- ps[ Label == qs$Label[x] & PC %in% pcs, Delta]
                        names(rr2) <- ps[ Label == qs$Label[x] & PC %in% pcs, PC] 
                        cv2 <- as.matrix(covar.m[[x]])
                        stopifnot(row.names(cv1) == names(rr2))
                        bd <- bhattacharyya.dist(mu1 = rr1, mu2 = rr2, Sigma1 = cv1, Sigma2 = cv2 )
                        bd
                    } )

                dd[, bhat.dist:=bhd]
                dd
})  %>% rbindlist


bh.list
fwrite(bh.list, "../data/bh_dist.tsv", sep="\t")
