#############################################################################
#### Looking for potential novel associations in Myositis driver SNPs   #####
#############################################################################

# Author: Guillermo Reales
# Date: 2022-02-11

# Background: We want to look for possible novel associations among driver SNPs of significant components. Given significance (FDR <1%) for a given component, are driver SNPs significant at FDR 5% amongst other driver SNPs?
# This script is meant to be run at the HPC, where the reduced summary statistics live
path_to_reduced = "../../../03-Bases/IMD_basis/reduced_datasets"

# Load libraries
library(data.table)
library(magrittr)


# Prepare data

manifest.translator <- fread("../data/Manifest_build_translator.tsv")
manifest.translator[,pid19:=paste(CHR19, BP19, sep = ":")][, pid38:=paste(CHR38, BP38, sep = ":")]
manifest.translator <- manifest.translator[,c(1,8:9)]
SNP.manifest <- merge(manifest.translator, copy(cupcake::SNP.manifest), by.x = "pid19", by.y = "pid")

rotmat = cupcake::rot.pca

# Quickly check significance for traits and PCs (see Myositis_on_IMDbasis_Report_202202.html for
# more details).
#
# Miller's Myositis and significant PCs (FDR 1%)
# 
# Myositis: PC1, 12, 13
# PM: PC13
# DM: PC1
# JDM: PC1
#
# Rothwell's myositis and significant PCs (FDR 1%)
#
# Myositis (IIM): PC1, 2, 3, 8, 9, 12
# PM: PC2, 3, 9
# Jo1+: PC3, 9
# DM: PC1, 12
# JDM: PC12

files = dir(path_to_reduced, "Miller|Rothwell")
files = files[ files != "IBM_Rothwell_up_1-ft.tsv"] # IBM not significant
files
# [1] "DMY_Miller_26291516_1-ft.tsv" "DMY_Rothwell_up_1-ft.tsv"    
# [3] "IIM_Rothwell_up_1-ft.tsv"     "JDM_Miller_26291516_1-ft.tsv"
# [5] "JDM_Rothwell_up_1-ft.tsv"     "JO1M_Rothwell_up_1-ft.tsv"   
# [7] "MYO_Miller_26291516_1-ft.tsv" "PM_Miller_26291516_1-ft.tsv" 
# [9] "PM_Rothwell_up_1-ft.tsv"     

# Assign PCs to files
pcs  <- list(1, c(1, 12), # DMY
	     c(1,2,3,8,9,12), #IIM
	     1, c(1,12), # JDM
	     c(3,9), # Jo1+
	     c(1,12,13), # Miller's myositis
     	     13, c(2,3,9)) # PM
files.pcs  <- mapply(function(X, Y) { list(list(file =X, pcs =Y))}, X = files, Y = pcs)

fdr.all.files <- lapply(files.pcs, function(j){

	tt  <- fread(paste0(path_to_reduced,"/", j$file))
	fts <- j$pcs

	fdr.per.pc <- lapply(fts, function(i){
	
		r1  <- data.table(pid19 = names(rotmat[rotmat[,i] !=0, i]), rot = rotmat[rotmat[, i] !=0, i])
		mt  <- merge(SNP.manifest, r1, by="pid19")
		tt2 <- merge(mt[, .(SNPID, pid38, rot)], tt, by="pid38")
		tt2[, FDR:=p.adjust(P, method="BH")][, PC:=paste0("PC", i)]
		tt2 <- tt2[ FDR < 0.05]
		tt2
	
		     })
	fdr.per.pc
	fdr.per.pc <- rbindlist(fdr.per.pc)
	fdr.per.pc[, file:= j$file ]

	     })
fdr.all.files  <- rbindlist(fdr.all.files)
setnames(fdr.all.files, "SNPID.x", "SNPID")
fdr.all.files[, SNPID.y:= NULL]

# Remove SNPs with opposing sign
fdr.all.files <- fdr.all.files[ sign(rot) == sign(BETA)] 

# Rename stuff
fdr.all.files[, file:= sapply(strsplit(file, "_u|_2"), `[`, 1) %>% gsub("_", " / ", .) %>% gsub("IIM|MYO", "Myositis",.) 

# Save
fwrite(fdr.all.files, "../tables/MainTable_driverSNP_validation.tsv", sep="\t")



