f###################################
##### Processing eqtl data    ######
####################################

# Author: Guillermo Reales
# Date: 2021-09-14


# We want to include eQTL data in specific blood cells (see below) from eQTLcatalogue. 
# However, downloading and extracting relevant SNPs from source has proved challenging. 
# Thus, I downloaded the full datasets of interest to our HPC storage, and will process them and add them to our causal SNP dataset.
# Note: This script is intended to be run in the HPC

library(data.table)
library(readr)
library(magrittr)
library(ggplot2)
library(cowplot)
# The following packages are Bioconductor's so they need to be installed differently
library(GenomicRanges)
library(biomaRt)
setDTthreads(10)


# Load basic data
ecat <- fread("../../../03-Bases/cell_basis_v3_varimax/Reports/Combined_eQTL_credible_paths.tsv") ## HPC path, might not work elsewhere
snp <- fread("../data/Myositis_IMD_driver_SNPs.tsv")
snp[, pid:=paste(CHR38,BP38, sep=":")] # Choose hg38 coordinates

#snp[,region:=paste0(CHR, ":", BP,"-",BP)][,pid:=paste(CHR, BP, sep=":")]
setnames(snp, "SNPID", "SNPID.basis")

# Filter ecat data by microarray
tissues <- c("B cell","CD16+ monocyte", "CD4+ T cell", "CD8+ T cell", "macrophage", "monocyte", "neutrophil", "NK cell", "platelet", "T cell", "Tfh cell", "Th1 cell", "Th1-17 cell", "Th17 cell", "Th2 cell", "Treg memory", "Treg naive")
ecat.filt <- ecat[tissue_label %in% tissues & (quant_method == "ge" | quant_method == "microarray") & condition_label == "naive",]
ecat.filt[,studylabel:= paste(study, quant_method, tissue_label, sep="_")] # To identify studies later

# Save the info on eQTL dataset used and their corresponding studies as a Suppl. Table
fwrite(ecat.filt, "../tables/SuppTable_eQTL_datasets.tsv", sep ="\t")


# Import local credible sets
crfilespath <- "~/rds/rds-cew54-wallace-share/Data/expr/eqtl-catalogue/credible_sets/"
crfiles <- sapply(strsplit(ecat.filt$ftp_paths_credible, split="\\/"), `[`, 9)
crfilespath  <- paste0(crfilespath, crfiles)

# Download all credible info for selected datasets
get.credible.info.local<- lapply(1:length(crfilespath), function(x){
  
	message("Working on ", crfilespath[x])
	ds <- fread(crfilespath[x])
	ds[,pid:=paste(chromosome, position, sep=":")]
	ds <- ds[pid %in% snp$pid]
	ds[, study:=ecat.filt$studylabel[x]]
})
crlocal <- rbindlist(get.credible.info.local) # 538 obs

# Import local datasets, and keep relevant SNPs
mafilespath <- "~/rds/rds-cew54-wallace-share/Data/expr/eqtl-catalogue/"
mafiles <- sapply(strsplit(ecat.filt$ftp_path, split="\\/"), `[`, 11)
mafilespath  <- paste0(mafilespath, mafiles)

get.eqtl.info.local <- lapply(1:length(mafilespath), function(x){
	message("Working on ", mafilespath[x])
	ds <- fread(mafilespath[x])
	ds[,pid:=paste(chromosome, position, sep=":")]
	ds <- ds[pid %in% snp$pid]
	ds[, study:=ecat.filt$studylabel[x]]
 
	})
malocal  <- rbindlist(get.eqtl.info.local, use.names=TRUE, fill=TRUE)
# Let's do it a bit differently than last time. Let's merge before looking for gene ID names

causal <- merge.data.table(malocal, crlocal, by.x = c("molecular_trait_id", "chromosome", "position", "ref", "alt","study", "pid", "variant"), by.y =  c("molecular_trait_id", "chromosome", "position", "ref", "alt","study", "pid", "variant")) %>% unique()
causal[,tissue_label:= sapply(strsplit(study, split="_"), `[`, 4)][is.na(tissue_label), tissue_label:=sapply(strsplit(study, split="_"), `[`, 3)] 

# Retrieve gene info info from biomart
# List EnsemblIDs
ids <- unique(causal$gene_id)
length(ids) # 103 different genes

# Use biomaRt for the task
mart <- useMart("ensembl")
mart <- useDataset(dataset = "hsapiens_gene_ensembl", mart)
genelist <- getBM(attributes = c("ensembl_gene_id", "external_gene_name", "description", "chromosome_name", "start_position", "end_position"), filters = "ensembl_gene_id", values = ids, mart = mart) %>% as.data.table()

names(genelist) <- c("gene_id","gene_name", "description", "gene_chr", "gene_start", "gene_end")
# Note: not all ensembl gene ids in info do still exist in Ensembl nor have names. Specifically, I found 4 deprecated IDs and 76 nameless genes. To amend this a bit, I downloaded the gene description and their coordinates
info.gn <- merge.data.table(causal, genelist, by= "gene_id", all = TRUE)
# Fix weird coercion and ending
info.gn[ref == "TRUE", ref:="T"][alt == "TRUE", alt:="T"][,pid:=paste(chromosome, position, sep=":")][,rsid:=gsub(pattern = "\\r", rsid, replacement = "")]

# Also, incorporate which SNP drives which component as an additional info
info.gn <- merge.data.table(info.gn, unique(snp[,c("pid", "SNPID.basis","ref_a1","ref_a2", "rot1" , "rot2", "rot3", "rot8","rot9", "rot12", "rot13")]), by = "pid")
info.gn[gene_name == "", gene_name:=gene_id][, snp_study:=paste(rsid,ref,alt, study, tissue_label, sep="_")][, study_tissue:=paste(study, tissue_label, sep="_")]
# Check alleles are ok, they should be!
all(info.gn$ref == info.gn$ref_a1)
all(info.gn$alt == info.gn$ref_a2)


# Finally, let's save the results 
fwrite(info.gn, "Causal_eQTLs_full_20210915.tsv", sep="\t")

