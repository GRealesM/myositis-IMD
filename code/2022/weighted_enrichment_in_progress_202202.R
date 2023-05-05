######################################
### Driver SNP rotation-aware GSEA
######################################

# Date: 2022-02-07
# Guillermo Reales

# Load libraries
library(data.table)


# Load data 
# SNP manifest
manifest.translator <- fread("../data/Manifest_build_translator.tsv")
manifest.translator[,pid19:=paste(CHR19, BP19, sep = ":")][, pid38:=paste(CHR38, BP38, sep = ":")]
manifest.translator <- manifest.translator[,c(1,8:9)]
SNP.manifest <- merge(manifest.translator, copy(cupcake::SNP.manifest), by.x = "pid19", by.y = "pid")
rotmat <- as.data.table(cupcake::rot.pca, keep.rownames = TRUE)
names(rotmat)[1] <- "pid19"
names(rotmat) <- gsub("PC", "rot", names(rotmat))
SNP.manifest <- merge(SNP.manifest, rotmat, by="pid19")

# Annotations
otg <- readRDS("../data/imdbasis-driver-snps-otg-annotations.RDS")
otg <- merge(SNP.manifest[,.(pid38, SNPID, ref_a1, ref_a2)], otg, by.x="pid38", by.y="pid")
# Check concordance in alleles
all(otg$ref_a1 == otg$ref_allele)
all(otg$ref_a2 == otg$alt_allele) # FALSE, let's see
nrow(otg[otg$ref_a2 != otg$alt_allele]) # 152 instances of alt allele discordance. We'll weed them out
otg <- otg[ref_a2 == alt_allele]
# Remove redundant columns
otg[, c("source_list", "source_score_list", "wh.vep") := NULL]

# Include rotation values in otg
otg <- merge(otg, manifest.translator[, c("pid19", "pid38")], by = "pid38")
otg <- merge(otg, rotmat, by="pid19")


# Add eQTL data
eqtl <- fread("../data/Causal_eQTLs_full_20210915.tsv") 
# Noticed a bug in the procesing of snp_study and study_tissue. Let's briefly fix it
eqtl[, snp_study:=gsub(pattern = "(_[A-Za-z0-9\\+ -]+)$", "", snp_study, perl = TRUE)][, study_tissue:=gsub(pattern = "(_[A-Za-z0-9\\+ -]+)$", "", study_tissue, perl = TRUE)][,study:=gsub(pattern = "(_[A-Za-z0-9\\+ -]+)$", "", study, perl = TRUE)][,quant_method:=gsub("([A-Za-z0-9_]+)_([a-z]+)", "\\2", study, perl = TRUE)][,study:=gsub(pattern = "(_[a-z]+)$", "", study, perl = TRUE)]

# Get Pathway data
pathways=scan("https://biit.cs.ut.ee/gprofiler/static/gprofiler_full_hsapiens.name.gmt", what="",sep="\n")
ss=strsplit(pathways,"\t")
pathway_ids=sapply(ss, "[[", 1)
pathway_term_names=sapply(ss, "[[", 2)
pathway_genes=lapply(ss, "[", -c(1:2))
names(pathway_genes)=pathway_term_names

## Restrict analysis to GO, REAC, and WP
idx <- grep("GO|WP|REAC", pathway_ids)
pathway_genes <- pathway_genes[idx]
pathway_term_names <- pathway_term_names[idx]
pathway_ids <- pathway_ids[idx]

## Get universe
universe_genes <- pathway_genes |> unlist() |> unique()



# Map genes and weights

weight.mapping <- lapply(c(2,8:9,12), function(pc){
  pc = paste0("rot", pc)
  rs <- SNP.manifest[get(pc) != 0, .(SNPID, get(pc))]
  rs$V2 <- abs(rs$V2)
  otg.pc <- otg[ SNPID %in% rs$SNPID & tss_score == 1, .(SNPID, hgnc_symbol)][ hgnc_symbol != ""] |> unique()  # We'll use gene names and remove nameless genes
  eqtl.pc <- eqtl[  SNPID.basis %in% rs$SNPID, .(SNPID.basis, gene_name)][ gene_name != ""] |> unique() # Genes are likely to be repeated in different tissues, so we remove duplicates
  annot.pc <- rbindlist(list(otg.pc, eqtl.pc), use.names = FALSE) |> unique() # Again remove duplicates in OTG and eQTL
  annot.pc <- merge(annot.pc, rs, by="SNPID")
  names(annot.pc)[3] <- "weight"
  annot.pc2 <- annot.pc[, .(sum_weights = sum(weight)), by = hgnc_symbol]
  list(SNP_gene = annot.pc, gene.weights = annot.pc2)
})
names(weight.mapping) <- paste0("PC", c(2,8:9,12))


# Genes mapped to driver SNPs with their respective sum_weights
pc_winfo <- weight.mapping$PC2$gene.weights

# Function to assign weights to genes in pathways
assign_pweights <- function(pathway_list, pc_winfo){
  pwlist <- lapply(pathway_list, function(pw){
            pw <- list(genes = pw, weights = ifelse(pw %in% pc_winfo$hgnc_symbol, pc_winfo$sum_weights, 0))
    })
}

# Assign "true" weights to genes in pathways
test.pw <- assign_pweights(pathway_genes, pc_winfo)

# Create a vector of length equal to all genes in all pathways, containing our weights for the number of genes mapped + zeroes for the rest.
# Then sample and assign the sampled weights to genes in pathways
null.sumweights <- c(pc_winfo$sum_weights, rep(0, length(universe_genes) - length(pc_winfo$sum_weights)))
null.sumweights <- sample(null.sumweights)
pc_winfo_null <- data.table(hgnc_symbol = universe_genes, sum_weights = null.sumweights)
null.pw <- assign_pweights(pathway_genes, pc_winfo_null)

# Compute t-bar
res <- sapply(null.pw, function(Ag){ sum(Ag$weights)/length(Ag$weights)})
size <- sapply(null.pw, function(Ag){ length(Ag$weights)})

# Plot!
plot(size,res, xlim = c(0, 1000))




