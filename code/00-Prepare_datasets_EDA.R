#########################################
##                                     ##
##     PREPARING DATASETS AND EDA      ##
##                                     ##
#########################################

# Author: Guillermo Reales 
# Date last updated: 2023/09/06

# Background: This script will prepare the datasets for downstream analyses and prepare some figures.

# This script will
# * Import projections and metadata files
# * Apply QC to projections, apply necessary relabelling removing redundant or non-informative datasets, and focus on IMD.
# * Apply FDR to projections.
# * Remove redundant datasets for improved visualisation
# * Create some figures and supplementary tables, and save datasets for downstream analyses
# * Prepare datasets for DPMUnc


##########################################



## Load required packages

library(data.table)
library(magrittr)
library(ggplot2)
library(cowplot)
library(pheatmap)
library(reshape2)
library(cupcake)

## Load datasets

# Internal NOTE: Before publication, consider doing a previous steps to prepare projection datasets including authorised
# datasets only, in case we need to publish the original projections before QC. This would involve removing PAPS and other
# private datasets, and maybe updating the References in the fields in the original q table (and maybe in the Metadata too?).

# We'll need the SNP manifest. This file contains the 566 SNPs in the features 
SNP.manifest <- cupcake::SNP.manifest

# Meta-data table, containing information about the datasets
m <- fread("../data/Metadata_20230906-v1.tsv")
m <- m[, .(Trait, First_Author, Reference, Trait_ID_2.0, Trait_long, Trait_class, N0, N1, N, Population, Public)]

# Projection QC table, containing information about the projections and their quality
q <- fread("../data/QC_IMD_basis_20230503-v1.tsv")
q <- merge(q, m, by="Trait") # Get metadata together

# Take the chance now to update some references in the q table
q[ First_Author == "Rothwell", Reference:="36580032"]
q[ First_Author == "Lessard", c("First_Author", "Reference"):=list("Khatri", "35896530")][ First_Author == "Wong", Reference:="doi:10.17863/CAM.51022"]


# Projection table
p <- fread("../data/Projection_IMD_basis_20230503-v1.tsv")
p[, Var.Delta:=as.numeric(Var.Delta)]
b <- p[1:169] # Basis traits
p <- p[179:nrow(p)] # projection table without basis traits
p[, z:=NULL]


##########################################

## Apply QC
#
# (1) Remove projections with low SNP match, Immunochip (ex. Rothwell), unauthorised (PAPS) and Neale. Keep IMD only
# (2) Prepare proper labels for selected datasets
# (3) Apply FDR procedure and remove datasets with FDR overall > 1%
# (4) Remove redundant datasets for visualisation and further analyses


### (1) Remove projections with low SNP match, Immunochip (ex. Rothwell), unauthorised (PAPS) and Neale. Keep IMD only.


# Remove datasets with missing overall p-values
summary(q)
q[is.na(q$overall_p)] 
# MDD_Wray_2, Major depression disorder, that failed to project by low SNP match.
# T2D_Gaulton 
q <- q[!is.na(q$overall_p)] 


# Remove datasets with <80% SNP match. This will include datasets using targeted arrays only (ie. ImmunoChip).
# Check SNP match
c(lessthan95 = nrow(q[q$nSNP < nrow(SNP.manifest)*.95,]), lessthan80 = nrow(q[q$nSNP < nrow(SNP.manifest)*.8,]), lessthan50 = nrow(q[q$nSNP < nrow(SNP.manifest)*.5,]))
c(lessthan95 = nrow(q[q$nSNP < nrow(SNP.manifest)*.95,])/nrow(q), lessthan80 = nrow(q[q$nSNP < nrow(SNP.manifest)*.8,])/nrow(q), lessthan50 = nrow(q[q$nSNP < nrow(SNP.manifest)*.5,])/nrow(q))
# 4% of datasets (248) have <80% SNP match

qf <- q[nSNP >= max(nSNP) * 0.8 ]


# Prepare list of datasets to remove, even prior to FDR

fbd <- c("PAPS_Casares_up_1") # Private. Cannot use
uninf <- c("D3_BLOOD_FinnGen_FinnGenR7_1", "D3_IMMUNEMECHANISMNAS_FinnGen_FinnGenR7_1", "ASTE_Ferreira_29083406_1", "RHEC_Ferreira_29083406_1",  "AUTOIMMUNE_FinnGen_FinnGenR7_1","AUTOIMMUNE_NONTHYROID_FinnGen_FinnGenR7_1", "AUTOIMMUNE_NONTHYROID_STRICT_FinnGen_FinnGenR7_1", "D3_IMMUNEMECHANISM_FinnGen_FinnGenR7_1", "ANAPH_SHOCK_ADVER_EFFECT_CORRE_DRUG_MEDICAM_PROPE_ADMINISTE_FinnGen_FinnGenR7_1") # Declared uninformative,  or not case-control
basis.datasets <- c("CD_DeLange_28067908_1", "PSC_Ji_27992413_1", "UC_DeLange_28067908_1", "SLE_Bentham_26502338_1", "PBC_Cordell_26394269_1", "IGAN_Kiryluk_25305756_1", "CEL_Dubois_20190752_1", "MS_IMSGC_21833088_1", "AST_Demenais_29273806_1", "VIT_Jin_27723757_1", "RA_Okada_24390342_1", "LADA_Cousminer_30254083_1", "T1D_Cooper_doi101101120022_1") # Traits used to train the basis, and thus overfitted
neale <- q[First_Author == "Neale", Trait] # Remove Neale, as we have PanUKBB
excl <- c(fbd, uninf, basis.datasets, neale)

qf <- qf[!Trait %in% excl]

# Keep IMD traits only
qf <- qf[Trait_class == "IMD"] # Interested in IMDs for now

# Check traits by class. In this case we have only one class
table(qf$Trait_class)
# IMD 
# 479


### (2) Prepare proper labels for selected datasets

qf[,Label:=Trait_long] %>% .[!grepl("myositis", Trait_long, ignore.case = TRUE), Label:=gsub(" \\(UKBB\\)", "", Label)] %>% .[!grepl("myositis", Trait_long, ignore.case = TRUE) , Label:=gsub(" \\(FinnGen\\)", "", Label)] %>% .[, Label:=gsub(" \\(FG\\)", "", Label)] # Remove UKBB/FinnGen stuff save for myositis datasets
# Relabel some Myositis datasets -- this is the place to modify the labels
qf[ First_Author == "Rothwell", Label:=gsub("Idiopathic Inflammatory Myopathies \\(IIM, Myositis\\)", "IIM", Label)]
qf[ First_Author == "Miller", Label:=gsub("Myositis", "IIM", Label, fixed = TRUE)]
qf[ First_Author %in% c("Rothwell", "Miller"), Label:=paste0(Label, " (", First_Author, ")")]
# Rename PBC and PSC
qf[grepl("sclerosing", Label), Label:= "Primary sclerosing cholangitis"]
qf[grepl("chirrosis", Label), Label:= "Primary biliary cholangitis"]
qf[grepl("MS-disease", Label), Label:="Multiple Sclerosis"] # Fix weird FinnGen MS label

# Create a filtered projection table with new labels too. This will keep selected datasets only

pf <- merge(p, qf[,c("First_Author","Trait", "Trait_ID_2.0", "Trait_long","Trait_class", "Population", "Label")], by = "Trait")


### (3) Apply FDR procedure and remove datasets with FDR overall > 1%

# Apply 1% FDR correction to overall p for all remaining datasets
qf[, FDR.overall := p.adjust(overall_p, method = "BH"), by="Trait_class"] # Only IMD, so trait class shouldn't matter

# Apply FDR.PC to projections
pf[, FDR.PC:=p.adjust(P, method = "BH"), by = c("PC", "Trait_class")][, stars:=ifelse(!is.na(FDR.PC) & FDR.PC<0.05,"○","")][ FDR.PC < 0.01 , stars:="●"]

##############################

# Create Supplementary table (Myositis)

## Now we have added the FDR information to the filtered (qf, pf) tables, we can extract them to create Supplementary tables

# At this point, prepare Supplementary Table with coloc info
qmyo <- qf[ First_Author %in% c("Miller", "Rothwell")]
qmyo <- qmyo[, .(Label, First_Author, Reference, N0, N1, N, FDR.overall )][order(First_Author, Label)]
qmyo
# fwrite(qmyo, "../tables/ST_Myo_info.tsv", sep="\t")

#############################


# Remove non-overall-significant datasets

# Keep significant datasets only in qs
qs <- qf[FDR.overall < 0.01]
nrow(qs) / nrow(qf)
# 57.5% significant

# Do the same for ps
ps <- pf[Trait %in% qs$Trait ]
ps[, PC:=factor(PC, levels = paste0("PC", 1:13))]


# Check how many are significant
tsig <- copy(qf)
tsig[, sig.overall:=ifelse(FDR.overall < 0.01, "Y", "N")]
table(tsig$sig.overall)
#   N   Y 
# 202 274 

# Which myositis datasets didn't make the cut?
qf[grepl("myositis", Trait_long, ignore.case = TRUE) & !Trait %in% qs$Trait, .(First_Author, Trait_ID_2.0, Trait_long)]
# 6 datasets. These include:
#    First_Author       Trait_ID_2.0                                 Trait_long 
# 1:      PanUKBB         20002_1322                   Myositis/Myopathy (UKBB) 
# 2:     Rothwell                IBM                    Inclusion Body Myositis 
# 3:      FinnGen   M13_DERMATOMYOTH            Other dermatomyositis (FinnGen) 
# 4:      FinnGen M13_DERMATOPOLYNAS Dermatopolymyositis, unspecified (FinnGen) 
# 5:       Sakaue                 PM                               Polymyositis 
# 6:      PanUKBB              ph770    Myalgia and myositis unspecified (UKBB) 


# Which features are relevant for myositis?
ps[grepl("myositis", Trait_long, ignore.case = TRUE) & FDR.PC < 0.01, .(PC, First_Author, Trait_long, FDR.PC, stars)][order(PC)]

# We have some significant Myositis for PC1, 2, 3, 8, 9, 12, and 13


### (4) Remove redundant datasets for visualisation and further analyses

# Let's try to remove redundant datasets in a less manual way
qs2 <- copy(qs)

# There are a number of duplicated diseases that we can revisit
# CD/UC/IBD
qs2[grepl("Crohn|Chron's|colitis|UC |Ulcerative|bowel|IBD", Label, ignore.case = TRUE), .(Trait, First_Author, Label, N0,N1,N, nSNP, overall_p, mscomp, FDR.overall)][order(N1, decreasing = TRUE)]
tr <- c(setdiff(qs2[grepl("Crohn|Chron's|colitis|UC |Ulcerative|bowel|IBD", Label, ignore.case = TRUE), Trait], c("K11_CROHN_FinnGen_FinnGenR7_1", "K11_ULCER_FinnGen_FinnGenR7_1", "K11_IBD_FinnGen_FinnGenR7_1"))) # We selected Liu, which has large N but we suspect overlaps with DeLange, so we'll choose 3 FinnGen representatives instead.

# T1D
qs2[grepl("Type 1 diabetes|Type1", Label, ignore.case = TRUE),  .(Trait, First_Author, Label, N0,N1,N, nSNP, overall_p, mscomp, FDR.overall)][order(N1, decreasing = TRUE)]
tr <- c(tr, setdiff(qs2[grepl("Type 1 diabetes|Type1", Label, ignore.case = TRUE), Trait], "E4_DM1_FinnGen_FinnGenR7_1")) 

# Asthma
qs2[grepl("Asthma", Label, ignore.case = TRUE), .(Trait, First_Author, Label, N0,N1,N, nSNP, overall_p, mscomp, FDR.overall)][order(N1, decreasing = TRUE)]
tr <- c(tr, setdiff(qs2[grepl("Asthma", Label, ignore.case = TRUE), Trait], c("ph495_PanUKBB_PanUKBBR2_1", "ASTAO_Ferreira_30929738_1", "ASTCO_Ferreira_30929738_1")))

# Celiac disease
qs2[grepl("celiac|coeliac", Label, ignore.case = TRUE)][order(N1, decreasing = TRUE)]
tr <- c(tr, setdiff(qs2[grepl("celiac|coeliac", Label, ignore.case = TRUE), Trait], c("K11_COELIAC_FinnGen_FinnGenR7_1")))

# (Rheumatoid) Arthritis
qs2[grepl("Arthritis", Label, ignore.case = TRUE),.(Trait, First_Author, Label, N0,N1,N, nSNP, overall_p, mscomp, FDR.overall, Population)][order(N1, decreasing = TRUE)]
tr <- c(tr, setdiff(qs2[grepl("Arthritis", Label, ignore.case = TRUE), Trait], c("M13_RHEUMA_FinnGen_FinnGenR7_1", "JIA_LopezIsac_33106285_1"))) # Take RA and JIA at this step
tr <- c(tr, "M13_JUVERHEU_FinnGen_FinnGenR7_1") #  Remove juvenile rheuma, as it's essentially the same as JIA

# Multiple sclerosis
qs2[grepl("multiple sclerosis", Label, ignore.case = TRUE),.(Trait, First_Author, Label, N0,N1,N, nSNP, overall_p, mscomp, FDR.overall, Population)][order(N1, decreasing = TRUE)]
tr <- c(tr, setdiff(qs2[grepl("multiple sclerosis", Label, ignore.case = TRUE), Trait], c("G6_MS_FinnGen_FinnGenR7_1"))) # largest, non-IMSGC dataset

# PBC
qs2[grepl("biliary", Label, ignore.case = TRUE)][order(N1, decreasing = TRUE)]
tr <- c(tr, setdiff(qs2[grepl("biliary", Label, ignore.case = TRUE), Trait], c("CHIRBIL_PRIM_FinnGen_FinnGenR7_1")))

# PSC
qs2[grepl("sclerosing", Label, ignore.case = TRUE)][order(N1, decreasing = TRUE)]
tr <- c(tr, setdiff(qs2[grepl("sclerosing", Label, ignore.case = TRUE), Trait], c("K11_CHOLANGI_FinnGen_FinnGenR7_1")))

# Lupus
qs2[grepl("Lupus|SLE", Label, ignore.case = TRUE)][order(N1, decreasing = TRUE)]
tr <- c(tr, setdiff(qs2[grepl("Lupus|SLE", Label, ignore.case = TRUE), Trait], c("M13_SLE_FinnGen_FinnGenR7_1"))) # Julia has Bentham individuals. Replace by FinnGen

# Vitiligo
qs2[grepl("vitiligo", Label, ignore.case = TRUE)][order(N1, decreasing = TRUE)]
tr <- c(tr, setdiff(qs2[grepl("vitiligo", Label, ignore.case = TRUE), Trait], c("L12_VITILIGO_FinnGen_FinnGenR7_1")))

# COPD
qs2[grepl("COPD|chronic obstructive", Label, ignore.case = TRUE)][order(N1, decreasing = TRUE)]
tr <- c(tr, setdiff(qs2[grepl("COPD|chronic obstructive", Label, ignore.case = TRUE), Trait], "J10_COPD"))

# Hypothyroidism
qs2[grepl("Hypothyroidism", Label, ignore.case = TRUE)][order(N1, decreasing = TRUE)]
tr <- c(tr, setdiff(qs2[grepl("Hypothyroidism", Label, ignore.case = TRUE), Trait], c("E4_HYTHY_AI_STRICT_FinnGen_FinnGenR7_1"))) # Take strict definition

# Hyperthyroidism
qs2[grepl("hyperthyroidism|graves", Label, ignore.case = TRUE)][order(N1, decreasing = TRUE)]
tr <- c(tr, setdiff(qs2[grepl("hyperthyroidism|graves", Label, ignore.case = TRUE), Trait], c("20002_1225_PanUKBB_PanUKBBR2_1","E4_GRAVES_STRICT_FinnGen_FinnGenR7_1"))) # UKBB has larger sample size. Take FinnGen Graves' too

# Thyroiditis
qs2[grepl("thyroiditis", Label, ignore.case = TRUE)][order(N1, decreasing = TRUE)]
tr <- c(tr, setdiff(qs2[grepl("thyroiditis", Label, ignore.case = TRUE) , Trait], c("E4_THYROIDITAUTOIM_FinnGen_FinnGenR7_1"))) # Choose, Autoimmune thyroiditis, even if it has smaller sample size than general thyroiditis, since we want to capture the immune component

# Psoriasis
qs2[grepl("Psoria", Label, ignore.case = TRUE)][order(N1, decreasing = TRUE)]
tr <- c(tr, setdiff(qs2[grepl("Psoria", Label, ignore.case = TRUE), Trait], c("L12_PSORIASIS_FinnGen_FinnGenR7_1", "L12_PSORI_ARTHRO_FinnGen_FinnGenR7_1")))


# Rhinitis
qs2[grepl("rhinitis|allerg", Label, ignore.case = TRUE)][order(N1, decreasing = TRUE)]
tr <- c(tr, setdiff(qs2[grepl("rhinitis|allerg", Label, ignore.case = TRUE), Trait], c("20002_1387_PanUKBB_PanUKBBR2_1")))

# Dermatitis
qs2[grepl("Dermatitis", Label, ignore.case = TRUE)][order(N1, decreasing = TRUE)]
tr <- c(tr, setdiff(qs2[grepl("Dermatitis", Label, ignore.case = TRUE), Trait], c("ATD_Paternoster_26482879_1")))

# Sjogren
qs2[grepl("Sjogren|Sjögren|Sicca", Label, ignore.case = TRUE)][order(N1, decreasing = TRUE)]
tr <- c(tr, setdiff(qs2[grepl("Sjogren|Sjögren|Sicca", Label, ignore.case = TRUE), Trait], c("SJOS_Lessard_up_1")))

# Nasal polyps
qs2[grepl("polyp", Label, ignore.case = TRUE)][order(N1, decreasing = TRUE)]
tr <- c(tr, setdiff(qs2[grepl("polyp", Label, ignore.case = TRUE), Trait], c("J10_NASALPOLYP_FinnGen_FinnGenR7_1")))

# Systemic sclerosis
qs2[grepl("systemic sclerosis", Label, ignore.case = TRUE)][order(N1, decreasing = TRUE)]
tr <- c(tr, setdiff(qs2[grepl("systemic sclerosis", Label, ignore.case = TRUE), Trait], c("SSC_LopezIsac_31672989_1")))

# Gout and ankylosing
qs2[grepl("gout|ankylosing", Label, ignore.case = TRUE)][order(N1, decreasing = TRUE)]
tr <- c(tr, setdiff(qs2[grepl("gout|ankylosing", Label, ignore.case = TRUE), Trait], c("GOUT_FinnGen_FinnGenR7_1", "M13_ANKYLOSPON_FinnGen_FinnGenR7_1")))

# Myastenia gravis
qs2[Label == "Myastenia Gravis"][order(N1, decreasing = TRUE)]
tr <- c(tr, setdiff(qs2[Label == "Myastenia Gravis" , Trait], c("MYG_Chia_35074870_1")))

# Rosacea
qs2[grepl("rosacea", Label, ignore.case = TRUE)][order(N1, decreasing = TRUE)]
tr <- c(tr, setdiff(qs2[grepl("rosacea", Label, ignore.case = TRUE) , Trait], c("L12_ROSACEA_FinnGen_FinnGenR7_1")))

# Addison's disease
qs2[grepl("addison", Label, ignore.case = TRUE)][order(N1, decreasing = TRUE)]
tr <- c(tr, setdiff(qs2[grepl("addison", Label, ignore.case = TRUE) , Trait], c("AAD_Eriksson_33574239_1")))

# We'll also remove all myositis from FinnGen, as they're not informative
tr <- c(tr, "DERMATOPOLY_FG_FinnGen_FinnGenR7_1", "M13_MYOSITIS_FinnGen_FinnGenR7_1", "M13_POLYMYO_FinnGen_FinnGenR7_1", "M13_DERMATOPOLY_FinnGen_FinnGenR7_1")



qs2 <- qs2[!Trait %in% tr] # First pass
length(unique(qs2$Label))
# 66  

# Check for duplicates, if there are still
dp <- qs2[duplicated(Label), unique(Label)]
qsdup <- qs2[Label %in% dp]
qsdup # No duplicates

# Apply this to projections as well
ttk2 <- qs2$Trait # Traits to keep
ps2 <- ps[Trait %in% ttk2]

ps2[, PC:=factor(PC, levels = paste0("PC", 1:13))]

fwrite(ps2, "../data/ps2.tsv", sep="\t")
fwrite(qs2, "../data/qs2.tsv", sep="\t")

# Prepare supplementary tables for ALL datasets and projections

# Save filtered datasets - not sure if we'll show these

aq <- copy(qf)
ap <- copy(pf)

aq[, sig.overall:=ifelse(Trait %in% qs$Trait, "Yes", "No")][, in.selection:=ifelse(Trait %in% qs2$Trait, "Yes", "No")]
ap[, sig.overall:=ifelse(Trait %in% qs$Trait, "Yes", "No")][, in.selection:=ifelse(Trait %in% qs2$Trait, "Yes", "No")]

aq <- aq[,.(Trait, Label, First_Author, Reference, N0, N1,N, Population, nSNP, mscomp, overall_p, FDR.overall, sig.overall, in.selection)]
ap <- ap[,.(Trait, Label, First_Author, Population, PC, Delta, Var.Delta, P, FDR.PC, stars, sig.overall, in.selection)]

fwrite(aq, "../tables/ST_all_datasets.tsv", sep ="\t")
fwrite(ap, "../tables/ST_all_projections.tsv", sep ="\t")

##########################################

### Prepare some basic figures

## Figure 1 -- A heatmap of myositis projections

pspm <- ps2[grepl("Miller|Rothwell", Label, ignore.case = TRUE)]

PCorder <- paste0("PC", 1:13)
hmcol <- rev(colorRampPalette(c("#67001F", "#B2182B", "#D6604D", "#F4A582", "#FDDBC7", "#F7F7F7", "#D1E5F0", "#92C5DE", "#4393C3", "#2166AC", "#053061"))(100))
Mmp <- acast(pspm[,c("PC", "Label", "Delta")], Label ~ PC) # PC, Trait, and Delta columns only
Mmp.stars <- acast(pspm[,c("PC","Label","stars")], Label ~ PC)
Mmp <- Mmp[,PCorder]
Mmp.stars <- Mmp.stars[,PCorder]
range <- max(abs(Mmp))

# Create heatmap
Mphm <- pheatmap(Mmp,  breaks = seq(-range, range, length.out = 100), 
                 cluster_cols = FALSE, cluster_rows = FALSE, display_numbers = Mmp.stars,
                 fontsize_row = 8.4, fontsize_number = 11, color = hmcol, 
                 annotation_names_row = FALSE, annotation_legend = TRUE)
Mphm

# Save Figure 1
ggsave("../figures/Myositis_allsources_heatmap.png", Mphm, width = 6, height = 2.5, bg="white")
ggsave("../figures/Myositis_allsources_heatmap.svg", Mphm, width = 6, height = 2.5, bg="white")

system("sed -i \"s/ textLength=\'[^\']*\'//\" ../figures/Myositis_allsources_heatmap.svg") # Trick to make the svg file text be more easily editable

###

## Internal figure -- Heatmap of all 66 projections

# Remove all  datasets without at least one FDR 1% significant PC
#a1s <- ps2[FDR.PC < 0.05 | grepl("myositis|IIM", Label, ignore.case = TRUE), Trait] %>% unique 

#ps2s <- ps2[Trait %in% a1s]

Map <- acast(ps2[,c("PC", "Label", "Delta")], Label ~ PC) # PC, Trait, and Delta columns only
Map.stars <- acast(ps2[,c("PC","Label","stars")], Label ~ PC)
Map <- Map[,PCorder]
Map.stars <- Map.stars[,PCorder]
range <- max(abs(Map))

# We have many datasets, so let's highlight myositis
# From https://github.com/raivokolde/pheatmap/issues/48
# use this function to make row or column names bold
# parameters:
#   mat: the matrix passed to pheatmap
#   rc_fun: either rownames or colnames
#   rc_names: vector of names that should appear in boldface
make_bold_names <- function(mat, rc_fun, rc_names) {
  bold_names <- rc_fun(mat)
  ids <- rc_names %>% match(rc_fun(mat))
  ids %>%
    purrr::walk(
      function(i)
        bold_names[i] <<-
        bquote(bold(.(rc_fun(mat)[i]))) %>%
        as.expression()
    )
  bold_names
}

papsb <- grep("myositis|IIM", rownames(Map), value = TRUE, ignore.case = TRUE)

# Create heatmap
Mahm <- pheatmap(Map,  breaks = seq(-range, range, length.out = 100), 
                 cluster_cols = FALSE, display_numbers = Map.stars,
                 fontsize_row = 8.4, fontsize_number = 11, color = hmcol, 
                 annotation_names_row = FALSE, annotation_legend = TRUE,
                 labels_row = make_bold_names(Map, rownames, papsb))

Mahm
ggsave("../figures/Myositis_IMD_heatmap.png", Mahm, width = 8, height = 14, bg="white")

##########################################


## Prepare inputs for DPMUnc

# At this point, we prepare the data for DPMUnc using significant PCs for myositis 
ps2[grepl("myositis|IIM", Label, ignore.case = TRUE) & FDR.PC < 0.01, .(PC, Label, stars)][, .N, by = PC]

dpmunc.ds <- ps2[PC %in% paste0("PC", c(1,2,3,8,9,12,13)), .(PC, Delta, Var.Delta, Label)]
dpmunc.delta <- reshape(dpmunc.ds[, .(PC, Delta, Label)], idvar="Label", timevar = "PC", direction = "wide")
dpmunc.var   <- reshape(dpmunc.ds[, .(PC, Var.Delta, Label)], idvar="Label", timevar = "PC", direction = "wide")

fwrite(dpmunc.delta, "../data/Myo_7PC_Delta.tsv", sep = "\t")
fwrite(dpmunc.var, "../data/Myo_7PC_Var.tsv", sep = "\t")



##########################################

sessionInfo()
# R version 4.3.1 (2023-06-16)
# Platform: x86_64-pc-linux-gnu (64-bit)
# Running under: Pop!_OS 22.04 LTS
# 
# Matrix products: default
# BLAS:   /usr/lib/x86_64-linux-gnu/blas/libblas.so.3.10.0 
# LAPACK: /usr/lib/x86_64-linux-gnu/lapack/liblapack.so.3.10.0
# 
# locale:
#   [1] LC_CTYPE=en_GB.UTF-8       LC_NUMERIC=C               LC_TIME=en_GB.UTF-8        LC_COLLATE=en_GB.UTF-8     LC_MONETARY=en_GB.UTF-8    LC_MESSAGES=en_GB.UTF-8   
# [7] LC_PAPER=en_GB.UTF-8       LC_NAME=C                  LC_ADDRESS=C               LC_TELEPHONE=C             LC_MEASUREMENT=en_GB.UTF-8 LC_IDENTIFICATION=C       
# 
# time zone: Europe/London
# tzcode source: system (glibc)
# 
# attached base packages:
#   [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] cupcake_0.1.0.0   reshape2_1.4.4    pheatmap_1.0.12   cowplot_1.1.1     ggplot2_3.4.2     magrittr_2.0.3    data.table_1.14.8
# 
# loaded via a namespace (and not attached):
#   [1] tidyr_1.3.0         utf8_1.2.3          generics_0.1.3      stringi_1.7.12      lattice_0.21-8      digest_0.6.33       evaluate_0.21       grid_4.3.1         
# [9] RColorBrewer_1.1-3  fastmap_1.1.1       plyr_1.8.8          Matrix_1.6-0        survival_3.5-5      purrr_1.0.1         fansi_1.0.4         scales_1.2.1       
# [17] snpStats_1.50.0     textshaping_0.3.6   cli_3.6.1           rlang_1.1.1         munsell_0.5.0       splines_4.3.1       withr_2.5.0         yaml_2.3.7         
# [25] tools_4.3.1         dplyr_1.1.2         colorspace_2.1-0    BiocGenerics_0.46.0 vctrs_0.6.3         R6_2.5.1            lifecycle_1.0.3     zlibbioc_1.46.0    
# [33] stringr_1.5.0       ragg_1.2.5          pkgconfig_2.0.3     pillar_1.9.0        gtable_0.3.3        glue_1.6.2          Rcpp_1.0.11         systemfonts_1.0.4  
# [41] xfun_0.39           tibble_3.2.1        tidyselect_1.2.0    rstudioapi_0.15.0   knitr_1.43          htmltools_0.5.5     svglite_2.1.1       rmarkdown_2.23     
# [49] compiler_4.3.1  