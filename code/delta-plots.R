# This script will make internal delta plots to allow us see the distribution of various traits across the features
#
# Unlike with the main analyses, we'll include all sort of traits, rather than just IMD, as our goal is to better interpret
# the features with the latest dataset collection.
#
# Author: Guillermo Reales
# Date: 2023-08-31

# Load packages
library(data.table)
library(magrittr)
library(ggplot2)
library(cowplot)
library(reshape2)


# Load files

# We'll need the SNP manifest
SNP.manifest <- cupcake::SNP.manifest

# Meta-data table, containing information about the datasets
m <- fread("../data/Metadata_20230608-v1.tsv")
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
b[, Trait:= paste0(Trait, " (Basis)")]
p <- p[179:nrow(p)] # projection table without basis traits
p[, z:=NULL]


#######################################################

# Now it's time to do our proper EDA on the projections. This will involve:
#
# (1) Remove projections with low SNP match, Immunochip (ex. Rothwell), unauthorised (PAPS) and Neale. Keep IMD only
# (2) Prepare proper labels for selected datasets
# (3) Apply FDR procedure and remove datasets with FDR overall > 1%
# (4) Remove redundant datasets for visualisation and further analyses


### (1) Remove projections with low SNP match, Immunochip (ex. Rothwell), unauthorised (PAPS) and Neale. Keep IMD only.


# Remove unauthorised projection (PAPS)
p <- p[Trait != "PAPS_Casares_up_1"]
q <- q[Trait != "PAPS_Casares_up_1"] # Remove PAPS


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

# Remove datasets used to build the IMD basis, which are overfitted
basis.datasets <- c("CD_DeLange_28067908_1", "PSC_Ji_27992413_1", "UC_DeLange_28067908_1", "SLE_Bentham_26502338_1", "PBC_Cordell_26394269_1", "IGAN_Kiryluk_25305756_1", "CEL_Dubois_20190752_1", "MS_IMSGC_21833088_1", "AST_Demenais_29273806_1", "VIT_Jin_27723757_1", "RA_Okada_24390342_1", "LADA_Cousminer_30254083_1", "T1D_Cooper_doi101101120022_1")

qf <- qf[!Trait %in% basis.datasets]

# Remove Neale, as we already have PanUKBBR2
qf <- qf[First_Author != "Neale"]

# Keep IMD traits only
#qf <- qf[Trait_class == "IMD"] # Interested in IMDs for now

# Check traits by class. In this case we have only one class
table(qf$Trait_class)
# IMD 
# 524 


########################################


# (2) Prepare proper labels for selected datasets

qf[,Label:=Trait_long] %>% .[!grepl("myositis", Trait_long, ignore.case = TRUE), Label:=gsub(" \\(UKBB\\)", "", Label)] %>% .[!grepl("myositis", Trait_long, ignore.case = TRUE) , Label:=gsub(" \\(FinnGen\\)", "", Label)] %>% .[, Label:=gsub(" \\(FG\\)", "", Label)] # Remove UKBB/FinnGen stuff save for myositis datasets
# Relabel some Myositis datasets -- this is the place to modify the labels
qf[ First_Author == "Rothwell", Label:=gsub("Idiopathic Inflammatory Myopathies \\(IIM, Myositis\\)", "IIM", Label)]
qf[ First_Author == "Miller", Label:=gsub("Myositis", "IIM", Label, fixed = TRUE)]
qf[ First_Author %in% c("Rothwell", "Miller"), Label:=paste0(Label, " (", First_Author, ")")]
# Rename PBC and PSC
qf[grepl("sclerosing", Label), Label:= "Primary sclerosing cholangitis"]
qf[grepl("chirrosis", Label), Label:= "Primary biliary cholangitis"]
qf[grepl("MS-disease", Label), Label:="Multiple Sclerosis"] # Fix weird FinnGen MS label


# Create a filtered projection table with new labels too

pf <- merge(p, qf[,c("First_Author","Trait", "Trait_ID_2.0", "Trait_long","Trait_class", "Population", "Label")], by = "Trait")


###########################################

### (3) Apply FDR procedure and remove datasets with FDR overall > 1%

# Apply 1% FDR correction to overall p for all remaining datasets
qf[, FDR.overall := p.adjust(overall_p, method = "BH"), by="Trait_class"] # Only IMD, so trait class shouldn't matter

# Apply FDR.PC to projections
pf[, FDR.PC:=p.adjust(P, method = "BH"), by = c("PC", "Trait_class")][, stars:=ifelse(!is.na(FDR.PC) & FDR.PC<0.05,"○","")][ FDR.PC < 0.01 , stars:="●"]
# Filter pf by overall significant traits


####

## Create significant-only datasets

# Keep significant datasets only in qs
qs <- qf[FDR.overall < 0.01]
nrow(qs) / nrow(qf)
# 58.6% significant

# Do the same for ps
ps <- pf[Trait %in% qs$Trait ]
ps[, PC:=factor(PC, levels = paste0("PC", 1:13))]


###########################################

# (4) Remove redundant datasets for visualisation and further analyses

# Let's try to remove redundant datasets in a less manual way
qs2 <- copy(qs)

# There are a number of duplicated diseases that we can revisit
# CD
qs2[grepl("Crohn|Chron", Label)][order(N1, decreasing = TRUE)]
tr <- c(setdiff(qs2[grepl("Crohn|Chron", Label), Trait], "CD_Liu_26192919_1")) # We select CD Liu, exclude the rest

# UC
qs2[grepl("colitis|UC |Ulcerative", Label, ignore.case = TRUE)][order(N1, decreasing = TRUE)]
tr <- c(tr, setdiff(qs2[grepl("colitis|UC |Ulcerative", Label), Trait], "UC_Liu_26192919_1")) # We select UC Liu, exclude the rest

# IBD
qs2[grepl("bowel|IBD ", Label, ignore.case = TRUE)][order(N1, decreasing = TRUE)]
tr <- c(tr, setdiff(qs2[grepl("bowel|IBD ", Label, ignore.case = TRUE), Trait], "IBD_Liu_26192919_1")) #  There's a larger, DeLange IBD, but we stick with Liu for consistency

# T1D
qs2[grepl("Type 1 diabetes|Type1", Label, ignore.case = TRUE)][order(N1, decreasing = TRUE)]
tr <- c(tr, setdiff(qs2[grepl("Type 1 diabetes|Type1", Label, ignore.case = TRUE), Trait], "T1D_Chiou_34012112_1"))

# Asthma
qs2[grepl("Asthma", Label, ignore.case = TRUE)][order(N1, decreasing = TRUE)]
tr <- c(tr, setdiff(qs2[grepl("Asthma", Label, ignore.case = TRUE), Trait], c("ph495_PanUKBB_PanUKBBR2_1", "ASTAO_Ferreira_30929738_1", "ASTCO_Ferreira_30929738_1")))

# COPD
qs2[grepl("COPD|chronic obstructive", Label, ignore.case = TRUE)][order(N1, decreasing = TRUE)]
tr <- c(tr, setdiff(qs2[grepl("COPD|chronic obstructive", Label, ignore.case = TRUE), Trait], "J10_COPD"))

# (Rheumatoid) Arthritis
qs2[grepl("Arthritis", Label, ignore.case = TRUE)][order(N1, decreasing = TRUE)]
tr <- c(tr, setdiff(qs2[grepl("Arthritis", Label, ignore.case = TRUE), Trait], c("RA_Ishigaki_36333501_3", "JIA_LopezIsac_33106285_1")))

# Hypothyroidism
qs2[grepl("Hypothyroidism", Label, ignore.case = TRUE)][order(N1, decreasing = TRUE)]
tr <- c(tr, setdiff(qs2[grepl("Hypothyroidism", Label, ignore.case = TRUE), Trait], c("HYPOTHYROIDISM_FinnGen_FinnGenR7_1")))

# Hyperthyroidism
qs2[grepl("hyperthyroidism|graves", Label, ignore.case = TRUE)][order(N1, decreasing = TRUE)]
tr <- c(tr, setdiff(qs2[grepl("hyperthyroidism|graves", Label, ignore.case = TRUE), Trait], c("20002_1225_PanUKBB_PanUKBBR2_1")))

# Thyroiditis
qs2[grepl("thyroiditis", Label, ignore.case = TRUE)][order(N1, decreasing = TRUE)]
tr <- c(tr, setdiff(qs2[grepl("thyroiditis", Label, ignore.case = TRUE) , Trait], c("E4_THYROIDITAUTOIM_FinnGen_FinnGenR7_1")))

# Psoriasis
qs2[grepl("Psoria", Label, ignore.case = TRUE)][order(N1, decreasing = TRUE)]
tr <- c(tr, setdiff(qs2[grepl("Psoria", Label, ignore.case = TRUE), Trait], c("L12_PSORIASIS_FinnGen_FinnGenR7_1", "L12_PSORI_ARTHRO_FinnGen_FinnGenR7_1")))

# Graves' disease
qs2[grepl("Graves", Label, ignore.case = TRUE)][order(N1, decreasing = TRUE)]
tr <- c(tr, setdiff(qs2[grepl("Graves", Label, ignore.case = TRUE), Trait], c("E4_GRAVES_STRICT_FinnGen_FinnGenR7_1")))

# Rhinitis
qs2[grepl("rhinitis", Label, ignore.case = TRUE)][order(N1, decreasing = TRUE)]
tr <- c(tr, setdiff(qs2[grepl("rhinitis", Label, ignore.case = TRUE), Trait], c("20002_1387_PanUKBB_PanUKBBR2_1")))

# Diabetic related stuff
qs2[grepl("Diabetic|Diabetes[, ]", Label, ignore.case = TRUE)][order(N1, decreasing = TRUE)]
tr <- c(tr, setdiff(qs2[grepl("Diabetic|Diabetes[, ]", Label, ignore.case = TRUE), Trait], c("DM_RETINOPATHY_FinnGen_FinnGenR7_1", "DM_NEPHROPATHY_FinnGen_FinnGenR7_1")))

# Dermatitis
qs2[grepl("Dermatitis", Label, ignore.case = TRUE)][order(N1, decreasing = TRUE)]
tr <- c(tr, setdiff(qs2[grepl("Dermatitis", Label, ignore.case = TRUE), Trait], c("ATD_Paternoster_26482879_1")))

# Sjogren
qs2[grepl("Sjogren|Sjögren|Sicca", Label, ignore.case = TRUE)][order(N1, decreasing = TRUE)]
tr <- c(tr, setdiff(qs2[grepl("Sjogren|Sjögren|Sicca", Label, ignore.case = TRUE), Trait], c("SJOS_Lessard_up_1")))

# Lupus
qs2[grepl("Lupus|SLE", Label, ignore.case = TRUE)][order(N1, decreasing = TRUE)]
tr <- c(tr, setdiff(qs2[grepl("Lupus|SLE", Label, ignore.case = TRUE), Trait], c("SLE_Julia_29848360_1")))

# Allergic
qs2[grepl("allerg", Label, ignore.case = TRUE)][order(N1, decreasing = TRUE)]
tr <- c(tr, setdiff(qs2[grepl("allerg", Label, ignore.case = TRUE), Trait], c("20002_1387_PanUKBB_PanUKBBR2_1")))

# Multiple sclerosis
qs2[grepl("multiple sclerosis", Label, ignore.case = TRUE)][order(N1, decreasing = TRUE)]
tr <- c(tr, setdiff(qs2[grepl("multiple sclerosis", Label, ignore.case = TRUE), Trait], c("G6_MS_FinnGen_FinnGenR7_1"))) # largest, non-IMSGC dataset


# Nasal polyps
qs2[grepl("polyp", Label, ignore.case = TRUE)][order(N1, decreasing = TRUE)]
tr <- c(tr, setdiff(qs2[grepl("polyp", Label, ignore.case = TRUE), Trait], c("J10_NASALPOLYP_FinnGen_FinnGenR7_1")))

# PBC
qs2[grepl("biliary", Label, ignore.case = TRUE)][order(N1, decreasing = TRUE)]
tr <- c(tr, setdiff(qs2[grepl("biliary", Label, ignore.case = TRUE), Trait], c("CHIRBIL_PRIM_FinnGen_FinnGenR7_1")))

# PSC
qs2[grepl("sclerosing", Label, ignore.case = TRUE)][order(N1, decreasing = TRUE)]
tr <- c(tr, setdiff(qs2[grepl("sclerosing", Label, ignore.case = TRUE), Trait], c("K11_CHOLANGI_FinnGen_FinnGenR7_1")))

# Systemic sclerosis
qs2[grepl("systemic sclerosis", Label, ignore.case = TRUE)][order(N1, decreasing = TRUE)]
tr <- c(tr, setdiff(qs2[grepl("systemic sclerosis", Label, ignore.case = TRUE), Trait], c("SSC_LopezIsac_31672989_1")))

# Autoimmune disorders
tr <- c(tr, qs2[grepl("autoimmune dis", Label, ignore.case = TRUE), Trait]) 

# Celiac disease
qs2[grepl("celiac|coeliac", Label, ignore.case = TRUE)][order(N1, decreasing = TRUE)]
tr <- c(tr, setdiff(qs2[grepl("celiac|coeliac", Label, ignore.case = TRUE), Trait], c("K11_COELIAC_FinnGen_FinnGenR7_1")))

# Gout and ankylosing
qs2[grepl("gout|ankylosing", Label, ignore.case = TRUE)][order(N1, decreasing = TRUE)]
tr <- c(tr, setdiff(qs2[grepl("gout|ankylosing", Label, ignore.case = TRUE), Trait], c("GOUT_FinnGen_FinnGenR7_1", "M13_ANKYLOSPON_FinnGen_FinnGenR7_1")))

# Vitiligo
qs2[grepl("vitiligo", Label, ignore.case = TRUE)][order(N1, decreasing = TRUE)]
tr <- c(tr, setdiff(qs2[grepl("vitiligo", Label, ignore.case = TRUE), Trait], c("L12_VITILIGO_FinnGen_FinnGenR7_1")))

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

tr <- c(tr, "M13_JUVERHEU_FinnGen_FinnGenR7_1") #  Remove juvenile rheuma, as it's essentially the same as JIA

qs2 <- qs2[!Trait %in% tr] # First pass
length(unique(qs2$Label))
# 73

# Check for duplicates, if there are still
dp <- qs2[duplicated(Label), unique(Label)]
qsdup <- qs2[Label %in% dp]
qsdup # Some duplicates

ps2 <- copy(ps)

qs2[ Trait %in% qsdup$Trait, Label:=paste0(Label, " / ", First_Author)]
ps2[ Trait %in% qsdup$Trait, Label:=paste0(Label, " / ", First_Author)]

tr2 <- qsdup[grepl("20002_|Chen_32888493_[1-4,6]|OTHER", Trait), Trait]
qs2 <- qs2[!Trait %in% tr2]

# Second pass
dp <- qs2[duplicated(Label), unique(Label)]
qsdup <- qs2[Label %in% dp]
qsdup # No duplicates

# Apply this to projections as well
ttk2 <- qs2$Trait # Traits to keep
ps2 <- ps2[Trait %in% ttk2]

ps2[, PC:=factor(PC, levels = paste0("PC", 1:13))]

# Prepare Delta plots with all available traits

make.custom.forest.plots <- function(ps, bs, feature, palette = tfp, fsize = NULL, threshold_basis=0, threshold_gen = 0, remove_sum = FALSE, remove_perc = FALSE, remove_side = "none"){
  pdelta <- ps[stars == "●"]
  pdelta <- rbind(pdelta, bs, fill = TRUE)
  pdelta[Var.Delta != 0, ci:=sqrt(Var.Delta) * 1.96][is.na(ci),ci:=0][is.na(Label), Label:=Trait]
  pdelta[ ci == 0, Trait_class:= "IMD"] # Basis traits are BC
  ds.col <- data.table(Trait_class = names(tfp), colours=tfp)
  # Need to update names in tfp for scale_colour_manual
  tfp.manual <- tfp
  names(tfp.manual) <- tfp.manual
  
  pdelta <- merge(pdelta, ds.col, by = "Trait_class", all.x=TRUE)
  
  ft <- paste0("PC", feature)
  dt <- pdelta[PC == ft][order(Delta, decreasing = TRUE)][!(Trait_class == "BC" & abs(Delta) < threshold_basis)][ abs(Delta) > threshold_gen]
  
  if(remove_sum) dt <- dt[!grepl("Sum ", Label, ignore.case = TRUE)]
  if(remove_perc) dt <- dt[!grepl("percentage of", Label, ignore.case = TRUE)]
  if(remove_side == "positive") dt <- dt[Delta <= 0]
  if(remove_side == "negative") dt <- dt[Delta >= 0]
  
  fplot <- ggplot(dt, aes(x = reorder(Label, -Delta), y = Delta, ymin=Delta-ci, ymax=Delta+ci, colour = colours))+
    geom_pointrange()+
    #scale_colour_manual(values = c("red" = "red", "#26547c" = "#26547c", "#049F76" = "#049F76", "#E09D00" = "#E09D00", "#F3B8A5" = "#F3B8A5", "black" = "black"))+
    scale_colour_manual(values = tfp.manual)+
    geom_hline(yintercept = 0, col="red", lty=2)+
    coord_flip()+
    #xlab("Traits")+
    ylab("Delta")+
    ggtitle(ft)+
    theme_minimal()+
    theme(axis.text.y = element_text(colour = dt$colours, size = fsize), legend.position = "none", axis.title.y = element_blank())
  fplot
}

tfp <- c(BC = "#CF000F", BMK = "#2E8856", IMD = "#1460AA", INF = "#B8860B", CAN = "#E65722", OTH ="#1C2833", PSD="#708090")


make.custom.forest.plots(ps2, b, feature = 2, palette = tfp, fsize = 8, threshold_basis=0, remove_sum = TRUE, remove_perc = TRUE, remove_side = "none")

make.custom.forest.plots(ps2, b, feature = 8, palette = tfp, fsize = 8, threshold_basis=0, remove_sum = TRUE, remove_perc = TRUE, remove_side = "none")

make.custom.forest.plots(ps2, b, feature = 9, palette = tfp, fsize = 8, threshold_basis=0, remove_sum = TRUE, remove_perc = TRUE, remove_side = "none")

make.custom.forest.plots(ps2, b, feature = 12, palette = tfp, fsize = 8, threshold_basis=0, remove_sum = TRUE, remove_perc = TRUE, remove_side = "none")

make.custom.forest.plots(ps2[Trait_class != "OTH"], b, feature = 3, palette = tfp, fsize = 8, threshold_basis=0, remove_sum = TRUE, remove_perc = TRUE, remove_side = "positive")
  
