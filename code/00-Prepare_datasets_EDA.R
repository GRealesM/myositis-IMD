#########################################
##                                     ##
##     PREPARING DATASETS AND EDA      ##
##                                     ##
#########################################

# Author: Guillermo Reales 
# Date last updated: 2023/04/10

# Background: This script will prepare the datasets for downstream analyses and prepare some figures.

# This script will
# * Apply FDR to projections.
# * Remove redundant datasets for improved visualisation
# * Create some figures and supplementary tables, and save datasets for downstream analyses
# * Prepare datasets for DPMUnc


##########################################

setwd("/home/gr440/rds/rds-cew54-basis/Projects/myositis-IMD/code")

## Load required packages

library(data.table)
library(magrittr)
library(ggplot2)
library(cowplot)
library(pheatmap)
library(reshape2)
library(cupcake)

## Load datasets

pf <- fread("../data/pf.tsv")
qf <- fread("../data/qf.tsv")

### (1) Apply FDR procedure and remove datasets with FDR overall > 1%

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


### (2) Remove redundant datasets for visualisation and further analyses

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
qs2[grepl("Hypothyroidism|thyroiditis|thyroid disease", Label, ignore.case = TRUE)][order(N1, decreasing = TRUE)]
tr <- c(tr, setdiff(qs2[grepl("Hypothyroidism|thyroiditis|thyroid disease", Label, ignore.case = TRUE), Trait], c("E4_HYTHY_AI_STRICT_FinnGen_FinnGenR7_1"))) # Take strict definition

# Hyperthyroidism
qs2[grepl("hyperthyroidism|graves", Label, ignore.case = TRUE)][order(N1, decreasing = TRUE)]
tr <- c(tr, setdiff(qs2[grepl("hyperthyroidism|graves", Label, ignore.case = TRUE), Trait], "20002_1225_PanUKBB_PanUKBBR2_1")) # UKBB has larger sample size. 

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

# We'll also remove all myositis from FinnGen, as they're not informative. Also biomedrhe, as we already have rheuma
tr <- c(tr, "DERMATOPOLY_FG_FinnGen_FinnGenR7_1", "M13_MYOSITIS_FinnGen_FinnGenR7_1", "M13_POLYMYO_FinnGen_FinnGenR7_1", "M13_DERMATOPOLY_FinnGen_FinnGenR7_1", "RX_RHEUMA_BIOLOGICAL_FinnGen_FinnGenR7_1")


qs2 <- qs2[!Trait %in% tr] # First pass
length(unique(qs2$Label))
# 62  

qs2[, .(Trait, Label, N0, N1, N, FDR.overall)][order(Label)]

# Check for duplicates, if there are still
dp <- qs2[duplicated(Label), unique(Label)]
qsdup <- qs2[Label %in% dp]
qsdup # No duplicates

# Apply this to projections as well
ttk2 <- qs2$Trait # Traits to keep
ps2 <- ps[Trait %in% ttk2]

ps2[, PC:=factor(PC, levels = paste0("PC", 1:13))]

# fwrite(ps2, "../data/ps2.tsv", sep="\t")
# fwrite(qs2, "../data/qs2.tsv", sep="\t")

# Prepare supplementary tables for ALL datasets and projections

# Save filtered datasets - not sure if we'll show these

aq <- copy(qf)
ap <- copy(pf)

aq[, sig.overall:=ifelse(Trait %in% qs$Trait, "Yes", "No")][, in.selection:=ifelse(Trait %in% qs2$Trait, "Yes", "No")]
ap[, sig.overall:=ifelse(Trait %in% qs$Trait, "Yes", "No")][, in.selection:=ifelse(Trait %in% qs2$Trait, "Yes", "No")]

aq <- aq[,.(Trait, Label, First_Author, Reference, N0, N1,N, Population, nSNP, mscomp, overall_p, FDR.overall, sig.overall, in.selection)]
ap <- ap[,.(Trait, Label, First_Author, Population, PC, Delta, Var.Delta, P, FDR.PC, stars, sig.overall, in.selection)]

# fwrite(aq, "../tables/ST_all_datasets.tsv", sep ="\t")
# fwrite(ap, "../tables/ST_all_projections.tsv", sep ="\t")

##########################################

### Prepare some basic figures

## Figure 1 -- A heatmap of myositis projections

pspm <- ps2[grepl("Miller|Rothwell", Label, ignore.case = TRUE)]

pspm[, Label:=gsub("Miller", "M", Label)][, Label:=gsub("Rothwell", "R", Label)][, Label:=gsub("Juvenile Dermatomyositis", "JDM", Label)][, Label:=gsub("Dermatomyositis", "DM", Label)][, Label:=gsub("Polymyositis", "PM", Label)][, Label:=gsub("Inclusion body myositis", "IBM", Label, ignore.case = T)][, Label:=gsub("Jo1\\+ Myositis", "Anti-Jo1+", Label, ignore.case = T)]

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
# ggsave("../figures/Fig1_Myositis_allsources_heatmap.png", Mphm, width = 4.5, height = 2.5, bg="white")
# ggsave("../figures/Fig1_Myositis_allsources_heatmap.svg", Mphm, width = 4.5, height = 2.5, bg="white")
# 
# system("sed -i \"s/ textLength=\'[^\']*\'//\" ../figures/Myositis_allsources_heatmap.svg") # Trick to make the svg file text be more easily editable

###

## Internal figure -- Heatmap of all 66 projections


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
# ggsave("../figures/Myositis_IMD_heatmap.png", Mahm, width = 8, height = 14, bg="white")


### 

## Figure SXX -- Delta plot of all myositis datasets across all 7 PCs

myoc <- c(`PM (R)` = "#CF000F", `PM (M)` = "#CF000F", `DM (R)` = "#2E8856", `DM (M)` = "#2E8856", `IIM (R)` = "#1460AA", `IIM (M)` = "#1460AA", `JDM (M)` = "#B8860B", `JDM (R)` = "#B8860B", `IBM (R)` = "#E65722", `Anti-Jo1+ (R)` ="#1C2833")

pmyo <- pf[(First_Author %in% c("Miller", "Rothwell") ) & PC %in% paste0("PC", c(1:3, 8:9, 12:13))][, PC:=factor(PC, levels = paste0("PC", c(1:3, 8:9, 12:13)))]
# | Trait_ID_2.0 %in% c("M13_DERMATOPOLY", "M13_POLYMYO")
pmyo[, Label:=gsub("Inclusion Body Myositis", "IBM", Label)][, Label:=gsub("Juvenile Dermatomyositis", "JDM", Label)][, Label:=gsub("Dermatomyositis", "DM", Label)][, Label:=gsub("Polymyositis", "PM", Label)][, Label:=gsub("Jo1\\+ Myositis", "Anti-Jo1+", Label)][, Label:=gsub("Miller", "M", Label)][, Label:=gsub("Rothwell", "R", Label)][, Label:=gsub("FinnGen", "FG", Label)][, Label:=gsub("Dermatopolymyositis", "DPM", Label)]
pmyo[, ci:=sqrt(Var.Delta) * 1.96]

dpm <- ggplot(pmyo, aes(x = Delta, y = Label, xmin=Delta-ci, xmax=Delta+ci, colour = Label))+
  geom_pointrange()+
  geom_vline(xintercept = 0, col="red", lty=2)+ 
  scale_colour_manual(values = myoc)+
  xlab("Delta")+
  facet_grid(PC~.,  scales = "free", space = "free", switch = "y")+
  theme_cowplot(11)+
  theme(legend.position = "none", strip.text.y.left = element_text(angle = 0), axis.title.y = element_blank())
dpm

#ggsave("../figures/FigSXX_deltaplot_myo.png", dpm, height = 10, width = 6, bg = "white")


###

# Figure SXX - Myositis + IMD across the 7 key PCs

dpdf <- ps2[ PC %in% paste0("PC", c(1:3, 8:9, 12:13)) & stars == "●"]
dpdf[, Label:=gsub("Miller", "M", Label)][, Label:=gsub("Rothwell", "R", Label)]
dpdf[, Label:=gsub("Juvenile Dermatomyositis", "JDM", Label)][, Label:=gsub("Dermatomyositis", "DM", Label)]
dpdf[, Label:=gsub("Polymyositis", "PM", Label)][, Label:=gsub("Jo1\\+ Myositis", "Anti-Jo1+", Label)]
dpdf[, colours:=ifelse(First_Author %in% c("Miller", "Rothwell"), "#CF000F", "#1460AA")]
dpdf[, ci:=sqrt(Var.Delta) * 1.96]

PCs <- paste0("PC", c(1:3, 8:9, 12:13))

kdp <- lapply(setNames(PCs, PCs), function(x){
       
        dt <- dpdf[PC==x][order(Delta, decreasing = TRUE),]
        dpm <- ggplot(dt, aes(x = Delta, y = reorder(Label, -Delta), xmin=Delta-ci, xmax=Delta+ci, colour = colours))+
          geom_pointrange()+
          geom_vline(xintercept = 0, col="red", lty=2)+ 
          scale_colour_manual(values = c("#CF000F" = "#CF000F", "#1460AA" =  "#1460AA"))+
          xlab("Delta")+
          ylab("Traits")+
          ggtitle(x)+
          theme_cowplot(11)+
          theme(axis.text.y = element_text(colour = dt$colours), legend.position = "none")
        dpm
})

# ggsave("../figures/FigSXX_deltaplots_PC1.png", kdp$PC1, bg = "white", height =5, width = 5)
# ggsave("../figures/FigSXX_deltaplots_PC2.png", kdp$PC2, bg = "white", height =4, width = 5)
# ggsave("../figures/FigSXX_deltaplots_PC3.png", kdp$PC3, bg = "white", height =6, width = 5)
# ggsave("../figures/FigSXX_deltaplots_PC8.png", kdp$PC8, bg = "white", height =4, width = 5)
# ggsave("../figures/FigSXX_deltaplots_PC9.png", kdp$PC9, bg = "white", height =3, width = 5)
# ggsave("../figures/FigSXX_deltaplots_PC12.png", kdp$PC12, bg = "white", height =5, width = 5)
# ggsave("../figures/FigSXX_deltaplots_PC13.png", kdp$PC13, bg = "white", height =5, width = 5)



##########################################


## Prepare inputs for DPMUnc

# At this point, we prepare the data for DPMUnc using significant PCs for myositis 
ps2[grepl("myositis|IIM", Label, ignore.case = TRUE) & FDR.PC < 0.01, .(PC, Label, stars)][, .N, by = PC]

dpmunc.ds <- ps2[PC %in% paste0("PC", c(1,2,3,8,9,12,13)), .(PC, Delta, Var.Delta, Label)]
dpmunc.delta <- reshape(dpmunc.ds[, .(PC, Delta, Label)], idvar="Label", timevar = "PC", direction = "wide")
dpmunc.var   <- reshape(dpmunc.ds[, .(PC, Var.Delta, Label)], idvar="Label", timevar = "PC", direction = "wide")

# Save
# fwrite(dpmunc.delta, "../data/Myo_7PC_Delta.tsv", sep = "\t")
# fwrite(dpmunc.var, "../data/Myo_7PC_Var.tsv", sep = "\t")


## Additionaly, we'll make a secondary selection, using all 13 PCs


dpmunc.13.ds <- ps2[, .(PC, Delta, Var.Delta, Label)]
dpmunc.13.delta <- reshape(dpmunc.13.ds[, .(PC, Delta, Label)], idvar="Label", timevar = "PC", direction = "wide")
dpmunc.13.var   <- reshape(dpmunc.13.ds[, .(PC, Var.Delta, Label)], idvar="Label", timevar = "PC", direction = "wide")

# Save
# fwrite(dpmunc.13.delta, "../data/Myo_13PC_Delta.tsv", sep = "\t")
# fwrite(dpmunc.13.var, "../data/Myo_13PC_Var.tsv", sep = "\t")


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
#   [1] cupcake_0.1.0.0   reshape2_1.4.4    pheatmap_1.0.12   cowplot_1.1.1     ggplot2_3.4.3     magrittr_2.0.3    data.table_1.14.8
# 
# loaded via a namespace (and not attached):
#   [1] Matrix_1.6-0        gtable_0.3.4        dplyr_1.1.3         compiler_4.3.1      tidyselect_1.2.0    Rcpp_1.0.11         stringr_1.5.0       splines_4.3.1      
# [9] scales_1.2.1        lattice_0.21-8      R6_2.5.1            plyr_1.8.8          labeling_0.4.3      generics_0.1.3      BiocGenerics_0.46.0 tibble_3.2.1       
# [17] snpStats_1.50.0     munsell_0.5.0       pillar_1.9.0        RColorBrewer_1.1-3  rlang_1.1.1         utf8_1.2.3          stringi_1.7.12      cli_3.6.1          
# [25] withr_2.5.0         zlibbioc_1.46.0     grid_4.3.1          rstudioapi_0.15.0   lifecycle_1.0.3     vctrs_0.6.3         glue_1.6.2          farver_2.1.1       
# [33] survival_3.5-5      fansi_1.0.4         colorspace_2.1-0    purrr_1.0.2         tools_4.3.1         pkgconfig_2.0.3  