##########################################
### Selecting diseases for clustering ####
##########################################

#
# Guillermo Reales
# 14-12-2021
#

# We want to select diseases to cluster with myositis datasets. To that end, we'll select only those that are overall significant and 
# significant for at least one of the 7 key components for myositis.

# Load libraries and data

library(cupcake)
library(data.table)
library(magrittr)
library(stringr)

proj.table <- fread("../data/Projection_IMD_basis_20211126-v1.tsv")
QC.table <- fread("../data/QC_IMD_basis_20211126-v1.tsv")
metadata <- fread("../data/Metadata_20211214-v1.tsv")
QC.table <- merge.data.table(QC.table, metadata, all.x = TRUE)

manifest.translator <- fread("../data/Manifest_build_translator.tsv")
manifest.translator[,pid19:=paste(CHR19, BP19, sep = ":")][, pid38:=paste(CHR38, BP38, sep = ":")]
manifest.translator <- manifest.translator[,c(1,8:9)]
SNP.manifest <- merge(manifest.translator, copy(cupcake::SNP.manifest), by.x = "pid19", by.y = "pid")


# Note HBC_Chen_6  have some issues that I couldn't resolve so far, so I'll exclude them
# We'll also exclude Lessard, since we don't have enough information about it yet -- and it's hence not in the Metadata table
## EOSC_Chen_3 used to have some problems in the blood cell basis, but we'll keep it
##### & Trait!="EOSC_Chen_32888493_3"

QC.table <- QC.table[!Trait %in% c("HBC_Chen_32888493_6" , "SJOS_Lessard_up_1")]
proj.table <- proj.table[!Trait %in% c("HBC_Chen_32888493_6" , "SJOS_Lessard_up_1")]

# Remove datasets with low (<80%) SNP coverage and Immunochip, except for Rothwell, which is myositis but is imputed.
# Remove datasets with <80% SNP match
# Remove remove datasets used to build the IMD basis, which are overfitted
basis.datasets <- c("CD_DeLange_28067908_1", "PSC_Ji_27992413_1", "UC_DeLange_28067908_1", "SLE_Bentham_26502338_1", "PBC_Cordell_26394269_1", "IGAN_Kiryluk_25305756_1", "CEL_Dubois_20190752_1", "MS_IMSGC_21833088_1", "AST_Demenais_29273806_1", "VIT_Jin_27723757_1", "RA_Okada_24390342_1", "LADA_Cousminer_30254083_1", "T1D_Cooper_doi101101120022_1")

QC.filt <- QC.table[nSNP >= max(nSNP) * 0.8 ][!Trait %in% basis.datasets][First_Author == "Rothwell"  | !grepl("ImmunoChip", Chip, ignore.case = TRUE)]

# Check number of datasets per category.
table(QC.filt$Trait_class)

# Create filtered projection table too. 
PT.filt <- merge(proj.table, QC.filt[,c("First_Author","Trait", "Trait_ID_2.0", "Trait_long","Trait_class", "Collection", "Population")], by = "Trait")

# Apply 1% FDR correction to overall p for all remaining datasets
QC.filt[, FDR.overall := p.adjust(overall_p, method = "BH"), by="Trait_class"]
QC.sig <- QC.filt[FDR.overall < 0.01,]
table(QC.sig$Trait_class)

# Apply 5% FDR correction by trait and PC to projections, then filter by overall significant traits
PT.filt[, FDR.PC:=p.adjust(P, method = "BH"), by = c("PC", "Trait_class")][, stars:=ifelse(!is.na(FDR.PC) & FDR.PC<0.05,"○","")][ FDR.PC < 0.01 , stars:="●"] ### Change here significance thresholds for the figure!
# ◑ Alternative
# • Original significant
# Filter PT.filt by overall significant traits
PT.sig <- PT.filt[Trait %in% QC.sig$Trait,]

## NEWS - We'll use PanUKBB instead of Neale's for having larger sample sizes. However Neale's might have some traits PanUKBB doesn't, so we'll keep those
nealenotpan <- setdiff(PT.sig[ First_Author == "Neale"]$Trait_ID_2.0, PT.sig[ First_Author == "PanUKBB"]$Trait_ID_2.0)
PT.sig <- PT.sig[ First_Author != "Neale" | Trait_ID_2.0 %in% nealenotpan]

# Restrain the selection to diseases significant for at least one of the key components, then remove Biomarkers
keypcs <- paste0("PC",c(1:3, 8:9, 12:13))
diseaselist <- PT.sig[PC %in% keypcs][stars == "●"][Trait_class != "BMK"]$Trait %>% unique()

# Create initial projection and QC table for processing
QC7 <- QC.sig[Trait %in% diseaselist]
PT7 <- PT.sig[Trait %in% diseaselist]


# Define custom Guille annotation function
make.labels <- function(x){
  x[,Label:=Trait_long] %>% 
    .[, Label:=gsub(" \\(UKBB\\)", "", Label)] %>% 
    .[, Label:=gsub(" \\(FinnGen\\)", "", Label)] %>% 
    .[, Label:=gsub(" \\(FG\\)", "", Label)] %>% 
    .[, Label:=gsub("Crohn disease \\(strict definition, all UC cases excluded\\)", "Crohn's, strict", Label)] %>% 
    .[, Label:=gsub("Eosinophilic Granulomatosis with Polyangiitis", "EGPA", Label)] %>% 
    .[ First_Author == "FinnGen", Label:=paste0(str_trunc(Label, width = 50), " / ", First_Author)] %>% 
    .[ First_Author != "FinnGen", Label:=paste0(str_trunc(paste0(Label, " (", Population,")"), width = 50), " / ", First_Author)] %>% 
    .[, Label:=gsub("Neale", "UKBB", Label)] %>% 
    .[, Label:=gsub("PanUKBB", "UKBB", Label)] # Addition to accommodate PanUKBB
  
}

QC7 <- make.labels(QC7)
# Remove FinnGen duplicates, keep the ones with largest sample size
dups <- QC7$Label[duplicated(QC7$Label)]
# This way we select for exclusion the ones with the lowest number of cases, or the ones with the lowest number of controls, if there's a draw.
dstoexclude <- QC7[Label %in% dups, c("Trait", "N0", "N1", "Label")][order(Label, N1, N0)][, .SD[1], by=Label]$Trait 

QC7 <- QC7[!Trait %in% dstoexclude]
PT7 <- PT7[Trait %in% QC7$Trait] # Update projection table
PT7 <- make.labels(PT7)

fwrite(QC7[, .(Trait, First_Author, Label)], "../data/raw_disease_for_clustering_list.tsv", sep = "\t")


