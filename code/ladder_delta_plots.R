# Ladder plots for myositis presentation

# We want to use the ladder (ex-forest) plots from the IMD basis for PC1, PC12, and PC13, using the traits used in clustering.

# Load libraries
library(cupcake)
library(data.table)
library(magrittr)
library(ggplot2)
library(cowplot)

# Define useful naming functions
# Define Chris' shorttrait function
shorten_names=function(x) {
  x %>%
    sub("_FinnGenR5_1","",.) %>%
    sub("_[0-9]+_1$","",.) %>%
    sub("_up_1","",.) %>%
    sub("_FinnGen"," FG",.) %>%
    sub("colitis\\/not","colitis UKBB",.) %>%
    sub("(UKBB)","UKBB",.) %>%
    sub(".*\\/","",.)
}

# Define custom Guille annotation function
make.labels <- function(x){
  x[,Label:=Trait_long] %>% 
    .[, Label:=gsub(" \\(UKBB\\)", "", Label)] %>% 
    .[, Label:=gsub(" \\(FinnGen\\)", "", Label)] %>% 
    .[, Label:=gsub(" \\(FG\\)", "", Label)] %>% 
    .[, Label:=paste0(Label, " / ", First_Author)] %>% 
    .[, Label:=gsub("Neale", "UKBB", Label)] %>% 
    .[, Label:=gsub("Crohn disease \\(strict definition, all UC cases excluded\\)", "Crohn's, strict", Label)] %>% 
    .[, Label:=gsub("Eosinophilic Granulomatosis with Polyangiitis", "EGPA", Label)]
}


# Load data
proj.table <- fread("../../Bases/IMD_basis/Projections/Projection_IMD_basis_20210907-v1.tsv")

if(is.character(proj.table$Var.Delta)){
  proj.table$Var.Delta <- as.numeric(proj.table$Var.Delta)
}

QC.table <- fread("../../Bases/IMD_basis/Projections/QC_IMD_basis_20210907-v1.tsv")
metadata <- fread("../data/Metadata_20210812-v1.tsv")
QC.table <- merge.data.table(QC.table, metadata, all.x = TRUE)
load("../data/bigpsm.RData") # To extract the labels and traits

# Filter Basis traits and compute FDR for significance

# Note HBC_Chen_6 and EOSC_Chen_3 have some issues that I couldn't resolve so far, so I'll exclude it
##### & Trait!="EOSC_Chen_32888493_3"
QC.table <- QC.table[Trait!="HBC_Chen_32888493_6" ,]
proj.table <- proj.table[Trait!="HBC_Chen_32888493_6",]

# Extract basis projections for later
PT.basis <- copy(cupcake::basis.trait.proj)
setnames(PT.basis, c("delta", "trait"), c("Delta", "Label"))
PT.basis[,proj:=NULL][, stars:=""]

# Remove Immunochip datasets and those with <80% SNP match
QC.filt <- QC.table[QC.table$nSNP >= max(QC.table$nSNP) * 0.8 & !grepl("ImmunoChip",QC.table$Chip, ignore.case = TRUE), ]
# Remove  FinnGenR2 and FinnGenR3, and UKBB
QC.filt <- QC.filt[!grepl("FinnGenR[23]", Trait, perl = TRUE)]
#Also, remove datasets used to build the IMD basis, which are overfitted
basis.datasets <- c("CD_DeLange_28067908_1", "PSC_Ji_27992413_1", "UC_DeLange_28067908_1", "SLE_Bentham_26502338_1", "PBC_Cordell_26394269_1", "IGAN_Kiryluk_25305756_1", "CEL_Dubois_20190752_1", "MS_IMSGC_21833088_1", "AST_Demenais_29273806_1", "VIT_Jin_27723757_1", "RA_Okada_24390342_1", "LADA_Cousminer_30254083_1", "T1D_Cooper_doi101101120022_1")
QC.filt <- QC.filt[!Trait %in% basis.datasets]
# Check number of datasets per category, and keep IMD, BMK, and INF only
#table(QC.filt$Trait_class)
QC.filt <- QC.filt[Trait_class %in% c("IMD", "INF", "BMK", "OTH")]
table(QC.filt$Trait_class)
# Create filtered projection table too. At this step we'll remove (1) Basis projections, (2) Projections from similar traits to those used to create the basis, and (3) All datasets not IMD, INF, BMK and OTH.
PT.filt <- merge(proj.table, QC.filt[,c("First_Author","Trait", "Trait_ID_2.0", "Trait_long","Trait_class", "Collection", "Population")], by = "Trait")

# Apply 1% FDR correction to overall p for all remaining datasets
QC.filt[, FDR.overall := p.adjust(overall_p, method = "BH"), by="Trait_class"]
QC.sig <- QC.filt[FDR.overall < 0.01,]
table(QC.sig$Trait_class)

# Apply 1% FDR correction by trait and PC to projections, then filter by overall significant traits
PT.filt[, FDR.PC:=p.adjust(P, method = "BH"), by = c("PC", "Trait_class")][, stars:=ifelse(!is.na(FDR.PC) & FDR.PC<0.01,"•","")] ### Change here significance thresholds for the figure!
# Filter PT.filt by overall significant traits
PT.sig <- PT.filt[Trait %in% QC.sig$Trait,]

# Now we'll prepare PT.sig to merge with cluster traits
PT.sig <- PT.sig[!grepl("PanUKBB", First_Author)]                                                                     # Remove PanUKBB, as we're using Neale UKBB only
PT.sig[, shorttrait:=ifelse(First_Author == "Neale", Trait_long, Trait)][, shorttrait:=shorten_names(shorttrait)] # Create a shorttrait column and apply function
PT.sig <- PT.sig[shorttrait %in% rownames(bigpsm)]   
# Create appropriate labels
PT.sig <- make.labels(PT.sig)

# Make ladder plots
PT.delta <- PT.sig[stars == "•",]
PT.delta <- rbind(PT.basis, PT.delta, fill = TRUE)
PT.delta[Var.Delta != 0, ci:=sqrt(Var.Delta) * 1.96][is.na(ci),ci:=0][is.na(Label), Label:=Trait]
PCs <- unique(PT.delta$PC)
forest.plots <- lapply(setNames(PCs, PCs), function(x){
  dt <- PT.delta[PC == x,][order(Delta, decreasing = TRUE),][,colours:=ifelse(ci != 0, "black", "red")][grepl("Miller",Label), colours:="darkgreen"]
  fplot <- ggplot(dt, aes(x = reorder(Label, -Delta), y = Delta, ymin=Delta-ci, ymax=Delta+ci, colour = colours))+
    geom_pointrange()+
    scale_colour_manual(values = c("red" = "red", "black" = "black", "darkgreen" = "darkgreen"))+
    geom_hline(yintercept = 0, col="red", lty=2)+
    coord_flip()+
    xlab("Traits")+
    ylab("Delta")+
    ggtitle(paste("Ladder Plot ", x, sep = ""))+
    theme_minimal()+
    theme(axis.text.y = element_text(colour = dt$colours), legend.position = "none")
  fplot
})

forest.plots[["PC1"]]
forest.plots[["PC12"]]
forest.plots[["PC13"]]

ggsave("../figures/Ladder_plot_PC1.png", forest.plots[["PC1"]], units = "mm", width = 180, height = 220, dpi = 300, bg="white")
ggsave("../figures/Ladder_plot_PC12.png", forest.plots[["PC12"]], units = "mm", width = 130, height = 159, dpi = 300, bg="white")
ggsave("../figures/Ladder_plot_PC13.png", forest.plots[["PC13"]], units = "mm", width = 180, height = 220, dpi = 300, bg="white")

