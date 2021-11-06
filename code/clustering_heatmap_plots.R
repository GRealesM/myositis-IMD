## Preparing Cluster Heatmap

# Load libraries
library(data.table)
library(magrittr)
library(pheatmap)
library(colorspace)
library(RColorBrewer)


## Load Data
load("../data/bigpsm.RData")
load("../data/delta_clustered_data.RData") # For Ann colours
rm(msub,msum)
md <- fread("../data/Metadata_20210812-v1.tsv")

# Fix labels
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

# Use Metadata to match shorttraits
md <- md[!grepl("PanUKBB", First_Author)]                                                                     # Remove PanUKBB, as we're using Neale UKBB only
md[, shorttrait:=ifelse(First_Author == "Neale", Trait_long, Trait)][, shorttrait:=shorten_names(shorttrait)] # Create a shorttrait column and apply function
md <- md[shorttrait %in% rownames(bigpsm)][match(rownames(bigpsm), shorttrait)]                               # Match and put in the same order as the matrix
md <- make.labels(md)                                                                                         # Create new labels
rownames(bigpsm) <- md$Label
rownames(ann$ann) <- md$Label

# Rebrand call to cluster
names(ann$ann)[1] <- "cluster"
names(ann$colors)[2] <- "cluster"

# Build Heatmap
# Palette for HeatMap
#pal1 <- sequential_hcl(100, palette = "Terrain2")
#pal1 <- sequential_hcl(100, palette = "Grays")
#pal1 <- colorRampPalette(c("#e3f2fd","#bbdefb","#90caf9","#64b5f6","#42a5f5","#2196f3","#1e88e5","#1976d2","#1565c0","#0d47a1"))(100) # Blue scales
#pal1 <- colorRampPalette(c("#dad7cd","#a3b18a","#588157","#3a5a40","#344e41"))(100) #  Green scale
#pal1 = rev(colorRampPalette(c("#10002b","#240046","#3c096c","#5a189a","#7b2cbf","#9d4edd","#c77dff","#e0aaff"))(100) ) # Purple scale
#pal1 = rev(colorRampPalette(c("#007f5f","#2b9348","#55a630","#80b918","#aacc00","#bfd200","#d4d700","#dddf00","#eeef20","#ffff3f"))(100) ) # Yellow and green scale
#pal1 = rev(colorRampPalette(c("#b7094c","#a01a58","#892b64","#723c70","#5c4d7d","#455e89","#2e6f95","#1780a1","#0091ad"))(100) ) # Blue and purple
#pal1 = colorRampPalette(c("#083d77","#ebebd3","#f4d35e","#ee964b","#f95738"))(100) # Blue orange
pal1 = rev(colorRampPalette(c("#0a1128","#001f54","#034078","#1282a2","#fefcfb"))(100) ) # Oxford blue & White
#pal1 <- colorRampPalette(c("red", "blue"))(100)
plot(rep(1,100),col=pal1,pch=15,cex=5)
# Palette for clusters
pal2 <- brewer.pal(12, name = "Paired")

hmp <- pheatmap(bigpsm, color = pal1, cluster_rows = TRUE, cluster_cols = TRUE, annotation_row=ann$ann[,"cluster", drop=FALSE], annotation_colors = list(cluster=ann$colors$cluster), fontsize_row = 8)
hmp
# Save the Heatmap
save_pheatmap_png <- function(x, filename, width=1200, height=1000, res = 150) {
  png(filename, width = width, height = height, res = res)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

save_pheatmap_png(hmp, filename = "../figures/clustering_heatmap_v2.png", width = 1500,height = 1500)

