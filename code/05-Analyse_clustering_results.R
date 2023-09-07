#########################################
##                                     ##
##    ANALYSING CLUSTERING RESULTS    ##
##                                     ##
#########################################

# Author: Guillermo Reales 
# Date last updated: 2023/09/07

# Background: After running DPMUnc and computing the Bhattacharyya distance, we want to look at the results together

# This script will
# * Import DPMUnc and Bhattacharyya results.
# * Call clusters from Bhattacharyya distances.
# * Create joint heatmaps to visualise the clustering structure across the two methods
# * Create an alluvial diagram to compare both structures and allow us to select IMD for coloc.

##########################################


# Load required packages
library(data.table)
library(stringr)
library(reshape2)
library(pheatmap)
library(ggplot2)
library(ggalluvial)

# Load helper function

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

# Load Bhattacharyya 

bhd <- fread("../data/bh_dist.tsv")
bhd[, logbhat:=log(1 + bhat.dist)]
bhd[, T1:=str_trunc(T1, 50, ellipsis = " [...]")][, T2:=str_trunc(T2, 50, ellipsis = " [...]")]

# Call clusters in Bhattacharyya

# We keep the initial matrices for make_bold_names(), since I don't have the time to adapt it now
bh.cl <- acast(bhd, T1~T2, value.var = "bhat.dist") %>% as.dist %>% hclust(., method = "average") # Cluster by raw
bhatta.d <-  acast(bhd, T1~T2, value.var = "logbhat") # Show log



bh.ph <- pheatmap(bhatta.d,
                  cluster_cols = bh.cl, 
                  cluster_rows = bh.cl,
                  annotation_names_row = FALSE,
                  show_colnames = FALSE,
                  annotation_legend = FALSE,
                  cutree_rows = 9,
                  treeheight_col=0,
                  fontsize_row=8)


tbh <- bh.ph$tree_row
tbhc <- cutree(tbh, k=9) # manually adjusted to capture visually-selected IMD
bhatta.sel <- names(tbhc[tbhc ==6]) # Myositis cluster is 9
bhatta.sel
bht <- data.table(Label = names(tbhc), Bhattacharyya = tbhc)

# Save bhattacharyya clustering
fwrite(bht, "../data/bhattacharyya_clustering.tsv", sep="\t")

# Load DPMUnc clustering results

# We'll show DPMUnc annotations as well. Note, this file is generated later, but we recover it now
resdpmunc <- readRDS("../data/Myo_7PC_psm_data.rds")
resdpmunc <- data.table(Label = names(resdpmunc$calls$cl), DPMUnc = resdpmunc$calls$cl)
resdpmunc[, Label:=stringr::str_trunc(Label, 50, ellipsis = " [...]")]

#########################################

### Prepare joint heatmap fwith Bhattacharyya distance with both DPMUnc and Bhatta annotations


### Full palette from DPMUnc code. We'll use an adapted one (below)
# palette <- c("#77AADD", "#000000", "#9E0142", "#D53E4F", "#FDAE61", "#66C2A5", "#ABDDA4", "#F46D43", "#3288BD", "#5E4FA2",
#              "#FEE08B", "#E6F598", "#771155", "#AA4488", "#CC99BB", "#114477", "#774411", "#EEEEEE", "#117777", "#117744", 
#              "#44AA77", "#88CCAA", "#777711", "#44AAAA", "#AAAA44", "#77CCCC", "#DDDD77", "#AA7744", "#DDAA77", "#771122", "#AA4455", "#DD7788") # Colours for clusters  

# Set colours, this palette is adapted from the DPMUnc code, but removing one colour

# Legacy palette
# palette1 <- c("#77AADD", "#000000", "#9E0142", "#D53E4F", "#F46D43", "#FEE08B", "#ABDDA4", "#66C2A5",            "#5E4FA2", 
#               "#FDAE61", "#E6F598", "#771155", "#AA4488", "#CC99BB", "#114477", "#774411", "#EEEEEE", "#117777", "#117744", 
#               "#44AA77", "#88CCAA", "#777711", "#44AAAA", "#AAAA44", "#77CCCC", "#DDDD77", "#AA7744", "#DDAA77", "#771122", "#AA4455", "#DD7788")

# An alternative palette with stronger colours
# palette1 <- c("#117744", "#114477", "#9E0142", "#D53E4F", "#F46D43", "#FEE08B", "#ABDDA4", "#5E4FA2", "#66C2A5",             
#               "#FDAE61", "#E6F598", "#771155", "#AA4488", "#CC99BB", "#114477", "#774411", "#EEEEEE", "#117777", "#117744", 
#               "#44AA77", "#88CCAA", "#777711", "#44AAAA", "#AAAA44", "#77CCCC", "#DDDD77", "#AA7744", "#DDAA77", "#771122", "#AA4455", "#DD7788")

palette1 <- c("#9E0142", "#77AADD", "#F46D43", "#ABDDA4", "#5E4FA2", "#000000",  "#D53E4F", "#117777",           "#FEE08B", 
              "#FDAE61", "#E6F598", "#771155", "#AA4488", "#CC99BB", "#114477", "#774411", "#EEEEEE", "#117777", "#117744", 
              "#44AA77", "#88CCAA", "#777711", "#44AAAA", "#AAAA44", "#77CCCC", "#DDDD77", "#AA7744", "#DDAA77", "#771122", "#AA4455", "#DD7788")

bhattacolp <- rev(colorRampPalette(RColorBrewer::brewer.pal(n = 9, name = "Spectral"))(100))

# Prepare annotation object

ann <- merge(bht, resdpmunc, by = "Label")
ann <- data.frame(ann[, 2:3], row.names = ann$Label) # Format

dpcol <- palette1[1:max(ann$DPMUnc)]
names(dpcol)  <- as.character(1:max(ann$DPMUnc))
bhcol <- palette1[1:max(ann$Bhattacharyya)]
names(bhcol)  <- as.character(1:max(ann$Bhattacharyya))

annotations <- list(ann = ann, colors = list( Bhattacharyya = bhcol,DPMUnc = dpcol))
myob <- grep("myositis|IIM", rownames(bhatta.d), value = TRUE, ignore.case = TRUE)

# Make heatmap with Bhatta contents and Bhata/DPMUnc clustering annotations

bh.ph <- pheatmap(bhatta.d,
                  cluster_cols = bh.cl, 
                  cluster_rows = bh.cl,
                  annotation_names_row = FALSE,
                  annotation_legend = FALSE,
                  treeheight_col=0,
                  fontsize_row=8,
                  annotation_row = annotations$ann,
                  annotation_col = annotations$ann,
                  annotation_colors = annotations$colors,
                  color = bhattacolp, 
                  show_colnames = FALSE, 
                  labels_row = make_bold_names(bhatta.d, rownames, myob))

# Save
# ggsave("../figures/Myositis_bhattacharyya_heatmap.svg", bh.ph, width = 9, height = 9.5, bg="white")
# system("sed -i \"s/ textLength=\'[^\']*\'//\" ../figures/Myositis_bhattacharyya_heatmap.svg") # Trick to make the svg file text be more easily editable

# Note: This figure will be edited to add a rectangle highlighting the Bhatta myositis group.


#########################################

### Prepare joint heatmap with DPMUnc PSM and both DPMUnc and Bhatta annotations

## To do on HPC

########################################

### Prepare alluvial plot
# Recover annotation object
rescl <- ann %>% data.table(keep.rownames = TRUE)

clsum=rescl[,.(y=.N),by=c("DPMUnc","Bhattacharyya")]


# Plot Relationship between disease clusterings
sdiag <- ggplot(clsum, aes(y=y, axis1=DPMUnc,axis2=Bhattacharyya)) + geom_alluvium(aes(fill=factor(DPMUnc))) +
  geom_stratum(width = 1/5,  fill="white", color = "black") +
  #geom_label(stat = "stratum", aes(label = after_stat(stratum))) +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)))+
  scale_x_discrete(limits = c("DPMUnc","Bhattacharyya"), expand = c(.05, .05)) +
  scale_fill_manual(values = palette1[1:10])+
  scale_y_continuous("Count") +
  theme_minimal()+
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), 
        axis.title.y = element_blank(), panel.grid = element_blank(), 
        legend.position = "none", axis.text.x = element_text(size = 13))
sdiag
# Save
# ggsave("../figures/sdiag_Bh_DP.png", sdiag, height = 9, width = 8, bg = "white")
# ggsave("../figures/sdiag_Bh_DP.svg", sdiag, height = 9, width = 8, bg = "white")
# system("sed -i \"s/ textLength=\'[^\']*\'//\" ../figures/sdiag_Bh_DP.svg")






