#########################################
##                                     ##
##    ANALYSING CLUSTERING RESULTS     ##
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
library(magrittr)
library(stringr)
library(reshape2)
library(pheatmap)
library(ggplot2)
library(ggalluvial)

# For PSM plot only
library(R.cache)
library(mcclust)
library(clue)


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
bh.cl <- acast(bhd, T1~T2, value.var = "bhat.dist") %>% as.dist %>% hclust(., method = "complete") # Cluster by raw
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
bhatta.sel <- names(tbhc[tbhc ==5]) # Myositis cluster is 9
bhatta.sel
bht <- data.table(Label = names(tbhc), Bhattacharyya = tbhc)

# Save bhattacharyya clustering
# fwrite(bht, "../data/bhattacharyya_clustering.tsv", sep="\t")

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

# Figure S11 - Heatmap with Bhatta contents and Bhata/DPMUnc clustering annotations

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
# ggsave("../figures/FigS11_Myositis_7PCs_BHDP_heatmap.svg", bh.ph, width = 9, height = 9.5, bg="white")
# ggsave("../figures/FigS11_Myositis_7PCs_BHDP_heatmap.png", bh.ph, width = 9, height = 9.5, bg="white")
# system("sed -i \"s/ textLength=\'[^\']*\'//\" ../figures/FigS11_Myositis_7PCs_BHDP_heatmap.svg") # Trick to make the svg file text be more easily editable

# Note: This figure will be edited to add a rectangle highlighting the Bhatta myositis group.


########################################

### Prepare alluvial plot
# Recover annotation object
rescl <- ann %>% data.table(keep.rownames = TRUE)

clsum=rescl[,.(y=.N),by=c("DPMUnc","Bhattacharyya")]


# Figure 2 -- Plot Relationship between disease clusterings

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
# ggsave("../figures/Fig2_Myositis_7PCs_sdiag_BHDP.png", sdiag, height = 9, width = 8, bg = "white")
# ggsave("../figures/Fig2_Myositis_7PCs_sdiag_BHDP.svg", sdiag, height = 9, width = 8, bg = "white")
# system("sed -i \"s/ textLength=\'[^\']*\'//\" ../figures/Fig2_Myositis_7PCs_sdiag_BHDP.svg") # Easier to edit in Inkscape.
# Note: Final Figure 2 was manually edited in Inkscape to modify font sizes and highlight diseases involved.  


#########################################

### Prepare joint heatmap with DPMUnc PSM and both DPMUnc and Bhatta annotations

## To do on HPC. This bit will be a bit longer.

# Some helper functions
# NOTE: next function was in scripts/utils.R. Modified from original to remove burnin
calc_psms <- function(datasets, burnin) {  
    allocs=lapply(paste0(datasets, "/clusterAllocations.csv"), fread) ## read the allocations
    # This line is essential for some reason
    allocs %<>% lapply(., function(x) as.matrix(x[1:nrow(x),]))
    message("Removing burnin.")
    allocs %<>% lapply(., function(x) x[-c(1:burnin),]) # Remove burnin at this stage

    bigalloc = do.call(rbind, allocs)    
    bigpsm=calc_psm(bigalloc,burn=0) ## make a psm, don't discard any burn in because already discarded
    psms = lapply(allocs, function(x) calc_psm(x, burn=0))
    return(list(bigpsm=bigpsm, psms=psms))
}

raw_calc_psm=function(x,burn=0.5) {
  n=nrow(x)
  if(burn>0)
    x=x[ (burn*n) : n, , drop=FALSE]
  if(any(is.na(x)))
    x=x[ apply(!is.na(x),1,all), ]
  unq=unique(as.vector(x))
  ## print(unq)
  m=matrix(0,ncol(x),ncol(x))
  for(k in unq) {
    xk=matrix(as.numeric(x==k),nrow(x),ncol(x))
    ## m=m + t(xk) %*% xk
    m=m + crossprod(xk)
  }
  psm=m/nrow(x)
  psm
}

calc_psm <- addMemoization(raw_calc_psm)

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

# Figure S11 - Heatmap with DPMUnc contents and Bhata/DPMUnc clustering annotations

# Process. This is just adapted from the PSM 

exp="Myo_7PC"
burnin = 250000
input_dir="../data/DPMUnc_results/"

datasets = dir(input_dir)[grepl(exp, dir(input_dir))] # Little adaptation to original function. From experiment name we retrieve all directories for different seeds.
datasets = paste0(input_dir, datasets)

obsData = read.table(paste0("../data/", exp, "_Delta.tsv"),
                     header=1, row.names=1, quote="", sep="\t")
obsVars = read.table(paste0("../data/", exp, "_Var.tsv"),
                     header=1, row.names=1, quote="", sep="\t")

# bht = fread("../data/bhattacharyya_clustering.tsv") # already loaded


# Some names are too long for display -- and won't match bhattacharyya, which are truncated already
rownames(obsData) = stringr::str_trunc(rownames(obsData), 50, ellipsis = " [...]")


message("Preparing PSMs.")
result = calc_psms(datasets, burnin)
bigpsm = result$bigpsm
psms = result$psms

print(isSymmetric(bigpsm))
print(max(bigpsm))
print(min(bigpsm))
print(diag(bigpsm))
rownames(bigpsm) = rownames(obsData)
colnames(bigpsm) = rownames(obsData)

# Call clusters
calls = minbinder(as.dist(1 - bigpsm), method = "comp") ## calls
hclust.comp <- hclust(as.dist(1 - bigpsm), method = "complete")

# Prepare annotations
ann <- data.table(Label = names(calls$cl), DPMUnc = calls$cl)
ann <- merge(ann, bht, by = "Label")

ann <- data.frame(ann[, 2:3], row.names = ann$Label)
dpcol <- palette1[1:max(ann$DPMUnc)]
names(dpcol)  <- as.character(1:max(ann$DPMUnc))
bhcol <- palette[1:max(ann$Bhattacharyya)]
names(bhcol)  <- as.character(1:max(ann$Bhattacharyya))

annotations <- list(ann = ann, colors = list(DPMUnc = dpcol, Bhattacharyya = bhcol))

message("Creating PSM heatmap.")
psm_heatmap = pheatmap(bigpsm,
                       show_rownames = TRUE,
                       show_colnames = FALSE,
                       cluster_rows = hclust.comp,
                       cluster_cols = hclust.comp,
                       annotation_names_row = FALSE,
                       annotation_legend = FALSE,
                       treeheight_col=0,
                       fontsize_row=8,
                       annotation_row = annotations$ann,
                       annotation_col = annotations$ann,
                       color=colorRampPalette((RColorBrewer::brewer.pal(n = 7,
                                                name = "Blues")))(100),
                       annotation_colors = annotations$colors,
                       labels_row = make_bold_names(bigpsm, rownames, myob))

# Save
# ggsave("../figures/FigS12_Myositis_7PCs_DPBH_heatmap.svg", psm_heatmap, width = 9, height = 9.5, bg="white")
# ggsave("../figures/FigS12_Myositis_7PCs_DPBH_heatmap.png", psm_heatmap, width = 9, height = 9.5, bg="white")
# system("sed -i \"s/ textLength=\'[^\']*\'//\" ../figures/FigS12_Myositis_7PCs_DPBH_heatmap.svg") # Trick to make the svg file text be more easily editable


#################################################

### Create a Supplementary table 4, including the diseases included in coloc

qs <- fread("../data/qs2.tsv")

selres <- rescl[DPMUnc == 2 & Bhattacharyya %in% c(5,7,9)] # 26 diseases -- 9 myositis + 17 IMD
stcd1 <- qs[Label %in% selres$rn & grepl("Miller|Rothwell", First_Author), .(Trait, Label, First_Author, Reference, N0, N1, N)][order(Label)]
stcd1[, coloc_Label:=c("DM (M)", "DM (R)", "IIM (M)", "IIM (R)", "Anti-Jo1+ (R)", "JDM (M)", "JDM (R)", "PM (M)", "PM (R)")]
stcd2 <- qs[Label %in% selres$rn & !grepl("Miller|Rothwell", First_Author), .(Trait, Label, First_Author, Reference, N0, N1, N)][order(Label)]
stcd2[, coloc_Label:=c("CR(E)ST", "EOMG", "Felty", "HyperThy", "HypoThy", "IgG+ NMO", "JIA", "LOMG", "MPO+ AAV", "MG", "PR", "PBC", "RA", "SjS", "SSc", "SLE", "GPA")]
stcd <- rbind(stcd1, stcd2)

# Save
# fwrite(stcd, "../tables/ST4_coloc_diseases.tsv", sep = "\t")
    
#################################################

# Save DPMUnc/Bhatta results as well

# dpbh.res <- merge( qs[,.(Trait, Label, First_Author, Reference, N0, N1, N)], rescl, by.y = "rn", by.x = "Label")

# Save
# fwrite(dpbh.res, "../tables/DPMUnc_BH_res.tsv", sep="\t")



sessionInfo()
# R version 4.3.1 (2023-06-16)
# Platform: x86_64-redhat-linux-gnu (64-bit)
# Running under: Rocky Linux 8.8 (Green Obsidian)

# Matrix products: default
# BLAS/LAPACK: /usr/lib64/libopenblaso-r0.3.15.so;  LAPACK version 3.9.0

# locale:
#  [1] LC_CTYPE=en_GB.UTF-8       LC_NUMERIC=C              
#  [3] LC_TIME=en_GB.UTF-8        LC_COLLATE=en_GB.UTF-8    
#  [5] LC_MONETARY=en_GB.UTF-8    LC_MESSAGES=en_GB.UTF-8   
#  [7] LC_PAPER=en_GB.UTF-8       LC_NAME=C                 
#  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
# [11] LC_MEASUREMENT=en_GB.UTF-8 LC_IDENTIFICATION=C       

# time zone: GB
# tzcode source: system (glibc)

# attached base packages:
# [1] stats     graphics  grDevices utils     datasets  methods   base     

# other attached packages:
#  [1] magrittr_2.0.3    clue_0.3-64       mcclust_1.0.1     lpSolve_5.6.18   
#  [5] R.cache_0.16.0    ggalluvial_0.12.5 reshape2_1.4.4    ggplot2_3.4.3    
#  [9] pheatmap_1.0.12   stringr_1.5.0     data.table_1.14.8

# loaded via a namespace (and not attached):
#  [1] gtable_0.3.4       dplyr_1.1.3        compiler_4.3.1     tidyselect_1.2.0  
#  [5] Rcpp_1.0.11        cluster_2.1.4      tidyr_1.3.0        systemfonts_1.0.4 
#  [9] scales_1.2.1       R6_2.5.1           plyr_1.8.8         labeling_0.4.3    
# [13] generics_0.1.3     tibble_3.2.1       munsell_0.5.0      svglite_2.1.1     
# [17] pillar_1.9.0       RColorBrewer_1.1-3 R.utils_2.12.2     rlang_1.1.1       
# [21] utf8_1.2.3         stringi_1.7.12     cli_3.6.1          withr_2.5.0       
# [25] digest_0.6.33      grid_4.3.1         lifecycle_1.0.3    R.oo_1.25.0       
# [29] R.methodsS3_1.8.2  vctrs_0.6.3        glue_1.6.2         farver_2.1.1      
# [33] fansi_1.0.4        colorspace_2.1-0   purrr_1.0.2        tools_4.3.1       
# [37] pkgconfig_2.0.3   