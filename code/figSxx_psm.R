##################################################
# Creating PSM with DPMUnc and Bhatta clustering #
##################################################

# Date: 04/08/2023
# Guillermo Reales

# NOTE: We needed to load r-4.0.2 since r-4.1.3 didn't have png/cairo capabilities

setwd("/home/gr440/rds/rds-cew54-basis/Projects/myositis-IMD/code")

# Load packages
library(data.table)
setDTthreads(15)
library(magrittr)
library(R.cache)
library(pheatmap)
library(mcclust)
library(clue)

# Here we adapted helper functions from Kath's repo (https://github.com/nichollskc/clusteringPublicGeneExpr/blob/26e1fad109009563c8a0dc1254a54f796176f6c5/scripts/psm_plots.R)
# Also from scripts/utils.R (https://github.com/nichollskc/clusteringPublicGeneExpr/blob/26e1fad109009563c8a0dc1254a54f796176f6c5/scripts/utils.R)


# Define some variables


palette <- c("#77AADD", "#000000", "#9E0142", "#D53E4F", "#F46D43", "#FEE08B", "#ABDDA4", "#66C2A5",             "#5E4FA2", 
             "#FDAE61", "#E6F598", "#771155", "#AA4488", "#CC99BB", "#114477", "#774411", "#EEEEEE", "#117777", "#117744", 
             "#44AA77", "#88CCAA", "#777711", "#44AAAA", "#AAAA44", "#77CCCC", "#DDDD77", "#AA7744", "#DDAA77", "#771122", "#AA4455", "#DD7788")
exp="Myo_7PC"
burnin = 250000
input_dir="../data/DPMUnc_results/"
outdir="../figures/"

#"#3288BD",

#####

# Helper functions


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

####################################################

# Process

    message("Loading data...")
    datasets = dir(input_dir)[grepl(exp, dir(input_dir))] # Little adaptation to original function. From experiment name we retrieve all directories for different seeds.
    datasets = paste0(input_dir, datasets)

    obsData = read.table(paste0("../data/", exp, "_Delta.tsv"),
                         header=1, row.names=1, quote="", sep="\t")
    obsVars = read.table(paste0("../data/", exp, "_Var.tsv"),
                         header=1, row.names=1, quote="", sep="\t")

    bht = fread("../data/bhattacharyya_clustering.tsv")
    
    
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

    calls = minbinder(bigpsm, method = "comp") ## calls
    
    # Make some names bold
    myobold <- grep("myositis|IIM", rownames(bigpsm), ignore.case = TRUE, value = TRUE)

    # # Get cluster colours automatically
    # annotations = get_ann_colors(calls$cl, bhv, obsData)

    # # Note that here we call it mclust for legacy reasons, but the calls are actually Bhattacharyya
    # names(annotations$ann)  <- c("DPMUnc", "Bhattacharyya")
    # annotations = list(ann = annotations$ann, colors = list(DPMUnc = annotations$colors$cluster, Bhattacharyya = annotations$colors$mclust ))

    # An alternative annotation strategy -- get_ann_colours is complicated and failed to properly match Bhattacharyya IMD to their corresponding clusters

    ann <- data.table(Label = names(calls$cl), DPMUnc = calls$cl)
    ann <- merge(ann, bht, by = "Label")
    
    ann <- data.frame(ann[, 2:3], row.names = ann$Label)
    dpcol <- palette[1:max(ann$DPMUnc)]
    names(dpcol)  <- as.character(1:max(ann$DPMUnc))
    bhcol <- palette[1:max(ann$Bhattacharyya)]
    names(bhcol)  <- as.character(1:max(ann$Bhattacharyya))

    annotations <- list(ann = ann, colors = list(DPMUnc = dpcol, Bhattacharyya = bhcol))

    # Same hclust calculation that minbinder does prior to choosing optimal cut point
    hclust.comp <- hclust(as.dist(1 - bigpsm), method = "complete")
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
                           width=9,
                           height=9.5,
                           labels_row = make_bold_names(bigpsm, rownames, myobold),
                           filename=paste0(outdir, exp, "_DPMUnc_Bh_psm_heatmap.png"))

svg(paste0(outdir, exp, "_DPMUnc_Bh_psm_heatmap.svg"), width = 9, height = 9.5)
  psm_heatmap
dev.off()

sessionInfo()
# R version 4.0.2 (2020-06-22)
# Platform: x86_64-pc-linux-gnu (64-bit)
# Running under: Rocky Linux 8.7 (Green Obsidian)

# Matrix products: default
# BLAS:   /usr/local/software/spack/spack-0.11.2/opt/spack/linux-rhel7-x86_64/gcc-5.4.0/r-4.0.2-xyx46xbuw2lmofomvrkwuty5rlez6to6/rlib/R/lib/libRblas.so
# LAPACK: /usr/local/software/spack/spack-0.11.2/opt/spack/linux-rhel7-x86_64/gcc-5.4.0/r-4.0.2-xyx46xbuw2lmofomvrkwuty5rlez6to6/rlib/R/lib/libRlapack.so

# locale:
#  [1] LC_CTYPE=en_GB.UTF-8       LC_NUMERIC=C              
#  [3] LC_TIME=en_GB.UTF-8        LC_COLLATE=en_GB.UTF-8    
#  [5] LC_MONETARY=en_GB.UTF-8    LC_MESSAGES=en_GB.UTF-8   
#  [7] LC_PAPER=en_GB.UTF-8       LC_NAME=C                 
#  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
# [11] LC_MEASUREMENT=en_GB.UTF-8 LC_IDENTIFICATION=C       

# attached base packages:
# [1] stats     graphics  grDevices utils     datasets  methods   base     

# other attached packages:
# [1] clue_0.3-64       mcclust_1.0.1     lpSolve_5.6.18    pheatmap_1.0.12  
# [5] R.cache_0.16.0    magrittr_2.0.3    data.table_1.14.8

# loaded via a namespace (and not attached):
#  [1] cluster_2.1.4      munsell_0.5.0      colorspace_2.1-0   R6_2.5.1          
#  [5] rlang_1.0.6        stringr_1.5.0      tools_4.0.2        grid_4.0.2        
#  [9] gtable_0.3.1       R.oo_1.25.0        cli_3.6.0          digest_0.6.31     
# [13] lifecycle_1.0.3    farver_2.1.1       purrr_1.0.1        RColorBrewer_1.1-3
# [17] vctrs_0.5.2        R.utils_2.12.2     glue_1.6.2         stringi_1.7.12    
# [21] compiler_4.0.2     scales_1.2.1       R.methodsS3_1.8.2