#########################################
##                                     ##
##     VISUALISING COLOC RESULTS       ##
##                                     ##
#########################################

# Author: Guillermo Reales 
# Date last updated: 2023/09/12

# Background: After running coloc, we can take a closer look at the regions where we found positive results.

# This script will
# * Load coloc results and the dense-SNP datasets.
# * Generate manhattan plots of the key regions.
# * Create H4 and H3 panels.
# * Extract coloc information for main table


##########################################


## Load packages and required
library(data.table)
library(magrittr)
library(coloc)
library(annotSnpStats)
library(ggplot2)
library(cowplot)
setDTthreads(20)
setwd("/home/gr440/rds/rds-cew54-basis/Projects/myositis-IMD/code")

#######

## Load key files

coloc <- fread("../data/coloc_results-v2.tsv")

# Import mapped genes 
mg <- fread("../data/mapped.genes-v2.tsv") %>% unique
mg <- mg[,.(SNPID, pid, nearestGene)]


# Add some info to the coloc table, and extract rsids to map

coloc <- merge(coloc, mg, by="pid") # First on driver SNPs
setnames(coloc, c("SNPID", "nearestGene"), c("driver.rsid", "driver.nearestgene"))
coloc <- merge(coloc, mg, by.x="bestsnp", by.y = "pid", all.x = TRUE) # Then on candidate snps. Bear in mind that we only mapped candidate SNPs with H4 > 0.5
setnames(coloc, c("SNPID", "nearestGene"), c("bestsnp.rsid", "bestsnp.nearestgene"))


# Fix labellings according to what we know
# To get a reference
keypids <- coloc[ H4 > 0.5, unique(pid)]
c2 <- copy(coloc)
c2 <- c2[pid %in% keypids]

c2[, .(driver.rsid, driver.nearestgene)] %>% unique # 11 SNPs

## This SNP had pairwise FDR > 0.5 for the only H4 > 0.5 SNP, so we'll remove it
# rs991817 - MYL2 (SH2B3)
# c2[ driver.rsid == "rs991817", driver.nearestgene:="SH2B3"] # Just this
c2[ driver.rsid == "rs991817" & H4 > 0.5 ] #  Remove this SNP because it had pairwise_fdr > 0.5
c2 <- c2[ driver.rsid != "rs991817"]

# 1. rs7073236 - IL2RA
c2[ driver.rsid == "rs7073236", c("bestsnp.rsid", "bestsnp.nearestgene"):= list("rs706778", "IL2RA")] # No changes needed, let's just update empty slots to use bestsnp

# 2. rs11066320 - RPL6.
c2[ driver.rsid == "rs11066320", c("bestsnp.rsid", "bestsnp.nearestgene"):= list("rs597808", "SH2B3")]

# 3. rs2476601 - RSBN1 (PTPN22)
c2[ driver.rsid == "rs2476601", driver.nearestgene:= "PTPN22"]
c2[ driver.rsid == "rs2476601", c("bestsnp.rsid", "bestsnp.nearestgene"):= list("rs2476601", "PTPN22")]

# 4. rs1160542 - AFF3
c2[ driver.rsid == "rs1160542", c("bestsnp.rsid", "bestsnp.nearestgene"):= list("rs11692867","AFF3")]

# 5. rs2286896 - NAB1 (STAT4)
c2[ driver.rsid == "rs2286896", driver.nearestgene:= "STAT4"]
c2[ driver.rsid == "rs2286896", c("bestsnp.rsid", "bestsnp.nearestgene"):= list("rs3821236", "STAT4")]

# 6.  rs669607 - CMC1 (EOMES)
c2[ driver.rsid == "rs669607", driver.nearestgene:= "EOMES"]
c2[ driver.rsid == "rs669607", c("bestsnp.rsid", "bestsnp.nearestgene"):= list("rs819991", "EOMES")]

# 7. rs394378 - GARIN3 (ITK)
c2[ driver.rsid == "rs394378", driver.nearestgene:= "ITK/HAVCR2"]
c2[ driver.rsid == "rs394378", c("bestsnp.rsid", "bestsnp.nearestgene"):= list("rs163315", "ITK/HAVCR2")]

# 8. rs1855025 - CCR6
c2[ driver.rsid == "rs1855025", c("bestsnp.rsid", "bestsnp.nearestgene"):= list("rs1571878", "CCR6")] # Update empty slots to use bestsnp

# 9. rs10488631 - IRF5
c2[ driver.rsid == "rs10488631", driver.nearestgene:= "IRF5/TNPO3"]
c2[ driver.rsid == "rs10488631", c("bestsnp.rsid", "bestsnp.nearestgene"):= list("rs13236009", "IRF5/TNPO3")]

# 10. rs13277113 - BLK
c2[ driver.rsid == "rs13277113", c("bestsnp.rsid", "bestsnp.nearestgene"):= list("rs2736340", "BLK")] # Update selected SNP

# Add gene/ driver SNP label
c2[ ,dlabel:=paste( bestsnp.nearestgene, bestsnp.rsid, sep=" - ")]

# Add novelty 

c2[H4 > 0.5 & pbest.myos.region > 5e-8, .(driver.rsid, bestsnp.rsid, bestsnp.nearestgene)] %>% unique
c2[ , novel.hit :=ifelse(bestsnp.nearestgene %in% c("STAT4", "PTPN22"), "No", "Yes")] # Rm STAT4, PTPN22.

######################################

# Prepare H4 panel

sumc <- c2[pid %in% keypids, .(dlabel, trait.myos,  trait.other, H4)]

# Remove IMD without sig coloc results
nocoloc <- setdiff(sumc[ H4 < 0.5, unique(trait.other)] , sumc[ H4 > 0.5, unique(trait.other)]) 
sumc <- sumc[!trait.other %in% nocoloc]


sumc <- sumc %>% tidyr::complete(dlabel, trait.myos, trait.other ) %>% as.data.table()
sumc[, dlabel:=factor(dlabel, levels = rev(unique(dlabel)))]
sumc[, trait.other:=factor(trait.other, levels = rev(unique(trait.other)))]
sumc[,flag:=ifelse(H4>.8, "High", ifelse(H4>.5, "Med", "Low"))]
sumc <- sumc[, .(dlabel, trait.myos, trait.other, H4, flag)]


gp1 <-  ggplot(sumc[ grepl(pattern = "^DM|JDM \\(M\\)", trait.myos) ], aes(x =  trait.other, y = dlabel, colour=flag, fill = H4)) +
              geom_tile( color = "black",
                        lwd = 0.2,
                        linetype = 1) +
              scale_fill_gradient(limits = c(0,1), na.value = "white")+
              geom_tile(aes( colour = flag),
                          lwd=1, linetype = 1, data = sumc[grepl(pattern = "^DM|JDM \\(M\\)", trait.myos) & H4 > 0.25]) + # Little trick to keep tiles in place
              scale_colour_manual(values=c(High="green",Med="yellow", Low ="#FFFFFF00"), guide=guide_none()) + # And make rectangles under 0.5 transparent
              
              #scale_x_discrete(position = "top")+
              theme_minimal() +
              theme(panel.grid.major = element_blank(),
                    axis.title = element_blank(),
                    axis.text.x = element_blank(),
                    legend.position = "none",
                    plot.margin = unit(c(0, 0.2, 0, 0.3), "cm")
                    )+
              facet_grid(cols = vars(trait.myos), scales = "free", space = "free",switch = "y")+
              labs(fill = "PP")

gp2 <-  ggplot(sumc[ grepl(pattern = "PM|JDM \\(R\\)", trait.myos) ], aes(x =  trait.other, y = dlabel, colour=flag, fill = H4)) +
              geom_tile( color = "black",
                        lwd = 0.2,
                        linetype = 1) +
              scale_fill_gradient(limits = c(0,1), na.value = "white")+
              geom_tile(aes( colour = flag),
                          lwd=1, linetype = 1, data = sumc[grepl(pattern = "PM|JDM \\(R\\)", trait.myos) & H4 > 0.25]) + # Little trick to keep tiles in place
              scale_colour_manual(values=c(High="green",Med="yellow", Low ="#FFFFFF00"), guide=guide_none()) + # And make rectangles under 0.5 transparent
              
              #scale_x_discrete(position = "top")+
              theme_minimal() +
              theme(panel.grid.major = element_blank(),
                    axis.title = element_blank(),
                    axis.text.x = element_blank(),
                    legend.position = "none",
                    plot.margin = unit(c(0, 0.2, 0, 0.3), "cm")
                    )+
              facet_grid(cols = vars(trait.myos), scales = "free", space = "free",switch = "y")+
              labs(fill = "PP")

gp3 <-  ggplot(sumc[ grepl(pattern = "IIM|Jo1", trait.myos) ], aes(x =  trait.other, y = dlabel, colour=flag, fill = H4)) +
              geom_tile( color = "black",
                        lwd = 0.2,
                        linetype = 1) +
              scale_fill_gradient(limits = c(0,1), na.value = "white")+
              geom_tile(aes( colour = flag),
                          lwd=1, linetype = 1, data = sumc[grepl(pattern = "IIM|Jo1", trait.myos) & H4 > .5]) +
              scale_colour_manual(values=c(High="green",Med="yellow"), guide=guide_none()) +
              #scale_x_discrete(position = "top")+
              theme_minimal() +
              theme(panel.grid.major = element_blank(),
                    axis.title = element_blank(),
                    axis.text.x = element_text(angle = 270, hjust=0, vjust = 0.5, size = 11),
                    legend.position = "bottom",
                    plot.margin = unit(c(0, 0.2, 0, 0.3), "cm")
                    )+
              guides(fill = guide_colorbar(title.vjust = 0.8))+
              facet_grid(cols = vars(trait.myos), scales = "free", space = "free",switch = "y")+
              labs(fill = "PP")

gp <- plot_grid(gp1, gp2, gp3, nrow = 3, rel_heights = c(0.75,0.75,1.4))
ggsave("../figures/driverSNP_H4_rnocoloc.png", gp, height =7 , width = 6.8, bg="white")

###############################


# Prepare H3 panel

sumc.h3 <- c2[pid %in% keypids, .(dlabel, trait.myos,  trait.other, H3)]
sumc.h3 <- sumc.h3 %>% tidyr::complete(dlabel, trait.myos, trait.other ) %>% as.data.table()
sumc.h3[, dlabel:=factor(dlabel, levels = rev(unique(dlabel)))]
sumc.h3[, trait.other:=factor(trait.other, levels = rev(unique(trait.other)))]
sumc.h3[,flag:=ifelse(H3>.8, "High", ifelse(H3>.5, "Med", "Low"))]
sumc.h3 <- sumc.h3[, .(dlabel, trait.myos, trait.other, H3, flag)]


gp1.h3 <-  ggplot(sumc.h3[ grepl(pattern = "^DM|JDM \\(M\\)", trait.myos) ], aes(x =  trait.other, y = dlabel, colour=flag, fill = H3)) +
              geom_tile( color = "black",
                        lwd = 0.2,
                        linetype = 1) +
              scale_fill_gradient(limits = c(0,1), na.value = "white")+
              geom_tile(aes( colour = flag),
                          lwd=1, linetype = 1, data = sumc.h3[grepl(pattern = "^DM|JDM \\(M\\)", trait.myos) & H3 > 0.25]) + # Little trick to keep tiles in place
              scale_colour_manual(values=c(High="green",Med="yellow", Low ="#FFFFFF00"), guide=guide_none()) + # And make rectangles under 0.5 transparent
              
              #scale_x_discrete(position = "top")+
              theme_minimal() +
              theme(panel.grid.major = element_blank(),
                    axis.title = element_blank(),
                    axis.text.x = element_blank(),
                    legend.position = "none",
                    plot.margin = unit(c(0, 0.2, 0, 0.3), "cm")
                    )+
              facet_grid(cols = vars(trait.myos), scales = "free", space = "free",switch = "y")+
              labs(fill = "PP")

gp2.h3 <-  ggplot(sumc.h3[ grepl(pattern = "PM|JDM \\(R\\)", trait.myos) ], aes(x =  trait.other, y = dlabel, colour=flag, fill = H3)) +
              geom_tile( color = "black",
                        lwd = 0.2,
                        linetype = 1) +
              scale_fill_gradient(limits = c(0,1), na.value = "white")+
              geom_tile(aes( colour = flag),
                          lwd=1, linetype = 1, data = sumc.h3[grepl(pattern = "PM|JDM \\(R\\)", trait.myos) & H3 > 0.25]) + # Little trick to keep tiles in place
              scale_colour_manual(values=c(High="green",Med="yellow", Low ="#FFFFFF00"), guide=guide_none()) + # And make rectangles under 0.5 transparent
              
              #scale_x_discrete(position = "top")+
              theme_minimal() +
              theme(panel.grid.major = element_blank(),
                    axis.title = element_blank(),
                    axis.text.x = element_blank(),
                    legend.position = "none",
                    plot.margin = unit(c(0, 0.2, 0, 0.3), "cm")
                    )+
              facet_grid(cols = vars(trait.myos), scales = "free", space = "free",switch = "y")+
              labs(fill = "PP")

gp3.h3 <-  ggplot(sumc.h3[ grepl(pattern = "IIM|Jo1", trait.myos) ], aes(x =  trait.other, y = dlabel, colour=flag, fill = H3)) +
              geom_tile( color = "black",
                        lwd = 0.2,
                        linetype = 1) +
              scale_fill_gradient(limits = c(0,1), na.value = "white")+
              geom_tile(aes( colour = flag),
                          lwd=1, linetype = 1, data = sumc.h3[grepl(pattern = "IIM|Jo1", trait.myos) & H3 > .25]) +
              scale_colour_manual(values=c(High="green",Med="yellow",Low ="#FFFFFF00"), guide=guide_none()) +
              #scale_x_discrete(position = "top")+
              theme_minimal() +
              theme(panel.grid.major = element_blank(),
                    axis.title = element_blank(),
                    axis.text.x = element_text(angle = 270, hjust=0, vjust = 0.5, size = 11),
                    legend.position = "bottom",
                    plot.margin = unit(c(0, 0.2, 0, 0.3), "cm")
                    )+
              guides(fill = guide_colorbar(title.vjust = 0.8))+
              facet_grid(cols = vars(trait.myos), scales = "free", space = "free",switch = "y")+
              labs(fill = "PP")

gp.h3 <- plot_grid(gp1.h3, gp2.h3, gp3.h3, nrow = 3, rel_heights = c(0.75,0.75,1.4))

# Save
ggsave("../figures/driverSNP_H3.png", gp.h3, height =7 , width = 8, bg="white")






###############################

## Prepare dataset for main table

## First we'll get the p-values for the table

files = list(   `DM (M)` = "DMY_Miller_26291516_1-hg38.tsv.gz",
                `DM (R)` = "DMY_Rothwell_up_1-hg38.tsv.gz",
                `JDM (M)` = "JDM_Miller_26291516_1-hg38.tsv.gz",
                `JDM (R)` = "JDM_Rothwell_up_1-hg38.tsv.gz",
                `Anti-Jo1+ (R)` = "JO1M_Rothwell_up_1-hg38.tsv.gz",
                `IIM (M)` = "MYO_Miller_26291516_1-hg38.tsv.gz",
                `IIM (R)` = "IIM_Rothwell_up_1-hg38.tsv.gz",
                `PM (M)` =   "PM_Miller_26291516_1-hg38.tsv.gz",
                `PM (R)` = "PM_Rothwell_up_1-hg38.tsv.gz") %>%
    lapply(.,  function(f) (file.path("~/rds/rds-cew54-basis/02-Processed",f)))


data=lapply(files, fread)
studies <- names(data)

mg.top <- mg[SNPID %in% c2[H4 > 0.5, unique(bestsnp.rsid)], .(SNPID, pid)]

# Extract required SNPs
d2 <- lapply(1:9, function(i){
        x <- data[[i]]
        x[, pid := paste(CHR38, BP38, sep=":")]
        x  <- x[ pid %in% mg.top$pid , .(pid, SNPID, REF, ALT, BETA, SE, P)]
        x[, study := studies[i]]
        x

}) %>% rbindlist

mg.top <- merge(mg.top, d2[,.(pid, P, study)], by = "pid")

c3 <- merge(c2, mg.top[, .(SNPID, P, study)], by.x=c("trait.myos", "bestsnp.rsid"), by.y=c("study", "SNPID"))

cts <- c3[H4 > 0.5, .(pid, driver.rsid, bestsnp.rsid, bestsnp.nearestgene, P, trait.myos, trait.other, pairwise_fdr, H4, novel.hit)][order(bestsnp.nearestgene, bestsnp.rsid, trait.myos, trait.other)]
names(cts) <- c("pid", "Driver SNP", "Top candidate SNP", "Open Targets candidate gene", "Top SNP P-value", "Myositis", "IMD", "Pairwise FDR", "H4", "Novel")

fwrite(cts, "../tables/MT_coloc_results.tsv", sep = "\t")


##########################################

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
# [1] cowplot_1.1.1      ggplot2_3.4.3      annotSnpStats_0.99 snpStats_1.50.0   
# [5] Matrix_1.5-4.1     survival_3.5-5     coloc_5.2.2        magrittr_2.0.3    
# [9] data.table_1.14.8 

# loaded via a namespace (and not attached):
#  [1] gtable_0.3.4        jsonlite_1.8.7      dplyr_1.1.3        
#  [4] compiler_4.3.1      crayon_1.5.2        tidyselect_1.2.0   
#  [7] Rcpp_1.0.11         gridExtra_2.3       tidyr_1.3.0        
# [10] susieR_0.12.35      splines_4.3.1       scales_1.2.1       
# [13] lattice_0.21-8      R6_2.5.1            plyr_1.8.8         
# [16] labeling_0.4.3      generics_0.1.3      mixsqp_0.3-48      
# [19] BiocGenerics_0.46.0 viridis_0.6.4       tibble_3.2.1       
# [22] munsell_0.5.0       pillar_1.9.0        rlang_1.1.1        
# [25] utf8_1.2.3          reshape_0.8.9       viridisLite_0.4.2  
# [28] cli_3.6.1           withr_2.5.0         zlibbioc_1.46.0    
# [31] grid_4.3.1          irlba_2.3.5.1       lifecycle_1.0.3    
# [34] vctrs_0.6.3         glue_1.6.2          farver_2.1.1       
# [37] fansi_1.0.4         colorspace_2.1-0    purrr_1.0.2        
# [40] tools_4.3.1         matrixStats_1.0.0   pkgconfig_2.0.3 




###############################

## Recreate data object for Manhattan plots from dense-SNP data.

## NOTE: We might not use the coloc Manhattan plots in the end, so we might remove part of this code.

# Load files, as before

# # Take the opportunity to create proper labels
# files = list(   `DM (M)` = "DMY_Miller_26291516_1-hg38.tsv.gz",
#                 `DM (R)`= "DMY_Rothwell_up_1-hg38.tsv.gz",
#                 `JDM (M)` = "JDM_Miller_26291516_1-hg38.tsv.gz",
#                 `JDM (R)` = "JDM_Rothwell_up_1-hg38.tsv.gz",
#                 `Anti-Jo1+ (R)` = "JO1M_Rothwell_up_1-hg38.tsv.gz",
#                 `IIM (M)` = "MYO_Miller_26291516_1-hg38.tsv.gz",
#                 `IIM (R)` = "IIM_Rothwell_up_1-hg38.tsv.gz",
#                 `PM (M)` =   "PM_Miller_26291516_1-hg38.tsv.gz",
#                 `PM (R)` = "PM_Rothwell_up_1-hg38.tsv.gz",

#                 `CR(E)ST` = "CREST_FinnGen_FinnGenR7_1-hg38.tsv.gz", 
#                 EOMG = "MYGEO_Renton_25643325_1-hg38.tsv.gz",
#                 Felty =  "FELTY_FinnGen_FinnGenR7_1-hg38.tsv.gz", 
#                 HyperThy = "20002_1225_PanUKBB_PanUKBBR2_1-hg38.tsv.gz",
#                 HypoThy = "E4_HYTHY_AI_STRICT_FinnGen_FinnGenR7_1-hg38.tsv.gz", 
#                 `IgG+ NMO` =  "NMOIGGp_Estrada_29769526_1-hg38.tsv.gz", 
#                 JIA = "JIA_LopezIsac_33106285_1-hg38.tsv.gz",
#                 LOMG = "MYGLO_Renton_25643325_1-hg38.tsv.gz",
#                 `MPO+ AAV` = "AAVMPO_Wong_up_1-hg38.tsv.gz", 
#                 MG = "MYG_Chia_35074870_1-hg38.tsv.gz",
#                 PR = "M13_PALINDROMIC_FinnGen_FinnGenR7_1-hg38.tsv.gz",
#                 PBC = "CHIRBIL_PRIM_FinnGen_FinnGenR7_1-hg38.tsv.gz", 
#                 RA = "M13_RHEUMA_FinnGen_FinnGenR7_1-hg38.tsv.gz", 
#                 SjS =   "SJOS_Lessard_up_1-hg38.tsv.gz",
#                 SSc  = "SSC_LopezIsac_31672989_1-hg38.tsv.gz",
#                 SLE = "M13_SLE_FinnGen_FinnGenR7_1-hg38.tsv.gz",
#                 GPA = "M13_WEGENER_FinnGen_FinnGenR7_1-hg38.tsv.gz") %>% 
#     lapply(.,  function(f) (file.path("~/rds/rds-cew54-basis/02-Processed",f)))
# stopifnot(all(file.exists(unlist(files))))

# data=lapply(files, fread)


# ## some datasets are not dense datasets. 
# names(data)[c(10, 12:14, 20:22, 25:26)] = c("CR(E)ST.local", "Felty.local","HyperThy.local", "HypoThy.local", "PR.local", "PBC.local", "RA.local" ,"SLE.local","GPA.local")

# dir.create("../data/fg_sumstats")

# if(!file.exists("../data/fg_sumstats/finngen_R7_CREST.gz"))
#     system("wget https://storage.googleapis.com/finngen-public-data-r7/summary_stats/finngen_R7_CREST.gz -O ../data/fg_sumstats/finngen_R7_CREST.gz")
# if(!file.exists("../data/fg_sumstats/finngen_R7_FELTY.gz"))
#     system("wget https://storage.googleapis.com/finngen-public-data-r7/summary_stats/finngen_R7_FELTY.gz -O ../data/fg_sumstats/finngen_R7_FELTY.gz")

# # Note: this PanUKBB file will need some work to reformat, see below
# if(!file.exists("../data/fg_sumstats/20002_1225_PanUKBB_1-hg38.tsv.gz")){
#     system("wget https://pan-ukb-us-east-1.s3.amazonaws.com/sumstats_flat_files/categorical-20002-both_sexes-1225.tsv.bgz -O ../data/fg_sumstats/20002_1225_PanUKBB_PanUKBBR2_1.bgz")
#     system("Rscript processing_panUKBB.R") # This script will prepare the PanUKBB file to be used by coloc
# }   
# if(!file.exists("../data/fg_sumstats/finngen_R7_E4_HYTHY_AI_STRICT.gz"))
#     system("wget https://storage.googleapis.com/finngen-public-data-r7/summary_stats/finngen_R7_E4_HYTHY_AI_STRICT.gz -O ../data/fg_sumstats/finngen_R7_E4_HYTHY_AI_STRICT.gz")
# if(!file.exists("../data/fg_sumstats/finngen_R7_M13_PALINDROMIC.gz"))
#     system("wget https://storage.googleapis.com/finngen-public-data-r7/summary_stats/finngen_R7_M13_PALINDROMIC.gz -O ../data/fg_sumstats/finngen_R7_M13_PALINDROMIC.gz")
# if(!file.exists("../data/fg_sumstats/finngen_R7_CHIRBIL_PRIM.gz"))
#     system("wget https://storage.googleapis.com/finngen-public-data-r7/summary_stats/finngen_R7_CHIRBIL_PRIM.gz -O ../data/fg_sumstats/finngen_R7_CHIRBIL_PRIM.gz")
# if(!file.exists("../data/fg_sumstats/finngen_R7_M13_RHEUMA.gz"))
#     system("wget https://storage.googleapis.com/finngen-public-data-r7/summary_stats/finngen_R7_M13_RHEUMA.gz -O ../data/fg_sumstats/finngen_R7_M13_RHEUMA.gz")
# if(!file.exists("../data/fg_sumstats/finngen_R7_M13_SLE.gz"))
#     system("wget https://storage.googleapis.com/finngen-public-data-r7/summary_stats/finngen_R7_M13_SLE.gz -O ../data/fg_sumstats/finngen_R7_M13_SLE.gz")
# if(!file.exists("../data/fg_sumstats/finngen_R7_M13_WEGENER.gz"))
#     system("wget https://storage.googleapis.com/finngen-public-data-r7/summary_stats/finngen_R7_M13_WEGENER.gz -O ../data/fg_sumstats/finngen_R7_M13_WEGENER.gz")


# # Incorporate new files into data


# crest=fread("../data/fg_sumstats/finngen_R7_CREST.gz")
# crest$pid=paste(crest[["#chrom"]],crest$pos,sep=":")
# table(crest$pid %in% snps$pid38)
# setnames(crest, c("#chrom","pos"), c("CHR38","BP38"))
# data$`CR(E)ST`=crest[,.(pid, CHR38, BP38, REF=ref, ALT=alt, BETA=beta, SE=sebeta, P=pval)]

# Felty=fread("../data/fg_sumstats/finngen_R7_FELTY.gz")
# Felty$pid=paste(Felty[["#chrom"]],Felty$pos,sep=":")
# table(Felty$pid %in% snps$pid38)
# setnames(Felty, c("#chrom","pos"), c("CHR38","BP38"))
# data$Felty=Felty[,.(pid, CHR38, BP38, REF=ref, ALT=alt, BETA=beta, SE=sebeta, P=pval)]

# hyperthy=fread("../data/fg_sumstats/20002_1225_PanUKBB_1-hg38.tsv.gz")
# hyperthy[, pid := paste(CHR38,BP38,sep=":")]
# table(hyperthy$pid %in% snps$pid38)
# data$HyperThy=hyperthy[,.(pid, CHR38, BP38, REF, ALT, BETA, SE, P)]

# hypothy=fread("../data/fg_sumstats/finngen_R7_E4_HYTHY_AI_STRICT.gz")
# hypothy$pid=paste(hypothy[["#chrom"]],hypothy$pos,sep=":")
# table(hypothy$pid %in% snps$pid38)
# setnames(hypothy, c("#chrom","pos"), c("CHR38","BP38"))
# data$HypoThy=hypothy[,.(pid, CHR38, BP38, REF=ref, ALT=alt, BETA=beta, SE=sebeta, P=pval)]

# PR=fread("../data/fg_sumstats/finngen_R7_M13_PALINDROMIC.gz")
# PR$pid=paste(PR[["#chrom"]],PR$pos,sep=":")
# table(PR$pid %in% snps$pid38)
# setnames(PR, c("#chrom","pos"), c("CHR38","BP38"))
# data$PR=PR[,.(pid, CHR38, BP38, REF=ref, ALT=alt, BETA=beta, SE=sebeta, P=pval)]

# PBC=fread("../data/fg_sumstats/finngen_R7_CHIRBIL_PRIM.gz")
# PBC$pid=paste(PBC[["#chrom"]],PBC$pos,sep=":")
# table(PBC$pid %in% snps$pid38)
# setnames(PBC, c("#chrom","pos"), c("CHR38","BP38"))
# data$PBC=PBC[,.(pid, CHR38, BP38, REF=ref, ALT=alt, BETA=beta, SE=sebeta, P=pval)]

# RA=fread("../data/fg_sumstats/finngen_R7_M13_RHEUMA.gz")
# RA$pid=paste(RA[["#chrom"]],RA$pos,sep=":")
# table(RA$pid %in% snps$pid38)
# setnames(RA, c("#chrom","pos"), c("CHR38","BP38"))
# data$RA=RA[,.(pid, CHR38, BP38, REF=ref, ALT=alt, BETA=beta, SE=sebeta, P=pval)]

# SLE=fread("../data/fg_sumstats/finngen_R7_M13_SLE.gz")
# SLE$pid=paste(SLE[["#chrom"]],SLE$pos,sep=":")
# table(SLE$pid %in% snps$pid38)
# setnames(SLE, c("#chrom","pos"), c("CHR38","BP38"))
# data$SLE=SLE[,.(pid, CHR38, BP38, REF=ref, ALT=alt, BETA=beta, SE=sebeta, P=pval)]

# GPA=fread("../data/fg_sumstats/finngen_R7_M13_WEGENER.gz")
# GPA$pid=paste(GPA[["#chrom"]],GPA$pos,sep=":")
# table(GPA$pid %in% snps$pid38)
# setnames(GPA, c("#chrom","pos"), c("CHR38","BP38"))
# data$GPA=GPA[,.(pid, CHR38, BP38, REF=ref, ALT=alt, BETA=beta, SE=sebeta, P=pval)]


# index <- coloc # we were using index instead of coloc. so let's simply use this other name

# # Now explore the SNPs in more detail
# index[ H4>.5, unique(pid)]
# # 11 SNPs
# index[ H4>.5, unique(trait.other)]
# # We have coloc associations with 11 (out of 17) IMDs
# # These are
# #  [1] "RA"       "HypoThy"  "JIA"      "HyperThy" "MPO+ AAV" "PBC"     
# #  [7] "MG"       "LOMG"     "SSc"      "SjS"      "SLE"     

# ## plot these signals
# plotter=function(pid,w=1e+6) {
#     chr=sub(":.*","",pid)  %>% as.numeric()
#     bp=sub(".*:","",pid)  %>% as.numeric()
#     st=bp-w
#     en=bp+w
#     wh=which(index$pid==pid & (index$H4>.5))
#     print(index[wh])
#     rsid = index[wh , unique(driver.rsid)]
#     traits=unique(c(index$trait.myos[wh],index$trait.other[wh]))
#     dp=lapply(data[traits], function(d)
#         d[CHR38==chr & BP38>st & BP38<en & !is.na(SE) & !is.na(BETA) & !duplicated(pid),
#           .(pid,CHR38,BP38,BETA,SE,P)])
#     for(i in seq_along(traits))
#         dp[[i]]$trait=traits[i]
#     dp %<>% rbindlist()
#     ggplot(dp, aes(x=BP38, y=-log10(P))) + 
#             geom_point() + 
#             facet_grid(trait~.,scales="free_y") + 
#             geom_vline(xintercept=bp,col="red")+
#             ggtitle(paste0(pid, " / ", rsid))+
#             theme_cowplot() + 
#             theme(axis.title.x = element_blank(),
#                   strip.background = element_rect(colour="black", fill = "white"),
#                   )
# }

# plots <- lapply(index[ H4>.5, unique(pid)], plotter)
# names(plots) <- index[ H4>.5, unique(pid)]

# index[ H4>.5, unique(pid)]

# index[ H4 > 0.5, .(bestsnp, bestsnp.novel)] %>% unique


# ggsave("../figures/coloc_chr1.png", plots$`1:113834946`, height = 12, width = 8, bg="white")
# ggsave("../figures/coloc_chr2_1.png", plots$`2:100215693`, height = 5, width = 8, bg="white")
# ggsave("../figures/coloc_chr2_2.png", plots$`2:190670850`, height = 5, width = 8, bg="white")
# ggsave("../figures/coloc_chr3.png", plots$`3:28029953`, height = 5, width = 8, bg="white")
# ggsave("../figures/coloc_chr5.png", plots$`5:157185077`,height = 5, width = 8, bg="white")
# ggsave("../figures/coloc_chr6.png", plots$`6:167124106`, height = 6, width = 8, bg="white")
# ggsave("../figures/coloc_chr7.png", plots$`7:128954129`, height = 14, width = 8, bg="white")
# ggsave("../figures/coloc_chr8.png", plots$`8:11491677`, height = 8, width = 8, bg="white")
# ggsave("../figures/coloc_chr8.png", plots$`10:6064589`, height = 5, width = 8, bg="white")
# ggsave("../figures/coloc_chr12_1.png", plots$`12:110972733`, height = 5, width = 8, bg="white")
# ggsave("../figures/coloc_chr12_2.png", plots$`12:112468611`, height = 6, width = 8, bg="white")

# # ggsave("../figures/coloc_chr11.png", plots$`11:64329761`, height = 8, width = 8, bg="white")
# # ggsave("../figures/coloc_chr17.png", plots$`17:75373341`, height = 5, width = 8, bg="white")
# # ggsave("../figures/coloc_chr1.png", plots$`1:113834946`, height = 11, width = 8, bg="white")
# # ggsave("../figures/coloc_chr4.png", plots$`4:122194347`, height = 6, width = 8, bg="white")
# # ggsave("../figures/coloc_chr6.png", plots$`6:167124106`, height = 6, width = 8, bg="white")
# # ggsave("../figures/coloc_chr7_2.png", plots$`7:37397251`, height = 5, width = 8, bg="white")