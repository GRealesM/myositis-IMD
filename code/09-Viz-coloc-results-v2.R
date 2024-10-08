#########################################
##                                     ##
##     VISUALISING COLOC RESULTS       ##
##                                     ##
#########################################

# Author: Guillermo Reales 
# Date last updated: 2023/04/30

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

coloc <- fread("../data/coloc_results-v3.tsv")

# Import mapped genes 
mg <- fread("../data/mapped.genes-v3.tsv") %>% unique
mg <- mg[,.(SNPID, pid, nearestGene)]


# Add some info to the coloc table, and extract rsids to map

coloc <- merge(coloc, mg, by="pid", all.x = TRUE) # First on driver SNPs
setnames(coloc, c("SNPID", "nearestGene"), c("driver.rsid", "driver.nearestgene"))
coloc <- merge(coloc, mg, by.x="bestsnp", by.y = "pid", all.x = TRUE) # Then on candidate snps. Bear in mind that we only mapped candidate SNPs with H4 > 0.5
setnames(coloc, c("SNPID", "nearestGene"), c("bestsnp.rsid", "bestsnp.nearestgene"))

# Prepare Supplementary table 5
st5 <- copy(coloc)
st5 <- st5[, .(pid,  chr, bp, driver.rsid, driver.nearestgene, 
                  trait.myos, trait.other, pdriver.myos, fdr.myos, fdr.other, 
                  pairwise_fdr, nsnps, driver, H0, H1, 
                  H2, H3, H4, bestsnp, bestsnp.rsid, 
                  bestsnp.nearestgene, bestsnp.pp, pbest.myos, pbest.other, pbest.myos.region, 
                  pbest.other.region)]
names(st5) <- c("Driver SNP coordinates","Chromosome","Base pair","Driver SNP rsid", "Driver SNP nearest gene (OpenTargets Genetics)", 
"Myositis trait", "IMD","Driver SNP P-value (Myositis dataset)", "Myositis FDR","IMD FDR",
"Pairwise FDR", "N SNPs in region","Features driven by driver SNP","H0","H1",
"H2","H3","H4","Candidate SNP coordinates","Candidate SNP rsid", 
"Candidate SNP nearest gene (OpenTargets Genetics)","Candidate SNP PP", "Candidate SNP P-value (Myositis dataset)","Candidate SNP P-value (IMD dataset)","Lowest P-value in the region (Myositis dataset)",
"Lowest P-value in the region (IMD dataset)")

# fwrite(st5, "../tables/ST5_full_coloc_results.tsv", sep ="\t")

# Fix labellings according to what we know
# To get a reference
keypids <- coloc[ H4 > 0.5, unique(pid)]
c2 <- copy(coloc)
c2 <- c2[pid %in% keypids]

c2[, .(driver.rsid, driver.nearestgene)] %>% unique # 8 SNPs
c2[ H4 > 0.5 & pbest.myos.region > 5e-8, .(pid, driver.rsid, bestsnp.rsid, driver.nearestgene, bestsnp.nearestgene)] %>% unique

## This SNP had pairwise FDR > 0.5 for the only H4 > 0.5 SNP, so we'll remove it
# 1. rs7072793 - IL2RA
c2[driver.rsid == "rs7072793"] # Nothing to do
c2[ driver.rsid == "rs7072793", c("bestsnp.rsid", "bestsnp.nearestgene"):= list("rs7072793", "IL2RA")]

# 2. rs2476601 - RSBN1 (PTPN22)
c2[ driver.rsid == "rs2476601", driver.nearestgene:= "PTPN22"]
c2[ driver.rsid == "rs2476601", c("bestsnp.rsid", "bestsnp.nearestgene"):= list("rs2476601", "PTPN22")]

# 3. rs5754217 - RIMBP3C (UBE2L3)
c2[ driver.rsid == "rs5754217", driver.nearestgene:= "UBE2L3"]
c2[ driver.rsid == "rs5754217", c("bestsnp.rsid", "bestsnp.nearestgene"):= list("rs5754217", "UBE2L3")]

# 4. rs10196612 - PLCL1
c2[driver.rsid == "rs10196612" & H4 > 0.5]
c2[ driver.rsid == "rs10196612", c("bestsnp.rsid", "bestsnp.nearestgene"):= list("rs6738825", "PLCL1")]

# 5.  rs669607 - CMC1 (EOMES)
c2[driver.rsid == "rs669607" & H4 > 0.5]
c2[ driver.rsid == "rs669607", driver.nearestgene:= "EOMES"]
c2[ driver.rsid == "rs669607", c("bestsnp.rsid", "bestsnp.nearestgene"):= list("rs819991", "EOMES")]

# 6. rs394378 - GARIN3 (ITK)
c2[driver.rsid == "rs394378" & H4 > 0.5]
c2[ driver.rsid == "rs394378", driver.nearestgene:= "ITK/HAVCR2"]
c2[ driver.rsid == "rs394378", c("bestsnp.rsid", "bestsnp.nearestgene"):= list("rs163315", "ITK/HAVCR2")] # Will use rs163315 

# 9. rs10488631 - IRF5
c2[driver.rsid == "rs10488631" & H4 > 0.5]
c2[ driver.rsid == "rs10488631", driver.nearestgene:= "IRF5/TNPO3"]
c2[ driver.rsid == "rs10488631", c("bestsnp.rsid", "bestsnp.nearestgene"):= list("rs13236009", "IRF5/TNPO3")] # will use rs13236009

# 10. rs13277113 - BLK
c2[driver.rsid == "rs13277113" & H4 > 0.5]
c2[ driver.rsid == "rs13277113", c("bestsnp.rsid", "bestsnp.nearestgene"):= list("rs2736340", "BLK")] # will use rs2736340

# Add gene/ driver SNP label
c2[ ,dlabel:=paste( bestsnp.nearestgene, bestsnp.rsid, sep=" - ")]

# Add novelty 

c2[H4 > 0.5 & pbest.myos.region > 5e-8, .(driver.rsid, bestsnp.rsid, bestsnp.nearestgene)] %>% unique
c2[ , novel.hit :=ifelse(bestsnp.nearestgene %in% c("PTPN22"), "No", "Yes")] # Rm PTPN22, as it was GW-sig in Miller 2015.

######################################

# Prepare H4 panel

sumc <- c2[pid %in% keypids , .(dlabel, trait.myos,  trait.other, H4)]

# Remove IMD without sig coloc results
nocoloc <- setdiff(sumc[ H4 < 0.5, unique(trait.other)] , sumc[ H4 > 0.5, unique(trait.other)]) 
sumc <- sumc[!trait.other %in% nocoloc]

sumc <- sumc %>% tidyr::complete(dlabel, trait.myos, trait.other ) %>% as.data.table()
sumc[, dlabel:=factor(dlabel, levels = rev(unique(dlabel)))]
sumc[, trait.other:=factor(trait.other, levels = rev(unique(trait.other)))]
sumc[,flag:=ifelse(H4>.8, "High", ifelse(H4>.5, "Med", "Low"))]
sumc <- sumc[, .(dlabel, trait.myos, trait.other, H4, flag)]
unique(sumc$trait.myos)

# Remove empty rows
rtr <- sumc[, all(is.na(H4)), by = c("trait.myos", "dlabel")][ V1 == TRUE]
rtr[, empty:=paste0(trait.myos, "/", dlabel)]
sumc[, ec:=paste0(trait.myos, "/", dlabel)]
sumc <- sumc[!ec %in% rtr$empty]
sumc[, ec:=NULL]

# Attempt of new figure

# gp1 <- ggplot(sumc[ trait.myos == "DM (R)" ], aes(x =  trait.other, y = dlabel, colour=flag, fill = H4)) +
#               geom_tile( color = "black", lwd = 0.2, linetype = 1) +
#               scale_fill_gradient(limits = c(0,1), na.value = "white")+
#               geom_tile(aes( colour = flag), lwd=1, linetype = 1, data = sumc[trait.myos == "DM (R)" & H4 > 0.25]) + # Little trick to keep tiles in place
#               scale_colour_manual(values=c(High="green",Med="yellow", Low ="#FFFFFF00"), guide=guide_none()) + # And make rectangles under 0.5 transparent
#               theme_minimal() +
#               theme(panel.grid.major = element_blank(),
#                     axis.title = element_blank(),
#                     axis.text.x = element_blank(),
#                     legend.position = "none",
#                     plot.margin = unit(c(0, 0.2, 0, 0.3), "cm")
#                     )+
#               facet_grid(cols = vars(trait.myos), scales = "free", space = "free",switch = "y")+
#               labs(fill = "PP")

# gp2 <- ggplot(sumc[ trait.myos == "PM (M)" ], aes(x =  trait.other, y = dlabel, colour=flag, fill = H4)) +
#               geom_tile( color = "black", lwd = 0.2, linetype = 1) +
#               scale_fill_gradient(limits = c(0,1), na.value = "white")+
#               geom_tile(aes( colour = flag), lwd=1, linetype = 1, data = sumc[trait.myos == "PM (M)" & H4 > 0.25]) + # Little trick to keep tiles in place
#               scale_colour_manual(values=c(High="green",Med="yellow", Low ="#FFFFFF00"), guide=guide_none()) + # And make rectangles under 0.5 transparent
#               theme_minimal() +
#               theme(panel.grid.major = element_blank(),
#                     axis.title = element_blank(),
#                     axis.text.x = element_blank(),
#                     legend.position = "none",
#                     plot.margin = unit(c(0, 0.2, 0, 0.3), "cm")
#                     )+
#               facet_grid(cols = vars(trait.myos), scales = "free", space = "free",switch = "y")+
#               labs(fill = "PP")

# gp3 <- ggplot(sumc[ trait.myos == "PM (R)" ], aes(x =  trait.other, y = dlabel, colour=flag, fill = H4)) +
#               geom_tile( color = "black", lwd = 0.2, linetype = 1) +
#               scale_fill_gradient(limits = c(0,1), na.value = "white")+
#               geom_tile(aes( colour = flag), lwd=1, linetype = 1, data = sumc[trait.myos == "PM (R)" & H4 > 0.25]) + # Little trick to keep tiles in place
#               scale_colour_manual(values=c(High="green",Med="yellow", Low ="#FFFFFF00"), guide=guide_none()) + # And make rectangles under 0.5 transparent
#               theme_minimal() +
#               theme(panel.grid.major = element_blank(),
#                     axis.title = element_blank(),
#                     axis.text.x = element_blank(),
#                     legend.position = "none",
#                     plot.margin = unit(c(0, 0.2, 0, 0.3), "cm")
#                     )+
#               facet_grid(cols = vars(trait.myos), scales = "free", space = "free",switch = "y")+
#               labs(fill = "PP")
# gp4 <- ggplot(sumc[ trait.myos == "IIM (M)" ], aes(x =  trait.other, y = dlabel, colour=flag, fill = H4)) +
#               geom_tile( color = "black", lwd = 0.2, linetype = 1) +
#               scale_fill_gradient(limits = c(0,1), na.value = "white")+
#               geom_tile(aes( colour = flag), lwd=1, linetype = 1, data = sumc[trait.myos == "IIM (M)" & H4 > 0.25]) + # Little trick to keep tiles in place
#               scale_colour_manual(values=c(High="green",Med="yellow", Low ="#FFFFFF00"), guide=guide_none()) + # And make rectangles under 0.5 transparent
#               theme_minimal() +
#               theme(panel.grid.major = element_blank(),
#                     axis.title = element_blank(),
#                     axis.text.x = element_blank(),
#                     legend.position = "none",
#                     plot.margin = unit(c(0, 0.2, 0, 0.3), "cm")
#                     )+
#               facet_grid(cols = vars(trait.myos), scales = "free", space = "free",switch = "y")+
#               labs(fill = "PP")

# gp5 <- ggplot(sumc[ trait.myos == "IIM (R)" ], aes(x =  trait.other, y = dlabel, colour=flag, fill = H4)) +
#               geom_tile( color = "black",
#                         lwd = 0.2,
#                         linetype = 1) +
#               scale_fill_gradient(limits = c(0,1), na.value = "white")+
#               geom_tile(aes( colour = flag),
#                           lwd=1, linetype = 1, data = sumc[trait.myos == "IIM (R)" & H4 > .5]) +
#               scale_colour_manual(values=c(High="green",Med="yellow"), guide=guide_none()) +
#               #scale_x_discrete(position = "top")+
#               theme_minimal() +
#               theme(panel.grid.major = element_blank(),
#                     axis.title = element_blank(),
#                     axis.text.x = element_text(angle = 270, hjust=0, vjust = 0.5, size = 11),
#                     legend.position = "bottom",
#                     plot.margin = unit(c(0, 0.2, 0, 0.3), "cm")
#                     )+
#               guides(fill = guide_colorbar(title.vjust = 0.8))+
#               facet_grid(cols = vars(trait.myos), scales = "free", space = "free",switch = "y")+
#               labs(fill = "PP")

# gp <- plot_grid(gp1, gp2, gp3, gp4, gp5, nrow = 5, rel_heights = c(0.13, 0.09, 0.09, 0.09, 0.6), align = "v")

# ggsave("../figures/Figure3_driverSNP_H4_rnocoloc-v2.png", gp, height =5.8 , width = 4.5, bg="white")

# Figure 3 - Coloc panel

sumc[, dlabelm := factor(paste0(dlabel," - ", trait.myos))]
gp.alt <- ggplot(sumc, aes(x =  trait.other, y = trait.myos, colour=flag, fill = H4)) +
              geom_tile( color = "black", lwd = 0.2, linetype = 1) +
              scale_fill_gradient(limits = c(0,1), na.value = "white")+
              geom_tile(aes( colour = flag),
                          lwd=1, linetype = 1, data = sumc[H4 > 0.45]) +
              scale_colour_manual(values=c(High="green",Med="yellow"), guide=guide_none()) +
              #scale_x_discrete(position = "top")+
              theme_minimal() +
              theme(panel.grid.major = element_blank(),
                    axis.title = element_blank(),
                    axis.text.x = element_text(angle = 270, hjust=0, vjust = 0.5, size = 11),
                    legend.position = "bottom",
                    plot.margin = unit(c(0.1, 0.2, 0, 0.3), "cm"),
                    strip.text.y.left = element_text(hjust = 1, angle = 0)
                    )+
              facet_grid(forcats::fct_rev(dlabel)~., scales = "free_y", space = "free_y", switch = "y")+
              scale_y_discrete(position = "right")+
              guides(fill = guide_colorbar(title.vjust = 0.8))+

              labs(fill = "PP H4")

            #   forcats::fct_rev(dlabel)

            #   ,scales="free_y", switch = "y"
gp.alt
ggsave("../figures/Fig3_driverSNP_H4_rnocoloc-v2.png", gp.alt, height =5.2 , width = 4.9, bg="white")

# Journal REQUIRES TIFF format. Dear journal, be better!
tiff("../figures/Fig3_driverSNP_H4_rnocoloc-v2.tiff", units="in", width=4.9, height=5.2, res=300)
gp.alt
dev.off()


###############################

# We didn't use this figure in the end


# # Prepare H3 panel

# sumc.h3 <- c2[pid %in% keypids, .(dlabel, trait.myos,  trait.other, H3)]
# sumc.h3 <- sumc.h3 %>% tidyr::complete(dlabel, trait.myos, trait.other ) %>% as.data.table()
# sumc.h3[, dlabel:=factor(dlabel, levels = rev(unique(dlabel)))]
# sumc.h3[, trait.other:=factor(trait.other, levels = rev(unique(trait.other)))]
# sumc.h3[,flag:=ifelse(H3>.8, "High", ifelse(H3>.5, "Med", "Low"))]
# sumc.h3 <- sumc.h3[, .(dlabel, trait.myos, trait.other, H3, flag)]

# # Remove empty rows
# rtr.h3 <- sumc.h3[, all(is.na(H3)), by = c("trait.myos", "dlabel")][ V1 == TRUE]
# rtr.h3[, empty:=paste0(trait.myos, "/", dlabel)]
# sumc.h3[, ec:=paste0(trait.myos, "/", dlabel)]
# sumc.h3 <- sumc.h3[!ec %in% rtr.h3$empty]
# sumc.h3[, ec:=NULL]

# # Figure S13 - H3 panel

# gp.h3.alt <- ggplot(sumc.h3, aes(x =  trait.other, y = trait.myos, colour=flag, fill = H3)) +
#               geom_tile( color = "black", lwd = 0.2, linetype = 1) +
#               geom_tile(aes( colour = flag),
#                           lwd=1, linetype = 1, data = sumc.h3[H3 > 0.49]) +
#               scale_colour_manual(values=c(High="green",Med="yellow"), guide=guide_none()) +
#               scale_fill_gradient(limits = c(0,1), na.value = "white")+
#               theme_minimal() +
#               theme(panel.grid.major = element_blank(),
#                     axis.title = element_blank(),
#                     axis.text.x = element_text(angle = 270, hjust=0, vjust = 0.5, size = 11),
#                     legend.position = "bottom",
#                     plot.margin = unit(c(0.1, 0.2, 0, 0.3), "cm"),
#                     strip.text.y.left = element_text(hjust = 1, angle = 0)
#                     )+
#               guides(fill = guide_colorbar(title.vjust = 0.8))+
#               facet_grid(forcats::fct_rev(dlabel)~., scales = "free_y", space = "free_y", switch = "y")+
#               scale_y_discrete(position = "right")+
#               labs(fill = "PP H3")
# gp.h3.alt
# ggsave("../figures/FigureS13_driverSNP_H3.png", gp.h3.alt, height =5.3 , width = 5.3, bg="white")


####


# gp1.h3 <-  ggplot(sumc.h3[ grepl(pattern = "^DM|JDM \\(M\\)", trait.myos) ], aes(x =  trait.other, y = dlabel, colour=flag, fill = H3)) +
#               geom_tile( color = "black",
#                         lwd = 0.2,
#                         linetype = 1) +
#               scale_fill_gradient(limits = c(0,1), na.value = "white")+
#               geom_tile(aes( colour = flag),
#                           lwd=1, linetype = 1, data = sumc.h3[grepl(pattern = "^DM|JDM \\(M\\)", trait.myos) & H3 > 0.25]) + # Little trick to keep tiles in place
#               scale_colour_manual(values=c(High="green",Med="yellow", Low ="#FFFFFF00"), guide=guide_none()) + # And make rectangles under 0.5 transparent
              
#               #scale_x_discrete(position = "top")+
#               theme_minimal() +
#               theme(panel.grid.major = element_blank(),
#                     axis.title = element_blank(),
#                     axis.text.x = element_blank(),
#                     legend.position = "none",
#                     plot.margin = unit(c(0, 0.2, 0, 0.3), "cm")
#                     )+
#               facet_grid(cols = vars(trait.myos), scales = "free", space = "free",switch = "y")+
#               labs(fill = "PP")

# gp2.h3 <-  ggplot(sumc.h3[ grepl(pattern = "PM|JDM \\(R\\)", trait.myos) ], aes(x =  trait.other, y = dlabel, colour=flag, fill = H3)) +
#               geom_tile( color = "black",
#                         lwd = 0.2,
#                         linetype = 1) +
#               scale_fill_gradient(limits = c(0,1), na.value = "white")+
#               geom_tile(aes( colour = flag),
#                           lwd=1, linetype = 1, data = sumc.h3[grepl(pattern = "PM|JDM \\(R\\)", trait.myos) & H3 > 0.25]) + # Little trick to keep tiles in place
#               scale_colour_manual(values=c(High="green",Med="yellow", Low ="#FFFFFF00"), guide=guide_none()) + # And make rectangles under 0.5 transparent
              
#               #scale_x_discrete(position = "top")+
#               theme_minimal() +
#               theme(panel.grid.major = element_blank(),
#                     axis.title = element_blank(),
#                     axis.text.x = element_blank(),
#                     legend.position = "none",
#                     plot.margin = unit(c(0, 0.2, 0, 0.3), "cm")
#                     )+
#               facet_grid(cols = vars(trait.myos), scales = "free", space = "free",switch = "y")+
#               labs(fill = "PP")

# gp3.h3 <-  ggplot(sumc.h3[ grepl(pattern = "IIM|Jo1", trait.myos) ], aes(x =  trait.other, y = dlabel, colour=flag, fill = H3)) +
#               geom_tile( color = "black",
#                         lwd = 0.2,
#                         linetype = 1) +
#               scale_fill_gradient(limits = c(0,1), na.value = "white")+
#               geom_tile(aes( colour = flag),
#                           lwd=1, linetype = 1, data = sumc.h3[grepl(pattern = "IIM|Jo1", trait.myos) & H3 > .25]) +
#               scale_colour_manual(values=c(High="green",Med="yellow",Low ="#FFFFFF00"), guide=guide_none()) +
#               #scale_x_discrete(position = "top")+
#               theme_minimal() +
#               theme(panel.grid.major = element_blank(),
#                     axis.title = element_blank(),
#                     axis.text.x = element_text(angle = 270, hjust=0, vjust = 0.5, size = 11),
#                     legend.position = "bottom",
#                     plot.margin = unit(c(0, 0.2, 0, 0.3), "cm")
#                     )+
#               guides(fill = guide_colorbar(title.vjust = 0.8))+
#               facet_grid(cols = vars(trait.myos), scales = "free", space = "free",switch = "y")+
#               labs(fill = "PP")

# gp.h3 <- plot_grid(gp1.h3, gp2.h3, gp3.h3, nrow = 3, rel_heights = c(0.75,0.75,1.4))

# # Save
# ggsave("../figures/driverSNP_H3.png", gp.h3, height =7 , width = 8, bg="white")


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

fwrite(cts, "../tables/MT_coloc_results-v2.tsv", sep = "\t")


##########################################

sessionInfo()
# R version 4.3.3 (2024-02-29)
# Platform: x86_64-redhat-linux-gnu (64-bit)
# Running under: Rocky Linux 8.9 (Green Obsidian)

# Matrix products: default
# BLAS/LAPACK: /usr/lib64/libopenblaso-r0.3.15.so;  LAPACK version 3.9.0

# locale:
#  [1] LC_CTYPE=en_GB.UTF-8       LC_NUMERIC=C               LC_TIME=en_GB.UTF-8        LC_COLLATE=en_GB.UTF-8     LC_MONETARY=en_GB.UTF-8    LC_MESSAGES=en_GB.UTF-8   
#  [7] LC_PAPER=en_GB.UTF-8       LC_NAME=C                  LC_ADDRESS=C               LC_TELEPHONE=C             LC_MEASUREMENT=en_GB.UTF-8 LC_IDENTIFICATION=C       

# time zone: GB
# tzcode source: system (glibc)

# attached base packages:
# [1] stats     graphics  grDevices utils     datasets  methods   base     

# other attached packages:
# [1] cowplot_1.1.3      ggplot2_3.5.0      annotSnpStats_0.99 snpStats_1.52.0    Matrix_1.6-5       survival_3.5-8     coloc_5.2.3        magrittr_2.0.3    
# [9] data.table_1.15.4 

# loaded via a namespace (and not attached):
#  [1] viridis_0.6.5       utf8_1.2.4          generics_0.1.3      tidyr_1.3.1         lattice_0.22-6      grid_4.3.3          R.oo_1.26.0         plyr_1.8.9         
#  [9] jsonlite_1.8.8      R.utils_2.12.3      reshape_0.8.9       mixsqp_0.3-54       gridExtra_2.3       purrr_1.0.2         fansi_1.0.6         viridisLite_0.4.2  
# [17] scales_1.3.0        textshaping_0.3.7   cli_3.6.2           rlang_1.1.3         crayon_1.5.2        R.methodsS3_1.8.2   munsell_0.5.1       splines_4.3.3      
# [25] susieR_0.12.35      withr_3.0.0         tools_4.3.3         dplyr_1.1.4         colorspace_2.1-0    BiocGenerics_0.48.1 vctrs_0.6.5         R6_2.5.1           
# [33] matrixStats_1.3.0   lifecycle_1.0.4     zlibbioc_1.48.0     ragg_1.3.0          httpgd_2.0.1        irlba_2.3.5.1       pkgconfig_2.0.3     pillar_1.9.0       
# [41] gtable_0.3.4        glue_1.7.0          Rcpp_1.0.12         systemfonts_1.0.6   tibble_3.2.1        tidyselect_1.2.1    farver_2.1.1        unigd_0.1.1        
# [49] labeling_0.4.3      compiler_4.3.3     
