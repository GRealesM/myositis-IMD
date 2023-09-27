#########################################
##                                     ##
##  EXTRACTING P-VALUES FOR KEY SNPS   ##
##                                     ##
#########################################

# Author: Guillermo Reales 
# Date last updated: 2023/09/14

# Background:  In some instances, the top candidate SNP differs among tests for the same driver SNPs. This is likely
# an artefact, as it's likely capturing a single signal. Thus, among the results, we chose the most likely top candidate SNP,
# based on H4 and pairwise FDR values. The p-values of the selected top candidate SNP might not come with coloc results when these mismatches occur,
# so we'll retrieve them here.
# We'll also take the chance to made some ladder plots of SNPs of interestAfter running coloc, we can take a closer look at the regions where we found positive results.

# This script will
# * Load coloc results and the dense-SNP datasets.
# * Extract the P-values of relevant candidate SNPs in the required myositis datasets.
# * Create ladder plots, if necessary


##########################################


setwd("/home/gr440/rds/rds-cew54-basis/Projects/myositis-IMD/code")

# Load pkgs

library(data.table)
setDTthreads(15)
library(magrittr)

# Set variables and load summary statistics


files = list(   dmy.m = "DMY_Miller_26291516_1-hg38.tsv.gz",
                dmy.r = "DMY_Rothwell_up_1-hg38.tsv.gz",
                jdm.m = "JDM_Miller_26291516_1-hg38.tsv.gz",
                jdm.r = "JDM_Rothwell_up_1-hg38.tsv.gz",
                jo1m.r = "JO1M_Rothwell_up_1-hg38.tsv.gz",
                myo.m = "MYO_Miller_26291516_1-hg38.tsv.gz",
                myo.r = "IIM_Rothwell_up_1-hg38.tsv.gz",
                pm.m =   "PM_Miller_26291516_1-hg38.tsv.gz",
                pm.r = "PM_Rothwell_up_1-hg38.tsv.gz") %>%
    lapply(.,  function(f) (file.path("~/rds/rds-cew54-basis/02-Processed",f)))


man=fread("~/rds/rds-cew54-basis/03-Bases/IMD_basis/SNP.manifest.38.tsv")
snps=fread("~/rds/rds-cew54-basis/03-Bases/IMD_basis/Manifest_build_translator.tsv")
snps[, pid38:=paste(CHR38, BP38, sep=":")][, pid19:=paste(CHR19, BP19, sep=":")]
snps <- merge(snps, man, by="pid38")

data=lapply(files, fread)



# Import tsv for main table
mt <- fread("../tables/MT_coloc_results.tsv")

mt[, .SD[which.max(H4)] , by=driver.rsid]
#    driver.rsid          pid driver.nearestgene bestsnp.rsid bestsnp.novel
#  1:   rs1160542  2:100215693               AFF3   rs11692867           Yes
#  2:  rs13277113   8:11491677                BLK    rs2736340           Yes
#  3:   rs1855025  6:167124106               CCR6    rs1571878           Yes
#  4:    rs669607   3:28029953               CMC1     rs819991           Yes
#  5:    rs394378  5:157185077             GARIN3     rs163315           Yes
#  6:   rs7073236   10:6064589              IL2RA     rs706778           Yes
#  7:  rs10488631  7:128954129               IRF5   rs13236009           Yes
#  8:   rs2286896  2:190670850               NAB1    rs3821236            No
#  9:   rs2476601  1:113834946             PTPN22    rs2476601           Yes
# 10:  rs11066320 12:112468611               RPL6     rs597808           Yes
# 11:    rs991817 12:110972733              SH2B3    rs3184504           Yes
#      pbest.myos trait.myos trait.other pairwise_fdr        H4
#  1: 3.99102e-04    IIM (R)         JIA 2.585052e-01 0.6252904
#  2: 1.52900e-04    IIM (M)         SSc 6.325387e-02 0.8963674
#  3: 9.10520e-04    IIM (R)          RA 1.779860e-01 0.5160895
#  4: 1.05722e-03    IIM (R)         SSc 2.992867e-02 0.5254404
#  5: 3.78260e-05     DM (R)    HyperThy 3.031138e-02 0.8198758
#  6: 5.83929e-04    IIM (R)          RA 1.800270e-02 0.5754436
#  7: 6.41047e-05    IIM (R)         SSc 3.022458e-03 0.9105339
#  8: 1.01225e-05    IIM (R)         SSc 5.656435e-02 0.9563165
#  9: 1.62873e-07    IIM (R)          MG 2.198465e-05 0.9984643
# 10: 2.97654e-04    IIM (R)     HypoThy 1.766611e-01 0.7404610
# 11: 2.56157e-04    IIM (R)     HypoThy 6.528374e-01 0.7754211

# List top SNPs
#top <- c("rs11692867","rs2736340","rs819991","rs13236009","rs3821236","rs2476601","rs597808","rs597808")
# top <- c("rs2736340", "rs13236009","rs2476601") 
top <- "rs13236009"

# Import mapped genes to relate rsid to their pid
mg <- fread("../data/mapped.genes.tsv")
mg <- mg[rsID %in% top, .(rsID, pid)]

studies <- names(data)

# Extract required SNPs
d2 <- lapply(1:9, function(i){
        x <- data[[i]]
        x[, pid := paste(CHR38, BP38, sep=":")]
        x  <- x[ pid %in% mg$pid , .(pid, SNPID, REF, ALT, BETA, SE, P)]
        x[, study := studies[i]]
        x

}) %>% rbindlist

d2 <- merge(d2, mg, by="pid")

d2[ , .(rsID, pid, REF,ALT, P, study)]
#          rsID         pid REF ALT           P  study
# 1: rs13236009 7:129023119   T   G 2.01914e-03  dmy.m
# 2: rs13236009 7:129023119   T   G 1.00165e-02  dmy.r
# 3: rs13236009 7:129023119   T   G 7.47704e-01  jdm.m
# 4: rs13236009 7:129023119   T   G 2.05076e-01  jdm.r
# 5: rs13236009 7:129023119   T   G 1.16643e-03 jo1m.r
# 6: rs13236009 7:129023119   T   G 1.34854e-03  myo.m
# 7: rs13236009 7:129023119   T   G 6.41047e-05  myo.r
# 8: rs13236009 7:129023119   T   G 2.29058e-02   pm.m
# 9: rs13236009 7:129023119   T   G 5.17526e-02   pm.r


#####################################

## OUTDATED BUT MIGHT STILL BE USEFUL
## Extract info on rs2476601 to create a ladder plot
d3  <- d2[bestsnp.rsid == "rs2476601"]

ml <- data.table(mlabel = c("DM (M)", "DM (R)", "JDM (M)","JDM (R)", "Anti-Jo1+ (R)", "IIM (M)", "IIM (R)", "PM (M)", "PM (R)"),
                 trait.myos = c("dmy.m", "dmy.r", "jdm.m", "jdm.r", "jo1m.r", "myo.m", "myo.r", "pm.m", "pm.r"))
d3  <- merge(d3, ml, by.x = "study", by.y = "trait.myos")
myoc <- c(`PM (R)` = "#CF000F", `PM (M)` = "#CF000F", `PM (FG)` = "#CF000F", `DM (R)` = "#2E8856", `DM (M)` = "#2E8856", `IIM (R)` = "#1460AA", `IIM (M)` = "#1460AA", `JDM (M)` = "#B8860B", `JDM (R)` = "#B8860B", `IBM (R)` = "#E65722", `Anti-Jo1+ (R)` ="#1C2833", `DPM (FG)` = "#053061")

# Remove IIM as it's composite
d3 <- d3[!grepl("IIM", mlabel)]

library(ggplot2)
library(cowplot)
dpm <- ggplot(d3, aes(x = BETA, y = mlabel, xmin=BETA-SE, xmax=BETA+SE, colour = mlabel))+
  geom_pointrange()+
  geom_vline(xintercept = 0, col="red", lty=2)+
  scale_colour_manual(values = myoc)+
  xlab("Beta")+
  theme_cowplot(15)+
  theme(legend.position = "none", axis.title.y = element_blank())
dpm
ggsave("../figures/ladder_plot_rs2476601.png", dpm, height = 3, width = 5, bg="white")


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
# [1] magrittr_2.0.3    data.table_1.14.8

# loaded via a namespace (and not attached):
# [1] compiler_4.3.1    cli_3.6.1         R.methodsS3_1.8.2 jsonlite_1.8.7   
# [5] R.utils_2.12.2    rlang_1.1.1       R.oo_1.25.0 