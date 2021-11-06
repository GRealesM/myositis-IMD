###############################################
#### Meta-analysis of key myositis SNPs #######
### Author: Guillermo Reales        ###########
### Date: 2021/11/06                ###########
###############################################



## NOTE: This is a preliminary version, including Miller and Rothwell. Since we wanted to search for replication in independent datasets, we
# did a meta-analysis excluding Miller in meta_key_SNPs_myo.R. This script is legacy now.

# Background: We want to meta-analyse three key SNPs for myositis across studies. 
# For it, we'll use rmeta package and the reduced datasets for IMD basis.

# Note: Run in HPC
getwd()
# [1] "/rds/project/rds-YRWZsjDGyaU/Projects/myositis-IMD/code"

# Load libraries
#install.packages("rmeta")
library(data.table)
library(rmeta)
library(magrittr)

# Path to files. Change if necessary
fpath = '~/rds/rds-cew54-basis/03-Bases/IMD_basis/reduced_datasets/'

# SNPs and coordinates for FinnGen
rsids = c("rs10196612","rs4728142", "rs2476601")
crds = c("2:198029415", "7:128933913", "1:113834946")


# Load datasets
# Myositis
myo1 <- fread(paste0(fpath, "MYO_Miller_26291516_1-ft.tsv")) %>% .[SNPID %in% rsids, c("SNPID","BETA","SE")] %>% .[, study:="MYO_Miller"]
myo2 <- fread(paste0(fpath, "20002_1322_Neale_UKBB_1-ft.tsv"))  %>% .[pid38 %in% crds, c("SNPID","BETA","SE")] %>% .[, study:="MYO_UKBB"] %>% .[, SNPID:=rsids[c(3,1,2)]] # Need little adjustment to match pids and rsids
myo3 <- fread(paste0(fpath, "M13_MYOSITIS_FinnGen_FinnGenR5_1-ft.tsv")) %>% .[SNPID %in% rsids, c("SNPID","BETA","SE")] %>% .[, study:="MYO_FinnGen"]
myo4 <- fread(paste0(fpath,"IIM_Rothwell_up_1-ft.tsv")) %>% .[pid38 %in% crds, c("SNPID", "BETA","SE")] %>% .[, study:="MYO_Rothwell"] %>% .[, SNPID:=rsids[c(3,1,2)]] 
myo <- rbindlist(list(myo1,myo2,myo3,myo4))
myo
#          SNPID       BETA         SE        study
#  1:  rs2476601 -0.1810966 0.06765641   MYO_Miller
#  2: rs10196612  0.1653285 0.04288243   MYO_Miller
#  3:  rs4728142  0.1523248 0.03977460   MYO_Miller
#  4:  rs2476601 -0.1969849 0.19350606     MYO_UKBB
#  5: rs10196612 -0.1652154 0.11776730     MYO_UKBB
#  6:  rs4728142  0.1444594 0.11830524     MYO_UKBB
#  7:  rs2476601 -0.2148000 0.11000000  MYO_FinnGen
#  8: rs10196612  0.1219000 0.07960000  MYO_FinnGen
#  9:  rs4728142  0.1148000 0.07930000  MYO_FinnGen
# 10:  rs2476601 -0.2600580 0.04884170 MYO_Rothwell
# 11: rs10196612  0.1307940 0.03141090 MYO_Rothwell
# 12:  rs4728142  0.1208990 0.03138790 MYO_Rothwell

# DM
dm1 <- fread(paste0(fpath, "DMY_Miller_26291516_1-ft.tsv")) %>% .[SNPID %in% rsids, c("SNPID","BETA","SE")] %>% .[, study:="DM_Miller"]
dm2 <- fread(paste0(fpath, "M13_DERMATOPOLY_FinnGen_FinnGenR5_1-ft.tsv")) %>% .[SNPID %in% rsids, c("SNPID","BETA","SE")] %>% .[, study:="DM_FinnGen"]
dm3 <- fread(paste0(fpath,"DMY_Rothwell_up_1-ft.tsv")) %>% .[pid38 %in% crds, c("SNPID", "BETA","SE")] %>% .[, study:="DM_Rothwell"] %>% .[, SNPID:=rsids[c(3,1,2)]] 
dm <- rbindlist(list(dm1,dm2, dm3))
dm
#         SNPID        BETA         SE       study
# 1:  rs2476601 -0.04672777 0.09762755   DM_Miller
# 2: rs10196612  0.17961666 0.06187896   DM_Miller
# 3:  rs4728142  0.17949114 0.05739437   DM_Miller
# 4:  rs2476601 -0.64550000 0.14330000  DM_FinnGen
# 5: rs10196612  0.03740000 0.10010000  DM_FinnGen
# 6:  rs4728142  0.36080000 0.10010000  DM_FinnGen
# 7:  rs2476601 -0.10920000 0.08925710 DM_Rothwell
# 8: rs10196612  0.14804700 0.05549940 DM_Rothwell
# 9:  rs4728142  0.11392600 0.05581850 DM_Rothwell

# PM
pm1 <- fread(paste0(fpath, "PM_Miller_26291516_1-ft.tsv")) %>% .[SNPID %in% rsids, c("SNPID","BETA","SE")] %>% .[, study:="PM_Miller"]
pm2 <- fread(paste0(fpath, "PM_Sakaue_doi1011012020102320213652_1-ft.tsv"))%>% .[SNPID %in% rsids, c("SNPID","BETA","SE")] %>% .[, study:="PM_Sakaue"]
pm3 <- fread(paste0(fpath,"PM_Rothwell_up_1-ft.tsv")) %>% .[pid38 %in% crds, c("SNPID", "BETA","SE")] %>% .[, study:="PM_Rothwell"] %>% .[, SNPID:=rsids[c(3,1,2)]] 
pm <- rbindlist(list(pm1,pm2, pm3))
pm
#         SNPID       BETA         SE       study
# 1:  rs2476601 -0.4767877 0.11332527   PM_Miller
# 2: rs10196612  0.1470304 0.07181533   PM_Miller
# 3:  rs4728142  0.1712436 0.06662292   PM_Miller
# 4: rs10196612 -0.1511835 0.20078347   PM_Sakaue
# 5:  rs4728142  0.4454922 0.29878545   PM_Sakaue
# 6:  rs2476601 -0.3855370 0.07789430 PM_Rothwell
# 7: rs10196612  0.1236040 0.05275360 PM_Rothwell
# 8:  rs4728142  0.1219290 0.05314870 PM_Rothwell


# Meta-analysis time

# rs10196612
rs10.myo  <- myo[SNPID == "rs10196612"] 
meta.myo.rs10 <- meta.summaries(rs10.myo$BETA, rs10.myo$SE, names= rs10.myo$study, method="fixed") 

rs10.dm  <- dm[SNPID == "rs10196612"] 
meta.dm.rs10 <- meta.summaries(rs10.dm$BETA, rs10.dm$SE, names= rs10.dm$study, method="fixed") 

rs10.pm  <- pm[SNPID == "rs10196612"] 
meta.pm.rs10 <- meta.summaries(rs10.pm$BETA, rs10.pm$SE, names= rs10.pm$study, method="fixed") 

res.rs10  <- data.table(SNPID ="rs10196612", Trait = c("MYO","DM","PM"), BETA = c(meta.myo.rs10$summary, meta.dm.rs10$summary, meta.pm.rs10$summary), SE = c(meta.myo.rs10$se.summary, meta.dm.rs10$se.summary, meta.pm.rs10$se.summary), chisq_het = c(meta.myo.rs10$het[1], meta.dm.rs10$het[1], meta.pm.rs10$het[1]) , df_het = c(meta.myo.rs10$het[2], meta.dm.rs10$het[2], meta.pm.rs10$het[2]), P_het = c(meta.myo.rs10$het[3], meta.dm.rs10$het[3], meta.pm.rs10$het[3]))
res.rs10
#         SNPID Trait      BETA         SE chisq_het df_het      P_het
# 1: rs10196612   MYO 0.1285745 0.02365405  6.969972      3 0.07286104
# 2: rs10196612    DM 0.1439664 0.03819074  1.470702      2 0.47933707
# 3: rs10196612    PM 0.1196701 0.04159337  1.970464      2 0.37335260


# rs4728142

rs47.myo  <- myo[SNPID == "rs4728142"] 
meta.myo.rs47 <- meta.summaries(rs47.myo$BETA, rs47.myo$SE, names= rs47.myo$study, method="fixed") 

rs47.dm  <- dm[SNPID == "rs4728142"] 
meta.dm.rs47 <- meta.summaries(rs47.dm$BETA, rs47.dm$SE, names= rs47.dm$study, method="fixed") 

rs47.pm  <- pm[SNPID == "rs4728142"] 
meta.pm.rs47 <- meta.summaries(rs47.pm$BETA, rs47.pm$SE, names= rs47.pm$study, method="fixed") 

res.rs47  <- data.table(SNPID ="rs4728142", Trait = c("MYO","DM","PM"), BETA = c(meta.myo.rs47$summary, meta.dm.rs47$summary, meta.pm.rs47$summary), SE = c(meta.myo.rs47$se.summary, meta.dm.rs47$se.summary, meta.pm.rs47$se.summary),chisq_het = c(meta.myo.rs47$het[1], meta.dm.rs47$het[1], meta.pm.rs47$het[1]),df_het = c(meta.myo.rs47$het[2], meta.dm.rs47$het[2], meta.pm.rs47$het[2]),  P_het = c(meta.myo.rs47$het[3], meta.dm.rs47$het[3], meta.pm.rs47$het[3]))

res.rs47
#        SNPID Trait      BETA         SE chisq_het df_het      P_het
# 1: rs4728142   MYO 0.1318587 0.02307805 0.4443023      3 0.93094602
# 2: rs4728142    DM 0.1754201 0.03715632 4.6484369      2 0.09785989
# 3: rs4728142    PM 0.1468818 0.04115169 1.3529615      2 0.50840306


# rs2476601

rs24.myo  <- myo[SNPID == "rs2476601"] 
meta.myo.rs24 <- meta.summaries(rs24.myo$BETA, rs24.myo$SE, names= rs24.myo$study, method="fixed") 

rs24.dm  <- dm[SNPID == "rs2476601"] # Note: two studies only 
meta.dm.rs24 <- meta.summaries(rs24.dm$BETA, rs24.dm$SE, names= rs24.dm$study, method="fixed") 

# No Sakaue for this SNP, so only Miller and Rothwell available for meta-analsysis
rs24.pm  <- pm[SNPID == "rs2476601"] 
meta.pm.rs24 <- meta.summaries(rs24.pm$BETA, rs24.pm$SE, names= rs24.pm$study, method="fixed") 

res.rs24  <- data.table(SNPID ="rs2476601", Trait = c("MYO","DM", "PM"), BETA = c(meta.myo.rs24$summary, meta.dm.rs24$summary, meta.pm.rs24$summary), SE = c(meta.myo.rs24$se.summary, meta.dm.rs24$se.summary, meta.pm.rs24$se.summary), chisq_het = c(meta.myo.rs24$het[1], meta.dm.rs24$het[1], meta.pm.rs24$het[1]), df_het = c(meta.myo.rs24$het[2], meta.dm.rs24$het[2], meta.pm.rs24$het[2])  ,  P_het = c(meta.myo.rs24$het[3], meta.dm.rs24$het[3], meta.pm.rs24$het[3]))

res.rs24
#        SNPID Trait       BETA         SE  chisq_het df_het      P_het
# 1: rs2476601   MYO -0.2297037 0.03658776  0.9493449      3 0.81350632
# 2: rs2476601    DM -0.1792802 0.05985371 13.0448627      2 0.00147009
# 3: rs2476601    PM -0.4148158 0.06419262  0.4403296      1 0.50696339


res <- rbindlist(list(res.rs10, res.rs47, res.rs24))
# Add significance for difference from zero
res[, Z:=BETA/SE][, P:=pnorm(-abs(Z))*2]

fwrite(res, "../tables/3SNPs_meta_results_Nov21.tsv", sep="\t")

## Forest plot of results.

# We're interested in following-up rs2476601, so we'll make a forest plot including all available data for this SNP

# Let's bring Miller

milmyo <- fread(paste0(fpath, "MYO_Miller_26291516_1-ft.tsv")) %>% .[SNPID == "rs2476601", c("SNPID","BETA","SE")] %>% .[, study:="Myositis / Miller"]
miljdm <- fread(paste0(fpath, "JDM_Miller_26291516_1-ft.tsv")) %>% .[SNPID == "rs2476601", c("SNPID","BETA","SE")] %>% .[, study:="JDM / Miller"]
mildm <- fread(paste0(fpath, "DMY_Miller_26291516_1-ft.tsv")) %>% .[SNPID == "rs2476601", c("SNPID","BETA","SE")] %>% .[, study:="DM / Miller"]
milpm <- fread(paste0(fpath, "PM_Miller_26291516_1-ft.tsv")) %>% .[SNPID == "rs2476601", c("SNPID","BETA","SE")] %>% .[, study:="PM / Miller"]

# Let's bring Rothwell
rotmyo <-  fread(paste0(fpath,"IIM_Rothwell_up_1-ft.tsv")) %>% .[pid38 %in% "1:113834946", c("SNPID","BETA","SE")] %>% .[,SNPID:="rs2476601"] %>% .[, study:="Myositis / Rothwell"] 
rotjdm <- fread(paste0(fpath,"JDM_Rothwell_up_1-ft.tsv")) %>% .[pid38 %in% "1:113834946", c("SNPID","BETA","SE")] %>% .[,SNPID:="rs2476601"] %>% .[, study:="JDM / Rothwell"]
rotdm <-fread(paste0(fpath,"DMY_Rothwell_up_1-ft.tsv")) %>% .[pid38 %in% "1:113834946", c("SNPID","BETA","SE")] %>% .[,SNPID:="rs2476601"] %>% .[, study:="DM / Rothwell"] 
rotpm <- fread(paste0(fpath,"PM_Rothwell_up_1-ft.tsv")) %>% .[pid38 %in% "1:113834946", c("SNPID","BETA","SE")] %>% .[,SNPID:="rs2476601"] %>% .[, study:="PM / Rothwell"] 

# Plus FinnGen PM 
finpm <- fread(paste0(fpath, "M13_POLYMYO_FinnGen_FinnGenR5_1-ft.tsv")) %>% .[SNPID == "rs2476601", c("SNPID","BETA","SE")] %>% .[, study:="PM / FinnGen"]

# FinnGen + UKBB Meta
ukbbmyo <- fread(paste0(fpath, "20002_1322_Neale_UKBB_1-ft.tsv"))  %>% .[pid38 %in% "1:113834946", c("SNPID","BETA","SE")]  %>% .[,SNPID:="rs2476601"] %>% .[, study:="Myositis / UKBB"]  
finmyo <- fread(paste0(fpath, "M13_MYOSITIS_FinnGen_FinnGenR5_1-ft.tsv")) %>% .[SNPID == "rs2476601", c("SNPID","BETA","SE")] %>% .[, study:="Myositis / FinnGen"]

rs24total <- rbindlist(list(milmyo, miljdm, mildm, milpm, rotmyo, rotjdm, rotdm, rotpm, finpm, ukbbmyo,finmyo))

rs24total[,Trait:=gsub(" \\/.*$","", study)][,study:=gsub(".+\\/ ","", study)]

rs24total[, Trait:=factor(Trait, levels=rev(c("Myositis", "PM", "DM", "JDM")))]

# Plot!

library(ggplot2)
library(cowplot)
fplot <- ggplot(rs24total, aes(y = study, x = BETA, xmin=BETA-SE, xmax=BETA+SE, colour = Trait))+
  geom_pointrange()+
  geom_vline(xintercept = 0, col="red", lty=2)+
#  coord_flip()+
  xlab("log(OR)")+
  ggtitle("Effect sizes for rs2476601 in Myositis and subtypes")+
  facet_grid(Trait~., scales = "free", space = "free", switch = "y")+
  theme_cowplot(8)+
 theme(legend.position = "none", strip.text.y.left = element_text(angle = 0), axis.title.y = element_blank(), plot.title = element_text(size = 8))
#  fplot

 ggsave("../figures/202111_Nov21/rs2476601_plot.png", fplot, bg = "white", units = "mm", width=125,height=60)

