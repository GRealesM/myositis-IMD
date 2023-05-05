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
myo1 <- fread(paste0(fpath, "MYO_Miller_26291516_1-ft.tsv")) %>% .[SNPID %in% rsids, c("SNPID","BETA","SE","P")] %>% .[, study:="Miller"] %>% .[, trait:="Myositis"]
myo2 <- fread(paste0(fpath, "20002_1322_Neale_UKBB_1-ft.tsv"))  %>% .[pid38 %in% crds, c("SNPID","BETA","SE","P")] %>% .[, study:="UKBB"] %>% .[,trait:="Myositis"] %>% .[, SNPID:=rsids[c(3,1,2)]] # Need little adjustment to match pids and rsids
myo3 <- fread(paste0(fpath, "M13_MYOSITIS_FinnGen_FinnGenR5_1-ft.tsv")) %>% .[SNPID %in% rsids, c("SNPID","BETA","SE","P")] %>% .[, study:="FinnGen"] %>% .[, trait:="Myositis"]
myo4 <- fread(paste0(fpath,"IIM_Rothwell_up_1-ft.tsv")) %>% .[pid38 %in% crds, c("SNPID", "BETA","SE","P")] %>% .[, study:="Rothwell"] %>% .[, trait:="Myositis"] %>% .[, SNPID:=rsids[c(3,1,2)]] 
myo <- rbindlist(list(myo1,myo2,myo3,myo4))
myo
#          SNPID       BETA         SE           P    study    trait
#  1:  rs2476601 -0.1810966 0.06765641 7.43500e-03   Miller Myositis
#  2: rs10196612  0.1653285 0.04288243 1.15543e-04   Miller Myositis
#  3:  rs4728142  0.1523248 0.03977460 1.28300e-04   Miller Myositis
#  4:  rs2476601 -0.1969849 0.19350606 3.08689e-01     UKBB Myositis
#  5: rs10196612 -0.1652154 0.11776730 1.60648e-01     UKBB Myositis
#  6:  rs4728142  0.1444594 0.11830524 2.22059e-01     UKBB Myositis
#  7:  rs2476601 -0.2148000 0.11000000 5.08900e-02  FinnGen Myositis
#  8: rs10196612  0.1219000 0.07960000 1.25700e-01  FinnGen Myositis
#  9:  rs4728142  0.1148000 0.07930000 1.47800e-01  FinnGen Myositis
# 10:  rs2476601 -0.2600580 0.04884170 1.62873e-07 Rothwell Myositis
# 11: rs10196612  0.1307940 0.03141090 3.04470e-05 Rothwell Myositis
# 12:  rs4728142  0.1208990 0.03138790 1.17873e-04 Rothwell Myositis

# DM
dm1 <- fread(paste0(fpath, "DMY_Miller_26291516_1-ft.tsv")) %>% .[SNPID %in% rsids, c("SNPID","BETA","SE","P")] %>% .[, study:="Miller"] %>% .[, trait:="DM"]
dm2 <- fread(paste0(fpath, "M13_DERMATOPOLY_FinnGen_FinnGenR5_1-ft.tsv")) %>% .[SNPID %in% rsids, c("SNPID","BETA","SE","P")] %>% .[, study:="FinnGen"] %>% .[, trait:="DM"]
dm3 <- fread(paste0(fpath,"DMY_Rothwell_up_1-ft.tsv")) %>% .[pid38 %in% crds, c("SNPID", "BETA","SE","P")] %>% .[, study:="Rothwell"] %>% .[, SNPID:=rsids[c(3,1,2)]] %>% .[, trait:="DM"]
dm <- rbindlist(list(dm1,dm2, dm3))
dm
#         SNPID        BETA         SE           P    study trait
# 1:  rs2476601 -0.04672777 0.09762755 0.632200000   Miller    DM
# 2: rs10196612  0.17961666 0.06187896 0.003699460   Miller    DM
# 3:  rs4728142  0.17949114 0.05739437 0.001764000   Miller    DM
# 4:  rs2476601 -0.64550000 0.14330000 0.000006666  FinnGen    DM
# 5: rs10196612  0.03740000 0.10010000 0.708300000  FinnGen    DM
# 6:  rs4728142  0.36080000 0.10010000 0.000315100  FinnGen    DM
# 7:  rs2476601 -0.10920000 0.08925710 0.224673000 Rothwell    DM
# 8: rs10196612  0.14804700 0.05549940 0.007530510 Rothwell    DM
# 9:  rs4728142  0.11392600 0.05581850 0.041302500 Rothwell    DM

# PM
pm1 <- fread(paste0(fpath, "PM_Miller_26291516_1-ft.tsv")) %>% .[SNPID %in% rsids, c("SNPID","BETA","SE","P")] %>% .[, study:="Miller"]  %>% .[, trait:="PM"]
pm2 <- fread(paste0(fpath, "PM_Sakaue_doi1011012020102320213652_1-ft.tsv"))%>% .[SNPID %in% rsids, c("SNPID","BETA","SE","P")] %>% .[, study:="Sakaue"]  %>% .[, trait:="PM"]
pm3 <- fread(paste0(fpath,"PM_Rothwell_up_1-ft.tsv")) %>% .[pid38 %in% crds, c("SNPID", "BETA","SE","P")] %>% .[, study:="Rothwell"] %>% .[, SNPID:=rsids[c(3,1,2)]]  %>% .[, trait:="PM"]
pm4 <- fread(paste0(fpath, "M13_POLYMYO_FinnGen_FinnGenR5_1-ft.tsv")) %>% .[SNPID %in% rsids , c("SNPID","BETA","SE","P")] %>% .[, study:="FinnGen"]  %>% .[, trait:="PM"]
pm <- rbindlist(list(pm1,pm2, pm3, pm4))
pm
#          SNPID       BETA         SE            P    study trait
#  1:  rs2476601 -0.4767877 0.11332527 2.585000e-05   Miller    PM
#  2: rs10196612  0.1470304 0.07181533 4.062420e-02   Miller    PM
#  3:  rs4728142  0.1712436 0.06662292 1.016000e-02   Miller    PM
#  4: rs10196612 -0.1511835 0.20078347 4.514693e-01   Sakaue    PM
#  5:  rs4728142  0.4454922 0.29878545 1.359588e-01   Sakaue    PM
#  6:  rs2476601 -0.3855370 0.07789430 1.247330e-06 Rothwell    PM
#  7: rs10196612  0.1236040 0.05275360 1.896910e-02 Rothwell    PM
#  8:  rs4728142  0.1219290 0.05314870 2.182920e-02 Rothwell    PM
#  9:  rs2476601 -0.7855000 0.18910000 3.273000e-05  FinnGen    PM
# 10: rs10196612 -0.0069000 0.13180000 9.580000e-01  FinnGen    PM
# 11:  rs4728142  0.3619000 0.13160000 5.937000e-03  FinnGen    PM


# Meta-analysis time

# Here we'll meta-analyse some specific SNPs for which we have external dataset information from UKBB and FinnGen

# rs10196612 in Myositis
rs10.myo  <- myo[SNPID == "rs10196612" & study %in% c("UKBB", "FinnGen")] # Note: different sign 
meta.myo.rs10 <- meta.summaries(rs10.myo$BETA, rs10.myo$SE, names= rs10.myo$study, method="fixed") 

res.rs10  <- data.table(SNPID ="rs10196612", Trait = "Myositis", BETA = meta.myo.rs10$summary, SE = meta.myo.rs10$se.summary, chisq_het = meta.myo.rs10$het[1] , df_het = meta.myo.rs10$het[2], P_het = meta.myo.rs10$het[3])
res.rs10[, Z:=BETA/SE][, P:=pnorm(-abs(Z))*2]
res.rs10
#         SNPID    Trait       BETA         SE chisq_het df_het      P_het
# 1: rs10196612 Myositis 0.03186376 0.06594854  4.079883      1 0.04339668
#           Z         P
# 1: 0.483161 0.6289815



# rs4728142 in myositis
rs47.myo  <- myo[SNPID == "rs4728142" & study %in% c("UKBB", "FinnGen")] 
meta.myo.rs47 <- meta.summaries(rs47.myo$BETA, rs47.myo$SE, names= rs47.myo$study, method="fixed") 

res.rs47  <- data.table(SNPID ="rs4728142", Trait = "Myositis", BETA = meta.myo.rs47$summary, SE = meta.myo.rs47$se.summary, chisq_het = meta.myo.rs47$het[1] , df_het = meta.myo.rs47$het[2], P_het = meta.myo.rs47$het[3])
res.rs47[, Z:=BETA/SE][, P:=pnorm(-abs(Z))*2]
res.rs47
#        SNPID    Trait      BETA         SE  chisq_het df_het     P_het
# 1: rs4728142 Myositis 0.1239948 0.06587092 0.04336673      1 0.8350363
#          Z          P
# 1: 1.88239 0.05978305


# rs2476601 in PM
# We don't have info on PM for UKBB so we'll simply look at the independent FinnGen result 
pm[SNPID == "rs2476601" & study == "FinnGen"]
#        SNPID    BETA     SE         P   study trait
# 1: rs2476601 -0.7855 0.1891 3.273e-05 FinnGen    PM


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
rotjo1 <-  fread(paste0(fpath,"JO1M_Rothwell_up_1-ft.tsv")) %>% .[pid38 %in% "1:113834946", c("SNPID","BETA","SE")] %>% .[,SNPID:="rs2476601"] %>% .[, study:="Jo1+ / Rothwell"] 
rotibm  <- fread(paste0(fpath,"IBM_Rothwell_up_1-ft.tsv")) %>% .[pid38 %in% "1:113834946", c("SNPID","BETA","SE")] %>% .[,SNPID:="rs2476601"] %>% .[, study:="IBM / Rothwell"] 
 

# Plus FinnGen PM 
finpm <- fread(paste0(fpath, "M13_POLYMYO_FinnGen_FinnGenR5_1-ft.tsv")) %>% .[SNPID == "rs2476601", c("SNPID","BETA","SE")] %>% .[, study:="PM / FinnGen"]

# FinnGen + UKBB 
ukbbmyo <- fread(paste0(fpath, "20002_1322_Neale_UKBB_1-ft.tsv"))  %>% .[pid38 %in% "1:113834946", c("SNPID","BETA","SE")]  %>% .[,SNPID:="rs2476601"] %>% .[, study:="Myositis / UKBB"]  
finmyo <- fread(paste0(fpath, "M13_MYOSITIS_FinnGen_FinnGenR5_1-ft.tsv")) %>% .[SNPID == "rs2476601", c("SNPID","BETA","SE")] %>% .[, study:="Myositis / FinnGen"]

rs24total <- rbindlist(list(milmyo, miljdm, mildm, milpm, rotmyo, rotjdm, rotdm, rotpm, rotibm, rotjo1, finpm, ukbbmyo,finmyo))

rs24total[,Trait:=gsub(" \\/.*$","", study)][,study:=gsub(".+\\/ ","", study)]

rs24total[, Trait:=factor(Trait, levels=rev(c("Myositis", "Jo1+", "PM", "DM", "JDM", "IBM")))]

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

 ggsave("../figures/202111_Nov21/rs2476601_plot.png", fplot, bg = "white", units = "mm", width=125,height=70)

