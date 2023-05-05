###############################################
#### Meta-analysis of key myositis SNPs #######
### Author: Guillermo Reales        ###########
### Date: 2021/09/23                ###########
###############################################


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
#myo1 <- fread(paste0(fpath, "MYO_Miller_26291516_1-ft.tsv")) %>% .[SNPID %in% rsids, c("SNPID","BETA","SE")] %>% .[, study:="MYO_Miller"]
myo1 <- fread(paste0(fpath, "20002_1322_Neale_UKBB_1-ft.tsv"))  %>% .[pid38 %in% crds, c("SNPID","BETA","SE")] %>% .[, study:="MYO_UKBB"] %>% .[, SNPID:=rsids[c(3,1,2)]] # Need little adjustment to match pids and rsids
myo2 <- fread(paste0(fpath, "M13_MYOSITIS_FinnGen_FinnGenR5_1-ft.tsv")) %>% .[SNPID %in% rsids, c("SNPID","BETA","SE")] %>% .[, study:="MYO_FinnGen"]
myo <- rbindlist(list(myo1,myo2))

# DM
#dm1 <- fread(paste0(fpath, "DMY_Miller_26291516_1-ft.tsv")) %>% .[SNPID %in% rsids, c("SNPID","BETA","SE")] %>% .[, study:="DM_Miller"]
#dm2 <- fread(paste0(fpath, "M13_DERMATOPOLY_FinnGen_FinnGenR5_1-ft.tsv")) %>% .[SNPID %in% rsids, c("SNPID","BETA","SE")] %>% .[, study:="DM_FinnGen"]
#dm <- rbindlist(list(dm1,dm2))

# PM
#pm1 <- fread(paste0(fpath, "PM_Miller_26291516_1-ft.tsv")) %>% .[SNPID %in% rsids, c("SNPID","BETA","SE")] %>% .[, study:="PM_Miller"]
#pm2 <- fread(paste0(fpath, "PM_Sakaue_doi1011012020102320213652_1-ft.tsv"))%>% .[SNPID %in% rsids, c("SNPID","BETA","SE")] %>% .[, study:="PM_Sakaue"]
#pm <- rbindlist(list(pm1,pm2))

# Meta-analysis time

# rs10196612
rs10.myo  <- myo[SNPID == "rs10196612"] 
meta.myo.rs10 <- meta.summaries(rs10.myo$BETA, rs10.myo$SE, names= rs10.myo$study, method="fixed") 

#rs10.dm  <- dm[SNPID == "rs10196612"] 
#meta.dm.rs10 <- meta.summaries(rs10.dm$BETA, rs10.dm$SE, names= rs10.dm$study, method="fixed") 

#rs10.pm  <- pm[SNPID == "rs10196612"] 
#meta.pm.rs10 <- meta.summaries(rs10.pm$BETA, rs10.pm$SE, names= rs10.pm$study, method="fixed") 

#res.rs10  <- data.table(SNPID ="rs10196612", Trait = c("MYO","DM","PM"), BETA = c(meta.myo.rs10$summary, meta.dm.rs10$summary, meta.pm.rs10$summary), SE = c(meta.myo.rs10$se.summary, meta.dm.rs10$se.summary, meta.pm.rs10$se.summary), chisq_het = c(meta.myo.rs10$het[1], meta.dm.rs10$het[1], meta.pm.rs10$het[1]) , df_het = c(meta.myo.rs10$het[2], meta.dm.rs10$het[2], meta.pm.rs10$het[2]), P_het = c(meta.myo.rs10$het[3], meta.dm.rs10$het[3], meta.pm.rs10$het[3]))

# Myo only, excluding Miller
res.rs10  <- data.table(SNPID ="rs10196612", Trait = "MYO", BETA = meta.myo.rs10$summary, SE = meta.myo.rs10$se.summary, chisq_het = meta.myo.rs10$het[1] , df_het = meta.myo.rs10$het[2], P_het = meta.myo.rs10$het[3])

# rs4728142

rs47.myo  <- myo[SNPID == "rs4728142"] 
meta.myo.rs47 <- meta.summaries(rs47.myo$BETA, rs47.myo$SE, names= rs47.myo$study, method="fixed") 

#rs47.dm  <- dm[SNPID == "rs4728142"] 
#meta.dm.rs47 <- meta.summaries(rs47.dm$BETA, rs47.dm$SE, names= rs47.dm$study, method="fixed") 

#rs47.pm  <- pm[SNPID == "rs4728142"] 
#meta.pm.rs47 <- meta.summaries(rs47.pm$BETA, rs47.pm$SE, names= rs47.pm$study, method="fixed") 

#res.rs47  <- data.table(SNPID ="rs4728142", Trait = c("MYO","DM","PM"), BETA = c(meta.myo.rs47$summary, meta.dm.rs47$summary, meta.pm.rs47$summary), SE = c(meta.myo.rs47$se.summary, meta.dm.rs47$se.summary, meta.pm.rs47$se.summary),chisq_het = c(meta.myo.rs47$het[1], meta.dm.rs47$het[1], meta.pm.rs47$het[1]),df_het = c(meta.myo.rs47$het[2], meta.dm.rs47$het[2], meta.pm.rs47$het[2]),  P_het = c(meta.myo.rs47$het[3], meta.dm.rs47$het[3], meta.pm.rs47$het[3]))

res.rs47  <- data.table(SNPID ="rs4728142", Trait = c("MYO"), BETA = meta.myo.rs47$summary, SE = meta.myo.rs47$se.summary, chisq_het = meta.myo.rs47$het[1], df_het = meta.myo.rs47$het[2],  P_het = meta.myo.rs47$het[3])

# rs2476601

rs24.myo  <- myo[SNPID == "rs2476601"] 
meta.myo.rs24 <- meta.summaries(rs24.myo$BETA, rs24.myo$SE, names= rs24.myo$study, method="fixed") 

#rs24.dm  <- dm[SNPID == "rs2476601"] # Note: two studies only 
#meta.dm.rs24 <- meta.summaries(rs24.dm$BETA, rs24.dm$SE, names= rs24.dm$study, method="fixed") 

# No Sakaue for this SNP, so only Miller available and meta-analsysis not possible
#rs24.pm  <- pm[SNPID == "rs2476601"] 
#meta.pm.rs24 <- meta.summaries(rs24.pm$BETA, rs24.pm$SE, names= rs24.pm$study, method="fixed") 

#res.rs24  <- data.table(SNPID ="rs2476601", Trait = c("MYO","DM"), BETA = c(meta.myo.rs24$summary, meta.dm.rs24$summary), SE = c(meta.myo.rs24$se.summary, meta.dm.rs24$se.summary), chisq_het = c(meta.myo.rs24$het[1], meta.dm.rs24$het[1]), df_het = c(meta.myo.rs24$het[2], meta.dm.rs24$het[2])  ,  P_het = c(meta.myo.rs24$het[3], meta.dm.rs24$het[3]))

res.rs24  <- data.table(SNPID ="rs2476601", Trait = c("MYO"), BETA = meta.myo.rs24$summary, SE = meta.myo.rs24$se.summary, chisq_het = meta.myo.rs24$het[1], df_het = meta.myo.rs24$het[2],  P_het = meta.myo.rs24$het[3])

res <- rbindlist(list(res.rs10, res.rs47, res.rs24))

# Add significance for difference from zero
res[, Z:=BETA/SE][, P:=pnorm(-abs(Z))*2]

fwrite(res, "../tables/3SNPs_meta_results.tsv", sep="\t")


## Forest plot of results.

# We're interested in following-up rs2476601, so we'll make a forest plot including all available data for this SNP

# Let's bring Miller

milmyo <- fread(paste0(fpath, "MYO_Miller_26291516_1-ft.tsv")) %>% .[SNPID == "rs2476601", c("SNPID","BETA","SE")] %>% .[, study:="Myositis / Miller"]
miljdm <- fread(paste0(fpath, "JDM_Miller_26291516_1-ft.tsv")) %>% .[SNPID == "rs2476601", c("SNPID","BETA","SE")] %>% .[, study:="JDM / Miller"]
mildm <- fread(paste0(fpath, "DMY_Miller_26291516_1-ft.tsv")) %>% .[SNPID == "rs2476601", c("SNPID","BETA","SE")] %>% .[, study:="DM / Miller"]
milpm <- fread(paste0(fpath, "PM_Miller_26291516_1-ft.tsv")) %>% .[SNPID == "rs2476601", c("SNPID","BETA","SE")] %>% .[, study:="PM / Miller"]

# Plus FinnGen PM and Meta elements
finpm <- fread(paste0(fpath, "M13_POLYMYO_FinnGen_FinnGenR5_1-ft.tsv")) %>% .[SNPID == "rs2476601", c("SNPID","BETA","SE")] %>% .[, study:="PM / FinnGen"]

rs24.myo[, study:=c("Myositis / UKBB", "Myositis / FinnGen")]

# Refurbish results
res.r24 <- res[SNPID=="rs2476601", c("SNPID","BETA","SE")][, study:="Meta Myositis / UKBB + FinnGen"]

rs24total <- rbindlist(list(rs24.myo, milmyo, miljdm, mildm, milpm, finpm))

rs24total[,Trait:=gsub(" \\/.*$","", study)][,study:=gsub(".+\\/ ","", study)]

rs24total[, Trait:=factor(Trait, levels=rev(c("Meta Myositis", "Myositis", "PM", "DM", "JDM")))]

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

 ggsave("../figures/rs2476601_plot.png", fplot, bg = "white", units = "mm", width=125,height=60)

