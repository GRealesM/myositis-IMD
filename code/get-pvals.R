### Getting p-values of top candidate SNPs in their respective myositis datasets

# Date: 2023/08/01
# Author: Guillermo Reales

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


# Import mapped genes 
mg <- fread("../data/mapped.genes.tsv") %>% unique
mg <- mg[,.(SNPID, pid, nearestGene)]

# Import coloc table and add some info to the coloc table, and extract rsids to map
coloc <- fread("../tables/coloc_results_dfilt.tsv")
coloc <- merge(coloc, mg, by="pid") # First on driver SNPs
setnames(coloc, c("SNPID", "nearestGene"), c("driver.rsid", "driver.nearestgene"))
coloc <- merge(coloc, mg, by.x="bestsnp", by.y = "pid", all.x = TRUE) # Then on candidate snps. Bear in mind that we only mapped candidate SNPs with H4 > 0.5
setnames(coloc, c("SNPID", "nearestGene"), c("bestsnp.rsid", "bestsnp.nearestgene"))


# List top SNPs
top <- c("rs11692867","rs2736340","rs819991","rs13236009","rs3821236","rs2476601","rs597808","rs597808")
topids <- coloc[ bestsnp.rsid %in% top, .(bestsnp, bestsnp.rsid)]  %>% unique

# Extract required SNPs
d2 <- lapply(data, function(x){
        x[, pid:=paste(CHR38, BP38, sep=":")]
        x <- x[ pid %in% topids$bestsnp]
        x
})
names(d2) = names(data)
d2 <- lapply(1:9, function(i){
        d2[[i]][, study:=names(d2)[i]]
        d2[[i]][, .(pid, REF,ALT, P, study)]

}) %>% rbindlist

d2 <- merge(d2, topids, by.x="pid", by.y="bestsnp")
d2
#              pid REF ALT           P  study bestsnp.rsid
#  1: 12:111535554   A   G 5.66558e-02  dmy.m     rs597808
#  2: 12:111535554   A   G 2.11556e-01  dmy.r     rs597808
#  3: 12:111535554   A   G 8.73370e-02  jdm.m     rs597808
#  4: 12:111535554   A   G 1.14514e-02  jdm.r     rs597808
#  5: 12:111535554   A   G 1.84146e-01 jo1m.r     rs597808
#  6: 12:111535554   A   G 1.22923e-02  myo.m     rs597808
#  7: 12:111535554   A   G 2.97654e-04  myo.r     rs597808
#  8: 12:111535554   A   G 1.10678e-01   pm.m     rs597808
#  9: 12:111535554   A   G 1.84235e-01   pm.r     rs597808
# 10:  1:113834946   A   G 6.32200e-01  dmy.m    rs2476601
# 11:  1:113834946   A   G 2.24673e-01  dmy.r    rs2476601
# 12:  1:113834946   A   G 7.74000e-01  jdm.m    rs2476601
# 13:  1:113834946   A   G 3.25488e-01  jdm.r    rs2476601
# 14:  1:113834946   A   G 2.72829e-03 jo1m.r    rs2476601
# 15:  1:113834946   A   G 7.43500e-03  myo.m    rs2476601
# 16:  1:113834946   A   G 1.62873e-07  myo.r    rs2476601
# 17:  1:113834946   A   G 2.58500e-05   pm.m    rs2476601
# 18:  1:113834946   A   G 1.24733e-06   pm.r    rs2476601
# 19:  2:100143015   G   A 1.71914e-02  dmy.m   rs11692867
# 20:  2:100143015   G   A 2.04719e-01  dmy.r   rs11692867
# 21:  2:100143015   G   A 6.40861e-01  jdm.m   rs11692867
# 22:  2:100143015   G   A 2.73739e-01  jdm.r   rs11692867
# 23:  2:100143015   G   A 3.99187e-01 jo1m.r   rs11692867
# 24:  2:100143015   G   A 2.90187e-03  myo.m   rs11692867
# 25:  2:100143015   G   A 3.99102e-04  myo.r   rs11692867
# 26:  2:100143015   G   A 4.13220e-03   pm.m   rs11692867
# 27:  2:100143015   G   A 3.65442e-02   pm.r   rs11692867
# 28:  2:191038032   G   A 3.32900e-02  dmy.m    rs3821236
# 29:  2:191038032   G   A 7.36387e-03  dmy.r    rs3821236
# 30:  2:191038032   G   A 7.43300e-02  jdm.m    rs3821236
# 31:  2:191038032   G   A 9.55881e-02  jdm.r    rs3821236
# 32:  2:191038032   G   A 4.35381e-02 jo1m.r    rs3821236
# 33:  2:191038032   G   A 5.06200e-03  myo.m    rs3821236
# 34:  2:191038032   G   A 1.01225e-05  myo.r    rs3821236
# 35:  2:191038032   G   A 1.39500e-01   pm.m    rs3821236
# 36:  2:191038032   G   A 6.10401e-03   pm.r    rs3821236
# 37:   3:28033679   G   T 1.13920e-01  dmy.m     rs819991
# 38:   3:28033679   G   T 2.91717e-02  dmy.r     rs819991
# 39:   3:28033679   G   T 4.39394e-01  jdm.m     rs819991
# 40:   3:28033679   G   T 2.91194e-01  jdm.r     rs819991
# 41:   3:28033679   G   T 1.64494e-01 jo1m.r     rs819991
# 42:   3:28033679   G   T 5.23297e-02  myo.m     rs819991
# 43:   3:28033679   G   T 1.05722e-03  myo.r     rs819991
# 44:   3:28033679   G   T 2.30833e-01   pm.m     rs819991
# 45:   3:28033679   G   T 1.16681e-01   pm.r     rs819991
# 46:  7:129023119   T   G 2.01914e-03  dmy.m   rs13236009
# 47:  7:129023119   T   G 1.00165e-02  dmy.r   rs13236009
# 48:  7:129023119   T   G 7.47704e-01  jdm.m   rs13236009
# 49:  7:129023119   T   G 2.05076e-01  jdm.r   rs13236009
# 50:  7:129023119   T   G 1.16643e-03 jo1m.r   rs13236009
# 51:  7:129023119   T   G 1.34854e-03  myo.m   rs13236009
# 52:  7:129023119   T   G 6.41047e-05  myo.r   rs13236009
# 53:  7:129023119   T   G 2.29058e-02   pm.m   rs13236009
# 54:  7:129023119   T   G 5.17526e-02   pm.r   rs13236009
# 55:   8:11486464   C   T 2.69300e-03  dmy.m    rs2736340
# 56:   8:11486464   C   T 1.43846e-02  dmy.r    rs2736340
# 57:   8:11486464   C   T 4.18800e-03  jdm.m    rs2736340
# 58:   8:11486464   C   T 1.85708e-03  jdm.r    rs2736340
# 59:   8:11486464   C   T 2.28270e-01 jo1m.r    rs2736340
# 60:   8:11486464   C   T 1.52900e-04  myo.m    rs2736340
# 61:   8:11486464   C   T 1.92214e-04  myo.r    rs2736340
# 62:   8:11486464   C   T 2.96900e-02   pm.m    rs2736340
# 63:   8:11486464   C   T 1.03374e-02   pm.r    rs2736340
#              pid REF ALT           P  study bestsnp.rsid

coloc[bestsnp.rsid == "rs3821236" & trait.myos == "myo.r"]
