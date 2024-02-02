# This is a short snippet to apply coloc to anti-Jo1+ and T1D on rs2476601.


## Set working directory
setwd("/home/gr440/rds/rds-cew54-basis/Projects/myositis-IMD/code")

## Load packages and required
library(data.table)
setDTthreads(20)
library(magrittr)
library(coloc)

# Load list of files used in this study
fl <- fread("../tables/ST_all_datasets.tsv")

fl[Label %in% c("Jo1+ Myositis (Rothwell)","Type 1 diabetes") & in.selection == "Yes"]
# JO1M_Rothwell_up_1, E4_DM1_FinnGen_FinnGenR7_1

# rs2476601 coordinates in hg38
chr = 1
bp = 113834946

w <- 1e6 # Same window choice as for coloc analyses
st <- bp - w
en <- bp + w

# Get sample sizes from Metadata
d  <- fread("../data/Metadata_20230503-v1.tsv") %>% .[ Trait %in% c("JO1M_Rothwell_up_1", "E4_DM1_FinnGen_FinnGenR7_1")]
d[, s:=N1/(N1+N0)] # Proportion of cases, required for coloc


# Get SNP-dense data 
jo1 <- fread("~/rds/rds-cew54-basis/02-Processed/JO1M_Rothwell_up_1-hg38.tsv.gz")  
jo1[, CHR38 :=as.numeric(CHR38)][ , pid:=paste(CHR38, BP38, sep=":")]
jo1 <- jo1[CHR38 == chr & BP38 > st & BP38 < en & !is.na(SE) & !is.na(BETA) & !duplicated(pid)][, .(pid, CHR38, BP38, REF, ALT, BETA, SE, P)]
jo1[CHR38 == chr & BP38 == bp]
#            pid CHR38      BP38 REF ALT     BETA       SE          P
# 1: 1:113834946     1 113834946   A   G -0.40666 0.132924 0.00272829

t1d <- fread("https://storage.googleapis.com/finngen-public-data-r7/summary_stats/finngen_R7_E4_DM1.gz") # We need to redownload the FinnGen dataset
setnames(t1d, c("#chrom","pos"), c("CHR38","BP38"))
t1d[ , pid:=paste(CHR38, BP38, sep=":")]
t1d <- t1d[,.(pid, CHR38, BP38, REF=ref, ALT=alt, BETA=beta, SE=sebeta, P=pval)]
t1d <- t1d[CHR38 == chr & BP38 > st & BP38 < en & !is.na(SE) & !is.na(BETA) & !duplicated(pid)]
t1d[CHR38 == chr & BP38 == bp]
#            pid CHR38      BP38 REF ALT      BETA        SE           P
# 1: 1:113834946     1 113834946   A   G -0.400244 0.0201905 1.86982e-87

jo1.l <- list(snp=jo1$pid, beta=jo1$BETA, varbeta=jo1$SE^2, type="cc", s=d$s[1])
t1d.l <- list(snp=t1d$pid, beta=t1d$BETA, varbeta=t1d$SE^2, type="cc", s=d$s[2])
result=coloc.abf(jo1.l,t1d.l, p12 = 5e-6)
