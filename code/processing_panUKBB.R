# PROCESSING UKBB FILE

# For coloc analyses, we needed dense SNP datasets, which we lacked for PanUKBB and FinnGen
# Most selected datasets lacking dense SNP data were FinnGen, which are in hg38 and are easily usable.
# However, there was one PanUKBB dataset that needed a special treatment.

setwd("/home/gr440/rds/rds-cew54-basis/Projects/myositis-IMD/data/fg_sumstats")

library(data.table)
setDTthreads(20)

# Transform into .gz so we can import it in fread
system("zcat 20002_1225_PanUKBB_PanUKBBR2_1.bgz > 20002_1225_PanUKBB_PanUKBBR2_1.tsv ; gzip 20002_1225_PanUKBB_PanUKBBR2_1.tsv.gz")

f <- fread("20002_1225_PanUKBB_PanUKBBR2_1.tsv.gz")
f <- f[, .(chr, pos, ref, alt, beta_meta, se_meta, neglog10_pval_meta)]
names(f) <- c("CHR", "BP", "REF", "ALT", "BETA", "SE", "P")
f <- f[complete.cases(f)]
f[, P:=10^-P] # P-values are stored as -log10(P), so we'll revert them -- even though we'll likely don't need the P-values

# Save it
fwrite(f, "20002_1225_PanUKBB_PanUKBBR2_1-hr.tsv.gz", sep = "\t")

# Pipeline it 
system("~/rds/rds-cew54-basis/GWAS_tools/01-Pipeline/pipeline_v5.3.2_beta.sh -f 20002_1225_PanUKBB_PanUKBBR2_1-hr.tsv.gz")
