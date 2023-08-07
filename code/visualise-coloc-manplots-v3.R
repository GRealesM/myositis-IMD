
## Load packages and required
library(data.table)
library(magrittr)
library(coloc)
library(annotSnpStats)
library(ggplot2)
library(cowplot)
setDTthreads(14)
setwd("/home/gr440/rds/rds-cew54-basis/Projects/myositis-IMD/code")

#######
# Recreate data object
# Load files, as before
files = list(   `DM (M)` = "DMY_Miller_26291516_1-hg38.tsv.gz",
                `DM (R)`= "DMY_Rothwell_up_1-hg38.tsv.gz",
                `JDM (M)` = "JDM_Miller_26291516_1-hg38.tsv.gz",
                `JDM (R)` = "JDM_Rothwell_up_1-hg38.tsv.gz",
                `Anti-Jo1+ (R)` = "JO1M_Rothwell_up_1-hg38.tsv.gz",
                `IIM (M)` = "MYO_Miller_26291516_1-hg38.tsv.gz",
                `IIM (R)` = "IIM_Rothwell_up_1-hg38.tsv.gz",
                `PM (M)` =   "PM_Miller_26291516_1-hg38.tsv.gz",
                `PM (R)` = "PM_Rothwell_up_1-hg38.tsv.gz",
                ITP =  "D3_ITP_FinnGen_FinnGenR7_1-hg38.tsv.gz", 
                `IgG+ NMO` =  "NMOIGGp_Estrada_29769526_1-hg38.tsv.gz",
                PBC = "CHIRBIL_PRIM_FinnGen_FinnGenR7_1-hg38.tsv.gz",
                JIA = "JIA_LopezIsac_33106285_1-hg38.tsv.gz",
                `CR(E)ST` = "CREST_FinnGen_FinnGenR7_1-hg38.tsv.gz",
                MG  = "MYG_Chia_35074870_1-hg38.tsv.gz",
                Felty =  "FELTY_FinnGen_FinnGenR7_1-hg38.tsv.gz",
                SjS =   "SJOS_Lessard_up_1-hg38.tsv.gz",
                SSc  = "SSC_LopezIsac_31672989_1-hg38.tsv.gz",
                GPA = "M13_WEGENER_FinnGen_FinnGenR7_1-hg38.tsv.gz",
                `MPO+ AAV` = "AAVMPO_Wong_up_1-hg38.tsv.gz", 
                `PR3+ AAV` =  "AAVPR3_Wong_up_1-hg38.tsv.gz",
                HyperThy = "20002_1225_PanUKBB_PanUKBBR2_1-hg38.tsv.gz",
                HypoThy = "HYPOTHYROIDISM_FinnGen_FinnGenR7_1-hg38.tsv.gz",
                PR = "M13_PALINDROMIC_FinnGen_FinnGenR7_1-hg38.tsv.gz",
                BioMedRhe = "RX_RHEUMA_BIOLOGICAL_FinnGen_FinnGenR7_1-hg38.tsv.gz") %>%
    lapply(.,  function(f) (file.path("~/rds/rds-cew54-basis/02-Processed",f)))
stopifnot(all(file.exists(unlist(files))))

data=lapply(files, fread)

## itp, pbc, crest, felty, wegen, hyperthy, and hypoty are not dense datasets. 
names(data)[c(10,12,14,16,19,22:25)]= c("ITP.local", "PBC.local", "CREST.local",   "Felty.local", "GPA.local", "HyperThy.local",  "HypoThy.local", "PR.local", "BioMedRhe.local")

dir.create("../data/fg_sumstats")
if(!file.exists("../data/fg_sumstats/finngen_R7_D3_ITP.gz"))
    system("wget ***REMOVED***finngen_R7_D3_ITP.gz -O ../data/fg_sumstats/finngen_R7_D3_ITP.gz")
if(!file.exists("../data/fg_sumstats/finngen_R7_CHIRBIL_PRIM.gz"))
    system("wget ***REMOVED***finngen_R7_CHIRBIL_PRIM.gz -O ../data/fg_sumstats/finngen_R7_CHIRBIL_PRIM.gz")
if(!file.exists("../data/fg_sumstats/finngen_R7_CREST.gz"))
    system("wget ***REMOVED***finngen_R7_CREST.gz -O ../data/fg_sumstats/finngen_R7_CREST.gz")
if(!file.exists("../data/fg_sumstats/finngen_R7_FELTY.gz"))
    system("wget ***REMOVED***finngen_R7_FELTY.gz -O ../data/fg_sumstats/finngen_R7_FELTY.gz")
if(!file.exists("../data/fg_sumstats/finngen_R7_M13_WEGENER.gz"))
    system("wget ***REMOVED***finngen_R7_M13_WEGENER.gz -O ../data/fg_sumstats/finngen_R7_M13_WEGENER.gz")
if(!file.exists("../data/fg_sumstats/finngen_R7_AUTOIMMUNE_HYPERTHYROIDISM.gz"))
    system("wget ***REMOVED***finngen_R7_AUTOIMMUNE_HYPERTHYROIDISM.gz -O ../data/fg_sumstats/finngen_R7_AUTOIMMUNE_HYPERTHYROIDISM.gz")
if(!file.exists("../data/fg_sumstats/finngen_R7_HYPOTHYROIDISM.gz"))
    system("wget ***REMOVED***finngen_R7_HYPOTHYROIDISM.gz -O ../data/fg_sumstats/finngen_R7_HYPOTHYROIDISM.gz")
if(!file.exists("../data/fg_sumstats/finngen_R7_M13_PALINDROMIC.gz"))
    system("wget ***REMOVED***finngen_R7_M13_PALINDROMIC.gz -O ../data/fg_sumstats/finngen_R7_M13_PALINDROMIC.gz")
if(!file.exists("../data/fg_sumstats/finngen_R7_RX_RHEUMA_BIOLOGICAL.gz"))
    system("wget ***REMOVED***finngen_R7_RX_RHEUMA_BIOLOGICAL.gz -O ../data/fg_sumstats/finngen_R7_RX_RHEUMA_BIOLOGICAL.gz")


ITP=fread("../data/fg_sumstats/finngen_R7_D3_ITP.gz")
ITP$pid=paste(ITP[["#chrom"]], ITP$pos,sep=":")
table(ITP$pid %in% snps$pid38)
setnames(ITP, c("#chrom","pos"), c("CHR38","BP38"))
data$ITP=ITP[,.(pid, CHR38, BP38, REF=ref, ALT=alt, BETA=beta, SE=sebeta, P=pval)]

PBC=fread("../data/fg_sumstats/finngen_R7_CHIRBIL_PRIM.gz")
PBC$pid=paste(PBC[["#chrom"]],PBC$pos,sep=":")
table(PBC$pid %in% snps$pid38)
setnames(PBC, c("#chrom","pos"), c("CHR38","BP38"))
data$PBC=PBC[,.(pid, CHR38, BP38, REF=ref, ALT=alt, BETA=beta, SE=sebeta, P=pval)]

crest=fread("../data/fg_sumstats/finngen_R7_CREST.gz")
crest$pid=paste(crest[["#chrom"]],crest$pos,sep=":")
table(crest$pid %in% snps$pid38)
setnames(crest, c("#chrom","pos"), c("CHR38","BP38"))
data$`CR(E)ST`=crest[,.(pid, CHR38, BP38, REF=ref, ALT=alt, BETA=beta, SE=sebeta, P=pval)]

Felty=fread("../data/fg_sumstats/finngen_R7_FELTY.gz")
Felty$pid=paste(Felty[["#chrom"]],Felty$pos,sep=":")
table(Felty$pid %in% snps$pid38)
setnames(Felty, c("#chrom","pos"), c("CHR38","BP38"))
data$Felty=Felty[,.(pid, CHR38, BP38, REF=ref, ALT=alt, BETA=beta, SE=sebeta, P=pval)]

GPA=fread("../data/fg_sumstats/finngen_R7_M13_WEGENER.gz")
GPA$pid=paste(GPA[["#chrom"]],GPA$pos,sep=":")
table(GPA$pid %in% snps$pid38)
setnames(GPA, c("#chrom","pos"), c("CHR38","BP38"))
data$GPA=GPA[,.(pid, CHR38, BP38, REF=ref, ALT=alt, BETA=beta, SE=sebeta, P=pval)]

hyperthy=fread("../data/fg_sumstats/finngen_R7_AUTOIMMUNE_HYPERTHYROIDISM.gz")
hyperthy$pid=paste(hyperthy[["#chrom"]],hyperthy$pos,sep=":")
table(hyperthy$pid %in% snps$pid38)
setnames(hyperthy, c("#chrom","pos"), c("CHR38","BP38"))
data$HyperThy=hyperthy[,.(pid, CHR38, BP38, REF=ref, ALT=alt, BETA=beta, SE=sebeta, P=pval)]

hypothy=fread("../data/fg_sumstats/finngen_R7_HYPOTHYROIDISM.gz")
hypothy$pid=paste(hypothy[["#chrom"]],hypothy$pos,sep=":")
table(hypothy$pid %in% snps$pid38)
setnames(hypothy, c("#chrom","pos"), c("CHR38","BP38"))
data$HypoThy=hypothy[,.(pid, CHR38, BP38, REF=ref, ALT=alt, BETA=beta, SE=sebeta, P=pval)]

PR=fread("../data/fg_sumstats/finngen_R7_M13_PALINDROMIC.gz")
PR$pid=paste(PR[["#chrom"]],PR$pos,sep=":")
table(PR$pid %in% snps$pid38)
setnames(PR, c("#chrom","pos"), c("CHR38","BP38"))
data$PR=PR[,.(pid, CHR38, BP38, REF=ref, ALT=alt, BETA=beta, SE=sebeta, P=pval)]

BMR=fread("../data/fg_sumstats/finngen_R7_M13_PALINDROMIC.gz")
BMR$pid=paste(BMR[["#chrom"]],BMR$pos,sep=":")
table(BMR$pid %in% snps$pid38)
setnames(BMR, c("#chrom","pos"), c("CHR38","BP38"))
data$BioMedRhe=BMR[,.(pid, CHR38, BP38, REF=ref, ALT=alt, BETA=beta, SE=sebeta, P=pval)]


###########

# Load index (aka coloc results)

coloc <- fread("../data/coloc_results_dfilt-v3.tsv")

# Import mapped genes 
mg <- fread("../data/mapped.genes_v2.tsv") %>% unique
mg <- mg[,.(SNPID, pid, nearestGene)]

# Add some info to the coloc table, and extract rsids to map

coloc <- merge(coloc, mg, by="pid") # First on driver SNPs
setnames(coloc, c("SNPID", "nearestGene"), c("driver.rsid", "driver.nearestgene"))
coloc <- merge(coloc, mg, by.x="bestsnp", by.y = "pid", all.x = TRUE) # Then on candidate snps. Bear in mind that we only mapped candidate SNPs with H4 > 0.5
setnames(coloc, c("SNPID", "nearestGene"), c("bestsnp.rsid", "bestsnp.nearestgene"))

# Fix inappropriate mappings
coloc[ driver.rsid == "rs2476601", driver.nearestgene:="PTPN22"]
coloc[ bestsnp.rsid == "rs2476601", bestsnp.nearestgene:="PTPN22"]
coloc[ driver.rsid == "rs991817", driver.nearestgene:="SH2B3"]
coloc[ bestsnp.rsid == "rs991817", bestsnp.nearestgene:="SH2B3"]

# Add gene/ driver SNP label
coloc[, dlabel:=paste( driver.nearestgene, driver.rsid, sep=" / ")]

# Add novelty
coloc[H4 > 0.5 & pbest.myos.region > 5e-8, bestsnp.novel:="Yes"][is.na(bestsnp.novel), bestsnp.novel:="No"] # Add novelty

index <- coloc # we were using index instead of coloc. so let's simply use this other name


#############

# Now explore the SNPs in more detail
index[ H4>.5, unique(pid)]
# 8 SNPs
index[ H4>.5, unique(trait.other)]
# We have coloc associations with 8 (out of 14) IMDs
# These are
# [1] "JIA"      "HypoThy"  "MPO+ AAV" "PBC"      "HyperThy" "MG"       "SSc"     
# [8] "SjS" 

## plot these signals
plotter=function(pid,w=1e+6) {
    chr=sub(":.*","",pid)  %>% as.numeric()
    bp=sub(".*:","",pid)  %>% as.numeric()
    st=bp-w
    en=bp+w
    wh=which(index$pid==pid & (index$H4>.5))
    print(index[wh])
    rsid = index[wh , unique(driver.rsid)]
    traits=unique(c(index$trait.myos[wh],index$trait.other[wh]))
    dp=lapply(data[traits], function(d)
        d[CHR38==chr & BP38>st & BP38<en & !is.na(SE) & !is.na(BETA) & !duplicated(pid),
          .(pid,CHR38,BP38,BETA,SE,P)])
    for(i in seq_along(traits))
        dp[[i]]$trait=traits[i]
    dp %<>% rbindlist()
    ggplot(dp, aes(x=BP38, y=-log10(P))) + 
            geom_point() + 
            facet_grid(trait~.,scales="free_y") + 
            geom_vline(xintercept=bp,col="red")+
            ggtitle(paste0(pid, " / ", rsid))+
            theme_cowplot() + 
            theme(axis.title.x = element_blank(),
                  strip.background = element_rect(colour="black", fill = "white"),
                  )
}

plots <- lapply(index[ H4>.5, unique(pid)], plotter)
names(plots) <- index[ H4>.5, unique(pid)]

index[ H4>.5, unique(pid)]

index[ H4 > 0.5, .(bestsnp, bestsnp.novel)] %>% unique


ggsave("../figures/coloc_chr1.png", plots$`1:113834946`, height = 11, width = 8, bg="white")
ggsave("../figures/coloc_chr2_1.png", plots$`2:100215693`, height = 5, width = 8, bg="white")
ggsave("../figures/coloc_chr2_2.png", plots$`2:190670850`, height = 5, width = 8, bg="white")
ggsave("../figures/coloc_chr3.png", plots$`3:28029953`, height = 5, width = 8, bg="white")
ggsave("../figures/coloc_chr7.png", plots$`7:128954129`, height = 14, width = 8, bg="white")
ggsave("../figures/coloc_chr8.png", plots$`8:11491677`, height = 8, width = 8, bg="white")
ggsave("../figures/coloc_chr12_1.png", plots$`12:110972733`, height = 5, width = 8, bg="white")
ggsave("../figures/coloc_chr12_2.png", plots$`12:112468611`, height = 6, width = 8, bg="white")

# ggsave("../figures/coloc_chr11.png", plots$`11:64329761`, height = 8, width = 8, bg="white")
# ggsave("../figures/coloc_chr17.png", plots$`17:75373341`, height = 5, width = 8, bg="white")
# ggsave("../figures/coloc_chr1.png", plots$`1:113834946`, height = 11, width = 8, bg="white")
# ggsave("../figures/coloc_chr4.png", plots$`4:122194347`, height = 6, width = 8, bg="white")
# ggsave("../figures/coloc_chr6.png", plots$`6:167124106`, height = 6, width = 8, bg="white")
# ggsave("../figures/coloc_chr7_2.png", plots$`7:37397251`, height = 5, width = 8, bg="white")
