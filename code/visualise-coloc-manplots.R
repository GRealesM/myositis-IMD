## I wonder if there is way to identify the specific loci that led to finding
# 1) similarity to NMO and Sjogren's and 
# 2) Similarity to Scleroderma/CREST and dermatomyositis? 
# These group of loci might give us a clue to the neurological and vascular components of PAPS, respectively.

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
files = list(   dmy.m = "DMY_Miller_26291516_1-hg38.tsv.gz",
                dmy.r = "DMY_Rothwell_up_1-hg38.tsv.gz",
                jdm.m = "JDM_Miller_26291516_1-hg38.tsv.gz",
                jdm.r = "JDM_Rothwell_up_1-hg38.tsv.gz",
                jo1m.r = "JO1M_Rothwell_up_1-hg38.tsv.gz",
                myo.m = "MYO_Miller_26291516_1-hg38.tsv.gz",
                myo.r = "IIM_Rothwell_up_1-hg38.tsv.gz",
                pm.m =   "PM_Miller_26291516_1-hg38.tsv.gz",
                pm.r = "PM_Rothwell_up_1-hg38.tsv.gz",
                itp =  "D3_ITP_FinnGen_FinnGenR7_1-hg38.tsv.gz", 
                nmo =  "NMOIGGp_Estrada_29769526_1-hg38.tsv.gz",
                pbc = "CHIRBIL_PRIM_FinnGen_FinnGenR7_1-hg38.tsv.gz",
                jia = "JIA_LopezIsac_33106285_1-hg38.tsv.gz",
                crest = "CREST_FinnGen_FinnGenR7_1-hg38.tsv.gz",
                myag  = "MYG_Chia_35074870_1-hg38.tsv.gz",
                felty =  "FELTY_FinnGen_FinnGenR7_1-hg38.tsv.gz",
                sjos =   "SJOS_Lessard_up_1-hg38.tsv.gz",
                ssc  = "SSC_LopezIsac_31672989_1-hg38.tsv.gz",
                wegen = "M13_WEGENER_FinnGen_FinnGenR7_1-hg38.tsv.gz",
                mpoaav = "AAVMPO_Wong_up_1-hg38.tsv.gz", 
                pr3aav =  "AAVPR3_Wong_up_1-hg38.tsv.gz",
                hyperthy = "20002_1225_PanUKBB_PanUKBBR2_1-hg38.tsv.gz",
                hypothy = "HYPOTHYROIDISM_FinnGen_FinnGenR7_1-hg38.tsv.gz") %>%
    lapply(.,  function(f) (file.path("~/rds/rds-cew54-basis/02-Processed",f)))
stopifnot(all(file.exists(unlist(files))))

data=lapply(files, fread)

## itp, pbc, crest, felty, wegen, hyperthy, and hypoty are not dense datasets. 
names(data)[c(10,12,14,16,19,22:23)]= c("itp.local", "pbc.local", "crest.local",   "felty.local", "wegen.local", "hyperthy.local",  "hypothy.local")

dir.create("../data/fg_sumstats")
if(!file.exists("../data/fg_sumstats/finngen_R7_D3_ITP.gz"))
    system("wget https://storage.googleapis.com/finngen-public-data-r7/summary_stats/finngen_R7_D3_ITP.gz -O ../data/fg_sumstats/finngen_R7_D3_ITP.gz")
if(!file.exists("../data/fg_sumstats/finngen_R7_CHIRBIL_PRIM.gz"))
    system("wget https://storage.googleapis.com/finngen-public-data-r7/summary_stats/finngen_R7_CHIRBIL_PRIM.gz -O ../data/fg_sumstats/finngen_R7_CHIRBIL_PRIM.gz")
if(!file.exists("../data/fg_sumstats/finngen_R7_CREST.gz"))
    system("wget https://storage.googleapis.com/finngen-public-data-r7/summary_stats/finngen_R7_CREST.gz -O ../data/fg_sumstats/finngen_R7_CREST.gz")
if(!file.exists("../data/fg_sumstats/finngen_R7_FELTY.gz"))
    system("wget https://storage.googleapis.com/finngen-public-data-r7/summary_stats/finngen_R7_FELTY.gz -O ../data/fg_sumstats/finngen_R7_FELTY.gz")
if(!file.exists("../data/fg_sumstats/finngen_R7_M13_WEGENER.gz"))
    system("wget https://storage.googleapis.com/finngen-public-data-r7/summary_stats/finngen_R7_M13_WEGENER.gz -O ../data/fg_sumstats/finngen_R7_M13_WEGENER.gz")
if(!file.exists("../data/fg_sumstats/finngen_R7_AUTOIMMUNE_HYPERTHYROIDISM.gz"))
    system("wget https://storage.googleapis.com/finngen-public-data-r7/summary_stats/finngen_R7_AUTOIMMUNE_HYPERTHYROIDISM.gz -O ../data/fg_sumstats/finngen_R7_AUTOIMMUNE_HYPERTHYROIDISM.gz")
if(!file.exists("../data/fg_sumstats/finngen_R7_HYPOTHYROIDISM.gz"))
    system("wget https://storage.googleapis.com/finngen-public-data-r7/summary_stats/finngen_R7_HYPOTHYROIDISM.gz -O ../data/fg_sumstats/finngen_R7_HYPOTHYROIDISM.gz")


itp=fread("../data/fg_sumstats/finngen_R7_D3_ITP.gz")
itp$pid=paste(itp[["#chrom"]], itp$pos,sep=":")
table(itp$pid %in% snps$pid38)
setnames(itp, c("#chrom","pos"), c("CHR38","BP38"))
data$itp=itp[,.(pid, CHR38, BP38, REF=ref, ALT=alt, BETA=beta, SE=sebeta, P=pval)]

pbc=fread("../data/fg_sumstats/finngen_R7_CHIRBIL_PRIM.gz")
pbc$pid=paste(pbc[["#chrom"]],pbc$pos,sep=":")
table(pbc$pid %in% snps$pid38)
setnames(pbc, c("#chrom","pos"), c("CHR38","BP38"))
data$pbc=pbc[,.(pid, CHR38, BP38, REF=ref, ALT=alt, BETA=beta, SE=sebeta, P=pval)]

crest=fread("../data/fg_sumstats/finngen_R7_CREST.gz")
crest$pid=paste(crest[["#chrom"]],crest$pos,sep=":")
table(crest$pid %in% snps$pid38)
setnames(crest, c("#chrom","pos"), c("CHR38","BP38"))
data$crest=crest[,.(pid, CHR38, BP38, REF=ref, ALT=alt, BETA=beta, SE=sebeta, P=pval)]

felty=fread("../data/fg_sumstats/finngen_R7_FELTY.gz")
felty$pid=paste(felty[["#chrom"]],felty$pos,sep=":")
table(felty$pid %in% snps$pid38)
setnames(felty, c("#chrom","pos"), c("CHR38","BP38"))
data$felty=felty[,.(pid, CHR38, BP38, REF=ref, ALT=alt, BETA=beta, SE=sebeta, P=pval)]

wegen=fread("../data/fg_sumstats/finngen_R7_M13_WEGENER.gz")
wegen$pid=paste(wegen[["#chrom"]],wegen$pos,sep=":")
table(wegen$pid %in% snps$pid38)
setnames(wegen, c("#chrom","pos"), c("CHR38","BP38"))
data$wegen=wegen[,.(pid, CHR38, BP38, REF=ref, ALT=alt, BETA=beta, SE=sebeta, P=pval)]

hyperthy=fread("../data/fg_sumstats/finngen_R7_AUTOIMMUNE_HYPERTHYROIDISM.gz")
hyperthy$pid=paste(hyperthy[["#chrom"]],hyperthy$pos,sep=":")
table(hyperthy$pid %in% snps$pid38)
setnames(hyperthy, c("#chrom","pos"), c("CHR38","BP38"))
data$hyperthy=hyperthy[,.(pid, CHR38, BP38, REF=ref, ALT=alt, BETA=beta, SE=sebeta, P=pval)]

hypothy=fread("../data/fg_sumstats/finngen_R7_HYPOTHYROIDISM.gz")
hypothy$pid=paste(hypothy[["#chrom"]],hypothy$pos,sep=":")
table(hypothy$pid %in% snps$pid38)
setnames(hypothy, c("#chrom","pos"), c("CHR38","BP38"))
data$hypothy=hypothy[,.(pid, CHR38, BP38, REF=ref, ALT=alt, BETA=beta, SE=sebeta, P=pval)]

###########

# Load index (aka coloc results)

coloc <- fread("../tables/coloc_results_dfilt.tsv")

# Import mapped genes 
mg <- fread("../data/mapped.genes.tsv") %>% unique
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

# Add novelty on bestsnps
coloc[H4 > 0.5 & pbest.myos > 5e-8 & pbest.myos.region > 5e-8, bestsnp.novel:="Yes"] # Add novelty
# Add gene/ driver SNP label
coloc[, dlabel:=paste( driver.nearestgene, driver.rsid, sep=" / ")]

# Add proper labels
ml <- data.table(mlabel = c("DM (Miller)", "DM (Rothwell)", "JDM (Miller)","JDM (Rothwell)", "Jo1+ (Rothwell)", "Myositis (Miller)", "Myositis (Rothwell)", "PM (Miller)", "PM (Rothwell)"),
                 trait.myos = c("dmy.m", "dmy.r", "jdm.m", "jdm.r", "jo1m.r", "myo.m", "myo.r", "pm.m", "pm.r"))
coloc <- merge(coloc, ml, by="trait.myos")

index <- coloc # we were using index instead of coloc. so let's simply use this other name


#############

# Now explore the SNPs in more detail
index[ H4>.5, unique(pid)]
# 8 SNPs
index[ H4>.5, unique(trait.other)]
# We have coloc associations with 8 (out of 14) IMDs
# These are
# [1] "ssc"      "sjos"     "jia"      "hypothy"  "mpoaav"   "pbc"      "hyperthy"
# [8] "myag"

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
          .(pid,CHR38,BP38,BETA,SE,P=2*pnorm(-abs(BETA)/SE))])
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



