#########################################
##                                     ##
##            RUNNING COLOC            ##
##                                     ##
#########################################

# Author: Guillermo Reales 
# Date last updated: 2024/04/25

# Background: After running DPMUnc and computing the Bhattacharyya distance, we selected a reduced number of IMD to run coloc on. Then we'll run the
# pairwise FDR procedure to select SNPs and run coloc for myositis-IMD pairs at selected SNPs.

# This script will
# * Import the list of IMD.
# * Run the pairwise FDR procedure.
# * Import and download dense-SNP datasets when not directly available.
# * Select SNPs and myositis-IMD to run coloc on
# * Run coloc

##########################################


## Set working directory
setwd("/home/gr440/rds/rds-cew54-basis/Projects/myositis-IMD/code")

## Load packages and required
library(data.table)
library(magrittr)
library(coloc)
#remotes::install_github("chr1swallace/annotSnpStats")
library(annotSnpStats)
library(ggplot2)
library(cowplot)
setDTthreads(15)

fg7url <- "***REMOVED***"

dpres <- fread("../tables/ST4_coloc_diseases.tsv")

dpres[,.(Trait, Label, coloc_Label)] 
#                                      Trait                               Label    coloc_Label
#  1:                  DMY_Miller_26291516_1            Dermatomyositis (Miller)         DM (M)
#  2:                      DMY_Rothwell_up_1          Dermatomyositis (Rothwell)         DM (R)
#  3:                  MYO_Miller_26291516_1                        IIM (Miller)        IIM (M)
#  4:                      IIM_Rothwell_up_1                      IIM (Rothwell)        IIM (R)
#  5:                     JO1M_Rothwell_up_1            Jo1+ Myositis (Rothwell)  Anti-Jo1+ (R)
#  6:                  JDM_Miller_26291516_1   Juvenile Dermatomyositis (Miller)        JDM (M)
#  7:                      JDM_Rothwell_up_1 Juvenile Dermatomyositis (Rothwell)        JDM (R)
#  8:                   PM_Miller_26291516_1               Polymyositis (Miller)         PM (M)
#  9:                       PM_Rothwell_up_1             Polymyositis (Rothwell)         PM (R)
# 10:              CREST_FinnGen_FinnGenR7_1                    CR(E)ST syndrome        CR(E)ST # 10
# 11:                MYGEO_Renton_25643325_1        Early-onset Myastenia Gravis       Early-MG 
# 12:              FELTY_FinnGen_FinnGenR7_1                      Felty syndrome          Felty #12
# 13:         20002_1225_PanUKBB_PanUKBBR2_1      Hyperthyroidism/Thyrotoxicosis       HyperThy #13
# 14: E4_HYTHY_AI_STRICT_FinnGen_FinnGenR7_1   Hypothyroidism, strict autoimmune        HypoThy #14
# 15:             NMOIGGp_Estrada_29769526_1           IgG+ Neuromyelitis Optica       IgG+ NMO
# 16:               JIA_LopezIsac_33106285_1       Juvenile Idiopathic Arthritis            JIA
# 17:                MYGLO_Renton_25643325_1         Late-onset Myastenia Gravis        Late-MG
# 18:                       AAVMPO_Wong_up_1                            MPO+ AAV       MPO+ AAV
# 19:                    MYG_Chia_35074870_1                    Myastenia Gravis             MG
# 20:    M13_PALINDROMIC_FinnGen_FinnGenR7_1              Palindromic rheumatism             PR #20
# 21:       CHIRBIL_PRIM_FinnGen_FinnGenR7_1         Primary biliary cholangitis            PBC #21
# 22:         M13_RHEUMA_FinnGen_FinnGenR7_1                Rheumatoid arthritis             RA #22
# 23:                      SJOS_Lessard_up_1                  Sjögren's syndrome            SjS
# 24:               SSC_LopezIsac_31672989_1                  Systemic Sclerosis            SSc
# 25:            M13_SLE_FinnGen_FinnGenR7_1        Systemic lupus erythematosus            SLE #25
# 26:        M13_WEGENER_FinnGen_FinnGenR7_1              Wegener granulomatosis            GPA #26
#                                      Trait                               Label    coloc_Label



# Take the opportunity to create proper labels
files = list(   `DM (M)` = "DMY_Miller_26291516_1-hg38.tsv.gz",
                `DM (R)`= "DMY_Rothwell_up_1-hg38.tsv.gz",
                `JDM (M)` = "JDM_Miller_26291516_1-hg38.tsv.gz",
                `JDM (R)` = "JDM_Rothwell_up_1-hg38.tsv.gz",
                `Anti-Jo1+ (R)` = "JO1M_Rothwell_up_1-hg38.tsv.gz",
                `IIM (M)` = "MYO_Miller_26291516_1-hg38.tsv.gz",
                `IIM (R)` = "IIM_Rothwell_up_1-hg38.tsv.gz",
                `PM (M)` =   "PM_Miller_26291516_1-hg38.tsv.gz",
                `PM (R)` = "PM_Rothwell_up_1-hg38.tsv.gz",

                `CR(E)ST` = "CREST_FinnGen_FinnGenR7_1-hg38.tsv.gz", 
                EOMG = "MYGEO_Renton_25643325_1-hg38.tsv.gz",
                Felty =  "FELTY_FinnGen_FinnGenR7_1-hg38.tsv.gz", 
                HyperThy = "20002_1225_PanUKBB_PanUKBBR2_1-hg38.tsv.gz",
                HypoThy = "E4_HYTHY_AI_STRICT_FinnGen_FinnGenR7_1-hg38.tsv.gz", 
                `IgG+ NMO` =  "NMOIGGp_Estrada_29769526_1-hg38.tsv.gz", 
                JIA = "JIA_LopezIsac_33106285_1-hg38.tsv.gz",
                LOMG = "MYGLO_Renton_25643325_1-hg38.tsv.gz",
                `MPO+ AAV` = "AAVMPO_Wong_up_1-hg38.tsv.gz", 
                MG = "MYG_Chia_35074870_1-hg38.tsv.gz",
                PR = "M13_PALINDROMIC_FinnGen_FinnGenR7_1-hg38.tsv.gz",
                PBC = "CHIRBIL_PRIM_FinnGen_FinnGenR7_1-hg38.tsv.gz", 
                RA = "M13_RHEUMA_FinnGen_FinnGenR7_1-hg38.tsv.gz", 
                SjS =   "SJOS_Lessard_up_1-hg38.tsv.gz",
                SSc  = "SSC_LopezIsac_31672989_1-hg38.tsv.gz",
                SLE = "M13_SLE_FinnGen_FinnGenR7_1-hg38.tsv.gz",
                GPA = "M13_WEGENER_FinnGen_FinnGenR7_1-hg38.tsv.gz") %>% 
    lapply(.,  function(f) (file.path("~/rds/rds-cew54-basis/02-Processed",f)))
stopifnot(all(file.exists(unlist(files))))

man=fread("~/rds/rds-cew54-basis/03-Bases/IMD_basis/SNP.manifest.38.tsv")
snps=fread("~/rds/rds-cew54-basis/03-Bases/IMD_basis/Manifest_build_translator.tsv")
snps[, pid38:=paste(CHR38, BP38, sep=":")][, pid19:=paste(CHR19, BP19, sep=":")]
snps <- merge(snps, man, by="pid38")

data=lapply(files, fread, tmpdir = "tmp/")


aligner=function(d) {
    d[,pid:=paste(CHR38,BP38,sep=":")]
    ## table(d$pid %in% snps$pid38)
    d=d[pid %in% snps$pid38]
    message("snps found: ",nrow(d), " / ",nrow(snps))
    geno=toupper(paste(d$REF,d$ALT,sep="/"))
    cl=g.class(geno, snps$alleles[ match(d$pid, snps$pid38) ])
    print(table(cl))
    d[ cl %in% c("rev","revcomp"), BETA:=-BETA]
    d=d[ (cl!="impossible") ]
    d
}

data2=lapply(data, aligner)
for(nm in names(files)) 
    data2[[nm]]$trait=nm
data2  %<>% rbindlist(., fill=TRUE)

## restrict to PCs of interest: PC1, PC2, PC3, PC8, PC9, PC12, and PC13
rot=cupcake::rot.pca[,paste0("PC",c(1,2,3,8,9,12,13))]
all(rownames(rot) == snps$pid19) # Check the SNPs are in the same order -- they aren't, let's reorder them to match
snps <- snps[order(pid19)]
all(rownames(rot) == snps$pid19) # now they are
rownames(rot) <- snps$pid38
snps.use=snps$pid38[ rowSums(rot != 0)>0 ] # This works because SNPs are in the same order, even though they belong to different builds.
rot.use <- rot[snps.use,]
pcs.use <- colnames(rot.use)
# Note which SNPs are driver for each 
dv <- apply(rot.use, 1, function(row) pcs.use[which(row != 0)])
dv <- sapply(dv, paste, collapse=", ")
dv <- data.table(pid = names(dv), driver = dv)

# Keep driver SNPs for selected PCs only
data2=data2[ pid %in% snps.use ]
data2[,P:=2 * pnorm( -abs(BETA)/SE )]
data2[,fdr:=p.adjust(P, method="BH"), by="trait"]
summary(data2$fdr)

## viz
data2=data2[ order(as.numeric(CHR38),BP38)]
data2[,x:=cumsum(c(0, pmax(diff(BP38), 0))),by="trait"]
axisdf <- data2[, .(center=(max(x) + min(x))/2), by = CHR38]

# There are 26 datasets, so instead of plotting everything in the same plot, we'll split it into three
traits  <- unique(data2$trait)

mhp1 <- ggplot(data2[ trait %in% traits[1:8]], aes(x=x, y=-log10(fdr), col=factor(as.numeric(CHR38) %% 2))) + 
                    geom_point() +  
                    scale_x_continuous( label = axisdf$CHR38, breaks= axisdf$center ) + 
                    scale_y_continuous(expand = c(0, 0) ) +
                    scale_color_manual(values = rep(c("#3939ff", "#009700"), 22 )) +
                    theme_cowplot(10)+
                    theme(legend.position="none",
                           panel.border = element_blank(),
                           panel.grid.major.x = element_blank(),
                           panel.grid.minor.x = element_blank(),
                           panel.grid.minor.y = element_blank(),
                           axis.title.x = element_blank() )+
                    facet_grid(trait~., scales="free_y", switch = "y")

mhp2 <- ggplot(data2[ trait %in% traits[9:16]], aes(x=x, y=-log10(fdr), col=factor(as.numeric(CHR38) %% 2))) + 
                    geom_point() +  
                    scale_x_continuous( label = axisdf$CHR38, breaks= axisdf$center ) + 
                    scale_y_continuous(expand = c(0, 0) ) +
                    scale_color_manual(values = rep(c("#3939ff", "#009700"), 22 )) +
                    theme_cowplot(10)+
                    theme(legend.position="none",
                           panel.border = element_blank(),
                           panel.grid.major.x = element_blank(),
                           panel.grid.minor.x = element_blank(),
                           panel.grid.minor.y = element_blank(),
                           axis.title.x = element_blank() )+
                    facet_grid(trait~., scales="free_y", switch = "y")
mhp3 <- ggplot(data2[ trait %in% traits[17:26]], aes(x=x, y=-log10(fdr), col=factor(as.numeric(CHR38) %% 2))) + 
                    geom_point() +  
                    scale_x_continuous( label = axisdf$CHR38, breaks= axisdf$center ) + 
                    scale_y_continuous(expand = c(0, 0) ) +
                    scale_color_manual(values = rep(c("#3939ff", "#009700"), 22 )) +
                    theme_cowplot(10)+
                    theme(legend.position="none",
                           panel.border = element_blank(),
                           panel.grid.major.x = element_blank(),
                           panel.grid.minor.x = element_blank(),
                           panel.grid.minor.y = element_blank(),
                           axis.title.x = element_blank() )+
                    facet_grid(trait~., scales="free_y", switch = "y")


mhp <- plot_grid(mhp1, mhp2, mhp3, ncol = 3) # This can be improved, some chromosomes don't seem to align properly with their data points.

ggsave("../figures/manhattan_26IMD.png", mhp, height = 8, width = 17, bg = "white")



###

myos=data2[trait %in% traits[1:9], .(trait, pid, fdr)] # I put myositis datasets at the beginning, so these should be the ones.
not=data2[ !trait %in% traits[1:9], .(trait,pid,fdr)]
myos=merge(myos, not, by=c("pid"), suffixes=c(".myos",".other"), allow.cartesian = TRUE)
myos[,pairwise_fdr:=1 - (1-fdr.myos) * (1-fdr.other)] # P(H0 for either disease) = 1 - P(H1 for PAPS & other | p1, p2). Small means high prob of sharing
summary(myos$pairwise_fdr)
myos[ pairwise_fdr < 0.05 ]
#              pid trait.myos    fdr.myos trait.other    fdr.other pairwise_fdr
#           <char>     <char>       <num>      <char>        <num>        <num>
#   1:  10:6064303    IIM (R) 0.022563360    HyperThy 3.623032e-04  0.022917489
#   2:  10:6064303    IIM (R) 0.022563360     HypoThy 4.583642e-11  0.022563360
#   3:  10:6064303    IIM (R) 0.022563360         JIA 7.159483e-04  0.023263154
#   4:  10:6064303    IIM (R) 0.022563360          RA 2.822467e-06  0.022566119
#   5:  10:6064589    IIM (R) 0.022563360    HyperThy 3.837692e-04  0.022938470
#  ---                                                                         
# 108: 7:128977412    IIM (R) 0.003918517         SLE 5.162179e-08  0.003918568
# 109:  8:11491677    IIM (R) 0.004448022          RA 1.681830e-04  0.004615457
# 110:  8:11491677    IIM (R) 0.004448022         SjS 3.363790e-05  0.004481510
# 111:  8:11491677    IIM (R) 0.004448022         SSc 9.515893e-11  0.004448022
# 112:  8:11491677    IIM (R) 0.004448022         SLE 1.010627e-04  0.004548635

# 112 pairs with pairwise_fdr < 0.05

################################################################################

myos[ pairwise_fdr < 0.05 , unique(pid)]
# 23 unique SNPs

fwrite(myos, "../data/raw_fdr_results.tsv", sep="\t")

# We'll simply consider pairs with pairwise FDR < 0.5

index=myos[ pairwise_fdr < 0.05 ][order(pairwise_fdr)] 

index[ , c("chr","bp"):=tstrsplit(pid,":")  %>% lapply(., as.numeric) ]

# However, some of these fdr-significant driver SNPs are located in the same genomic region, 
# possibly capturing the same signals.
# To avoid having multiple SNPs in the same region, we'll cluster SNPs using a distance approach,
# and select one SNP per cluster, such that we have only one significant driver SNP per region
ss <- split(index[, .(pid, chr, bp)], index[,.(chr)]) # split by chr
ss  <- lapply(ss, unique) 

# This function will compute the distances, call clusters and select clusters with more than one SNP

f=function(d) {
    if(nrow(d)==1)
        return(NULL)
    dist <- as.dist(abs(outer(d$bp, d$bp, "-")))
    dc <- hclust(dist)
    calls <- cutree(dc, h = 1e+6)
    d[, cl:=calls]
    d <- d[, if(.N > 1) .SD, by = cl]
    if(nrow(d) == 0)
        return(NULL)
    d[, cl:=paste0("chr", unique(chr), ".", cl)]
    return(d)
}

cl.snps  <- ss %>% lapply(., f)  %>% rbindlist()

# To select which SNPs to keep, we'll look at pairwise fdr
withfdr <- merge(cl.snps, index[, .(pid, pairwise_fdr)])
tokeep <- withfdr[  , .SD[which.min(pairwise_fdr)] , by=cl][, pid] # For each cluster, keep the one with lowest pairwise_fdr
drop <- withfdr[!pid %in% tokeep, unique(pid)]
message("The following SNPs were dropped: ", paste0(drop, collapse = ", "))

nrow(index)
# 112
length(unique(index$pid))
# 23 SNPs

# Remove SNPs to drop
index <- index[!pid %in% drop]

nrow(index)
# 61
length(unique(index$pid))
# 13

## Next step, bring dense SNP datasets to run coloc

## some datasets are not dense datasets. 
names(data)[c(10, 12:14, 20:22, 25:26)] = c("CR(E)ST.local", "Felty.local","HyperThy.local", "HypoThy.local", "PR.local", "PBC.local", "RA.local" ,"SLE.local","GPA.local")

dir.create("../data/fg_sumstats")

if(!file.exists("../data/fg_sumstats/finngen_R7_CREST.gz"))
    system(paste0("wget ", fg7url, "finngen_R7_CREST.gz -O ../data/fg_sumstats/finngen_R7_CREST.gz"))
if(!file.exists("../data/fg_sumstats/finngen_R7_FELTY.gz"))
    system(paste0("wget ", fg7url, "finngen_R7_FELTY.gz -O ../data/fg_sumstats/finngen_R7_FELTY.gz"))

# Note: this PanUKBB file will need some work to reformat, see below
if(!file.exists("../data/fg_sumstats/20002_1225_PanUKBB_PanUKBBR2_1-hg38.tsv.gz")){
    system("wget https://pan-ukb-us-east-1.s3.amazonaws.com/sumstats_flat_files/categorical-20002-both_sexes-1225.tsv.bgz -O ../data/fg_sumstats/20002_1225_PanUKBB_PanUKBBR2_1.bgz")
    system("Rscript processing_panUKBB.R") # This script will prepare the PanUKBB file to be used by coloc
}   
if(!file.exists("../data/fg_sumstats/finngen_R7_E4_HYTHY_AI_STRICT.gz"))
    system(paste0("wget ", fg7url, "finngen_R7_E4_HYTHY_AI_STRICT.gz -O ../data/fg_sumstats/finngen_R7_E4_HYTHY_AI_STRICT.gz"))
if(!file.exists("../data/fg_sumstats/finngen_R7_M13_PALINDROMIC.gz"))
    system(paste0("wget ", fg7url, "finngen_R7_M13_PALINDROMIC.gz -O ../data/fg_sumstats/finngen_R7_M13_PALINDROMIC.gz"))
if(!file.exists("../data/fg_sumstats/finngen_R7_CHIRBIL_PRIM.gz"))
    system(paste0("wget ", fg7url, "finngen_R7_CHIRBIL_PRIM.gz -O ../data/fg_sumstats/finngen_R7_CHIRBIL_PRIM.gz"))
if(!file.exists("../data/fg_sumstats/finngen_R7_M13_RHEUMA.gz"))
    system(paste0("wget ", fg7url, "finngen_R7_M13_RHEUMA.gz -O ../data/fg_sumstats/finngen_R7_M13_RHEUMA.gz"))
if(!file.exists("../data/fg_sumstats/finngen_R7_M13_SLE.gz"))
    system(paste0("wget ", fg7url, "finngen_R7_M13_SLE.gz -O ../data/fg_sumstats/finngen_R7_M13_SLE.gz"))
if(!file.exists("../data/fg_sumstats/finngen_R7_M13_WEGENER.gz"))
    system(paste0("wget ", fg7url, "finngen_R7_M13_WEGENER.gz -O ../data/fg_sumstats/finngen_R7_M13_WEGENER.gz"))


# Incorporate new files into data

crest=fread("../data/fg_sumstats/finngen_R7_CREST.gz", tmpdir = "tmp/")
crest$pid=paste(crest[["#chrom"]],crest$pos,sep=":")
table(crest$pid %in% snps$pid38)
setnames(crest, c("#chrom","pos"), c("CHR38","BP38"))
data$`CR(E)ST`=crest[,.(pid, CHR38, BP38, REF=ref, ALT=alt, BETA=beta, SE=sebeta, P=pval)]

Felty=fread("../data/fg_sumstats/finngen_R7_FELTY.gz", tmpdir = "tmp/")
Felty$pid=paste(Felty[["#chrom"]],Felty$pos,sep=":")
table(Felty$pid %in% snps$pid38)
setnames(Felty, c("#chrom","pos"), c("CHR38","BP38"))
data$Felty=Felty[,.(pid, CHR38, BP38, REF=ref, ALT=alt, BETA=beta, SE=sebeta, P=pval)]

hyperthy=fread("../data/fg_sumstats/20002_1225_PanUKBB_PanUKBBR2_1-hg38.tsv.gz", tmpdir = "tmp/")
hyperthy[, pid := paste(CHR38,BP38,sep=":")]
table(hyperthy$pid %in% snps$pid38)
data$HyperThy=hyperthy[,.(pid, CHR38, BP38, REF, ALT, BETA, SE, P)]

hypothy=fread("../data/fg_sumstats/finngen_R7_E4_HYTHY_AI_STRICT.gz", tmpdir = "tmp/")
hypothy$pid=paste(hypothy[["#chrom"]],hypothy$pos,sep=":")
table(hypothy$pid %in% snps$pid38)
setnames(hypothy, c("#chrom","pos"), c("CHR38","BP38"))
data$HypoThy=hypothy[,.(pid, CHR38, BP38, REF=ref, ALT=alt, BETA=beta, SE=sebeta, P=pval)]

PR=fread("../data/fg_sumstats/finngen_R7_M13_PALINDROMIC.gz", tmpdir = "tmp/")
PR$pid=paste(PR[["#chrom"]],PR$pos,sep=":")
table(PR$pid %in% snps$pid38)
setnames(PR, c("#chrom","pos"), c("CHR38","BP38"))
data$PR=PR[,.(pid, CHR38, BP38, REF=ref, ALT=alt, BETA=beta, SE=sebeta, P=pval)]

PBC=fread("../data/fg_sumstats/finngen_R7_CHIRBIL_PRIM.gz", tmpdir = "tmp/")
PBC$pid=paste(PBC[["#chrom"]],PBC$pos,sep=":")
table(PBC$pid %in% snps$pid38)
setnames(PBC, c("#chrom","pos"), c("CHR38","BP38"))
data$PBC=PBC[,.(pid, CHR38, BP38, REF=ref, ALT=alt, BETA=beta, SE=sebeta, P=pval)]

RA=fread("../data/fg_sumstats/finngen_R7_M13_RHEUMA.gz", tmpdir = "tmp/")
RA$pid=paste(RA[["#chrom"]],RA$pos,sep=":")
table(RA$pid %in% snps$pid38)
setnames(RA, c("#chrom","pos"), c("CHR38","BP38"))
data$RA=RA[,.(pid, CHR38, BP38, REF=ref, ALT=alt, BETA=beta, SE=sebeta, P=pval)]

SLE=fread("../data/fg_sumstats/finngen_R7_M13_SLE.gz", tmpdir = "tmp/")
SLE$pid=paste(SLE[["#chrom"]],SLE$pos,sep=":")
table(SLE$pid %in% snps$pid38)
setnames(SLE, c("#chrom","pos"), c("CHR38","BP38"))
data$SLE=SLE[,.(pid, CHR38, BP38, REF=ref, ALT=alt, BETA=beta, SE=sebeta, P=pval)]

GPA=fread("../data/fg_sumstats/finngen_R7_M13_WEGENER.gz", tmpdir = "tmp/")
GPA$pid=paste(GPA[["#chrom"]],GPA$pos,sep=":")
table(GPA$pid %in% snps$pid38)
setnames(GPA, c("#chrom","pos"), c("CHR38","BP38"))
data$GPA=GPA[,.(pid, CHR38, BP38, REF=ref, ALT=alt, BETA=beta, SE=sebeta, P=pval)]


# Now we'll get the sample sizes, using the metadata file
md  <- fread("../data/Metadata_20230906-v1.tsv")
fd <- data.table(File_ID = files, name = names(files))
fd[, File_ID:=gsub("~/rds/rds-cew54-basis/02-Processed/", "", File_ID, fixed = TRUE)]
fd <- merge(fd, md[,.(File_ID, N0, N1)], by="File_ID")

n1=fd$N1
names(n1) <- fd$name
n0=fd$N0
names(n0) <- fd$name


## run coloc

d2l=function(d, trait) {
    list(snp=d$pid,
         beta=d$BETA,
         varbeta=d$SE^2,
         p=d$P, # Not necessary for coloc, but we'll use p-values to call novelty
         type="cc",
         s=n1[trait]/(n1[trait]+n0[trait]))
    }
w=1e+6 # large window choice
# This code will call the best SNP and find their P-values in both the myositis and the IMD dataset, so we don't need to run coloc twice.
# 1:nrow(index)
for(i in 1:nrow(index)) {
    message("Applying coloc on ", index$trait.myos[i], " and ", index$trait.other[i], " at SNP ", index$chr[i], ":", index$bp[i], ".")
    st=index$bp[i]-w
    en=index$bp[i]+w
    chr=index$chr[i]
    d1=d2l( data[[ index$trait.myos[i]]][CHR38==chr & BP38>st & BP38<en & !is.na(SE) & !is.na(BETA) & !duplicated(pid)],
           trait=index$trait.myos[[i]])
    d2=d2l( data[[ index$trait.other[i] ]][CHR38==chr & BP38>st & BP38<en & !is.na(SE) & !is.na(BETA) & !duplicated(pid)],
           trait=index$trait.other[[i]])
    
    result=coloc.abf(d1,d2, p12 = 5e-6) # Note: we modified the prior to be more conservative
    index[i ,c("nsnps","H0","H1","H2","H3","H4"):=as.list(result$summary)]
    best=result$results$snp[ which.max(result$results$SNP.PP.H4) ]
    w1=which(d1$snp==best)
    w2=which(d2$snp==best)
    index[i ,c("bestsnp","bestsnp.pp","pbest.myos","pbest.other", "pbest.myos.region", "pbest.other.region"):=
                 list(best, max(result$results$SNP.PP.H4),
                      d1$p[w1],
                      d2$p[w2],
                      min(d1$p),
                      min(d2$p))]
}

index[H4>.5]


# Add info on which PCs the SNPs are driver for
index <- merge(index, dv, by="pid", all.x=TRUE)

# Add P-value for the driver SNPs
p.myos <- unique(index[ , .(pid, trait.myos)], all.x=TRUE)
pex <- data2[trait %in% unique(p.myos$trait.myos) & pid %in% unique(p.myos$pid), .(pid, trait, P)]
p.myos <- merge(p.myos, pex, by.x=c("pid", "trait.myos"), by.y=c("pid", "trait"))
names(p.myos)[3]  <- "pdriver.myos"
index <- merge(index, p.myos, by=c("pid", "trait.myos"), all.x=TRUE)


index[ H4>.5 , .(trait.myos, trait.other,  fdr.myos, fdr.other, pairwise_fdr, H4, pid, bestsnp, bestsnp.pp, pdriver.myos)]

# Save results
fwrite(index, "../data/coloc_results-v3.tsv", sep="\t") 


####
# Some checks
####

ccs <- fread("../data/coloc_results.tsv")  %>% .[ H4 > 0.5]
ccs2 <- fread("../data/coloc_results-v3.tsv") %>% .[ H4 > 0.5]

comp <- paste(ccs$pid, ccs$trait.myos, ccs$trait.other, sep = "_")
comp2 <- paste(ccs2$pid, ccs2$trait.myos, ccs2$trait.other, sep = "_")

comp %in% comp2 %>% table
# FALSE   TRUE 
#     15    37 
# 15 lost (28.8%)
comp2 %in% comp %>% table
# FALSE  TRUE 
#    5    37 
# 5 new 

cs <- unique(ccs$pid)
cs2 <- unique(ccs2$pid)

cs %in% cs2 %>% table
# FALSE  TRUE 
#     6     5 
# 6 lost in new 
cs[!cs %in% cs2]
# "10:6064589"   "12:110972733" "12:112468611" "2:100215693"  "2:190670850"  "6:167124106"
# 10:6064589 corresponds to rs7073236, associated with IL2RA in the main table
# 12:110972733 corresponds to rs991817. We excluded it from the final main table, as its FDR was 0.65.
# 12:112468611 corresponds to rs11066320, associated with SH2B3 in the main table. A novel finding, but weak. 
# 2:100215693 
# 2:190670850
# 6:167124106
cs2 %in% cs %>% table
# FALSE  TRUE 
#     3     5 
# 3 new hits!
cs2[!cs2 %in% cs]
# "10:6064303" "22:21585386"  "2:198029415"
# 10:6064303 is rs7072793. Close to IL2RA. Coloc with RA (H4 = 0.57, pwFDR = 0.02), but associated with T1D in OpenTargets.
# 22:21585386 is rs5754217, an intron variant of UBE2L3. Coloc with SSc (H4 = 0.82, pwFDR = 0.005). Associated with some blood traits in OTG.
# 2:198029415 is rs10196612, an intron variant of PLCL1. Coloc with SjS (H4 = 0.64, pwFDR = 0.004). Associated with hay fever/rhinitis in OTG.

###########


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
#  [1] viridis_0.6.5       utf8_1.2.4          generics_0.1.3      lattice_0.22-6      grid_4.3.3          R.oo_1.26.0         plyr_1.8.9          jsonlite_1.8.8     
#  [9] R.utils_2.12.3      reshape_0.8.9       mixsqp_0.3-54       gridExtra_2.3       fansi_1.0.6         viridisLite_0.4.2   scales_1.3.0        textshaping_0.3.7  
# [17] cli_3.6.2           rlang_1.1.3         crayon_1.5.2        R.methodsS3_1.8.2   munsell_0.5.1       splines_4.3.3       susieR_0.12.35      withr_3.0.0        
# [25] tools_4.3.3         dplyr_1.1.4         colorspace_2.1-0    BiocGenerics_0.48.1 vctrs_0.6.5         R6_2.5.1            matrixStats_1.3.0   lifecycle_1.0.4    
# [33] zlibbioc_1.48.0     ragg_1.3.0          irlba_2.3.5.1       pkgconfig_2.0.3     pillar_1.9.0        gtable_0.3.4        glue_1.7.0          Rcpp_1.0.12        
# [41] systemfonts_1.0.6   tibble_3.2.1        tidyselect_1.2.1    farver_2.1.1        labeling_0.4.3      compiler_4.3.3      cupcake_0.1.0.0  