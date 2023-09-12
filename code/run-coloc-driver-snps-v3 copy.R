# RUNNING COLOC

# This script will run the coloc step.

# 20230904 - Changes vs previous versions: We replaced some basis traits datasets to avoid overfitting. This changed the selection of IMD 
# for coloc. Namely, we removed hypothy, MG, ITP and others, and included EOMG, LOMG, SLE, and RA

## Load packages and required
library(data.table)
library(magrittr)
library(coloc)
library(annotSnpStats)
library(ggplot2)
library(cowplot)
setDTthreads(15)
setwd("/home/gr440/rds/rds-cew54-basis/Projects/myositis-IMD/code")


dpres <- fread("../data/DPMUnc_res_v3.tsv") 
dpres[DPMUnc == 1, .(Trait,Label)] 
# Note: JDM (R) was allocated to a different cluster by DPMUnc, so it doesn't appear in this list. 
#                                        Trait                             Label
#  1: RX_RHEUMA_BIOLOGICAL_FinnGen_FinnGenR7_1  Biological medication for rheuma
#  2:                CREST_FinnGen_FinnGenR7_1                  CR(E)ST syndrome
#  3:                    DMY_Miller_26291516_1          Dermatomyositis (Miller)
#  4:                        DMY_Rothwell_up_1        Dermatomyositis (Rothwell)
#  5:                  MYGEO_Renton_25643325_1      Early-onset Myastenia Gravis
#  6:                FELTY_FinnGen_FinnGenR7_1                    Felty syndrome
#  7:           20002_1225_PanUKBB_PanUKBBR2_1    Hyperthyroidism/Thyrotoxicosis
#  8:                    MYO_Miller_26291516_1                      IIM (Miller)
#  9:                        IIM_Rothwell_up_1                    IIM (Rothwell)
# 10:               NMOIGGp_Estrada_29769526_1         IgG+ Neuromyelitis Optica
# 11:                       JO1M_Rothwell_up_1          Jo1+ Myositis (Rothwell)
# 12:                    JDM_Miller_26291516_1 Juvenile Dermatomyositis (Miller)
# 13:                 JIA_LopezIsac_33106285_1     Juvenile Idiopathic Arthritis
# 14:                  MYGLO_Renton_25643325_1       Late-onset Myastenia Gravis
# 15:                         AAVMPO_Wong_up_1                          MPO+ AAV
# 16:                         AAVPR3_Wong_up_1                          PR3+ AAV
# 17:      M13_PALINDROMIC_FinnGen_FinnGenR7_1            Palindromic rheumatism
# 18:                     PM_Miller_26291516_1             Polymyositis (Miller)
# 19:                         PM_Rothwell_up_1           Polymyositis (Rothwell)
# 20:         CHIRBIL_PRIM_FinnGen_FinnGenR7_1       Primary biliary cholangitis
# 21:           M13_RHEUMA_FinnGen_FinnGenR7_1              Rheumatoid arthritis
# 22:                        SJOS_Lessard_up_1                Sjögren's syndrome
# 23:                 SSC_LopezIsac_31672989_1                Systemic Sclerosis
# 24:              M13_SLE_FinnGen_FinnGenR7_1      Systemic lupus erythematosus
# 25:          M13_WEGENER_FinnGen_FinnGenR7_1            Wegener granulomatosis
#                                        Trait                             Label


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

                `IgG+ NMO` =  "NMOIGGp_Estrada_29769526_1-hg38.tsv.gz", 
                PBC = "CHIRBIL_PRIM_FinnGen_FinnGenR7_1-hg38.tsv.gz", # 11
                JIA = "JIA_LopezIsac_33106285_1-hg38.tsv.gz",
                `CR(E)ST` = "CREST_FinnGen_FinnGenR7_1-hg38.tsv.gz", # 13
                EOMG = "MYGEO_Renton_25643325_1-hg38.tsv.gz",
                LOMG = "MYGLO_Renton_25643325_1-hg38.tsv.gz",
                Felty =  "FELTY_FinnGen_FinnGenR7_1-hg38.tsv.gz", # 16
                SjS =   "SJOS_Lessard_up_1-hg38.tsv.gz",
                SSc  = "SSC_LopezIsac_31672989_1-hg38.tsv.gz",
                GPA = "M13_WEGENER_FinnGen_FinnGenR7_1-hg38.tsv.gz", #19
                `MPO+ AAV` = "AAVMPO_Wong_up_1-hg38.tsv.gz", 
                `PR3+ AAV` =  "AAVPR3_Wong_up_1-hg38.tsv.gz",
                HyperThy = "20002_1225_PanUKBB_PanUKBBR2_1-hg38.tsv.gz", #22
                PR = "M13_PALINDROMIC_FinnGen_FinnGenR7_1-hg38.tsv.gz", #23
                BioMedRhe = "RX_RHEUMA_BIOLOGICAL_FinnGen_FinnGenR7_1-hg38.tsv.gz", #24
                RA = "M13_RHEUMA_FinnGen_FinnGenR7_1-hg38.tsv.gz", # 25 
                SLE = "M13_SLE_FinnGen_FinnGenR7_1-hg38.tsv.gz") %>% #26
    lapply(.,  function(f) (file.path("~/rds/rds-cew54-basis/02-Processed",f)))
stopifnot(all(file.exists(unlist(files))))

man=fread("~/rds/rds-cew54-basis/03-Bases/IMD_basis/SNP.manifest.38.tsv")
snps=fread("~/rds/rds-cew54-basis/03-Bases/IMD_basis/Manifest_build_translator.tsv")
snps[, pid38:=paste(CHR38, BP38, sep=":")][, pid19:=paste(CHR19, BP19, sep=":")]
snps <- merge(snps, man, by="pid38")

data=lapply(files, fread)


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
snps.use=snps$pid38[ rowSums(rot)>0 ] # This works because SNPs are in the same order, even though they belong to different builds.
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

ggsave("../figures/manhattan_26IMD_v3.png", mhp, height = 8, width = 17, bg = "white")


myos=data2[trait %in% traits[1:9], .(trait, pid, fdr)] # I put myositis datasets at the beginning, so these should be the ones.
not=data2[ !trait %in% traits[1:9], .(trait,pid,fdr)]
myos=merge(myos, not, by=c("pid"), suffixes=c(".myos",".other"), allow.cartesian = TRUE)
myos[,pairwise_fdr:=1 - (1-fdr.myos) * (1-fdr.other)] # P(H0 for either disease) = 1 - P(H1 for PAPS & other | p1, p2). Small means high prob of sharing
summary(myos$pairwise_fdr)
myos[ pairwise_fdr < 0.05 ]

#              pid    trait.myos     fdr.myos trait.other    fdr.other  pairwise_fdr
#  1:   10:6064589       IIM (R) 1.799960e-02         JIA 8.090470e-04  1.879409e-02
#  2:   10:6064589       IIM (R) 1.799960e-02    HyperThy 3.620361e-04  1.835512e-02
#  3:   10:6064589       IIM (R) 1.799960e-02   BioMedRhe 8.117894e-03  2.597138e-02
#  4:   10:6064589       IIM (R) 1.799960e-02          RA 3.154588e-06  1.800270e-02
#  5: 11:118871133       IIM (R) 3.414525e-02         SjS 1.450420e-04  3.428534e-02
#  6: 11:118871133       IIM (R) 3.414525e-02         SSc 1.299060e-05  3.415780e-02
#  7: 11:118871133       IIM (R) 3.414525e-02          RA 1.027001e-02  4.406459e-02
#  8:  11:64329761       IIM (R) 5.992934e-03         JIA 1.908621e-02  2.496476e-02
#  9:  11:64329761       IIM (R) 5.992934e-03         SjS 2.443384e-02  3.028034e-02
# 10:  11:64362250       IIM (R) 8.432189e-03         SjS 3.515450e-03  1.191800e-02
# 11:  11:64362250       IIM (R) 8.432189e-03         SSc 2.458050e-02  3.280542e-02
# 12:  1:113834946       IIM (R) 2.196624e-05         PBC 4.542891e-02  4.544988e-02
# 13:  1:113834946       IIM (R) 2.196624e-05         JIA 5.734167e-11  2.196630e-05
# 14:  1:113834946       IIM (R) 2.196624e-05        LOMG 3.320395e-02  3.322519e-02
# 15:  1:113834946       IIM (R) 2.196624e-05    MPO+ AAV 3.569593e-04  3.789177e-04
# 16:  1:113834946       IIM (R) 2.196624e-05    HyperThy 4.349570e-15  2.196624e-05
# 17:  1:113834946       IIM (R) 2.196624e-05   BioMedRhe 4.953324e-13  2.196624e-05
# 18:  1:113834946       IIM (R) 2.196624e-05          RA 5.062804e-70  2.196624e-05
# 19:  1:113834946        PM (M) 6.539977e-03         JIA 5.734167e-11  6.539977e-03
# 20:  1:113834946        PM (M) 6.539977e-03        LOMG 3.320395e-02  3.952678e-02
# 21:  1:113834946        PM (M) 6.539977e-03    MPO+ AAV 3.569593e-04  6.894602e-03
# 22:  1:113834946        PM (M) 6.539977e-03    HyperThy 4.349570e-15  6.539977e-03
# 23:  1:113834946        PM (M) 6.539977e-03   BioMedRhe 4.953324e-13  6.539977e-03
# 24:  1:113834946        PM (M) 6.539977e-03          RA 5.062804e-70  6.539977e-03
# 25:  1:113834946        PM (R) 8.073321e-05         PBC 4.542891e-02  4.550597e-02
# 26:  1:113834946        PM (R) 8.073321e-05         JIA 5.734167e-11  8.073327e-05
# 27:  1:113834946        PM (R) 8.073321e-05        LOMG 3.320395e-02  3.328200e-02
# 28:  1:113834946        PM (R) 8.073321e-05    MPO+ AAV 3.569593e-04  4.376637e-04
# 29:  1:113834946        PM (R) 8.073321e-05    HyperThy 4.349570e-15  8.073321e-05
# 30:  1:113834946        PM (R) 8.073321e-05   BioMedRhe 4.953324e-13  8.073321e-05
# 31:  1:113834946        PM (R) 8.073321e-05          RA 5.062804e-70  8.073321e-05
# 32:  2:190670850        PM (R) 2.299677e-05         JIA 4.701598e-02  4.703789e-02
# 33:  2:190670850        PM (R) 2.299677e-05         SjS 3.745633e-04  3.975515e-04
# 34:  2:190670850        PM (R) 2.299677e-05         SSc 3.065160e-06  2.606186e-05
# 35:  2:191071078       IIM (R) 7.728768e-04         PBC 2.374151e-02  2.449604e-02
# 36:  2:191071078       IIM (R) 7.728768e-04         JIA 1.188817e-03  1.960775e-03
# 37:  2:191071078       IIM (R) 7.728768e-04         SjS 2.319724e-10  7.728770e-04
# 38:  2:191071078       IIM (R) 7.728768e-04         SSc 6.452880e-11  7.728769e-04
# 39:  2:191071078       IIM (R) 7.728768e-04          RA 4.829483e-06  7.777026e-04
# 40:  2:191071078       IIM (R) 7.728768e-04         SLE 5.241505e-05  8.252513e-04
# 41:   3:28029953       IIM (R) 1.799960e-02         SjS 9.040184e-04  1.888735e-02
# 42:   3:28029953       IIM (R) 1.799960e-02         SSc 1.214773e-02  2.992867e-02
# 43:  4:122194347       IIM (R) 3.182470e-02         JIA 1.116514e-02  4.263451e-02
# 44:  4:122194347       IIM (R) 3.182470e-02    PR3+ AAV 4.777696e-04  3.228726e-02
# 45:  5:157185077        DM (R) 3.005199e-02    HyperThy 2.674321e-04  3.031138e-02
# 46:  5:157185077       IIM (R) 1.198567e-02    HyperThy 2.674321e-04  1.224989e-02
# 47:  7:128933913       IIM (M) 3.245979e-02         JIA 1.169413e-02  4.377433e-02
# 48:  7:128933913       IIM (M) 3.245979e-02         SjS 1.394114e-31  3.245979e-02
# 49:  7:128933913       IIM (M) 3.245979e-02         SSc 1.413720e-06  3.246116e-02
# 50:  7:128933913       IIM (M) 3.245979e-02          RA 1.761528e-04  3.263023e-02
# 51:  7:128933913       IIM (M) 3.245979e-02         SLE 3.850379e-07  3.246016e-02
# 52:  7:128933913       IIM (R) 4.241146e-03         JIA 1.169413e-02  1.588568e-02
# 53:  7:128933913       IIM (R) 4.241146e-03         SjS 1.394114e-31  4.241146e-03
# 54:  7:128933913       IIM (R) 4.241146e-03         SSc 1.413720e-06  4.242554e-03
# 55:  7:128933913       IIM (R) 4.241146e-03    MPO+ AAV 1.814119e-02  2.230540e-02
# 56:  7:128933913       IIM (R) 4.241146e-03   BioMedRhe 3.148800e-02  3.559560e-02
# 57:  7:128933913       IIM (R) 4.241146e-03          RA 1.761528e-04  4.416552e-03
# 58:  7:128933913       IIM (R) 4.241146e-03         SLE 3.850379e-07  4.241529e-03
# 59:  7:128954129 Anti-Jo1+ (R) 4.905755e-02         SjS 8.469175e-32  4.905755e-02
# 60:  7:128954129 Anti-Jo1+ (R) 4.905755e-02         SSc 1.157436e-17  4.905755e-02
# 61:  7:128954129 Anti-Jo1+ (R) 4.905755e-02    MPO+ AAV 2.540366e-04  4.929912e-02
# 62:  7:128954129 Anti-Jo1+ (R) 4.905755e-02          RA 1.550458e-05  4.907229e-02
# 63:  7:128954129 Anti-Jo1+ (R) 4.905755e-02         SLE 2.622553e-08  4.905757e-02
# 64:  7:128954129       IIM (R) 3.022458e-03         JIA 1.467235e-02  1.765046e-02
# 65:  7:128954129       IIM (R) 3.022458e-03         SjS 8.469175e-32  3.022458e-03
# 66:  7:128954129       IIM (R) 3.022458e-03         SSc 1.157436e-17  3.022458e-03
# 67:  7:128954129       IIM (R) 3.022458e-03    MPO+ AAV 2.540366e-04  3.275727e-03
# 68:  7:128954129       IIM (R) 3.022458e-03          RA 1.550458e-05  3.037916e-03
# 69:  7:128954129       IIM (R) 3.022458e-03         SLE 2.622553e-08  3.022484e-03
# 70:  7:128977412 Anti-Jo1+ (R) 4.905755e-02    MPO+ AAV 2.540366e-04  4.929912e-02
# 71:  7:128977412 Anti-Jo1+ (R) 4.905755e-02          RA 1.550458e-05  4.907229e-02
# 72:  7:128977412 Anti-Jo1+ (R) 4.905755e-02         SLE 2.622553e-08  4.905757e-02
# 73:  7:128977412       IIM (R) 3.022458e-03         JIA 1.372170e-02  1.670269e-02
# 74:  7:128977412       IIM (R) 3.022458e-03    MPO+ AAV 2.540366e-04  3.275727e-03
# 75:  7:128977412       IIM (R) 3.022458e-03          RA 1.550458e-05  3.037916e-03
# 76:  7:128977412       IIM (R) 3.022458e-03         SLE 2.622553e-08  3.022484e-03
# 77:   8:11491677       IIM (R) 4.117055e-03         SjS 2.136693e-05  4.138334e-03
# 78:   8:11491677       IIM (R) 4.117055e-03         SSc 6.452880e-11  4.117055e-03
# 79:   8:11491677       IIM (R) 4.117055e-03          RA 1.511672e-04  4.267600e-03
# 80:   8:11491677       IIM (R) 4.117055e-03         SLE 6.161174e-05  4.178413e-03
#              pid    trait.myos     fdr.myos trait.other    fdr.other 

# 80 pairs with pairwise_fdr < 0.05

################################################################################

myos[ pairwise_fdr < 0.05 , unique(pid)]
# 14 SNPs, 10 independent genomic region

fwrite(myos, "../data/raw_fdr_results_v3.tsv", sep="\t")

# Next, to increase our chances of discovery, we'll relax our pairwise_fdr limits to 0.5.

# We want to ensure that we can compare index SNPs and trait.other, even if they're not pairwise_fdr < 0.5 for all myositis,
# so we'll select those SNPs with at pairwise_fdr < 0.5 and then coloc those for all combinations.
myos[, trait_snp:=paste0(trait.other, "_", pid)] # auxiliary variable for trait.other and index_snp pairs

# This way, in index_tspairs, we'll capture IMD-SNP pairs that have at least one pairwise_fdr with a myositis trait.
# Then, by selecting those IMD-SNP combinations in index_tspairs in index, we'll capture
# all myositis-IMD combinations at a given SNP where at least one myositis-IMD pair with pairwise_fdr < 0.5 exists.
index_tspairs <- myos[ pairwise_fdr < 0.5 , unique(trait_snp) ]
length(index_tspairs)
# 292 unique IMD-indexSNPs with pairwise_fdr < 0.5 pairs

index=myos[ pairwise_fdr < 0.5 | trait_snp %in% index_tspairs ][order(pairwise_fdr)] 

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
# 2517
length(unique(index$pid))
# 81 SNPs

# Remove SNPs to drop
index <- index[!pid %in% drop]

nrow(index)
# 1817
length(unique(index$pid))
# 62

## Next step, bring dense SNP datasets to run coloc

## itp, pbc, crest, felty, wegen, hyperthy, are not dense datasets. 
names(data)[c(11,13,16,19,22,23,24,25,26)] = c("PBC.local", "CREST.local",  "Felty.local", "GPA.local", "HyperThy.local",  "PR.local", "BioMedRhe.local", "RA.local", "SLE.local")

dir.create("../data/fg_sumstats")
if(!file.exists("../data/fg_sumstats/finngen_R7_CHIRBIL_PRIM.gz"))
    system("wget https://storage.googleapis.com/finngen-public-data-r7/summary_stats/finngen_R7_CHIRBIL_PRIM.gz -O ../data/fg_sumstats/finngen_R7_CHIRBIL_PRIM.gz")
if(!file.exists("../data/fg_sumstats/finngen_R7_CREST.gz"))
    system("wget https://storage.googleapis.com/finngen-public-data-r7/summary_stats/finngen_R7_CREST.gz -O ../data/fg_sumstats/finngen_R7_CREST.gz")
if(!file.exists("../data/fg_sumstats/finngen_R7_FELTY.gz"))
    system("wget https://storage.googleapis.com/finngen-public-data-r7/summary_stats/finngen_R7_FELTY.gz -O ../data/fg_sumstats/finngen_R7_FELTY.gz")
if(!file.exists("../data/fg_sumstats/finngen_R7_M13_WEGENER.gz"))
    system("wget https://storage.googleapis.com/finngen-public-data-r7/summary_stats/finngen_R7_M13_WEGENER.gz -O ../data/fg_sumstats/finngen_R7_M13_WEGENER.gz")

# Note: this PanUKBB file will need some work to reformat, see below
if(!file.exists("../data/fg_sumstats/20002_1225_PanUKBB_1-hg38.tsv.gz")){
    system("wget https://pan-ukb-us-east-1.s3.amazonaws.com/sumstats_flat_files/categorical-20002-both_sexes-1225.tsv.bgz -O ../data/fg_sumstats/20002_1225_PanUKBB_PanUKBBR2_1.bgz")
    system("Rscript processing_panUKBB.R") # This script will prepare the PanUKBB file to be used by coloc
}   

if(!file.exists("../data/fg_sumstats/finngen_R7_HYPOTHYROIDISM.gz"))
    system("wget https://storage.googleapis.com/finngen-public-data-r7/summary_stats/finngen_R7_HYPOTHYROIDISM.gz -O ../data/fg_sumstats/finngen_R7_HYPOTHYROIDISM.gz")
if(!file.exists("../data/fg_sumstats/finngen_R7_M13_PALINDROMIC.gz"))
    system("wget https://storage.googleapis.com/finngen-public-data-r7/summary_stats/finngen_R7_M13_PALINDROMIC.gz -O ../data/fg_sumstats/finngen_R7_M13_PALINDROMIC.gz")
if(!file.exists("../data/fg_sumstats/finngen_R7_RX_RHEUMA_BIOLOGICAL.gz"))
    system("wget https://storage.googleapis.com/finngen-public-data-r7/summary_stats/finngen_R7_RX_RHEUMA_BIOLOGICAL.gz -O ../data/fg_sumstats/finngen_R7_RX_RHEUMA_BIOLOGICAL.gz")
if(!file.exists("../data/fg_sumstats/finngen_R7_M13_RHEUMA.gz"))
    system("wget https://storage.googleapis.com/finngen-public-data-r7/summary_stats/finngen_R7_M13_RHEUMA.gz -O ../data/fg_sumstats/finngen_R7_M13_RHEUMA.gz")
if(!file.exists("../data/fg_sumstats/finngen_R7_M13_SLE.gz"))
    system("wget https://storage.googleapis.com/finngen-public-data-r7/summary_stats/finngen_R7_M13_SLE.gz -O ../data/fg_sumstats/finngen_R7_M13_SLE.gz")

# Incorporate new files into data

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

hyperthy=fread("../data/fg_sumstats/20002_1225_PanUKBB_1-hg38.tsv.gz")
hyperthy[, pid := paste(CHR38,BP38,sep=":")]
table(hyperthy$pid %in% snps$pid38)
data$HyperThy=hyperthy[,.(pid, CHR38, BP38, REF, ALT, BETA, SE, P)]

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

RA=fread("../data/fg_sumstats/finngen_R7_M13_RHEUMA.gz")
RA$pid=paste(RA[["#chrom"]],RA$pos,sep=":")
table(RA$pid %in% snps$pid38)
setnames(RA, c("#chrom","pos"), c("CHR38","BP38"))
data$RA=RA[,.(pid, CHR38, BP38, REF=ref, ALT=alt, BETA=beta, SE=sebeta, P=pval)]

SLE=fread("../data/fg_sumstats/finngen_R7_M13_SLE.gz")
SLE$pid=paste(SLE[["#chrom"]],SLE$pos,sep=":")
table(SLE$pid %in% snps$pid38)
setnames(SLE, c("#chrom","pos"), c("CHR38","BP38"))
data$SLE=SLE[,.(pid, CHR38, BP38, REF=ref, ALT=alt, BETA=beta, SE=sebeta, P=pval)]


# Now we'll get the sample sizes, using the metadata file
md  <- fread("../data/Metadata_20230503-v1.tsv")
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
fwrite(index, "../data/coloc_results_dfilt-v3.tsv", sep="\t") 

# Compared to earlier versions, this v3 should include:
# * RA, SLE, EOMG, and LOMG. ITP, MG, and HypoThy removed
# * Distance-filtered driver SNPs
# * Original p-values, rather than our calculated ones.

###########


sessionInfo()
# R version 4.1.3 (2022-03-10)
# Platform: x86_64-pc-linux-gnu (64-bit)
# Running under: Rocky Linux 8.8 (Green Obsidian)

# Matrix products: default
# BLAS/LAPACK: /usr/local/software/spack/spack-rhel8-20210927/opt/spack/linux-centos8-icelake/gcc-11.2.0/intel-oneapi-mkl-2021.4.0-s2cksi33smowj5zlqvmew37cufvztdkc/mkl/2021.4.0/lib/intel64/libmkl_gf_lp64.so.1

# locale:
#  [1] LC_CTYPE=en_GB.UTF-8       LC_NUMERIC=C              
#  [3] LC_TIME=en_GB.UTF-8        LC_COLLATE=en_GB.UTF-8    
#  [5] LC_MONETARY=en_GB.UTF-8    LC_MESSAGES=en_GB.UTF-8   
#  [7] LC_PAPER=en_GB.UTF-8       LC_NAME=C                 
#  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
# [11] LC_MEASUREMENT=en_GB.UTF-8 LC_IDENTIFICATION=C       

# attached base packages:
# [1] stats     graphics  grDevices utils     datasets  methods   base     

# other attached packages:
# [1] cowplot_1.1.1     ggplot2_3.4.1     annotSnpStats_1.1 snpStats_1.40.0  
# [5] survival_3.5-3    Matrix_1.5-3      coloc_5.1.0.1     magrittr_2.0.3   
# [9] data.table_1.14.8

# loaded via a namespace (and not attached):
#  [1] Rcpp_1.0.10         pillar_1.8.1        compiler_4.1.3     
#  [4] susieR_0.12.35      plyr_1.8.8          tools_4.1.3        
#  [7] R.methodsS3_1.8.2   R.utils_2.12.2      viridis_0.6.2      
# [10] zlibbioc_1.36.0     jsonlite_1.8.4      lifecycle_1.0.3    
# [13] tibble_3.1.8        gtable_0.3.1        lattice_0.20-45    
# [16] viridisLite_0.4.1   pkgconfig_2.0.3     rlang_1.0.6        
# [19] cli_3.6.0           parallel_4.1.3      gridExtra_2.3      
# [22] withr_2.5.0         dplyr_1.1.0         systemfonts_1.0.4  
# [25] generics_0.1.3      vctrs_0.5.2         grid_4.1.3         
# [28] tidyselect_1.2.0    cupcake_0.1.0.0     reshape_0.8.9      
# [31] glue_1.6.2          R6_2.5.1            textshaping_0.3.6  
# [34] fansi_1.0.4         mixsqp_0.3-48       irlba_2.3.5.1      
# [37] farver_2.1.1        BiocGenerics_0.36.1 scales_1.2.1       
# [40] matrixStats_0.63.0  splines_4.1.3       colorspace_2.1-0   
# [43] ragg_1.2.5          labeling_0.4.2      utf8_1.2.3         
# [46] munsell_0.5.0       crayon_1.5.2        R.oo_1.25.0   