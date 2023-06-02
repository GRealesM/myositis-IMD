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


dpres <- fread("../data/DPMUnc_res_v1.tsv") 
dpres[DPMUnc.cl == 1 | grepl("jia|felty|scleroderma|dermatopoly", Trait, ignore.case = TRUE), .(Trait,Label)] 
#                                             Trait                                    Label
#  1:                  AITD_Saevarsdottir_32581359_1              Autoimmune Thyroid disease
#  2: AUTOIMMUNE_HYPERTHYROIDISM_FinnGen_FinnGenR7_1              Autoimmune hyperthyroidism
#  3:         E4_THYROIDITAUTOIM_FinnGen_FinnGenR7_1                  Autoimmune thyroiditis
#  4:               CHIRBIL_PRIM_FinnGen_FinnGenR7_1              Biliary chirrosis, primary
#  5:       RX_RHEUMA_BIOLOGICAL_FinnGen_FinnGenR7_1        Biological medication for rheuma
#  6:                      CREST_FinnGen_FinnGenR7_1                        CR(E)ST syndrome
#  7:                          DMY_Miller_26291516_1                Dermatomyositis (Miller)
#  8:                              DMY_Rothwell_up_1              Dermatomyositis (Rothwell)
#  9:            M13_DERMATOPOLY_FinnGen_FinnGenR7_1           Dermatopolymyositis (FinnGen)
# 10:                        MYGEO_Renton_25643325_1            Early-onset Myastenia Gravis
# 11:                      FELTY_FinnGen_FinnGenR7_1                          Felty syndrome
# 12:                 20002_1225_PanUKBB_PanUKBBR2_1          Hyperthyroidism/Thyrotoxicosis
# 13:             HYPOTHYROIDISM_FinnGen_FinnGenR7_1 Hypothyroidism (congenital or acquired)
# 14:                     D3_ITP_FinnGen_FinnGenR7_1     Idiopathic thrombocytopenic purpura
# 15:                     NMOIGGp_Estrada_29769526_1               IgG+ Neuromyelitis Optica
# 16:                             JO1M_Rothwell_up_1                Jo1+ Myositis (Rothwell)
# 17:                          JDM_Miller_26291516_1       Juvenile Dermatomyositis (Miller)
# 18:                              JDM_Rothwell_up_1     Juvenile Dermatomyositis (Rothwell)
# 19:                       JIA_LopezIsac_33106285_1           Juvenile Idiopathic Arthritis
# 20:               M13_JUVERHEU_FinnGen_FinnGenR7_1                         Juvenile rheuma
# 21:                        MYGLO_Renton_25643325_1             Late-onset Myastenia Gravis
# 22:       L12_LOCALSCLERODERMA_FinnGen_FinnGenR7_1         Localized scleroderma [morphea]
# 23:                               AAVMPO_Wong_up_1                                MPO+ AAV
# 24:                            MYG_Chia_35074870_1                        Myastenia Gravis
# 25:                          MYO_Miller_26291516_1                       Myositis (Miller)
# 26:                              IIM_Rothwell_up_1                     Myositis (Rothwell)
# 27:                               AAVPR3_Wong_up_1                                PR3+ AAV
# 28:            M13_PALINDROMIC_FinnGen_FinnGenR7_1                  Palindromic rheumatism
# 29:                M13_POLYMYO_FinnGen_FinnGenR7_1                  Polymyositis (FinnGen)
# 30:                           PM_Miller_26291516_1                   Polymyositis (Miller)
# 31:                               PM_Rothwell_up_1                 Polymyositis (Rothwell)
# 32:                              PAPS_Casares_up_1       Primary antiphospholipid syndrome
# 33:                              SJOS_Lessard_up_1                      Sjögren's syndrome
# 34:                       SSC_LopezIsac_31672989_1                      Systemic Sclerosis
# 35:                M13_WEGENER_FinnGen_FinnGenR7_1                  Wegener granulomatosis
#                                              Trait                                   Label

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
                paps =   "PAPS_Casares_up_1-hg38.tsv.gz",
                sjos =   "SJOS_Lessard_up_1-hg38.tsv.gz",
                ssc  = "SSC_LopezIsac_31672989_1-hg38.tsv.gz",
                wegen = "M13_WEGENER_FinnGen_FinnGenR7_1-hg38.tsv.gz",
                mpoaav = "AAVMPO_Wong_up_1-hg38.tsv.gz", 
                pr3aav =  "AAVPR3_Wong_up_1-hg38.tsv.gz",
                hyperthy = "AUTOIMMUNE_HYPERTHYROIDISM_FinnGen_FinnGenR7_1-hg38.tsv.gz",
                hypothy = "HYPOTHYROIDISM_FinnGen_FinnGenR7_1-hg38.tsv.gz") %>%
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

# There are 24 datasets, so instead of plotting everything in the same plot, we'll split it into two
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
mhp3 <- ggplot(data2[ trait %in% traits[17:24]], aes(x=x, y=-log10(fdr), col=factor(as.numeric(CHR38) %% 2))) + 
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

ggsave("../figures/manhattan_24IMD.png", mhp, height = 8, width = 17, bg = "white")


myos=data2[trait %in% traits[1:9], .(trait, pid, fdr)] # I put myositis datasets at the beginning, so these should be the ones.
not=data2[ !trait %in% traits[1:9], .(trait,pid,fdr)]
myos=merge(myos, not, by=c("pid"), suffixes=c(".myos",".other"), allow.cartesian = TRUE)
myos[,pairwise_fdr:=1 - (1-fdr.myos) * (1-fdr.other)] # P(H0 for either disease) = 1 - P(H1 for PAPS & other | p1, p2). Small means high prob of sharing
summary(myos$pairwise_fdr)
myos[ pairwise_fdr < 0.05 ]

#              pid trait.myos     fdr.myos trait.other    fdr.other pairwise_fdr
#  1:   10:6064589      myo.r 1.799960e-02         jia 8.090470e-04 1.879409e-02
#  2:   10:6064589      myo.r 1.799960e-02     hypothy 2.875370e-04 1.828196e-02
#  3: 11:118871133      myo.r 3.414525e-02        sjos 1.450420e-04 3.428534e-02
#  4: 11:118871133      myo.r 3.414525e-02         ssc 1.299060e-05 3.415780e-02
#  5:  11:35245397      myo.r 1.170098e-02     hypothy 6.148507e-06 1.170706e-02
#  6:  11:64329761      myo.r 5.992934e-03         jia 1.908621e-02 2.496476e-02
#  7:  11:64329761      myo.r 5.992934e-03        sjos 2.443384e-02 3.028034e-02
#  8:  11:64329761      myo.r 5.992934e-03     hypothy 9.713330e-04 6.958446e-03
#  9:  11:64362250      myo.r 8.432189e-03        sjos 3.515450e-03 1.191800e-02
# 10:  11:64362250      myo.r 8.432189e-03         ssc 2.458050e-02 3.280542e-02
# 11:  11:64362250      myo.r 8.432189e-03     hypothy 8.776556e-03 1.713474e-02
# 12:  1:113834946      myo.r 2.196624e-05         pbc 4.542891e-02 4.544988e-02
# 13:  1:113834946      myo.r 2.196624e-05         jia 5.734167e-11 2.196630e-05
# 14:  1:113834946      myo.r 2.196624e-05        myag 1.841060e-08 2.198465e-05
# 15:  1:113834946      myo.r 2.196624e-05      mpoaav 3.569593e-04 3.789177e-04
# 16:  1:113834946      myo.r 2.196624e-05    hyperthy 5.367583e-12 2.196624e-05
# 17:  1:113834946      myo.r 2.196624e-05     hypothy 6.540492e-84 2.196624e-05
# 18:  1:113834946       pm.m 6.539977e-03         jia 5.734167e-11 6.539977e-03
# 19:  1:113834946       pm.m 6.539977e-03        myag 1.841060e-08 6.539995e-03
# 20:  1:113834946       pm.m 6.539977e-03      mpoaav 3.569593e-04 6.894602e-03
# 21:  1:113834946       pm.m 6.539977e-03    hyperthy 5.367583e-12 6.539977e-03
# 22:  1:113834946       pm.m 6.539977e-03     hypothy 6.540492e-84 6.539977e-03
# 23:  1:113834946       pm.r 8.073321e-05         pbc 4.542891e-02 4.550597e-02
# 24:  1:113834946       pm.r 8.073321e-05         jia 5.734167e-11 8.073327e-05
# 25:  1:113834946       pm.r 8.073321e-05        myag 1.841060e-08 8.075162e-05
# 26:  1:113834946       pm.r 8.073321e-05      mpoaav 3.569593e-04 4.376637e-04
# 27:  1:113834946       pm.r 8.073321e-05    hyperthy 5.367583e-12 8.073322e-05
# 28:  1:113834946       pm.r 8.073321e-05     hypothy 6.540492e-84 8.073321e-05
# 29:  2:190670850       pm.r 2.299677e-05         jia 4.701598e-02 4.703789e-02
# 30:  2:190670850       pm.r 2.299677e-05        sjos 3.745633e-04 3.975515e-04
# 31:  2:190670850       pm.r 2.299677e-05         ssc 3.065160e-06 2.606186e-05
# 32:  2:190670850       pm.r 2.299677e-05    hyperthy 4.022179e-02 4.024386e-02
# 33:  2:190670850       pm.r 2.299677e-05     hypothy 1.787724e-02 1.789983e-02
# 34:  2:191071078      myo.r 7.728768e-04         pbc 2.374151e-02 2.449604e-02
# 35:  2:191071078      myo.r 7.728768e-04         jia 1.188817e-03 1.960775e-03
# 36:  2:191071078      myo.r 7.728768e-04        paps 8.867650e-04 1.658956e-03
# 37:  2:191071078      myo.r 7.728768e-04        sjos 2.319724e-10 7.728770e-04
# 38:  2:191071078      myo.r 7.728768e-04         ssc 6.452880e-11 7.728769e-04
# 39:  2:191071078      myo.r 7.728768e-04    hyperthy 5.344977e-04 1.306961e-03
# 40:  2:191071078      myo.r 7.728768e-04     hypothy 2.181396e-17 7.728768e-04
# 41:   3:28029953      myo.r 1.799960e-02        sjos 9.040184e-04 1.888735e-02
# 42:   3:28029953      myo.r 1.799960e-02         ssc 1.214773e-02 2.992867e-02
# 43:  4:122194347      myo.r 3.182470e-02         jia 1.116514e-02 4.263451e-02
# 44:  4:122194347      myo.r 3.182470e-02      pr3aav 4.777696e-04 3.228726e-02
# 45:  5:157185077      dmy.r 3.005199e-02     hypothy 4.411007e-03 3.433044e-02
# 46:  5:157185077      myo.r 1.198567e-02    hyperthy 2.129210e-02 3.302257e-02
# 47:  5:157185077      myo.r 1.198567e-02     hypothy 4.411007e-03 1.634380e-02
# 48:  7:128933913      myo.m 3.245979e-02         jia 1.169413e-02 4.377433e-02
# 49:  7:128933913      myo.m 3.245979e-02        sjos 1.394114e-31 3.245979e-02
# 50:  7:128933913      myo.m 3.245979e-02         ssc 1.413720e-06 3.246116e-02
# 51:  7:128933913      myo.m 3.245979e-02     hypothy 7.132875e-03 3.936114e-02
# 52:  7:128933913      myo.r 4.241146e-03         jia 1.169413e-02 1.588568e-02
# 53:  7:128933913      myo.r 4.241146e-03        sjos 1.394114e-31 4.241146e-03
# 54:  7:128933913      myo.r 4.241146e-03         ssc 1.413720e-06 4.242554e-03
# 55:  7:128933913      myo.r 4.241146e-03      mpoaav 1.814119e-02 2.230540e-02
# 56:  7:128933913      myo.r 4.241146e-03     hypothy 7.132875e-03 1.134377e-02
# 57:  7:128954129     jo1m.r 4.905755e-02        sjos 8.469175e-32 4.905755e-02
# 58:  7:128954129     jo1m.r 4.905755e-02         ssc 1.157436e-17 4.905755e-02
# 59:  7:128954129     jo1m.r 4.905755e-02      mpoaav 2.540366e-04 4.929912e-02
# 60:  7:128954129      myo.r 3.022458e-03         jia 1.467235e-02 1.765046e-02
# 61:  7:128954129      myo.r 3.022458e-03        paps 1.559745e-02 1.857277e-02
# 62:  7:128954129      myo.r 3.022458e-03        sjos 8.469175e-32 3.022458e-03
# 63:  7:128954129      myo.r 3.022458e-03         ssc 1.157436e-17 3.022458e-03
# 64:  7:128954129      myo.r 3.022458e-03      mpoaav 2.540366e-04 3.275727e-03
# 65:  7:128954129      myo.r 3.022458e-03     hypothy 4.651927e-02 4.940113e-02
# 66:  7:128977412     jo1m.r 4.905755e-02      mpoaav 2.540366e-04 4.929912e-02
# 67:  7:128977412      myo.r 3.022458e-03         jia 1.372170e-02 1.670269e-02
# 68:  7:128977412      myo.r 3.022458e-03        paps 1.991110e-03 5.007550e-03
# 69:  7:128977412      myo.r 3.022458e-03      mpoaav 2.540366e-04 3.275727e-03
# 70:   8:11491677      myo.r 4.117055e-03        paps 1.888645e-03 5.997925e-03
# 71:   8:11491677      myo.r 4.117055e-03        sjos 2.136693e-05 4.138334e-03
# 72:   8:11491677      myo.r 4.117055e-03         ssc 6.452880e-11 4.117055e-03
#              pid trait.myos     fdr.myos trait.other    fdr.other pairwise_fdr


################################################################################

## We'll now select a broader set of SNPs to apply coloc to (Pairwise FDR < 0.5)

index=myos[ pairwise_fdr < 0.5 ][order(pairwise_fdr)]
index[ , c("chr","bp"):=tstrsplit(pid,":")  %>% lapply(., as.numeric) ]

# Remove SNPs in proximity (less than 1Mb apart)
f=function(d) {
    if(nrow(d)==1)
        return(d)
    dist=abs(outer(d$bp, d$bp, "-"))
    dist[ upper.tri(dist) ] = Inf
    diag(dist)=Inf
    drop=apply(dist, 1, min) < 1e+6
    return(d[ !drop ])
}

## Now apply this to dense SNP data
## can safely ignore warnings about no non-missing arguments to min
index=split(index, index[,.(trait.myos, trait.other ,chr)])  %>% lapply(., f)  %>% rbindlist()
nrow(index) ## 477


## itp, pbc, crest, felty, wegen, hyperthy, and hypoty are not dense datasets. 
names(data)[c(10,12,14,16,20,23:24)]= c("itp.local", "pbc.local", "crest.local",   "felty.local", "wegen.local", "hyperthy.local",  "hypothy.local")

dir.create("../data")
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

# Now we'll get the sample sizes, using the metadata file
md  <- fread("../data/Metadata_20230503-v1.tsv")
fd <- data.table(File_ID = files, name = names(files))
fd[, File_ID:=gsub("~/rds/rds-cew54-basis/02-Processed/", "", File_ID, fixed = TRUE)]
fd <- merge(fd, md[,.(File_ID, N0, N1)], by="File_ID")

n1=fd$N1
names(n1) <- fd$name
n0=fd$N0
names(n0) <- fd$name


## run 477 colocs

d2l=function(d, trait) {
    list(snp=d$pid,
         beta=d$BETA,
         varbeta=d$SE^2,
         type="cc",
         s=n1[trait]/(n1[trait]+n0[trait]))
    }
w=1e+6 # large window choice
for(i in 1:nrow(index)) {
    st=index$bp[i]-w
    en=index$bp[i]+w
    chr=index$chr[i]
    d1=d2l( data[[ index$trait.myos[i]]][CHR38==chr & BP38>st & BP38<en & !is.na(SE) & !is.na(BETA) & !duplicated(pid)],
           trait=index$trait.myos[[i]])
    d2=d2l( data[[ index$trait.other[i] ]][CHR38==chr & BP38>st & BP38<en & !is.na(SE) & !is.na(BETA) & !duplicated(pid)],
           trait=index$trait.other[[i]])
    index[i ,c("nsnps","H0","H1","H2","H3","H4"):=as.list(coloc.abf(d1,d2)$summary)]
}

index[pairwise_fdr < 0.05 | H4>.5]



# We're we interested on the best SNP for H4 >.5 
for(i in which(index$H4>.5)) {
    st=index$bp[i]-w
    en=index$bp[i]+w
    chr=index$chr[i]
    d1=d2l( data[[ index$trait.myos[i]]][CHR38==chr & BP38>st & BP38<en & !is.na(SE) & !is.na(BETA) & !duplicated(pid)],
           trait= index$trait.myos[i])
    d2=d2l( data[[ index$trait.other[i] ]][CHR38==chr & BP38>st & BP38<en & !is.na(SE) & !is.na(BETA) & !duplicated(pid)],
           trait=index$trait.other[i])
    result=coloc.abf(d1,d2)
    best=result$results$snp[ which.max(result$results$SNP.PP.H4) ]
    w1=which(d1$snp==best)
    w2=which(d2$snp==best)
    index[i ,c("bestsnp","bestsnp.pp","pbest.myos","pbest.other"):=
                 list(best, max(result$results$SNP.PP.H4),
                      2*pnorm( -abs(d1$beta)/sqrt(d1$varbeta))[w1],
                      2*pnorm( -abs(d2$beta)/sqrt(d2$varbeta))[w2])]
}

# Add info on which PCs the SNPs are driver for
index <- merge(index, dv, by="pid")

# Add P-value for the driver SNPs
p.myos <- unique(index[ , .(pid, trait.myos)])
pex <- data2[trait %in% unique(p.myos$trait.myos) & pid %in% unique(p.myos$pid), .(pid, trait, P)]
p.myos <- merge(p.myos, pex, by.x=c("pid", "trait.myos"), by.y=c("pid", "trait"))
names(p.myos)[3]  <- "pdriver.myos"
index <- 

merge(index, p.myos, by=c("pid", "trait.myos"))


index[ H4>.5 , .(trait.myos, trait.other,  fdr.myos, fdr.other, pairwise_fdr, H4, pid, bestsnp, bestsnp.pp, pdriver.myos)]

# Save results
fwrite(index, "../tables/coloc_results.tsv", sep="\t")

# Now explore the SNPs in more detail
index[ H4>.5, unique(pid)]
# 15 SNPs, comprising 12 regions
index[ H4>.5, unique(trait.other)]
# We have coloc associations with 12 (out of 15) IMDs
# These are
#  [1] "hypothy"  "sjos"     "jia"      "ssc"      "hyperthy" "myag"    
#  [7] "felty"    "mpoaav"   "pbc"      "pr3aav"   "paps"    

## plot these signals
plotter=function(pid,w=1e+6) {
    chr=sub(":.*","",pid)  %>% as.numeric()
    bp=sub(".*:","",pid)  %>% as.numeric()
    st=bp-w
    en=bp+w
    wh=which(index$pid==pid & (index$H4>.5|index$pairwise_fdr<.05))
    print(index[wh])
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
            ggtitle(pid)+
            theme_cowplot() + 
            theme(axis.title.x = element_blank(),
                  strip.background = element_rect(colour="black", fill = "white"),
                  )
}



plots <- lapply(index[ H4>.5, unique(pid)], plotter)
names(plots) <- index[ H4>.5, unique(pid)]

index[ H4>.5, unique(pid)]

p11.1  <- plots$`11:64329761`
p11.2 <- plots$`11:64362250`
p11 <- plot_grid(p11.1,p11.2, ncol =1, labels = NULL)
ggsave("../figures/coloc_chr11.png", p11, height = 8, width = 8, bg="white")

ggsave("../figures/coloc_chr12.png", plots$`12:112468611`, height = 8, width = 8, bg="white")
ggsave("../figures/coloc_chr17_1.png", plots$`17:39913696`, height = 8, width = 8, bg="white")
ggsave("../figures/coloc_chr17_2.png", plots$`17:75373341`, height = 8, width = 8, bg="white")
ggsave("../figures/coloc_chr2_1.png", plots$`2:100215693`, height = 8, width = 8, bg="white")
ggsave("../figures/coloc_chr2_2.p7:128954129
7:128977412ng", plots$`2:191071078`, height = 8, width = 8, bg="white")
ggsave("../figures/coloc_chr3.png", plots$`3:28029953`, height = 8, width = 8, bg="white")
ggsave("../figures/coloc_chr4.png", plots$`4:122194347`, height = 8, width = 8, bg="white")

p7.1  <- plots$`7:128933913`
p7.2  <- plots$`7:128954129`
p7.3  <- plots$`7:128977412`
p7 <- plot_grid(p7.1,p7.2,p7.3, ncol =1, labels = NULL)
ggsave("../figures/coloc_chr7_1.png", p7, height = 16, width = 8, bg="white")

ggsave("../figures/coloc_chr7_2.png", plots$`7:37397251`, height = 8, width = 8, bg="white")
ggsave("../figures/coloc_chr8.png", plots$`8:11491677`, height = 8, width = 8, bg="white")


## paps, crest, sjorgrens, jia, pbc - already genomewide sig in paps alone - STAT4 (intronic), STAT1
## https://genetics.opentargets.org/Variant/2_191079016_C_T/associations
## paps, sjorgrens, ssc - novel for paps - FAM167A, BLK
## https://genetics.opentargets.org/Variant/8_11486464_C_T/associations
## paps, ssc - novel for paps - TNPO3 (intronic), IRF5
## https://genetics.opentargets.org/Variant/7_129011468_C_T/associations


data$paps[ pid %in% index[H4>.5]$bestsnp ]$SNPID






