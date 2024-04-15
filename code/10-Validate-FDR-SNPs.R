#########################################
##                                     ##
##     VALIDATING SNP SELECTION        ##
##                                     ##
#########################################

# Author: Guillermo Reales 
# Date last updated: 2024/04/09

# Background: We wanted to validate our method to narrow down SNPs. To do so, we will select pairs of IMD from large GWAS and similarly sized as IIM GWASs.
# We'll 

# This script will
# * Project selected R5 IMD
# * Apply FDR together with filtered projections
# * Select FDR significant R5/pso and check which PCs they're significant for
# * Compute pairwise Bhattacharyya distance for each R5/pso trait with the set of selected IMD used with myositis
# * Select genetically close IMD to each of the R5/pso traits
# * Apply pairwise FDR and to find significant SNPs in the R5/pso traits
# * Check which SNPs replicate in the validation datasets (ie. R7/Tsoi)


##########################################


## Load packages and required data sets
library(data.table)
library(magrittr)
library(IMDtools)
library(cupcake)
library(annotSnpStats)
library(fpc)

setDTthreads(20)
setwd("/home/gr440/rds/rds-cew54-basis/Projects/myositis-IMD/code")
bpath <- "/home/gr440/rds/rds-cew54-basis/03-Bases/IMD_basis/"

###############################################

################################
##  Project selected R5 IMD  ###
################################


# The bits below are similar to 00-Prepare_datasets_EDA.R. Our goal is to make the projection datasets 
# similar to those we used for the main analysis, but including the list of desired traits for validation.
# These validation traits should 
# (1) NOT be in the original basis (eg. T1D, RA, MS, etc.)
# (2) Have two instances per trait, one smaller (comparable in size to myositis), and one larger, for validation.
# Based on this we chose to try the following IMD:
#  * Ankylosing spondylitis (FinnGen, PanUKBB,...).
#  * Biological medication for Rheuma (FinnGen)
#  * COPD (FinnGen)
#  * Gout (FinnGen)
#  * Hyperparathyroidism (FinnGen)
#  * Psoriasis (PanUKBB, Tsoi)

# For some of these, we'll use FinnGen R5/R7 for validation/discovery. We have the R7 projections already, but we need to reproject the R5 ones

# Meta-data table, containing information about the datasets
m <- fread("../data/Metadata_20230906-v1.tsv")
m <- m[, .(Trait, First_Author, Reference, Trait_ID_2.0, Trait_long, Trait_class, N0, N1, N, Population, Public)]
m <- m[ Trait_class == "IMD" & Public == "Y"]
m[grepl("Hyperparat", Trait_long) & grepl("FinnGen", First_Author)]

trp <- c("M13_ANKYLOSPON_FinnGen_FinnGenR5_1", "RX_RHEUMA_BIOLOGICAL_FinnGen_FinnGenR5_1", "J10_COPD_FinnGen_FinnGenR5_1", "GOUT_FinnGen_FinnGenR5_1", "E4_HYPERPARA_FinnGen_FinnGenR5_1")

# Now, we need to project these. The process will be analogous to that applied for the rest of projections

# Prepare manifest to filter
m38 <- fread("~/rds/rds-cew54-basis/GWAS_tools/04-Reduction/SNP.manifest.38.tsv")
m38[, c("CHR38", "BP38") := tstrsplit(pid38, split = ":")]
setnames(m38, c("ref_a1", "ref_a2"), c("REF", "ALT"))
m38 <- m38[, .(CHR38, BP38, REF, ALT)]

redf <- lapply(trp, function(x){
    input <- fread(paste0("~/rds/rds-cew54-basis/02-Processed/", x, "-hg38.tsv.gz"), tmpdir = "tmp")
    mincold <- c("CHR38", "BP38", "REF", "ALT", "BETA", "SE", "P")
    input <- input[, ..mincold]
    M <- g.align(input, m38)
    dups <- M[duplicated(pid), pid]
    if(length(dups) > 0) M <- M[!pid %in% dups] # Remove dups
    setnames(M, "pid", "pid38") # We need pid column to be named pid38
    return(M)
})
names(redf) <- trp

# Now for the projections

# Since files are in hg38 and the basis is in hg19, we use a dirty fix to liftover back to hg19 (easier to do at this stage than in the previous one, due to the huge size of files).
build_dict <- fread("~/rds/rds-cew54-basis/03-Bases/IMD_basis/Manifest_build_translator.tsv")
build_dict[, `:=`(pid38=paste(CHR38,BP38,sep=":"), pid=paste(CHR19,BP19,sep=":"))]
build_dict <- build_dict[,c("pid38", "pid")] 

# We load the traits that were used to create the basis first, then we'll append the rest of traits to that table
projected.basis <- copy(cupcake::basis.trait.proj)
projected.basis[, `:=` (Var.Delta=0, z=0, P=NA)]
setnames(projected.basis, c("delta", "trait"), c("Delta", "Trait"))

nfiles <- length(trp)
Trait <- rep(NA,nfiles)
nSNP <- rep(NA, nfiles)
overall_p <- rep(NA, nfiles)
mscomp <- rep(NA,nfiles)

projected.table <- lapply(trp, function(file){
	message("Projecting ", file)
    index <- which(trp == file)
	sm <- redf[[file]]
	# Some checks
	sm <- unique(sm)
	sm[sm == ""] <- NA # Some missing data might pass as empty string. This will fix that	
	sm <- na.omit(sm, cols = c("pid38", "BETA", "SE", "P"))
	dups <- sm$pid38[duplicated(sm$pid38)]

	if(length(dups) > 0){
		dupmessage= "This file has duplicated pids. I removed them prior to projection. You might want to check it."
		message(dupmessage)
		# Write to log
		write(paste0(ss.file,". ",dupmessage), logname, append=TRUE)
		sm <- sm[!pid38 %in% dups] # Remove all duplicated instances, to be safe
	}
	
	sm <- merge(sm, build_dict)

	# A bit of QC 
	Trait[index] <<- file
	nSNP[index] <<- nrow(sm)

  	projected.userdata <- tryCatch(expr = project_sparse(beta=sm$BETA, seb=sm$SE, pid=sm$pid)[,trait:=file][],
  	                      error =function(e) {
				      failmessage <- "Projection for this file had non-zero exit status, please check. Jumping to next file..."
				      message(failmessage)
				    #   write(paste0(ss.file,". ",failmessage), logname, append=TRUE)
				      return(NULL)})
  	if(is.null(projected.userdata)) {
		return(NULL)
	} else{
	projected.userdata[, proj:=NULL]
  	setnames(projected.userdata, c("var.proj", "delta", "p", "trait"), c("Var.Delta", "Delta", "P", "Trait")) # Changed var.proj to var.delta
  	setcolorder(projected.userdata, c( "Trait", "PC", "Delta", "Var.Delta", "p.overall", "z", "P"))

  	# More QC
  	overall_p[index] <<- projected.userdata$p.overall[1]
  	minc.idx <- which.min(projected.userdata$P)
  	mscomp[index] <<- sprintf("%s (%1.0e)",projected.userdata$PC[minc.idx],projected.userdata$P[minc.idx])
	}
	projected.userdata
}
) 

projected.table[sapply(projected.table, is.null)]  <- NULL 
projected.table <- rbindlist(projected.table)
projected.table[, Var.Delta:=as.numeric(Var.Delta)][, z:=NULL][, p.overall:=NULL]


# projected.full  <- rbind(projected.basis[,.(PC, Delta, Var.Delta, z, P, Trait)],projected.table[,.(PC, Delta, Var.Delta, z, P, Trait)])
QC.table <- data.table(Trait, nSNP, overall_p, mscomp)
QC.table <- merge(QC.table, m, by = "Trait")
QC.table[, Label:=gsub(" \\(FinnGen\\)", "", Trait_long)]

projected.table <- merge(projected.table, QC.table[,c("First_Author","Trait", "Trait_ID_2.0", "Trait_long","Trait_class", "Population", "Label")], by = "Trait")

fwrite(projected.table, "../data/Projection_FGR5_IMD_basis.tsv", sep = "\t")
fwrite(QC.table, "../data/QC_FGR5_IMD_basis.tsv", sep = "\t")

###############################################

# Put projections together

pr5 <- fread("../data/Projection_FGR5_IMD_basis.tsv")
qr5 <- fread("../data/QC_FGR5_IMD_basis.tsv")

pf <- fread("../data/pf.tsv")
qf <- fread("../data/qf.tsv")

pf <- rbind(pf,pr5)
qf <- rbind(qf,qr5)

# We've incorporated the new datasets. Now we can apply the FDR procedure as in 00-Prepar_datasets_EDA.R


############################################################################################################
##  Apply FDR to new projections together with old ones and check which PCs they're FDR significant for  ###
############################################################################################################


# Apply 1% FDR correction to overall p for all remaining datasets
qf[, FDR.overall := p.adjust(overall_p, method = "BH"), by="Trait_class"] # Only IMD, so trait class shouldn't matter

# Apply FDR.PC to projections
pf[, FDR.PC:=p.adjust(P, method = "BH"), by = c("PC", "Trait_class")][, stars:=ifelse(!is.na(FDR.PC) & FDR.PC<0.05,"○","")][ FDR.PC < 0.01 , stars:="●"]

# Focus on our traits of interest

trp2 <- c("M13_ANKYLOSPON_FinnGen_FinnGenR7_1", "RX_RHEUMA_BIOLOGICAL_FinnGen_FinnGenR7_1", "J10_COPD_FinnGen_FinnGenR7_1", "GOUT_FinnGen_FinnGenR7_1", "E4_HYPERPARA_FinnGen_FinnGenR7_1", "ph696.4_PanUKBB_PanUKBBR2_1")
toi <- c(trp,trp2)

qf[Trait %in% toi, .(Trait, Trait_long, N0, N1, FDR.overall)][order(Trait_long, Trait)]
#                                        Trait                                   Trait_long     N0    N1  FDR.overall
#                                       <char>                                       <char>  <int> <int>        <num>
#  1:       M13_ANKYLOSPON_FinnGen_FinnGenR5_1             Ankylosing spondylitis (FinnGen) 164682  1462 2.744777e-05
#  2:       M13_ANKYLOSPON_FinnGen_FinnGenR7_1             Ankylosing spondylitis (FinnGen) 227388  2252 3.055808e-12
#  3: RX_RHEUMA_BIOLOGICAL_FinnGen_FinnGenR5_1   Biological medication for rheuma (FinnGen) 216990  1802 1.278198e-18
#  4: RX_RHEUMA_BIOLOGICAL_FinnGen_FinnGenR7_1   Biological medication for rheuma (FinnGen) 306015  3139 1.117447e-36
#  5:             J10_COPD_FinnGen_FinnGenR5_1                               COPD (FinnGen) 186723  6915 5.909220e-05
#  6:             J10_COPD_FinnGen_FinnGenR7_1                               COPD (FinnGen) 259128 13255 5.800457e-12
#  7:                 GOUT_FinnGen_FinnGenR5_1                                         Gout 203546  3576 7.897608e-03
#  8:                 GOUT_FinnGen_FinnGenR7_1                                         Gout 284844  6031 3.064421e-04
#  9:         E4_HYPERPARA_FinnGen_FinnGenR5_1                Hyperparathyroidism (FinnGen) 211123  2928 3.440704e-01
# 10:         E4_HYPERPARA_FinnGen_FinnGenR7_1                Hyperparathyroidism (FinnGen) 297590  4386 2.657683e-02
# 11:              ph696.4_PanUKBB_PanUKBBR2_1                             Psoriasis (UKBB) 415468  2981 1.015866e-04

# Only hyperparathyroidism isn't significant. We'll use the rest of R5 + Psoriasis datasets

ttid <- qf[grepl("FinnGenR5|ph696.4_", Trait) & FDR.overall < 0.01, Trait] # List of focus IMD


## Check which PC is significant for each focus IMD
sigpcs <- lapply(ttid, function(i){
	    p.i <- pf[ Trait == i & FDR.PC < 0.01, PC]
		p.i
})
names(sigpcs) <- ttid

######################################################################################################################
### Compute pairwise Bhattacharyya distance for each R5/pso trait with the set of selected IMD used with myositis  ###
######################################################################################################################

# Load selected IMD set
ps <- fread("../data/ps2.tsv")
qs <- fread("../data/qs2.tsv") # From here, we'll take the list of IMD used with myositis. These are significant proj with redundant removed.

imdset <- qs[!grepl("Miller|Rothwell", Trait)][!Trait %in% trp2, Trait]
# Remove Myositis datasets and validation R7 datasets from the IMD pool
# Note: ph696.4_PanUKBB_PanUKBBR2_1 is a focus (not validation) dataset and it's in trp2, but not in qs, so it has no effect to use trp2 to filter/

# We modified the following function wrt 01-compute-bhatta... to
# * Take an R5/Pso dataset
# * Look up the number of PCs to use, significant at FDR 1% for each R5/pso dataset.
# * Generate the covar matrices from reduced_datasets/ but take the R5 from redf instead
# * Use Trait instead of Label, to facilitate things
# * Compute Bhattacharyya from our focus IMD to the rest only.

compute_bhatta <- function(tt){

		pcs <- sigpcs[[tt]]

		fs <- c(tt, imdset)


        covar.m <- lapply(fs, function(x){

                if(x %in% names(redf)){
					fss <- redf[[x]]
				} else{
					fss <- fread(paste0(bpath, "reduced_datasets/", x, "-ft.tsv"))
				}
				
                # Some checks
                fss <- unique(fss)
                fss[fss == ""] <- NA # Some missing data might pass as empty string. This will fix that	
                fss <- na.omit(fss, cols = c("pid38", "BETA", "SE", "P"))
                dups <- fss$pid38[duplicated(fss$pid38)]

                if(length(dups) > 0){
                    dupmessage= "This file has duplicated pids. I removed them prior to projection. You might want to check it."
                    message(dupmessage)
                    fss <- fss[!pid38 %in% dups] # Remove all duplicated instances, to be safe
                }
                fss <- merge(fss, build_dict)

                pids <- fss$pid
                seb <- fss$SE

                v <- seb * shrinkage[pids] * rot.pca[pids,]
                var.proj  <- t(v) %*% LD[pids,pids] %*% v

                var.proj <- var.proj[pcs, pcs]
                var.proj
        })

        dd <- data.table(T1 = tt, T2 = fs) # Create a data.table to store combinations. Remove own hits

		# Data for focus trait (ie. T1)
		rr1  <- pf[ Trait == tt & PC %in% pcs, Delta] # We can use pf to look for the projections
		names(rr1) <- pf[ Trait == tt & PC %in% pcs, PC]

		cv1 <- as.matrix(covar.m[[1]]) # We put the focus trait the first in fs, so this is fine
		
		# Then, take each trait in turn, extract necessary data, and compute the Bhattacharyya distance
        bhd <- sapply(seq_along(fs), function(x){
                rr2 <- pf[ Trait == fs[x] & PC %in% pcs, Delta]
                names(rr2) <- pf[ Trait == fs[x] & PC %in% pcs, PC] 
                cv2 <- as.matrix(covar.m[[x]]) # Since we used fs to create covar.m, they should be in the same order
                bd <- bhattacharyya.dist(mu1 = rr1, mu2 = rr2, Sigma1 = cv1, Sigma2 = cv2 )
                bd
                })

		dd[, bhat.dist:=bhd]
		dd <- dd[order(bhat.dist)]
		dd <- dd[-1] # Remove self-distances

}

bhdist <- lapply(ttid, compute_bhatta)
bhsel <- lapply(bhdist, function(x){
		x  <- x[1:17] # Take the closest 17 IMD, as we did with Myositis
}) 
bhsel %<>% rbindlist
unique(bhsel$T2)
# 44 IMD to consider for FDR

###########################################################################
## Apply pairwise FDR and to find significant SNPs in the R5/pso traits ###
###########################################################################

# We have the significant PCs and the list of genetically close IMDs.
# Now we'll extract the driver SNPs for each of these PCs

# Rotation matrix is in hg19, but SNPs are in hg38. We'll need to translate
man=fread("~/rds/rds-cew54-basis/03-Bases/IMD_basis/SNP.manifest.38.tsv")
snps=fread("~/rds/rds-cew54-basis/03-Bases/IMD_basis/Manifest_build_translator.tsv")
snps[, pid38:=paste(CHR38, BP38, sep=":")][, pid19:=paste(CHR19, BP19, sep=":")]
snps <- merge(snps, man, by="pid38")

# Update SNPs in rotation matrix
rot=cupcake::rot.pca
all(rownames(rot) == snps$pid19) # Check the SNPs are in the same order -- they aren't, let's reorder them to match
snps <- snps[order(pid19)]
all(rownames(rot) == snps$pid19) # now they are
rownames(rot) <- snps$pid38

# Sig driver SNPs
sigdriver <- lapply(sigpcs, function(i){
		if(length(i) == 1){
			rpc <- rot[, i]
			snps.use <- names(rpc[ rpc != 0])
		} else{
			rpc <- rot[, i]
			snps.use <- rownames(rpc[rowSums(rpc != 0) > 0, ])
		}
		
})
names(sigdriver) <- ttid

# Load reduced files for all genetically close IMD

redimd.l <- c("ph696.4_PanUKBB_PanUKBBR2_1",unique(bhsel$T2)) # psoriasis missing, so we add it here

redimd <- lapply(redimd.l, function(x){
		i <- fread(paste0(bpath, "reduced_datasets/", x, "-ft.tsv"))
		i <- i[order(CHR38, BP38)]
})
names(redimd) <- redimd.l

redf %<>% lapply(. , function(x){ x[order(CHR38,BP38)]})

redall <- c(redf, redimd)

# Next, take each focus dataset in turn, and filter it and its close IMD by their driver SNPs.

pw.fdr <- lapply(ttid, function(x){
		# Driver SNPs
		dsnp <- sigdriver[[x]]
		# Genetically close IMD
		cimd <- bhsel[T1 == x, T2]

		d1 <- redall[[x]][ pid38 %in% dsnp]
		d2 <- lapply(cimd, function(i){
				y <- redall[[i]][ pid38 %in% dsnp]
		})
		d <- c(list(d1),d2)
		nm <- c(x, cimd)
		names(d) <- nm
		for(xs in nm)
			d[[xs]][, IMD:=xs]

		d %<>% rbindlist(., fill=TRUE)
		d[,P2:=2 * pnorm( -abs(BETA)/SE )] # Compute new P from BETA and SE
		d[,fdr:=p.adjust(P2, method="BH"), by="IMD"]

		f1 <- d[IMD == x, .(IMD, pid38, fdr)]
		f2 <- d[IMD != x, .(IMD, pid38, fdr)]
		f <- merge(f1, f2, by=c("pid38"), suffixes=c(".focus",".other"), allow.cartesian = TRUE)
		f[,pairwise_fdr:=1 - (1-fdr.focus) * (1-fdr.other)]
		f <- f[order(pairwise_fdr)]
})
pw.fdr.sig <- lapply(pw.fdr, function(x){
			y <- x[ pairwise_fdr < 0.05]
			y
})
pw.snps <- lapply(pw.fdr.sig, function(x){
			y <- x[, unique(pid38)]
			y
})


# Now for the validation datasets!

val.l <- c("PSO_Tsoi_23143594_1","J10_COPD_FinnGen_FinnGenR7_1", "RX_RHEUMA_BIOLOGICAL_FinnGen_FinnGenR7_1")

redval <- lapply(val.l, function(x){
		i <- fread(paste0(bpath, "reduced_datasets/", x, "-ft.tsv"))
		i <- i[order(CHR38, BP38)]
})
names(redval) <- val.l

# Psoriasis
t1 <- redall$ph696.4_PanUKBB_PanUKBBR2_1[ pid38 %in% pw.snps[[1]]][,P:=2 * pnorm( -abs(BETA)/SE )] # Update P, as UKBB comes as log10
v1 <- redval$PSO_Tsoi_23143594_1[ pid38 %in% pw.snps[[1]]]
tv1 <- merge(t1[, .(pid38, P)], v1[, .(pid38, SNPID, P)], suffixes=c(".focus",".val"), by = "pid38")
tv1[, .(pid38, SNPID, P.focus, P.val)]

# COPD
t2 <- redall$J10_COPD_FinnGen_FinnGenR5_1[ pid38 %in% pw.snps[[3]]]
v2 <- redval$J10_COPD_FinnGen_FinnGenR7_1[ pid38 %in% pw.snps[[3]]]
tv2 <- merge(t2[, .(pid38, P)], v2[, .(pid38, P)], suffixes=c(".focus",".val"), by = "pid38")


# Bio med rheum
t3 <- redall$RX_RHEUMA_BIOLOGICAL_FinnGen_FinnGenR5_1[ pid38 %in% pw.snps[[5]]]
v3 <- redval$RX_RHEUMA_BIOLOGICAL_FinnGen_FinnGenR7_1[ pid38 %in% pw.snps[[5]]]
tv3 <- merge(t3[, .(pid38, P)], v3[, .(pid38, P)], suffixes=c(".focus",".val"), by = "pid38")

cat("Psoriasis.\n")
tv1
cat("COPD.\n")
tv2
cat("Bio Med Rheum.\n")
tv3

# Try validation with FinnGen R10

# Whilst it isn't ideal to validate using FinnGen on UKBB, we'll use one of the psoriasis datasets
f10url <- c("https://storage.googleapis.com/finngen-public-data-r10/summary_stats/finngen_R10_L12_PSORIASIS.gz",
			"https://storage.googleapis.com/finngen-public-data-r10/summary_stats/finngen_R10_J10_COPD.gz",
			"https://storage.googleapis.com/finngen-public-data-r10/summary_stats/finngen_R10_RX_RHEUMA_BIOLOGICAL.gz")

d10 <- lapply(f10url, fread, tmpdir = "tmp")
d10 <- lapply(d10, function(x){
			x[ , pid38:=paste(`#chrom`, pos, sep = ":")]
})

# Psoriasis (R10_L12_PSORIASIS, Cases = 10312, Controls = 397564)
v1.10 <- d10[[1]][ pid38 %in% pw.snps[[1]]]
setnames(v1.10, "pval", "P")
tv1.10 <- merge(t1[, .(pid38, P)], v1.10[, .(pid38, P)], suffixes=c(".focus",".val"), by = "pid38")
tv1.10[, .(pid38, P.focus, P.val)]

#          pid38      P.focus       P.val
#         <char>        <num>       <num>
# 1: 19:10467167 6.529498e-07 1.74201e-03
# 2: 1:197642552 7.541812e-04 9.16579e-06
# 3: 1:197662011 8.388510e-04 7.84784e-06
# 4: 2:162254026 8.845999e-06 1.60081e-06
# 5: 5:159412209 5.314136e-07 1.86939e-17
# 6: 6:159044945 1.189902e-04 3.13437e-02

# COPD (R10_J10_COPD, Cases = 20066, Controls = 338303)
v2.10 <- d10[[2]][ pid38 %in% pw.snps[[3]]]
setnames(v2.10, "pval", "P")
tv2.10 <- merge(t2[, .(pid38, P)], v3.10[, .(pid38, P)], suffixes=c(".focus",".val"), by = "pid38")
tv2.10[, .(pid38, P.focus, P.val)]

#          pid38   P.focus       P.val
#         <char>     <num>       <num>
# 1: 5:132435113 0.0003623 1.88969e-06
# 2: 5:132448701 0.0004056 2.38144e-06
# 3: 5:132461230 0.0001080 7.59889e-06

# Bio Med Rheum (R10_RX_RHEUMA_BIOLOGICAL, Cases = 5117, Controls = 407064)
v3.10 <- d10[[3]][ pid38 %in% pw.snps[[5]]]
setnames(v3.10, "pval", "P")
tv3.10 <- merge(t3[, .(pid38, P)], v3.10[, .(pid38, P)], suffixes=c(".focus",".val"), by = "pid38")
tv3.10[, .(pid38, P.focus, P.val)]

#          pid38   P.focus       P.val
#         <char>     <num>       <num>
# 1:  10:6064303 1.719e-05 1.22360e-04
# 2:  10:6064589 1.743e-05 1.22823e-04
# 3:  10:6066476 3.849e-07 2.83949e-04
# 4: 1:113834946 3.640e-11 3.94912e-21






# qf2 <- qf[Trait %in% toi]
# pf2 <- pf[Trait %in% toi]

# qf2[grepl("FinnGenR5", Trait), Label:=paste0(Label, " (R5)")][grepl("FinnGenR7", Trait), Label:=paste0(Label, " (R7)")]
# pf2[grepl("FinnGenR5", Trait), Label:=paste0(Label, " (R5)")][grepl("FinnGenR7", Trait), Label:=paste0(Label, " (R7)")]

# qf2 <- qf2[order(Label)]
# pf2 <- pf2[order(Label, PC)]

# qs <- qf2[FDR.overall < 0.01]
# ps <- pf2[Trait %in% qs$Trait]


# # We'll now extract the significant PCs and their associated driver SNPs for each "training" trait (ie. R5 + Psoriasis)
# ttid <- qs[grepl("FinnGenR5|ph696.4", Trait), Trait]

# # Rotation matrix is in hg19, but SNPs are in hg38. We'll need to translate
# man=fread("~/rds/rds-cew54-basis/03-Bases/IMD_basis/SNP.manifest.38.tsv")
# snps=fread("~/rds/rds-cew54-basis/03-Bases/IMD_basis/Manifest_build_translator.tsv")
# snps[, pid38:=paste(CHR38, BP38, sep=":")][, pid19:=paste(CHR19, BP19, sep=":")]
# snps <- merge(snps, man, by="pid38")

# # Update SNPs in rotation matrix
# rot=cupcake::rot.pca
# all(rownames(rot) == snps$pid19) # Check the SNPs are in the same order -- they aren't, let's reorder them to match
# snps <- snps[order(pid19)]
# all(rownames(rot) == snps$pid19) # now they are
# rownames(rot) <- snps$pid38

# # Sig PCs
# sigpcs <- lapply(ttid, function(i){
# 	    p.i <- ps[ Trait == i & FDR.PC < 0.01, PC]
# 		p.i
# })
# names(sigpcs) <- ttid

# # Sig driver SNPs
# sigdriver <- lapply(sigpcs, function(i){
# 		if(length(i) == 1){
# 			rpc <- rot[, i]
# 			snps.use <- names(rpc[ rpc != 0])
# 		} else{
# 			rpc <- rot[, i]
# 			snps.use <- rownames(rpc[rowSums(rpc != 0) > 0, ])
# 		}
		
# })
# ttil <- qs[Trait %in% ttid, Label] %>% gsub(" \\(R5\\)", "", .)
# names(sigdriver) <- ttil

# sigdriver <- rbindlist(lapply(seq_along(sigdriver), function(i) {
#   data.table(IMD = names(sigdriver)[i], pid = sigdriver[[i]])
# }))


# # Load all datasets
# dsn <- c(qs$Trait, "PSO_Tsoi_23143594_1")
# dsp <- paste0("~/rds/rds-cew54-basis/02-Processed/", dsn, "-hg38.tsv.gz")

# data <- lapply(dsp, fread)

# aligner=function(d) {
#     d[,pid:=paste(CHR38,BP38,sep=":")]
#     ## table(d$pid %in% snps$pid38)
#     d=d[pid %in% snps$pid38]
#     message("snps found: ",nrow(d), " / ",nrow(snps))
#     geno=toupper(paste(d$REF,d$ALT,sep="/"))
#     cl=g.class(geno, snps$alleles[ match(d$pid, snps$pid38) ])
#     print(table(cl))
#     d[ cl %in% c("rev","revcomp"), BETA:=-BETA]
#     d=d[ (cl!="impossible") ]
#     d
# }

# data=lapply(data, aligner)
# names(data) <- dsn
# for(nm in dsn) 
#     data[[nm]]$Trait=nm

# imdv <- rep(ttil, each = 2)
# for(i in 1:10)
# 	data[[i]]$IMD=imdv[i]

# data <- rbindlist(data, fill = TRUE)
# data <- data[, .(CHR38,BP38, pid, SNPID, REF, ALT, BETA, SE, P, Trait, IMD)]

# data <- merge(data, sigdriver)
# data[,P2:=2 * pnorm( -abs(BETA)/SE )]
# data2[,fdr:=p.adjust(P2, method="BH"), by="trait"]





# p <- fread("../tables/ST_all_projections.tsv")
# q <- fread("../tables/ST_all_datasets.tsv")

# q[grepl("Ankylosing", Label) & sig.overall == "Yes" ] #
# # Ankylosing (N1 = 495 and 2252, respectively)
# # c("ph715.2_PanUKBB_PanUKBBR2_1","M13_ANKYLOSPON_FinnGen_FinnGenR7_1")
# # Bio Med for rheuma (No R5!)
# # COPD (FinnGen only, smaller COPD are a subtype)
# # Gout (Only FinnGen is sig overall. There are smaller Gout (PanUKBB, Koettgen), but not sig overall)
# # Hyperparathyroidism (Only R7, and not significant overall)
# # Psoriasis (N1 = 2981 and 6995, respectively). Tsoi was ImmunoChip and was therefore not included
# # c("ph696.4_PanUKBB_PanUKBBR2_1", "L12_PSORIASIS_FinnGen_FinnGenR7_1")

# q[sig.overall == "Yes"][ order(Label, N1)]
