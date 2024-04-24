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
library(coloc)

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
tv2.10 <- merge(t2[, .(pid38, P)], v2.10[, .(pid38, P)], suffixes=c(".focus",".val"), by = "pid38")
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


#######################################
###   Run coloc in pairwise pairs   ###
#######################################

pfc <- copy(pw.fdr.sig) %>% rbindlist
pfc <- pfc[, .(pid38, IMD.focus, IMD.other, pairwise_fdr)]
pfc[, c("CHR38", "BP38"):=tstrsplit(pid38, split = ":", fixed = TRUE)][, c("CHR38", "BP38"):=list(as.numeric(CHR38), as.numeric(BP38))]
pfc <- pfc[order(IMD.focus, CHR38, BP38)]
#           pid38                                IMD.focus                            IMD.other pairwise_fdr CHR38      BP38
#          <char>                                   <char>                               <char>        <num> <num>     <num>
#  1: 5:132435113             J10_COPD_FinnGen_FinnGenR5_1          K11_IBD_FinnGen_FinnGenR7_1 4.978316e-02     5 132435113
#  2: 5:132448701             J10_COPD_FinnGen_FinnGenR5_1          K11_IBD_FinnGen_FinnGenR7_1 4.978316e-02     5 132448701
#  3: 5:132461230             J10_COPD_FinnGen_FinnGenR5_1          K11_IBD_FinnGen_FinnGenR7_1 3.891900e-02     5 132461230
#  4: 5:132461230             J10_COPD_FinnGen_FinnGenR5_1        K11_ULCER_FinnGen_FinnGenR7_1 3.919840e-02     5 132461230
#  5: 5:132461230             J10_COPD_FinnGen_FinnGenR5_1 L12_PSORI_ARTHRO_FinnGen_FinnGenR7_1 4.405671e-02     5 132461230
#  6: 1:113834946 RX_RHEUMA_BIOLOGICAL_FinnGen_FinnGenR5_1       M13_RHEUMA_FinnGen_FinnGenR7_1 1.777463e-08     1 113834946
#  7: 1:113834946 RX_RHEUMA_BIOLOGICAL_FinnGen_FinnGenR5_1       20002_1225_PanUKBB_PanUKBBR2_1 1.777464e-08     1 113834946
#  8: 1:113834946 RX_RHEUMA_BIOLOGICAL_FinnGen_FinnGenR5_1                     AAVMPO_Wong_up_1 6.705413e-04     1 113834946
#  9:  10:6064303 RX_RHEUMA_BIOLOGICAL_FinnGen_FinnGenR5_1       M13_RHEUMA_FinnGen_FinnGenR7_1 2.069803e-03    10   6064303
# 10:  10:6064303 RX_RHEUMA_BIOLOGICAL_FinnGen_FinnGenR5_1       20002_1225_PanUKBB_PanUKBBR2_1 2.446605e-03    10   6064303
# 11:  10:6064589 RX_RHEUMA_BIOLOGICAL_FinnGen_FinnGenR5_1       M13_RHEUMA_FinnGen_FinnGenR7_1 2.069803e-03    10   6064589
# 12:  10:6064589 RX_RHEUMA_BIOLOGICAL_FinnGen_FinnGenR5_1       20002_1225_PanUKBB_PanUKBBR2_1 2.462998e-03    10   6064589
# 13:  10:6066476 RX_RHEUMA_BIOLOGICAL_FinnGen_FinnGenR5_1       M13_RHEUMA_FinnGen_FinnGenR7_1 2.609654e-04    10   6066476
# 14:  10:6066476 RX_RHEUMA_BIOLOGICAL_FinnGen_FinnGenR5_1       20002_1225_PanUKBB_PanUKBBR2_1 5.598710e-03    10   6066476
# 15:  10:6066476 RX_RHEUMA_BIOLOGICAL_FinnGen_FinnGenR5_1                     AAVMPO_Wong_up_1 2.394969e-02    10   6066476
# 16: 1:197642552              ph696.4_PanUKBB_PanUKBBR2_1        K11_CROHN_FinnGen_FinnGenR7_1 2.618519e-02     1 197642552
# 17: 1:197642552              ph696.4_PanUKBB_PanUKBBR2_1 L12_LICHENPLANUS_FinnGen_FinnGenR7_1 4.418516e-02     1 197642552
# 18: 1:197662011              ph696.4_PanUKBB_PanUKBBR2_1        K11_CROHN_FinnGen_FinnGenR7_1 2.618519e-02     1 197662011
# 19: 1:197662011              ph696.4_PanUKBB_PanUKBBR2_1 L12_LICHENPLANUS_FinnGen_FinnGenR7_1 4.418516e-02     1 197662011
# 20: 2:162254026              ph696.4_PanUKBB_PanUKBBR2_1           E4_DM1_FinnGen_FinnGenR7_1 5.555605e-04     2 162254026
# 21: 2:162254026              ph696.4_PanUKBB_PanUKBBR2_1        K11_CROHN_FinnGen_FinnGenR7_1 2.409490e-03     2 162254026
# 22: 2:162254026              ph696.4_PanUKBB_PanUKBBR2_1 L12_PSORI_ARTHRO_FinnGen_FinnGenR7_1 3.821361e-03     2 162254026
# 23: 2:162254026              ph696.4_PanUKBB_PanUKBBR2_1 L12_LICHENPLANUS_FinnGen_FinnGenR7_1 1.124813e-02     2 162254026
# 24: 5:159412209              ph696.4_PanUKBB_PanUKBBR2_1 L12_PSORI_ARTHRO_FinnGen_FinnGenR7_1 3.038197e-04     5 159412209
# 25: 5:159412209              ph696.4_PanUKBB_PanUKBBR2_1        K11_CROHN_FinnGen_FinnGenR7_1 1.202625e-02     5 159412209
# 26: 6:159044945              ph696.4_PanUKBB_PanUKBBR2_1      K11_COELIAC_FinnGen_FinnGenR7_1 5.602305e-03     6 159044945
# 27: 6:159044945              ph696.4_PanUKBB_PanUKBBR2_1           E4_DM1_FinnGen_FinnGenR7_1 1.645586e-02     6 159044945
# 28: 6:159044945              ph696.4_PanUKBB_PanUKBBR2_1   D3_SARCOIDOSIS_FinnGen_FinnGenR7_1 4.026911e-02     6 159044945
# 29: 19:10467167              ph696.4_PanUKBB_PanUKBBR2_1 L12_PSORI_ARTHRO_FinnGen_FinnGenR7_1 3.038197e-04    19  10467167
# 30: 19:10467167              ph696.4_PanUKBB_PanUKBBR2_1             JIA_LopezIsac_33106285_1 5.810514e-03    19  10467167
#           pid38                                IMD.focus                            IMD.other pairwise_fdr CHR38      BP38

# Some of these may be redundant signals. We'll check which SNPs are in close proximity and keep the ones with lower pairwise FDR
pfc[, BP38[1] - BP38[3]] # 5:132435113
pfc[, BP38[2] - BP38[3]] # 5:132448701
# These two SNPs are very close to 5:132461230, we'll remove them
pfc[, BP38[9] - BP38[13]] # 10:6064303
pfc[, BP38[11] - BP38[13]] # 10:6064589

# 1:197642552 and 1:197642552 are very close. They have similar pairwise FDR. 
pfc[, BP38[16] - BP38[18]] # 1:197642552
redall$ph696.4_PanUKBB_PanUKBBR2_1[ BP38 %in% c(pfc$BP38[16], pfc$BP38[18]), .(pid38, P2 = 2 * pnorm( -abs(BETA)/SE ))]
# 1:197642552 has slightly lower P-value (7.4e-4) in psoriasis
redall$K11_CROHN_FinnGen_FinnGenR7_1[ BP38 %in% c(pfc$BP38[16], pfc$BP38[18])]
# 1:197662011 has slightly lower P-value in Chron's
redall$L12_LICHENPLANUS_FinnGen_FinnGenR7_1[ BP38 %in% c(pfc$BP38[16], pfc$BP38[18])]
# Virtually same P-value
# We'll arbitrarily keep 1:197642552
# Remove unwanted SNPs
pfc <- pfc[ -c(1,2, 9:12, 16:17 )]

iof <- c( unique(pfc$IMD.focus),unique(pfc$IMD.other))

# Load dense-SNP data
ffc <- lapply(iof, function(x){
		if(grepl("FinnGenR5", x)){
			fn <- gsub("_FinnGen_FinnGenR5_1", "", x)
			if(!file.exists(paste0("../data/fg_sumstats/finngen_R5_", fn,".gz"))) system(paste0("wget ***REMOVED***finngen_R5_", fn, ".gz -O ../data/fg_sumstats/finngen_R5_",fn,".gz"))
			y = fread(paste0("../data/fg_sumstats/finngen_R5_", fn,".gz"), tmpdir = "tmp")
			y[ , pid38:=paste(`#chrom`, pos, sep = ":")]
			y <- y[, .(`#chrom`, pos, ref, alt, beta, sebeta, pval)]
			setnames(y, c("#chrom", "pos", "ref", "alt", "beta", "sebeta", "pval"), c("CHR38", "BP38", "REF", "ALT", "BETA", "SE", "P"))
		}else if(grepl("FinnGenR7", x)){
			fn <- gsub("_FinnGen_FinnGenR7_1", "", x)
			if(!file.exists(paste0("../data/fg_sumstats/finngen_R7_", fn,".gz"))) system(paste0("wget ***REMOVED***finngen_R7_", fn, ".gz -O ../data/fg_sumstats/finngen_R7_",fn,".gz"))
			y = fread(paste0("../data/fg_sumstats/finngen_R7_", fn,".gz"), tmpdir = "tmp")
			y[ , pid38:=paste(`#chrom`, pos, sep = ":")]
			y <- y[, .(`#chrom`, pos, ref, alt, beta, sebeta, pval)]
			setnames(y, c("#chrom", "pos", "ref", "alt", "beta", "sebeta", "pval"), c("CHR38", "BP38", "REF", "ALT", "BETA", "SE", "P"))
		} else if(grepl("PanUKBBR2", x)){
			if(!file.exists(paste0("../data/fg_sumstats/",x,"-hg38.tsv.gz"))){
				# PanUKBB requires a bit more processing
				# Transform into .gz so we can import it in fread
				mf <- fread("../data/Metadata_20230906-v1.tsv")
				url <- mf[Trait == x, URL]
				system(paste0("wget  ", url, " -O ../data/fg_sumstats/", x, ".bgz"))
				system(paste0("zcat ../data/fg_sumstats/",x,".bgz > ../data/fg_sumstats/",x,".tsv ; gzip ../data/fg_sumstats/", x, ".tsv"))
				f <- fread(paste0("../data/fg_sumstats/",x, ".tsv.gz"), tmpdir = "tmp")
				f <- f[, .(chr, pos, ref, alt, beta_meta, se_meta, neglog10_pval_meta)]
				names(f) <- c("CHR", "BP", "REF", "ALT", "BETA", "SE", "P")
				f <- f[complete.cases(f)]
				f[, P:=10^-P] # P-values are stored as -log10(P), so we'll revert them -- even though we'll likely don't need the P-values

				# Save it
				fwrite(f, paste0("../data/fg_sumstats/",x,"-hr.tsv.gz"), sep = "\t")

				# Pipeline it 
				system(paste0("cd ../data/fg_sumstats/; ~/rds/rds-cew54-basis/GWAS_tools/01-Pipeline/pipeline_v5.3.2_beta.sh -f ", x,"-hr.tsv.gz"))
				unlink(paste0("../data/fg_sumstats/", x, ".bgz"))
				unlink(paste0("../data/fg_sumstats/", x, ".tsv.gz"))
				unlink(paste0("../data/fg_sumstats/", x, "-hr.tsv.gz"))
			}
				y <- fread(paste0("../data/fg_sumstats/",x,"-hg38.tsv.gz"), tmpdir = "tmp")
		} else{
			y = fread(paste0("~/rds/rds-cew54-basis/02-Processed/", x, "-hg38.tsv.gz"), tmpdir = "tmp")
			y <- y[, .(CHR38, BP38, REF, ALT, BETA, SE, P)]
		}

})
names(ffc) <- iof

mf <- fread("../data/Metadata_20230906-v1.tsv")
mf <- mf[Trait %in% iof, .(Trait, N0, N1)]
n1=mf$N1
names(n1) <- mf$Trait
n0=mf$N0
names(n0) <- mf$Trait

# Prepare data for coloc
d2l=function(d, trait, chr, st, en) {
	df <- d[CHR38==chr & BP38>st & BP38<en & !is.na(SE) & !is.na(BETA)]
	df[, pid:=paste(CHR38, BP38, sep = ":")]
	df <- df[!duplicated(pid)]
    list(snp=df$pid,
         beta=df$BETA,
         varbeta=df$SE^2,
         p=df$P, # Not necessary for coloc, but we'll use p-values to call novelty
         type="cc",
         s=n1[trait]/(n1[trait]+n0[trait]))
    }
w=1e+6 # large window choice

index = pfc # Let's rebrand pfc as index
# This code will call the best SNP and find their P-values in both the myositis and the IMD dataset, so we don't need to run coloc twice.
# 1:nrow(index)
for(i in 1:nrow(index)) {
    message("Applying coloc on ", index$IMD.focus[i], " and ", index$IMD.other[i], " at SNP ", index$CHR38[i], ":", index$BP38[i], ".")
    st=index$BP38[i]-w
    en=index$BP38[i]+w
    chr=index$CHR38[i]
    d1=d2l( ffc[[ index$IMD.focus[i]]], trait=index$IMD.focus[[i]], chr = chr, st = st, en = en)
    d2=d2l( ffc[[ index$IMD.other[i]]], trait=index$IMD.other[[i]], chr = chr, st = st, en = en)
    
    result=coloc.abf(d1,d2, p12 = 5e-6) # Note: we modified the prior to be more conservative
    index[i ,c("nsnps","H0","H1","H2","H3","H4"):=as.list(result$summary)]
    best=result$results$snp[ which.max(result$results$SNP.PP.H4) ]
    w1=which(d1$snp==best)
    w2=which(d2$snp==best)
    index[i ,c("bestsnp","bestsnp.pp","pbest.focus","pbest.other", "pbest.focus.region", "pbest.other.region"):=
                 list(best, max(result$results$SNP.PP.H4),
                      d1$p[w1],
                      d2$p[w2],
                      min(d1$p),
                      min(d2$p))]
}
index[, .(pid38, IMD.focus, IMD.other, H3, H4, bestsnp)]

# Save it 
fwrite(index, "../data/coloc_val_results.tsv", sep = "\t")

# index <- fread("../data/coloc_val_results.tsv")

##################################
# Compare to previous validation #
##################################

# Create unified validation files
vv1 <- merge(tv1, tv1.10[ ,.(pid38, P.val.R10 = P.val)])
vv1[, IMD.focus:="ph696.4_PanUKBB_PanUKBBR2_1"]
vv2 <- merge(tv2, tv2.10[ ,.(pid38, P.val.R10 = P.val)])
vv2[, IMD.focus:="J10_COPD_FinnGen_FinnGenR5_1"]
vv3 <- merge(tv3, tv3.10[ ,.(pid38, P.val.R10 = P.val)])
vv3[, IMD.focus:="RX_RHEUMA_BIOLOGICAL_FinnGen_FinnGenR5_1"]
vv <- rbindlist(list(vv1,vv2,vv3), fill = TRUE)
vv[, SNPID:=NULL]
vv[, b.val:=ifelse(P.focus > P.val, "g", "b")][, b.val.10:=ifelse(P.focus > P.val.R10, "g", "b")]
vv
#           pid38      P.focus        P.val   P.val.R10                                IMD.focus  b.val b.val.10
#          <char>        <num>        <num>       <num>                                   <char> <char>   <char>
#  1: 19:10467167 6.529498e-07 2.087094e-04 1.74201e-03              ph696.4_PanUKBB_PanUKBBR2_1      b        b
#  2: 1:197642552 7.541812e-04 5.449465e-01 9.16579e-06              ph696.4_PanUKBB_PanUKBBR2_1      b        g
#  3: 1:197662011 8.388510e-04 5.895143e-01 7.84784e-06              ph696.4_PanUKBB_PanUKBBR2_1      b        g
#  4: 2:162254026 8.845999e-06 2.370537e-08 1.60081e-06              ph696.4_PanUKBB_PanUKBBR2_1      g        g
#  5: 5:159412209 5.314136e-07 3.367099e-04 1.86939e-17              ph696.4_PanUKBB_PanUKBBR2_1      b        g
#  6: 6:159044945 1.189902e-04 2.851924e-01 3.13437e-02              ph696.4_PanUKBB_PanUKBBR2_1      b        b
#  7: 5:132435113 3.623000e-04 4.580470e-05 1.88969e-06             J10_COPD_FinnGen_FinnGenR5_1      g        g
#  8: 5:132448701 4.056000e-04 5.245780e-05 2.38144e-06             J10_COPD_FinnGen_FinnGenR5_1      g        g
#  9: 5:132461230 1.080000e-04 4.184370e-08 7.59889e-06             J10_COPD_FinnGen_FinnGenR5_1      g        g
# 10:  10:6064303 1.719000e-05 1.579940e-04 1.22360e-04 RX_RHEUMA_BIOLOGICAL_FinnGen_FinnGenR5_1      b        b
# 11:  10:6064589 1.743000e-05 1.604280e-04 1.22823e-04 RX_RHEUMA_BIOLOGICAL_FinnGen_FinnGenR5_1      b        b
# 12:  10:6066476 3.849000e-07 9.567310e-06 2.83949e-04 RX_RHEUMA_BIOLOGICAL_FinnGen_FinnGenR5_1      b        b
# 13: 1:113834946 3.640000e-11 1.957940e-15 3.94912e-21 RX_RHEUMA_BIOLOGICAL_FinnGen_FinnGenR5_1      g        g
# Note that for coloc we removed some SNPs and kept others. This didn't have an effect in validation behavior, because 
# both removed and included had the same behaviour.
# We kept 1:197662011 and removed 1:197642552. Both were bad for Tsoi validation, but good for R10 validation.
# We kept 5:132461230 and removed 5:132435113 and 5:132448701. All three are good.
# We kept 10:6066476 and removed 10:6064303 and 10:6064589. All three were bad.


tindex <- merge(index, vv, by=c("pid38", "IMD.focus"))
tindex[, .(pid38, IMD.focus, IMD.other, H4, P.focus, P.val, P.val.R10, b.val, b.val.10)]
fwrite(tindex, "../data/coloc_val_results.tsv", sep = "\t")

ss <- copy(tindex)
ss <- ss[, .(IMD.focus, maxH4 = max(H4), P.focus, b.val, b.val.10), by = pid38]  %>% unique
ss[, coloc.sig:=ifelse(maxH4 > 0.5, "y", "n")]
table(ss$b.val, ss$coloc.sig)
table(ss$b.val.10, ss$coloc.sig)

tindex[pid38 == "10:6066476"]
