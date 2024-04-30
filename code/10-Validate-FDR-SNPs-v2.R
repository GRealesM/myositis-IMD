#########################################
##                                     ##
##     VALIDATING SNP SELECTION        ##
##                                     ##
#########################################

# Author: Guillermo Reales 
# Date last updated: 2024/04/26

# Background: We wanted to validate our method to narrow down SNPs. To do so, we will select pairs of IMD from large GWAS and similarly sized as IIM GWASs.
# This script will
#
# 1. Validation via pairwiseFDR and coloc
# * Project selected R5 IMD
# * Apply FDR together with filtered projections
# * Select FDR significant R5/pso and check which PCs they're significant for
# * Compute pairwise Bhattacharyya distance for each R5/pso trait with the set of selected IMD used with myositis
# * Select genetically close IMD to each of the R5/pso traits
# * Apply pairwise FDR and to find significant SNPs in the R5/pso traits
# * Check which SNPs replicate in the validation datasets (ie. R10)
# * Run coloc in pairwise pairs
# 2. Validation via coloc
# 3. Validation via genome-wide threshold only.
# 4. Summary tables and plots.

##########################################


## Load packages and required data sets
library(data.table)
library(magrittr)
library(IMDtools)
library(cupcake)
library(annotSnpStats)
library(fpc)
library(coloc)
library(ggplot2)
library(cowplot)

setDTthreads(20)
setwd("/home/gr440/rds/rds-cew54-basis/Projects/myositis-IMD/code")
bpath <- "/home/gr440/rds/rds-cew54-basis/03-Bases/IMD_basis/"
# FinnGen Manifest requires registation to be shared, so DELETE after.
url10man <- "***REMOVED***"


#############################################################################

## Helper functions

repzeros <- function(data){
		for (j in seq_len(ncol(data)))
    	set(data,which(is.na(data[[j]])),j,0)
}

# This function will compute the distances, call clusters and select clusters with more than one SNP
f=function(d) {
	if(nrow(d)==1){
		d[, cl:=paste0("chr", unique(CHR38), ".1")]
		return(d)
	} # This is another fix. We don't want to remove unique SNPs
	dist <- as.dist(abs(outer(d$BP38, d$BP38, "-")))
	dc <- hclust(dist)
	calls <- cutree(dc, h = 1e+6)
	d[, cl:=calls]
	# d <- d[, if(.N > 1) .SD, by = cl]
	# if(nrow(d) == 0)
	# 	return(NULL)
	d[, cl:=paste0("chr", unique(CHR38), ".", cl)]
	return(d)
}

###########################################################################

###############################################

# 1. Validation via pairwiseFDR and coloc

###############################################

################################
##  Project selected R5 IMD  ###
################################


# The bits below are similar to 00-Prepare_datasets_EDA.R. Our goal is to make the projection datasets 
# similar to those we used for the main analysis, but including the list of desired traits for validation.
# These validation traits should 
# (1) NOT be in the original basis (eg. T1D, RA, MS, etc.)
# (2) Have two instances per trait, one smaller (comparable in size to myositis), and one larger, for validation.
# We'll choose a number of FinnGenR5 with 1000+ cases, with the expectation to use their R7 counterparts as validation datasets.

# For some of these, we'll use FinnGen R5/R7 for validation/discovery. We have the R7 projections already, but we need to reproject the R5 ones

# Meta-data table, containing information about the datasets
m <- fread("../data/Metadata_20230906-v1.tsv")
m <- m[, .(Trait, First_Author, Reference, Trait_ID_2.0, Trait_long, Trait_class, N0, N1, N, Population, Public)]
m <- m[ Trait_class == "IMD" & Public == "Y"]

# We'll focus on R5 with more than 1000 cases, and remove all basis traits.
mr5 <- m[Reference == "FinnGenR5" & N1 > 1000 & !grepl("Asthma|Ulcerative colitis|UC |Chron|Crohn|coeliac disease|Type 1 diabetes|Rheumatoid arthritis|lupus|multiple sclerosis|Inflammatory bowel disease|IBD|Autoimmune diseases|Certain disorders|Diseases of the blood", Trait_long, ignore.case = TRUE)]
mr5[, .(Trait, Trait_long, N1)] # Some redundant diseases, so we'll remove them
mr5 <- mr5[-c(1, 7,10:12,14,19,21,24,28,33,36:37,42:44,47,51,53:54)] # Remove redundant traits. In this case we use Hypothy strict autoimmune, since it has an equivalent in R10
trp <- mr5$Trait

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

# We've incorporated the new datasets. Now we can apply the FDR procedure as in 00-Prepare_datasets_EDA.R


############################################################################################################
##  Apply FDR to new projections together with old ones and check which PCs they're FDR significant for  ###
############################################################################################################


# Apply 1% FDR correction to overall p for all remaining datasets
qf[, FDR.overall := p.adjust(overall_p, method = "BH"), by="Trait_class"] # Only IMD, so trait class shouldn't matter

# Apply FDR.PC to projections
pf[, FDR.PC:=p.adjust(P, method = "BH"), by = c("PC", "Trait_class")][, stars:=ifelse(!is.na(FDR.PC) & FDR.PC<0.05,"○","")][ FDR.PC < 0.01 , stars:="●"]

# Identify FinnGenR7 equivalent datasets. In this case we'll simply remove them
tnv <- gsub("_FinnGen_FinnGenR5_1", "_FinnGen_FinnGenR7_1", trp) # Find validation datasets, too
qf[Trait %in% tnv]
# One of them is coded slightly differently in R7, so we'll rename it
tnv[10] <- "N14_GLOMER_NEPHRITIS_FinnGen_FinnGenR7_1"
qf[Trait %in% tnv]

qf[Trait %in% trp, .(Trait, Trait_long, N0, N1, FDR.overall)][order(Trait_long, Trait)]
#                                             Trait                                              Trait_long     N0    N1   FDR.overall
#                                            <char>                                                  <char>  <int> <int>         <num>
#  1:               STILL_ADULT_FinnGen_FinnGenR5_1                     Adult-onset Still disease (FinnGen) 198544  3403  1.312205e-40
#  2:        D3_AGRANULOCYTOSIS_FinnGen_FinnGenR5_1                               Agranulocytosis (FinnGen) 215755  1527  6.147468e-02
#  3:       L12_ALLERGICCONTACT_FinnGen_FinnGenR5_1                    Alergic contact dermatitis (FinnGen) 198740  2204  1.234425e-04
#  4: H7_ALLERGICCONJUNCTIVITIS_FinnGen_FinnGenR5_1                       Allergic conjunctivitis (FinnGen) 208959  9833  1.060750e-14
#  5:           ALLERG_RHINITIS_FinnGen_FinnGenR5_1                             Allergic rhinitis (FinnGen) 212387  5527  2.824028e-22
#  6:         L12_URTICA_ALLERG_FinnGen_FinnGenR5_1                            Allergic urticaria (FinnGen) 212464  1169  7.882002e-01
#  7:            M13_ANKYLOSPON_FinnGen_FinnGenR5_1                        Ankylosing spondylitis (FinnGen) 164682  1462  2.732807e-05
#  8:          L12_PSORI_ARTHRO_FinnGen_FinnGenR5_1                        Arthropathic psoriasis (FinnGen) 212242  1637  4.524920e-03
#  9:                L12_ATOPIC_FinnGen_FinnGenR5_1                             Atopic dermatitis (FinnGen) 198740  7024  7.032371e-25
# 10:      RX_RHEUMA_BIOLOGICAL_FinnGen_FinnGenR5_1              Biological medication for rheuma (FinnGen) 216990  1802  1.292849e-18
# 11:                  J10_COPD_FinnGen_FinnGenR5_1                                          COPD (FinnGen) 186723  6915  5.894021e-05
# 12:         CHILDHOOD_ALLERGY_FinnGen_FinnGenR5_1                  Childhood allergy (age < 16) (FinnGen) 216044  2748  1.853316e-10
# 13:      L12_DERMATITISECZEMA_FinnGen_FinnGenR5_1                         Dermatitis and eczema (FinnGen) 198740 20052  2.126305e-22
# 14:          GLOMER_NEPHRITIS_FinnGen_FinnGenR5_1                            Glomerulonephritis (FinnGen) 214179  4613  1.046384e-02
# 15:                  M13_GOUT_FinnGen_FinnGenR5_1                                          Gout (FinnGen) 147221  3576  1.218630e-02
# 16:              E4_HYPERPARA_FinnGen_FinnGenR5_1                           Hyperparathyroidism (FinnGen) 211123  2928  3.427849e-01
# 17:            HYPOTHYROIDISM_FinnGen_FinnGenR5_1       Hypothyroidism (congenital or acquired) (FinnGen)  59827 26342 2.736205e-108
# 18:             N14_INFCERVIX_FinnGen_FinnGenR5_1          Inflammatory disease of cervix uteri (FinnGen) 111858  1093  7.986116e-01
# 19:               N14_INFLUTH_FinnGen_FinnGenR5_1                Inflammatory disease of uterus (FinnGen) 111858  2177  5.480857e-01
# 20:          N14_FEMALEGENINF_FinnGen_FinnGenR5_1 Inflammatory diseases of female pelvic organs (FinnGen) 111858 11721  2.596041e-01
# 21:                   K11_IBS_FinnGen_FinnGenR5_1                      Irritable bowel syndrome (FinnGen) 182423  4605  5.401270e-01
# 22:          L12_LICHENPLANUS_FinnGen_FinnGenR5_1                                 Lichen planus (FinnGen) 212242  1865  4.150854e-12
# 23:     L12_LICHENSCLERATROPH_FinnGen_FinnGenR5_1                Lichen sclerosus et atrophicus (FinnGen) 207482  1021  2.318569e-04
# 24:            J10_NASALPOLYP_FinnGen_FinnGenR5_1                                   Nasal polyp (FinnGen) 167849  3236  5.173694e-35
# 25:         K11_ENERCOLNONINF_FinnGen_FinnGenR5_1            Noninfective enteritis and colitis (FinnGen) 210300  8492  1.530951e-88
# 26:        RHEU_ARTHRITIS_OTH_FinnGen_FinnGenR5_1                          Other arthritis (FG) (FinnGen) 172834  5109  6.147825e-08
# 27:              POLLENALLERG_FinnGen_FinnGenR5_1                                Pollen allergy (FinnGen) 214464  2972  7.882412e-09
# 28:           M13_POLYMYALGIA_FinnGen_FinnGenR5_1                        Polymyalgia rheumatica (FinnGen) 213145  1523  5.319368e-03
# 29:             L12_PSORIASIS_FinnGen_FinnGenR5_1                                     Psoriasis (FinnGen) 212242  4510  5.876136e-07
# 30:               L12_ROSACEA_FinnGen_FinnGenR5_1                                       Rosacea (FinnGen) 211139  1195  1.876695e-02
# 31:            D3_SARCOIDOSIS_FinnGen_FinnGenR5_1                                   Sarcoidosis (FinnGen) 215712  2046  4.126903e-08
# 32:           L12_SEBORRHOEIC_FinnGen_FinnGenR5_1                        Seborrhoeic dermatitis (FinnGen) 198740  1253  5.763488e-02
# 33:               M13_SJOGREN_FinnGen_FinnGenR5_1                      Sicca syndrome [Sjogren] (FinnGen) 213145  1290  1.076818e-12
# 34:         SPONDYLOARTHRITIS_FinnGen_FinnGenR5_1                             Spondyloarthritis (FinnGen) 198544  3037  1.161757e-04
#                                             Trait                                              Trait_long     N0    N1   FDR.overall

# Take a look at the significant primary IMD
qf[Trait %in% trp & FDR.overall < 0.01, .(Trait, Trait_long, N0, N1, FDR.overall)][order(Trait_long, Trait)]
# Only 23/34
ttid <- qf[Trait %in% trp & FDR.overall < 0.01, Trait] # List of focus IMD

#############################################################################
# Select FDR significant R5/pso and check which PCs they're significant for #
#############################################################################


## Check which PC is significant for each focus IMD
sigpcs <- lapply(ttid, function(i){
	    p.i <- pf[ Trait == i & FDR.PC < 0.01, PC]
		p.i
})
names(sigpcs) <- ttid

# L12_PSORI_ARTHRO_FinnGen_FinnGenR5_1 has no significant PCs. We'll remove it
sigpcs$L12_PSORI_ARTHRO_FinnGen_FinnGenR5_1  <- NULL
ttid <- ttid[ttid != "L12_PSORI_ARTHRO_FinnGen_FinnGenR5_1"]
# 22/34 now

sigpcs.dt <- data.table(IMD.focus = names(sigpcs), sig.PCs = lapply(sigpcs, paste, collapse=", ") %>% unlist)
fwrite(sigpcs.dt, "../data/R5_sig_PCs.tsv", sep="\t")

######################################################################################################################
### Compute pairwise Bhattacharyya distance for each R5/pso trait with the set of selected IMD used with myositis  ###
######################################################################################################################

# Load selected IMD set
ps <- fread("../data/ps2.tsv")
qs <- fread("../data/qs2.tsv") # From here, we'll take the list of IMD used with myositis. These are significant proj with redundant removed.

# Update tnv to remove equivalent R7 datasets
tnv2 <- tnv[tnv %in% gsub("_FinnGen_FinnGenR5_1","_FinnGen_FinnGenR7_1", ttid)] # In this case, no update of ID needed

imdset <- qs[!grepl("Miller|Rothwell", Trait)][!Trait %in% tnv2, Trait]
# Remove Myositis datasets and R7 datasets from the IMD pool

# We modified the following function wrt 01-compute-bhatta... to
# * Take an R5 dataset
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

################################################################
# Select genetically close IMD to each of the R5/pso traits    #
################################################################

bhsel <- lapply(bhdist, function(x){
		x  <- x[1:17] # Take the closest 17 IMD, as we did with Myositis
}) 
bhsel %<>% rbindlist
unique(bhsel$T2)
# 45 IMD to consider for FDR

# Some of the close diseases correspond to some of the initial R7 set (but won't be equivalents to any validation dataset)
bhsel[T2 %in% tnv]


#######################################################################
## Apply pairwise FDR and to find significant SNPs in the R5 traits ###
#######################################################################

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

redimd.l <- bhsel$T2

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
names(pw.fdr.sig) <- ttid
lapply(pw.fdr.sig,nrow)

# Some have no significant pairwise FDR SNP
zero.fdr <- names(Filter(function(x) nrow(x) == 0, pw.fdr.sig)) 

ttid <- ttid[!ttid %in% zero.fdr]

pw.fdr.sig[zero.fdr]  <- NULL
names(pw.fdr.sig)

# Extract pairwise-FDR sig SNPs
pw.snps <- lapply(pw.fdr.sig, function(x){
			y <- x[, unique(pid38)]
			y
})

#####################################################################
# Check which SNPs replicate in the validation datasets (ie. R10)   #
#####################################################################

# Now for the validation datasets! We'll use FinnGen R10 datasets for validation

r10man <- fread(url10man)
fid <- gsub("_FinnGen_FinnGenR5_1", "", ttid)
r10man <- r10man[phenocode %in% fid] # All present

# Let's take a quick look at what we have
r10man[,Trait:=paste0(phenocode, "_FinnGen_FinnGenR10_1")]
setnames(r10man, c("phenotype", "num_controls", "num_cases"), c("Trait_long", "N0", "N1"))

tv.man <- rbind(qf[Trait %in% ttid, .(Trait, Trait_long, N0, N1)], r10man[ ,.(Trait, Trait_long, N0, N1)])  %>% .[order(Trait, decreasing = TRUE)]
tv.man[, Trait_long:=gsub(" \\(FinnGen\\)", "", Trait_long)][, Trait_long:=gsub("Sjogren", "Sjögren", Trait_long)] 
tv.man <- tv.man %>% .[order(Trait_long)]
#                                              Trait                         Trait_long     N0    N1
#                                             <char>                             <char>  <int> <int>
#  1:                STILL_ADULT_FinnGen_FinnGenR5_1          Adult-onset Still disease 198544  3403
#  2:               STILL_ADULT_FinnGen_FinnGenR10_1          Adult-onset Still disease 372620   125
#  3:  H7_ALLERGICCONJUNCTIVITIS_FinnGen_FinnGenR5_1            Allergic conjunctivitis 208959  9833
#  4: H7_ALLERGICCONJUNCTIVITIS_FinnGen_FinnGenR10_1            Allergic conjunctivitis 388516 23665
#  5:            ALLERG_RHINITIS_FinnGen_FinnGenR5_1                  Allergic rhinitis 212387  5527
#  6:           ALLERG_RHINITIS_FinnGen_FinnGenR10_1                  Allergic rhinitis 392069 12240
#  7:                 L12_ATOPIC_FinnGen_FinnGenR5_1                  Atopic dermatitis 198740  7024
#  8:                L12_ATOPIC_FinnGen_FinnGenR10_1                  Atopic dermatitis 367046 15208
#  9:       RX_RHEUMA_BIOLOGICAL_FinnGen_FinnGenR5_1   Biological medication for rheuma 216990  1802
# 10:      RX_RHEUMA_BIOLOGICAL_FinnGen_FinnGenR10_1   Biological medication for rheuma 407064  5117
# 11:                   J10_COPD_FinnGen_FinnGenR5_1                               COPD 186723  6915
# 12:                  J10_COPD_FinnGen_FinnGenR10_1                               COPD 338303 20066
# 13:          CHILDHOOD_ALLERGY_FinnGen_FinnGenR5_1       Childhood allergy (age < 16) 216044  2748
# 14:         CHILDHOOD_ALLERGY_FinnGen_FinnGenR10_1       Childhood allergy (age < 16) 406071  6110
# 15:       L12_DERMATITISECZEMA_FinnGen_FinnGenR5_1              Dermatitis and eczema 198740 20052
# 16:      L12_DERMATITISECZEMA_FinnGen_FinnGenR10_1              Dermatitis and eczema 367046 45135
# 17:         E4_HYTHY_AI_STRICT_FinnGen_FinnGenR5_1  Hypothyroidism, strict autoimmune 175475 22997
# 18:        E4_HYTHY_AI_STRICT_FinnGen_FinnGenR10_1  Hypothyroidism, strict autoimmune 298847 45321
# 19:             J10_NASALPOLYP_FinnGen_FinnGenR5_1                        Nasal polyp 167849  3236
# 20:            J10_NASALPOLYP_FinnGen_FinnGenR10_1                        Nasal polyp 308457  6841
# 21:          K11_ENERCOLNONINF_FinnGen_FinnGenR5_1 Noninfective enteritis and colitis 210300  8492
# 22:         K11_ENERCOLNONINF_FinnGen_FinnGenR10_1 Noninfective enteritis and colitis 392974 19207
# 23:         RHEU_ARTHRITIS_OTH_FinnGen_FinnGenR5_1               Other arthritis (FG) 172834  5109
# 24:        RHEU_ARTHRITIS_OTH_FinnGen_FinnGenR10_1               Other arthritis (FG) 311210 12079
# 25:               POLLENALLERG_FinnGen_FinnGenR5_1                     Pollen allergy 214464  2972
# 26:              POLLENALLERG_FinnGen_FinnGenR10_1                     Pollen allergy 402286  6851
# 27:            M13_POLYMYALGIA_FinnGen_FinnGenR5_1             Polymyalgia rheumatica 213145  1523
# 28:           M13_POLYMYALGIA_FinnGen_FinnGenR10_1             Polymyalgia rheumatica 399355  3871
# 29:              L12_PSORIASIS_FinnGen_FinnGenR5_1                          Psoriasis 212242  4510
# 30:             L12_PSORIASIS_FinnGen_FinnGenR10_1                          Psoriasis 397564 10312
# 31:             D3_SARCOIDOSIS_FinnGen_FinnGenR5_1                        Sarcoidosis 215712  2046
# 32:            D3_SARCOIDOSIS_FinnGen_FinnGenR10_1                        Sarcoidosis 405620  4399
# 33:                M13_SJOGREN_FinnGen_FinnGenR5_1           Sicca syndrome [Sjögren] 213145  1290
# 34:               M13_SJOGREN_FinnGen_FinnGenR10_1           Sicca syndrome [Sjögren] 399355  2735
# 35:          SPONDYLOARTHRITIS_FinnGen_FinnGenR5_1                  Spondyloarthritis 198544  3037
# 36:         SPONDYLOARTHRITIS_FinnGen_FinnGenR10_1                  Spondyloarthritis 372620  6991
#                                              Trait                         Trait_long     N0    N1
# Check if validation datasets are larger in cases than primary trait, and remove if not
tv.man[, rel:=ifelse(grepl("R5", Trait), "R5", "R10")]
st <- dcast(tv.man, Trait_long ~ rel, value.var = "N1")
st[, largerval:= R5 < R10]
# Only Still disease has larger cases in R5 than R10, so we'll remove it

tv.man <- tv.man[Trait_long != "Adult-onset Still disease"]

ttid <- tv.man[grepl("R5", Trait), Trait]  %>% .[order(.)] # Update ttid
val.l <- tv.man[grepl("R10", Trait), Trait]  %>% .[order(.)] # Create list of validation datasets in the same order
# Checked they're the same

# Now we'll incorporate the P-values of the primary and the validation dataset for each IMD

valp <- lapply(seq_along(ttid), function(x){
	pn <- ttid[x]
	vn <- val.l[x]

	fd <- pw.fdr.sig[[pn]]
	psnp <- pw.snps[[pn]]
	pr <- redall[[pn]]
	if(!file.exists(paste0("../data/fg_sumstats/", vn, "-hg38.tsv.gz"))){
		urlv <- r10man[ Trait == vn, path_https]
		system(paste0("wget ", urlv, " -O ../data/fg_sumstats/", vn, "-hg38.tsv.gz"))
	}

	vr <- fread(paste0("../data/fg_sumstats/", vn, "-hg38.tsv.gz"), tmpdir = "tmp/")
	setnames(vr, c("#chrom", "pos", "ref", "alt", "pval"), c("CHR38", "BP38", "REF", "ALT", "P"))
	vr[, pid38:=paste(CHR38, BP38, sep = ":")]
	vr <- vr[, .(pid38, REF, ALT, P)]

	# Filter SNPs
	pr <- pr[pid38 %in% psnp]
	vr <- vr[pid38 %in% psnp]
	if(nrow(pr) != nrow(vr)) message(pn, " and ", vn, " don't have the same rows. Check.")
	
	ptg <- merge(pr[,.(pid38, REF, ALT, P)], vr[,.(pid38, REF, ALT,P)], suffixes=c(".focus",".val"), by = c("pid38", "REF", "ALT"), all.x = TRUE)

	ff <- merge(fd, ptg, by = "pid38", all.x = TRUE)
	ff

})
valp %<>% rbindlist

# Remove SNPs impossible to validate
valp[is.na(P.val)] # Only 1 SNP
valp <- valp[!is.na(P.val)]
valp[, c("CHR38", "BP38"):=tstrsplit(pid38, split = ":", fixed = TRUE)][, c("CHR38", "BP38"):=list(as.numeric(CHR38), as.numeric(BP38))]

valp[, b.val:=ifelse(P.focus > P.val, "good", "bad")]

v2 <- valp[, .(pid38, IMD.focus, P.focus, P.val, b.val)]  %>% unique
table(v2$IMD.focus, v2$b.val)
#                                                 bad good
#   ALLERG_RHINITIS_FinnGen_FinnGenR5_1             2    6
#   CHILDHOOD_ALLERGY_FinnGen_FinnGenR5_1           1    0
#   D3_SARCOIDOSIS_FinnGen_FinnGenR5_1              0    3
#   E4_HYTHY_AI_STRICT_FinnGen_FinnGenR5_1          2   47
#   H7_ALLERGICCONJUNCTIVITIS_FinnGen_FinnGenR5_1   4   16
#   J10_COPD_FinnGen_FinnGenR5_1                    0    3
#   J10_NASALPOLYP_FinnGen_FinnGenR5_1              1   45
#   K11_ENERCOLNONINF_FinnGen_FinnGenR5_1           9   49
#   L12_ATOPIC_FinnGen_FinnGenR5_1                  1   19
#   L12_DERMATITISECZEMA_FinnGen_FinnGenR5_1        1   24
#   L12_PSORIASIS_FinnGen_FinnGenR5_1               4   13
#   M13_POLYMYALGIA_FinnGen_FinnGenR5_1             1    0
#   M13_SJOGREN_FinnGen_FinnGenR5_1                 0    5
#   POLLENALLERG_FinnGen_FinnGenR5_1                1    4
#   RHEU_ARTHRITIS_OTH_FinnGen_FinnGenR5_1          0    1
#   RX_RHEUMA_BIOLOGICAL_FinnGen_FinnGenR5_1        3    1
#   SPONDYLOARTHRITIS_FinnGen_FinnGenR5_1           2    0

# fwrite(valp, "../data/R5R10_validation.tsv", sep = "\t")
valp <- fread("../data/R5R10_validation.tsv")

v3 <- v2[, .(count = .N), by=c("IMD.focus", "b.val")] %>% dcast(. , IMD.focus ~ b.val, value.var = "count")

# Replace NA by zeros in both columns
repzeros(v3)
v3[, pbad:=bad/(good+bad)]
sum(v3$bad)/(sum(v3$good)+sum(v3$bad))
# 0.119 (11.9%) are bad.


#######################################
###   Run coloc in pairwise pairs   ###
#######################################


pfc <- valp[, .(pid38, CHR38, BP38, IMD.focus, IMD.other, pairwise_fdr)][order(IMD.focus, CHR38, BP38)]

# In the next step, we'll remove SNPs in close proximity, if any
pfl <- split(pfc, pfc[,.(IMD.focus)])

cl.snps <- lapply(pfl, function(x){
	y <- x[,.(pid38, CHR38,BP38)]  %>% unique

	ss <- split(y[, .(pid38, CHR38, BP38)], y[,.(CHR38)]) # split by chr
	cl.snp  <- ss %>% lapply(., f)  %>% rbindlist(., use.names = TRUE)
	cl.snp
})

# Now use pairwise fdr to select SNPs
fpfc <- lapply(seq_along(pfl), function(x){
	pfx <- pfl[[x]]
	clx <- cl.snps[[x]]
	if(is.null(clx)) return(pfx)

	pfx <- merge(pfx, clx[,.(pid38,cl)], by = "pid38")
	tokeep <- pfx[  , .SD[which.min(pairwise_fdr)] , by=cl][, pid38] # For each cluster, keep the one with lowest pairwise_fdr
	pfx <- pfx[pid38 %in% tokeep]
	pfx
})  %>% rbindlist


# Load dense-SNP data
iof <- c( unique(fpfc$IMD.focus),unique(fpfc$IMD.other))

# This will take a while, if there are many datasets and they're not available in the data directory
ffc <- lapply(iof, function(x){
		message("Working on ", x, ".")
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

index = fpfc # Let's rebrand pfc as index
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

# Save coloc results
# fwrite(index, "../data/coloc_valR10_results.tsv", sep = "\t")

index <- fread("../data/coloc_valR10_results.tsv")

# Compare to previous validation
tindex <- merge(index, unique(valp[,.(pid38, IMD.focus, P.focus, P.val, b.val)]))

ss <- copy(tindex)
ss <- ss[, .(maxH4 = max(H4), b.val), by = c("pid38", "IMD.focus")]  %>% unique
ss[ b.val == "good" & maxH4 > 0.5, b.val05:="TP"][ b.val == "bad" & maxH4 <= 0.5,  b.val05:="TN"][ b.val == "good" & maxH4 <= 0.5,  b.val05:="FN"][ b.val == "bad" & maxH4 > 0.5,  b.val05:="FP"]
ss[ b.val == "good" & maxH4 > 0.8, b.val08:="TP"][ b.val == "bad" & maxH4 <= 0.8,  b.val08:="TN"][ b.val == "good" & maxH4 <= 0.8,  b.val08:="FN"][ b.val == "bad" & maxH4 > 0.8,  b.val08:="FP"]

ss05 <- ss[, .(count = .N), by=c("IMD.focus", "b.val05")] %>% dcast(. , IMD.focus ~ b.val05, value.var = "count")  
repzeros(ss05)
ss05[, sens:=TP/(TP + FN)][, spec:=TN/(TN + FP)][, FPR:=FP/(TP+FP)]
ss05
#                                         IMD.focus    FN    FP    TN    TP      sens  spec        FPR
#                                            <char> <int> <int> <int> <int>     <num> <num>      <num>
#  1:           ALLERG_RHINITIS_FinnGen_FinnGenR5_1     0     1     0     5 1.0000000   0.0 0.16666667
#  2:         CHILDHOOD_ALLERGY_FinnGen_FinnGenR5_1     0     1     0     0       NaN   0.0 1.00000000
#  3:            D3_SARCOIDOSIS_FinnGen_FinnGenR5_1     2     0     0     0 0.0000000   NaN        NaN
#  4:        E4_HYTHY_AI_STRICT_FinnGen_FinnGenR5_1    13     0     2    16 0.5517241   1.0 0.00000000
#  5: H7_ALLERGICCONJUNCTIVITIS_FinnGen_FinnGenR5_1     2     1     1     6 0.7500000   0.5 0.14285714
#  6:                  J10_COPD_FinnGen_FinnGenR5_1     1     0     0     0 0.0000000   NaN        NaN
#  7:            J10_NASALPOLYP_FinnGen_FinnGenR5_1     8     0     1    11 0.5789474   1.0 0.00000000
#  8:         K11_ENERCOLNONINF_FinnGen_FinnGenR5_1     8     1     4    22 0.7333333   0.8 0.04347826
#  9:                L12_ATOPIC_FinnGen_FinnGenR5_1     4     0     1     5 0.5555556   1.0 0.00000000
# 10:      L12_DERMATITISECZEMA_FinnGen_FinnGenR5_1     7     1     0     5 0.4166667   0.0 0.16666667
# 11:             L12_PSORIASIS_FinnGen_FinnGenR5_1     5     0     3     3 0.3750000   1.0 0.00000000
# 12:           M13_POLYMYALGIA_FinnGen_FinnGenR5_1     0     1     0     0       NaN   0.0 1.00000000
# 13:               M13_SJOGREN_FinnGen_FinnGenR5_1     0     0     0     2 1.0000000   NaN 0.00000000
# 14:              POLLENALLERG_FinnGen_FinnGenR5_1     1     1     0     2 0.6666667   0.0 0.33333333
# 15:        RHEU_ARTHRITIS_OTH_FinnGen_FinnGenR5_1     0     0     0     1 1.0000000   NaN 0.00000000
# 16:      RX_RHEUMA_BIOLOGICAL_FinnGen_FinnGenR5_1     0     1     0     1 1.0000000   0.0 0.50000000
# 17:         SPONDYLOARTHRITIS_FinnGen_FinnGenR5_1     0     0     2     0       NaN   1.0        NaN
ss05[,sum(FP)/(sum(FP) + sum(TP))]
# 0.09195402


ss08 <- ss[, .(count = .N), by=c("IMD.focus", "b.val08")] %>% dcast(. , IMD.focus ~ b.val08, value.var = "count")  
repzeros(ss08)
ss08[, sens:=TP/(TP + FN)][, spec:=TN/(TN + FP)][, FPR:=FP/(TP+FP)]
ss08
#                                         IMD.focus    FN    FP    TN    TP      sens  spec   FPR
#                                            <char> <int> <int> <int> <int>     <num> <num> <num>
#  1:           ALLERG_RHINITIS_FinnGen_FinnGenR5_1     0     0     1     5 1.0000000     1   0.0
#  2:         CHILDHOOD_ALLERGY_FinnGen_FinnGenR5_1     0     1     0     0       NaN     0   1.0
#  3:            D3_SARCOIDOSIS_FinnGen_FinnGenR5_1     2     0     0     0 0.0000000   NaN   NaN
#  4:        E4_HYTHY_AI_STRICT_FinnGen_FinnGenR5_1    16     0     2    13 0.4482759     1   0.0
#  5: H7_ALLERGICCONJUNCTIVITIS_FinnGen_FinnGenR5_1     3     0     2     5 0.6250000     1   0.0
#  6:                  J10_COPD_FinnGen_FinnGenR5_1     1     0     0     0 0.0000000   NaN   NaN
#  7:            J10_NASALPOLYP_FinnGen_FinnGenR5_1    12     0     1     7 0.3684211     1   0.0
#  8:         K11_ENERCOLNONINF_FinnGen_FinnGenR5_1    13     0     5    17 0.5666667     1   0.0
#  9:                L12_ATOPIC_FinnGen_FinnGenR5_1     7     0     1     2 0.2222222     1   0.0
# 10:      L12_DERMATITISECZEMA_FinnGen_FinnGenR5_1     8     0     1     4 0.3333333     1   0.0
# 11:             L12_PSORIASIS_FinnGen_FinnGenR5_1     7     0     3     1 0.1250000     1   0.0
# 12:           M13_POLYMYALGIA_FinnGen_FinnGenR5_1     0     1     0     0       NaN     0   1.0
# 13:               M13_SJOGREN_FinnGen_FinnGenR5_1     0     0     0     2 1.0000000   NaN   0.0
# 14:              POLLENALLERG_FinnGen_FinnGenR5_1     1     0     1     2 0.6666667     1   0.0
# 15:        RHEU_ARTHRITIS_OTH_FinnGen_FinnGenR5_1     0     0     0     1 1.0000000   NaN   0.0
# 16:      RX_RHEUMA_BIOLOGICAL_FinnGen_FinnGenR5_1     0     1     0     1 1.0000000     0   0.5
# 17:         SPONDYLOARTHRITIS_FinnGen_FinnGenR5_1     0     0     2     0       NaN     1   NaN
ss08[,sum(FP)/(sum(FP) + sum(TP))]
# 0.04761905

###############################################

# 2. Validation via coloc

###############################################

# A couple more things. Let's compare myositis coloc and validation's false discovery rate using 1 - H4

cc <- fread("../data/coloc_results-v3.tsv")

vh <- ss[ maxH4 > 0.5]
mh <- cc[ H4 > 0.5][, .(maxH4 = max(H4), pid, trait.myos), by = c("pid", "trait.myos")]  %>% unique

vh[, mean(1 - maxH4)]
# [1] 0.1524762
mh[, mean(1 - maxH4)]
# [1] 0.2230366

# What about H4 > 0.8 threshold?

vh[ maxH4 > 0.8, mean(1 - maxH4)]
# [1] 0.07311549
mh[ maxH4 > 0.8, mean(1 - maxH4)]
# [1] 0.1023415

###############################################

# 3. Validation via genome-wide threshold only.

###############################################

# Lastly, we'd like to investigate the proportion of lead SNPs in the R5 datasets that have smaller p-values in the R10 set.
# This would give an estimated fpr for the gw sig threshold calculated by the same manner

# To avoid having to re-run a lot of code, we can recreate the file list from tindex
ttid  <- unique(tindex$IMD.focus) %>% .[order(.)]
val.l <- gsub("R5", "R10", ttid)
fn <- gsub("_FinnGen_FinnGenR5_1", "", ttid) # R5 have slightly different naming system
lapply(fn, function(i){	file.exists(paste0("../data/fg_sumstats/finngen_R5_", i, ".gz"))}) %>% unlist
lapply(val.l, function(i){ file.exists(paste0("../data/fg_sumstats/", i, "-hg38.tsv.gz"))}) %>% unlist

lead.R5 <- lapply(fn, function(x){
		# Load and change names, just as before
		message("Working on ", x)
		y <- fread(paste0("../data/fg_sumstats/finngen_R5_", x, ".gz"), tmpdir = "tmp/") %>% .[ pval < 5e-8]
		if(nrow(y) == 0){ message(x, " doesn't have any gw-sig SNPs."); return(NULL)}
		y <- y[, .(`#chrom`, pos, ref, alt, beta, sebeta, pval)]
		setnames(y, c("#chrom", "pos", "ref", "alt", "beta", "sebeta", "pval"), c("CHR38", "BP38", "REF", "ALT", "BETA", "SE", "P"))
		y[ , pid38:=paste(CHR38, BP38, sep = ":")]
		# As before, split and assign clusters
		ss <- split(y[, .(pid38, CHR38, BP38)], y[,.(CHR38)]) # split by chr
		cl.snp  <- ss %>% lapply(., f)  %>% rbindlist(., use.names = TRUE)
		y <- merge(y, cl.snp[,.(pid38,cl)], by = "pid38")
		y <- y[, .SD[which.min(P)] , by=cl ]
		# Lastly, remove MHC, as we don't use it in our method
		mhc.pids <- y[as.numeric(CHR38) == 6 & BP38 > 20e6 & BP38 < 40e6, pid38]
		if(length(mhc.pids) > 0) y <- y[!pid38 %in% mhc.pids]
		if(nrow(y) == 0){ message(x, " doesn't have SNPs after MHC filtering."); return(NULL)}
		y
})
names(lead.R5) <- ttid
nsig.snps <- lapply(lead.R5, nrow)
nsig.snps <- data.table(IMD.focus = names(unlist(nsig.snps)), sig.snps = unlist(nsig.snps))
fwrite(nsig.snps, "../data/R5_sig_snps.tsv", sep = "\t")

wR10 <- lapply(seq_along(ttid), function(x){
		pn <- ttid[[x]]
		message("Working on ", pn)
		pf <- lead.R5[[pn]]
		if(is.null(pf)){ message("No SNPs to validate in ", pn); return(NULL)}
		pf <- pf[, .(pid38, CHR38, BP38, REF, ALT, P, cl)]
		vn <- val.l[x]
		message("Using ", vn, " as a validation dataset")
		vf <- fread(paste0("../data/fg_sumstats/", vn, "-hg38.tsv.gz"), tmpdir = "tmp/")
		vf <- vf[, .(`#chrom`, pos, ref, alt, beta, sebeta, pval)]
		setnames(vf, c("#chrom", "pos", "ref", "alt", "pval"), c("CHR38", "BP38", "REF", "ALT", "P"))
		vf[ , pid38:=paste(CHR38, BP38, sep = ":")]
		vf <- vf[pid38 %in% pf$pid38]
		ff <- merge(pf, vf[, .(pid38, REF, ALT, P)], by=c("pid38", "REF", "ALT"), suffixes=c(".focus",".val"), all.x = TRUE)
		ff[, IMD.focus := pn]
}) %>% rbindlist
wR10[, b.val:=ifelse(P.focus > P.val, "good", "bad")][, CHR38:=as.numeric(CHR38)]
wR10 <- wR10[!is.na(b.val)] # Remove SNPs that couldn't be validated
wR10 <- wR10[order(IMD.focus, CHR38, BP38 )]
# Save wR10
fwrite(wR10, "../data/R5R10_gwsig_validation.tsv", sep = "\t")
# Load it, if you want to save the previous step
wR10 <- fread("../data/R5R10_gwsig_validation.tsv")

ssr10 <- dcast(wR10[, .(pid38, IMD.focus, b.val)], IMD.focus ~ b.val)
sum(ssr10$bad)/(sum(ssr10$bad)+sum(ssr10$good))
# [1] 0.1533333

###############################################

# 4. Summary tables and plots.

###############################################

# To make the following tables, we can load results we've been saving throughout our process.

valp <- fread("../data/R5R10_validation.tsv")
index <- fread("../data/coloc_valR10_results.tsv")
tindex <- merge(index, unique(valp[,.(pid38, IMD.focus, P.focus, P.val, b.val)]))
cc <- fread("../data/coloc_results-v3.tsv")
wR10 <- fread("../data/R5R10_gwsig_validation.tsv")
nsig.snps <- fread("../data/R5_sig_snps.tsv")
sigpcs.dt <- fread("../data/R5_sig_PCs.tsv")
r10man <- fread(url10man)
m <- fread("../data/Metadata_20230906-v1.tsv")


# Summary table of the validation datasets
mimd <- m[Trait %in% unique(index$IMD.focus), .(IMD.focus = Trait, N1.R5 = N1)][order(IMD.focus)]
mimd <- merge(mimd, nsig.snps)
mimd <- merge(mimd, sigpcs.dt)
r10n <- r10man[phenocode %in% gsub("_FinnGen_FinnGenR5_1", "", mimd$IMD.focus), .(phenocode, num_cases, phenotype)]
r10n[, IMD.focus:=paste0(phenocode,"_FinnGen_FinnGenR5_1")]
mimd <- merge(mimd, r10n)
mimd <- mimd[, .(phenocode, phenotype, sig.PCs, sig.snps,N1.R5, num_cases)]
setnames(mimd, c("sig.PCs", "sig.snps", "num_cases"), c("sig.PC.R5", "gwsig.SNP.R5", "N1.R10"))
fwrite(mimd, "../tables/ST6_validation_R5_info.tsv", sep ="\t")

# Let's collect the numbers

# For this bit, you can scroll up and run the relevant bits of code to generate the relevant objects, 
# which are quick to create from the saved files we loaded right above.



sxt <- data.table(validation = c("pwFDR + coloc (0.5)",
						  "pwFDR + coloc (0.8)", 
						  "coloc (0.5)", 
						  "coloc (0.8)", 
						  "Empirical"), 						  
		  FPR = c(ss05[,sum(FP)/(sum(FP) + sum(TP))],
		  		  ss08[,sum(FP)/(sum(FP) + sum(TP))],
				  vh[, mean(1 - maxH4)],
				  vh[ maxH4 > 0.8, mean(1 - maxH4)],
				  sum(ssr10$bad)/(sum(ssr10$bad)+sum(ssr10$good))),
		`Validating SNPs` = c(ss05[,sum(TP)], 
							  ss08[,sum(TP)],
							  NA,
							  NA,
							  sum(ssr10$good)),
							  
		`Non-validating SNPs` = c(ss05[,sum(FP)], 
							  ss08[,sum(FP)],
							  NA,
							  NA,
							  sum(ssr10$bad))
							  )

fwrite(sxt, "../tables/ST7_validation_R5_results.tsv", sep ="\t")


#################################################


ssv <- copy(tindex)
ssv <- ssv[, .(maxH4 = max(H4), b.val, P.focus, P.val), by = c("pid38", "IMD.focus")]  %>% unique
ssv <- ssv[maxH4 > 0.5]

sp1 <- ggplot(ssv, aes(x = -log10(P.focus), y = -log10(P.val)))+
			geom_point()+
			# xlab("-log10(P) FinnGen R5")+
			ylab(bquote(-log[10](P)~R10))+
			annotate("text", x = 95, y = 95, colour = "darkgreen", size = 7, label = "79\n\n8")+
			geom_abline( colour = "red")+
			geom_vline(xintercept = -log10(5e-8), colour = "blue")+
			ggtitle("pwFDR + coloc (0.5)")+
			theme_cowplot()+
			theme(axis.title.x = element_blank())

sp1z <- ggplot(ssv, aes(x = -log10(P.focus), y = -log10(P.val)))+
			geom_point()+
			xlim(c(0,25))+
			ylim(c(0,50))+
			# xlab("-log10(P) FinnGen R5")+
			# ylab("-log10(P) FinnGen R10")+
			geom_abline( colour = "red")+
			geom_vline(xintercept = -log10(5e-8), colour = "blue")+
			ggtitle("zoomed in")+
			theme_cowplot()+
			theme(axis.title.x = element_blank(),
				axis.title.y = element_blank())

sp2 <- ggplot(ssv[ maxH4 > 0.8], aes(x = -log10(P.focus), y = -log10(P.val)))+
			geom_point()+
			annotate("text", x = 95, y = 95, colour = "darkgreen", size = 7, label = "60\n\n3")+
			ylab(bquote(-log[10](P)~R10))+
			geom_abline( colour = "red")+
			geom_vline(xintercept = -log10(5e-8), colour = "blue")+
			ggtitle("pwFDR + coloc (0.8)")+
			theme_cowplot()+
			theme(axis.title.x = element_blank())

sp2z <- ggplot(ssv[ maxH4 > 0.8], aes(x = -log10(P.focus), y = -log10(P.val)))+
			geom_point()+
			xlim(c(0,25))+
			ylim(c(0,50))+
			# xlab("-log10(P) FinnGen R5")+
			# ylab("-log10(P) FinnGen R10")+
			geom_abline( colour = "red")+
			geom_vline(xintercept = -log10(5e-8), colour = "blue")+
			ggtitle("zoomed in")+
			theme_cowplot()+
			theme(axis.title.x = element_blank(),
				axis.title.y = element_blank())

sp3 <- ggplot(wR10, aes(x = -log10(P.focus), y = -log10(P.val)))+
			geom_point()+
			annotate("text", x = 95, y = 95, colour = "darkgreen", size = 7, label = "127\n\n23")+
			xlab(bquote(-log[10](P)~R5))+
			ylab(bquote(-log[10](P)~R10))+
			ggtitle("Empirical")+
			geom_abline( colour = "red")+
			geom_vline(xintercept = -log10(5e-8), colour = "blue")+
			theme_cowplot()

sp3z <- ggplot(wR10, aes(x = -log10(P.focus), y = -log10(P.val)))+
			geom_point()+
			xlim(c(0,25))+
			ylim(c(0,50))+
			xlab(bquote(-log[10](P)~R5))+
			# ylab("-log10(P) FinnGen R10")+
			ggtitle("zoomed in")+
			geom_abline( colour = "red")+
			geom_vline(xintercept = -log10(5e-8), colour = "blue")+
			theme_cowplot()+
			theme(axis.title.y = element_blank())

sp <- plot_grid(sp1, sp1z, sp2, sp2z, sp3, sp3z, nrow = 3)
sp

ggsave("../figures/FigSXX_validation_SNPs_P_plot.png", sp, bg = "white", height = 8, width = 7)

# Now the same but for our validation snps







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
#  [1] coloc_5.2.3        fpc_2.2-11         annotSnpStats_0.99 snpStats_1.52.0    Matrix_1.6-5       survival_3.5-8     cupcake_0.1.0.0    IMDtools_1.0.0    
#  [9] magrittr_2.0.3     data.table_1.15.4 

# loaded via a namespace (and not attached):
#  [1] viridis_0.6.5       utf8_1.2.4          generics_0.1.3      class_7.3-22        robustbase_0.99-2   lattice_0.22-6      grid_4.3.3          R.oo_1.26.0        
#  [9] plyr_1.8.9          jsonlite_1.8.8      R.utils_2.12.3      nnet_7.3-19         reshape_0.8.9       mclust_6.1          gridExtra_2.3       mixsqp_0.3-54      
# [17] fansi_1.0.6         kernlab_0.9-32      viridisLite_0.4.2   scales_1.3.0        modeltools_0.2-23   cli_3.6.2           rlang_1.1.3         crayon_1.5.2       
# [25] R.methodsS3_1.8.2   munsell_0.5.1       splines_4.3.3       susieR_0.12.35      tools_4.3.3         flexmix_2.3-19      parallel_4.3.3      dplyr_1.1.4        
# [33] colorspace_2.1-0    ggplot2_3.5.0       BiocGenerics_0.48.1 vctrs_0.6.5         R6_2.5.1            matrixStats_1.3.0   stats4_4.3.3        lifecycle_1.0.4    
# [41] zlibbioc_1.48.0     MASS_7.3-60.0.1     irlba_2.3.5.1       cluster_2.1.6       pkgconfig_2.0.3     pillar_1.9.0        gtable_0.3.4        glue_1.7.0         
# [49] Rcpp_1.0.12         DEoptimR_1.1-3      tibble_3.2.1        tidyselect_1.2.1    prabclus_2.3-3      compiler_4.3.3      diptest_0.77-1     