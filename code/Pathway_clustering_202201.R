library(data.table)
library(magrittr)
library(FastJaccard)
library(pheatmap)
library(gprofiler2)
library(mcclust)
library(mclust)

# Get data. We'll use the results from our enrichment analysis (see SNP_OTG_followup_202201.Rmd)
resgp <- readRDS("../data/enrichment_results_OTG_eQTL_GOBP.WP.KEGG.REAC.RDS")
resgp <- resgp[ source != "KEGG"] # Remove KEGG, as it's not publicly available from gprofiler's pathway GMT file, as data source licenses restrict public access to the data
RES=split(resgp, resgp$PC)



# Extract all available pathways from Gprofiler sources (except for KEGG), this will help us compare clustering when all genes associated are considered vs genes associated with our driver SNPs
pathways=scan("https://biit.cs.ut.ee/gprofiler/static/gprofiler_full_hsapiens.name.gmt", what="",sep="\n")
ss=strsplit(pathways,"\t")
pathway_ids=sapply(ss, "[[", 1)
pathway_term_names=sapply(ss, "[[", 2)
pathway_genes=lapply(ss, "[", -c(1:2))
names(pathway_genes)=pathway_term_names

# Apply Jaccard to significant pathways in each PC, using their associated genes contained in OTG and eQTL annotations
JACC=lapply(RES, function(res) {
  genes=strsplit(res$intersection,", ")
  jacc=jaccard_symlist(genes)
})

# Apply Jaccard to significant pathways in each PC, using *all* associated genes in the Gprofiler database
JACC_PATHWAY=lapply(RES, function(res) {
  jacc_pathway=jaccard_symlist( pathway_genes[ res$term_name ])
})

# Apply maxpear across clustering in JACC and JACC_PATHWAY
CL=lapply(JACC, function(jacc) maxpear(jacc,method="avg")$cl)
CL_PATHWAY=lapply(JACC_PATHWAY, function(jacc) {
  jacc[ jacc>1 ]=1 ## why???
  maxpear(jacc,method="avg")$cl })

# Apply resulting clusters to original pathways
for(i in seq_along(JACC_PATHWAY)) {
  message(i)
  RES[[i]]$cl=CL[[i]]
  RES[[i]]$cl_pathway=CL_PATHWAY[[i]]
}

# Check by how many components are pathways shared
result=rbindlist(RES)
result[,ncopies:=.N,by=c("source","term_name")]
table(result$ncopies) 



# Compare general pathway clustering vs our pathway clustering
result[,.(ari=adjustedRandIndex(cl,cl_pathway)),by="PC"] # all in 0.3 -- 0.56
# PC        ari
# 1:  PC1 0.03815554
# 2: PC12 0.03004273
# 3: PC13 0.01870101
# 4:  PC2 0.03129530
# 5:  PC3 0.01067189
# 6:  PC8 0.11885560
# 7:  PC9 0.04094321

# Order clusters by number of copies and associated p-values
result=result[order(cl, ncopies, p_value)]
options(width=200)
result[term_size >=5 & ncopies==1 & !duplicated(cl),.(cl,source,term_name,term_size,intersection_size,p_value),by="PC"][order(PC,p_value)]


result1=result[p_value < 0.01]
result1[PC %in% c("PC2","PC3","PC8","PC9","PC12"),ncopies:=.N,by=c("source","term_name")]

## most significant exclusive pathway per cluster
result1=result1[order(cl_pathway, ncopies, p_value)]
result1[!(source %in% c("HP")) & term_size>=5 & ncopies==1 & !duplicated(cl_pathway),.(cl_pathway,source,term_name,term_size,intersection_size,p_value),by="PC"][order(PC,p_value)]

## most significant pathway per cluster
result1[,pc_cl:=paste(PC,cl_pathway,sep="_")]
result1=result1[order(pc_cl, p_value)]
result2=result1[!(source %in% c("HP")) & ncopies<=2 & term_size>=5 & p_value < 0.01][!duplicated(pc_cl),.(cl_pathway,source,term_name,term_size,intersection_size,p_value,ncopies),by="PC"][order(PC,cl_pathway,p_value)][PC %in% c("PC2","PC3","PC8","PC9","PC12")] %>% as.data.frame()

split(result2, result2$PC)