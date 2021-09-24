### Making Forest plot from cluster for myositis paper

library(ggplot2)
library(cowplot)
library(magrittr)
library(data.table)
guide <- fread("../data/decode-trait-shorttrait.csv")
metadata <- fread("../data/Metadata_20210812-v1.tsv") %>% .[ First_Author != "PanUKBB"] # We're considering only Neale's UKBB
load("../data/delta_clustered_data.RData")

# Data in delta_clustered_data.RData comes with 'shorttrait', which is a custom naming system. To start from a common ground, I will restore the original names.
# Here we're simply merging UKBB traits with metadata by their "Trait_long", and non-UKBB by the Trait, so we have a unified file with the original Trait and
# Trait_long columns.

tmsub <- msub[grepl("UKBB",shorttrait)] %>% merge(., guide) %>% merge(., metadata[,c("Trait", "First_Author", "Trait_long")], by.x ="Trait", by.y = "Trait_long", all.x = TRUE) %>% setnames(. , c("Trait", "Trait.y"), c("Trait_long", "Trait"))
tmsub2 <- msub[!grepl("UKBB",shorttrait)] %>% merge(., guide) %>% merge(., metadata[,c("Trait", "First_Author", "Trait_long")], by = "Trait", all.x = TRUE)
msub <- rbindlist(list(tmsub, tmsub2), use.names = TRUE)

# We can now create Label from Trait
make.labels <- function(x){
  x[,Label:=Trait_long] %>% 
  .[, Label:=gsub(" \\(UKBB\\)", "", Label)] %>% 
  .[, Label:=gsub(" \\(FinnGen\\)", "", Label)] %>% 
  .[, Label:=gsub(" \\(FG\\)", "", Label)] %>% 
  .[, Label:=paste0(Label, " / ", First_Author)] %>% 
  .[, Label:=gsub("Neale", "UKBB", Label)] %>% 
  .[, Label:=gsub("Crohn disease \\(strict definition, all UC cases excluded\\)", "Crohn's, strict", Label)] %>% 
  .[, Label:=gsub("Eosinophilic Granulomatosis with Polyangiitis", "EGPA", Label)]
    
  
}
make.labels(msub)


PCorder <- c("PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10","PC11","PC12","PC13")

cs <- ggplot() +
  ## shape=shorttrait=="overall",
  ## fill=shorttrait=="overall")) +
  geom_vline(xintercept=0,col="grey") +
  geom_vline(aes(xintercept=delta),data=msum) +
  geom_rect(aes(xmin=lci,xmax=uci,ymin=-Inf,ymax=Inf),colour=NA,data=msum,alpha=0.5,fill="grey") +
  geom_pointrange(aes(y=Label,x=delta,xmin=lci,xmax=uci,col=factor(cl)),shape=23,data=msub) +
  ## scale_shape_manual(values=c("TRUE"=23,"FALSE"=21)) +
  scale_fill_manual(values=c("TRUE"="black","FALSE"="white")) +
  scale_colour_manual(values=ann$colors$call) +
  labs(y="") +
  background_grid(major="y",size.major=0.1) +
  theme_minimal()+
  theme(legend.position="none", axis.text = element_text(size =8), axis.title.x =  element_text(size =11), strip.text.y = element_text(size = 11, angle = 0), strip.text.x = element_text(size = 11), panel.grid = element_blank()) +
  facet_grid(cl ~ factor(PC, levels = PCorder),space="free",scales="free")
cs

ggsave(filename = "../figures/clustering_forest_v1.png", plot = cs, units = "mm", width = 480, height = 300, dpi = 300, bg="white")
