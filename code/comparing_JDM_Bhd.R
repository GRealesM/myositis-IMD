# This is a short snippet to compare Bh distances results across JDM/DM and other IMD.

## Set working directory
setwd("/home/gr440/rds/rds-cew54-basis/Projects/myositis-IMD/code")

## Load packages and required
library(data.table)
setDTthreads(20)
library(magrittr)


bhd <- fread("../data/bh_dist.tsv")

bhd[grepl("Juvenile Dermatomyositis", T1) & T2 %in% c("Type 1 diabetes", "Juvenile Idiopathic Arthritis", "Rheumatoid arthritis") ]
#                                     T1                            T2 bhat.dist
# 1:   Juvenile Dermatomyositis (Miller)               Type 1 diabetes  7.396167
# 2:   Juvenile Dermatomyositis (Miller) Juvenile Idiopathic Arthritis  3.543378
# 3:   Juvenile Dermatomyositis (Miller)          Rheumatoid arthritis  4.851606
# 4: Juvenile Dermatomyositis (Rothwell)               Type 1 diabetes 14.617242
# 5: Juvenile Dermatomyositis (Rothwell) Juvenile Idiopathic Arthritis  6.583247
# 6: Juvenile Dermatomyositis (Rothwell)          Rheumatoid arthritis  8.318788
bhd[grepl("^Dermatomyositis", T1) & T2 %in% c("Type 1 diabetes", "Juvenile Idiopathic Arthritis", "Rheumatoid arthritis") ]
#                            T1                            T2 bhat.dist
# 1:   Dermatomyositis (Miller)               Type 1 diabetes 11.117949
# 2:   Dermatomyositis (Miller) Juvenile Idiopathic Arthritis  4.341990
# 3:   Dermatomyositis (Miller)          Rheumatoid arthritis  4.197110
# 4: Dermatomyositis (Rothwell)               Type 1 diabetes 13.780748
# 5: Dermatomyositis (Rothwell) Juvenile Idiopathic Arthritis  4.690899
# 6: Dermatomyositis (Rothwell)          Rheumatoid arthritis  4.364159