library(magrittr)
library(microbenchmark)
library(ggplot2)
library(data.table)
library(cowplot)

setDTthreads(20)
setwd("/home/gr440/rds/rds-cew54-basis/Projects/myositis-IMD/code")


# We can either run the following

## get some genotypes to sample from - quicker than rbinom multiple times
geno=list(rbinom(10000,2,0.3),
          rbinom(10000,2,0.25))
names(geno)=as.character(c(0.3,0.25))


microbenchmark(rbinom(20000,2,0.3),
               sample(geno[[1]], 20000, replace=TRUE))
## Unit: microseconds
##                                      expr      min        lq      mean   median
##                     rbinom(20000, 2, 0.3)  628.746  691.9415  774.2519  790.339
##  sample(geno[[1]], 20000, replace = TRUE) 1030.256 1119.2150 1243.1151 1265.715
##        uq      max neval cld
##   832.432  853.787   100  a 
##  1334.623 1378.719   100   b

y1=rep(c(1,0),c(500,20000))
y2=rep(c(1,0),c(500,20000))  %>% c(y1, .)

f=function(n1=500, m1=20000, # number of samples in round 1
           n2=1000, m2=40000, # number of samples in round 2
           fcase=0.3, fctl=0.25) { # allele freq in case & ctl
    g1=c(rbinom(n1,2,fcase), rbinom(m1,2,fctl)) # genotype in round 1
    g2=c(g1, rbinom(n2-n1,2,fcase), rbinom(m2-m1,2,fctl)) # genotype in round 2
    ## g1=c(sample(geno[[ as.character(fcase) ]], n1, replace=TRUE),
    ##      sample(geno[[ as.character(fctl) ]], m1, replace=TRUE)) # genotype in round 1
    ## g2=c(g1,
    ##      sample(geno[[ as.character(fcase) ]], n2-n1, replace=TRUE),
    ##      sample(geno[[ as.character(fctl) ]], m2-m1, replace=TRUE)) # genotype in round 2
    p1=t.test(g1 ~ y1)$p.value
    p2=t.test(g2 ~ y2)$p.value
    c(p1=p1,p2=p2)
}

set.seed(42)
alts=replicate(1000,f())  %>% t() # standard
alts2=replicate(2000,f(fcase=0.27))  %>% t() # lower power
nulls=replicate(10000,f(fcase=0.25))  %>% t()

dt=as.data.table(rbind(alts,alts2,nulls))
dt[,condition:=rep(c("high power","low power","null"), c(nrow(alts), nrow(alts2), nrow(nulls)))]
dt[,direction:=ifelse(p2 < p1, "down","up")]
dt[,change:=ifelse(p2 < p1, "more significant","less significant")]
dt[,sig1:=ifelse(p1 < 5e-2, "p < 0.05", "p >= 0.05")]
with(dt, table(condition, direction, sig1))
save(dt, file="../data/Simulation_data.RData")

# Or simply load the results! :)
load("../data/Simulation_data.RData")

dt[,.(n=sum(p1 < 0.05)),by="condition"]

theme_set(theme_cowplot(font_size=16))

## fig 1
ndt=dt[,.(n=.N,x=median(-log10(p1)),y=median(-log10(p2))),by=c("change","condition")]
ndt[,pc:=sprintf("%4.1f%%",100*n/sum(n)),by="condition"]
ndt[change=="less significant",x:=6]
ndt[change=="less significant",y:=3]
ndt[change=="more significant",x:=3]
ndt[change=="more significant",y:=8]

s1 <- ggplot(dt, aes(x=-log10(p1),y=-log10(p2),col=change)) +
        geom_point() +
        geom_label(aes(x=x,y=y,label=pc),data=ndt,size=7,show.legend=FALSE) +
        labs(x="Round 1: -log10(p)",y="Round 1 + 2: -log10(p)") +
        background_grid() +
        geom_abline(width=2) +
        theme(legend.position="bottom",
              legend.text = element_text(margin = margin(l = 4, r = 10, unit = "pt")),
              strip.background = element_rect(color="black", fill="#ffffff", size=1.5, linetype="solid")) +
        facet_grid(. ~ condition)
s1

## fig 2
ndt=dt[,.(n=.N,x=median(-log10(p1)),y=median(-log10(p2))),by=c("change","condition","sig1")]
ndt[,pc:=sprintf("%4.1f%%",100*n/sum(n)),by=c("sig1","condition")]
ndt[change=="less significant",x:=6]
ndt[change=="less significant",y:=3]
ndt[change=="more significant",x:=4]
ndt[change=="more significant",y:=8]

s2 <-  ggplot(dt[sig1=="p < 0.05"], aes(x=-log10(p1),y=-log10(p2),col=change)) +
          geom_point() +
          geom_label(aes(x=x,y=y,label=pc),data=ndt[sig1 == "p < 0.05"],size=7,show.legend=FALSE) +
          labs(x="Round 1: -log10(p)",y="Round 1 + 2: -log10(p)") +
          background_grid() +
          geom_abline(width=2) +
          theme(legend.position="none",
              strip.background = element_rect(color="black", fill="#ffffff", size=1.5, linetype="solid")) +
          facet_grid(sig1 ~ condition)
s2

## fig 3
sdt=dt[p1 < 0.05]
probs0=as.vector(table(sdt$condition))
s=function(high=.5,low=.25,null=.25) {
    probs1=c("high power"=high,"low power"=low,"null"=null)
    ## divide by observed freq
    probs=probs1/probs0
    probs=probs/sum(probs) # normalise
    idx=sample(1:nrow(sdt), prob=probs[sdt$condition], replace=TRUE)
    c(probs1,m=with(sdt[idx], mean(p2 < p1)))
}

fp=c(0.01,(1:9)/10,0.99)
m11=sapply(fp, function(i) replicate(100,s((1-i)/2,(1-i)/2,i),simplify=FALSE))  %>% do.call("rbind",.)  %>% as.data.table()
m21=sapply(fp, function(i) replicate(100,s(2*(1-i)/3,(1-i)/3,i),simplify=FALSE))  %>% do.call("rbind",.)  %>% as.data.table()
m12=sapply(fp, function(i) replicate(100,s((1-i)/3,2*(1-i)/3,i),simplify=FALSE))  %>% do.call("rbind",.)  %>% as.data.table()

m11[,ratio:="1:1"]
m21[,ratio:="2:1"]
m12[,ratio:="1:2"]
m=rbind(m11,m21,m12)

## figure 3
s3 <- ggplot(m, aes(x=null,y=m,col=ratio)) + geom_point() + geom_smooth(method="lm",se=FALSE) +
          background_grid() +
        #   ggtitle("High:low power dataset ratio")+
          scale_colour_discrete("Hi:lo power dataset ratio") +
          labs(x="False Positive Rate", y="p2 < p1 proportion") +
          theme(legend.position = "bottom",
          legend.justification = "center")+
          guides(colour = guide_legend(title.position="top", title.hjust = 0.5))
s3

# plot_grid(
#   plot_grid(s1, s2, nrow = 1, ncol = 2, labels = "AUTO"),
#   plot_grid(NULL, s3, NULL, nrow = 1, rel_widths = c(0.5, 1, 0.5), labels = "C"),
#   nrow = 2
# )

## Safe Figure S3

sxx <- plot_grid(s1,s2,s3, nrow = 3, labels = "AUTO")
ggsave2("../figures/FigS3_simulation.png", sxx, height = 12, width = 8, bg = "white")
