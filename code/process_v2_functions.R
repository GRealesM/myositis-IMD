
calc_psm=function(x,burn=0.5) {
  n=nrow(x)
  if(burn>0)
    x=x[ (burn*n) : n, , drop=FALSE]
  if(any(is.na(x)))
    x=x[ apply(!is.na(x),1,all), ]
  unq=unique(as.vector(x))
  ## print(unq)
  m=matrix(0,ncol(x),ncol(x))
  for(k in unq) {
    xk=matrix(as.numeric(x==k),nrow(x),ncol(x))
    ## m=m + t(xk) %*% xk
    m=m + crossprod(xk)
  }
  psm=m/nrow(x)
  rownames(psm)=rownames(obsData)
  psm
}

library(mcclust)
library(pheatmap)
show_plots=interactive()

library(bayesplot); color_scheme_set("viridis")
library(mcmcplots)
library(data.table)
library(abind)
library(magrittr)

read_alpha=function(dirs) {
  files=file.path(dirs,"alpha.csv")
  use=file.exists(files)
  alpha=lapply(which(use), function(i) {
    message(dirs[i])
    scan(files[i],what="") %>% as.numeric() %>% setdiff(., NA)
  })
  nmin=sapply(alpha,length) %>% min()
  alpha %<>% lapply(., "[", 1:nmin) %>% do.call("cbind",.)
  colnames(alpha)=dirs[use]
  alpha
}

read_allocs=function(dirs) {
  files=file.path(dirs,"clusterAllocations.csv")
  use=file.exists(files)
  allocs=lapply(which(use), function(i) {
    message(i)
    paste("ls -lh",files[i]) %>% system(., intern=TRUE) %>% print()
    fread(files[i], fill=TRUE)
  })
  nmin=sapply(allocs,nrow) #%>% min()
  allocs %<>% lapply(., function(x) as.matrix(x[1:nrow(x),])) #[1:nmin,]))
  names(allocs)=dirs[use]
  allocs
  ## use=sapply(allocs,ncol)==nrow(obsData)
  ## allocs[use]
}

extract_nallocs=function(allocs) {
  nallocs = lapply(allocs, function(x) apply(x, 1, function(y) length(unique(y))))
  n=sapply(nallocs, length) %>% min()
  lapply(nallocs, function(x) x[1:n]) %>% do.call("cbind",.)
}

get_calls=function(psm0) {
  ## print(str(psm))
  maxpear(psm=psm0,max.k=ceiling(nrow(psm0)/2),method="avg")
}

get_ann_colors=function(calls,verbose=TRUE) {
  palette= c("#771155", "#AA4488", "#CC99BB", "#114477", "#4477AA", "#77AADD", "#117777", "#44AAAA", "#77CCCC", "#117744", "#44AA77", "#88CCAA", "#777711", "#AAAA44", "#DDDD77", "#774411", "#AA7744", "#DDAA77", "#771122", "#AA4455", "#DD7788") %>%
    matrix(.,7,3,byrow=TRUE) %>%
    as.vector()
  ann=data.frame(row.names=rownames(obsData),
                 cluster=factor(calls),
                 myositis=factor(grepl("myogen|Miller|ph770|DERMATOPOLY|MYO",rownames(obsData))))
  if(verbose)
      print(subset(ann, ann$call %in% ann$call[ann$myositis=="TRUE"] ))
  ann_colors=list(myositis=c("TRUE"="cornflowerblue","FALSE"="white"))
  ncalls=length(unique(calls))
  if(ncalls <= length(palette))
    ann_colors$call=structure(palette[1:ncalls], names=unique(calls))
  list(ann=ann,colors=ann_colors)
}

shorten_names=function(x) {
  x %>%
    sub("_FinnGenR5_1","",.) %>%
    sub("_[0-9]+_1$","",.) %>%
    sub("_up_1","",.) %>%
    sub("_FinnGen"," FG",.) %>%
    sub("colitis\\/not","colitis UKBB",.) %>%
    sub("(UKBB)","UKBB",.) %>%
    sub(".*\\/","",.)
}

ladder_plot=function(.psm,.m,expandgroups=NULL,subgroups=NULL,savefile=NULL) {
  calls=get_calls(.psm)$cl
  ## m2=rbind(m,msum,fill=TRUE)
  ## lev=sort(unique(m$shorttrait)) %>% setdiff(., "overall") %>% c(., "overall") %>% rev()
  ## m2[,shorttrait:=factor(shorttrait,levels=lev)]
  msub=copy(m)
  msub[,cl:=calls[shorttrait]]
  if(!is.null(subgroups)) {
    msub=msub[cl %in% subgroups]
  }
  msum=msub[,.(shorttrait="overall",
            n=.N,
            delta=sum(delta / var)/sum(1/var),
            var=1/sum(1/var)),
            by=c("PC","cl")]
  msum[,lci:=delta-1.96*sqrt(var)]
  msum[,uci:=delta+1.96*sqrt(var)]
  msub$n=NA
  if(!is.null(expandgroups)) {
    msub=rbind(msub[ cl %in% expandgroups ],
               msum[ !(cl %in% expandgroups)])
    msub[shorttrait=="overall",shorttrait:=paste(n,"traits")]
  }
  ann=get_ann_colors(calls,verbose=FALSE)
  if(!is.null(savefile))
    save(msub, msum, ann, file=savefile)
  ggplot() +
                 ## shape=shorttrait=="overall",
                 ## fill=shorttrait=="overall")) +
    geom_vline(xintercept=0,col="grey") +
    geom_vline(aes(xintercept=delta),data=msum) +
    geom_rect(aes(xmin=lci,xmax=uci,ymin=-Inf,ymax=Inf),colour=NA,data=msum,alpha=0.5,fill="grey") +
    geom_pointrange(aes(y=shorttrait,x=delta,xmin=lci,xmax=uci,col=factor(cl)),shape=23,data=msub) +
    ## scale_shape_manual(values=c("TRUE"=23,"FALSE"=21)) +
    scale_fill_manual(values=c("TRUE"="black","FALSE"="white")) +
    scale_colour_manual(values=ann$colors$call) +
    labs(y="") +
    background_grid(major="y",size.major=0.1) +
    theme(legend.position="none") +
    facet_grid(cl ~ PC,space="free",scales="free")
}



plot_psm=function(PSM,n=1,what=c("psm","data"),nm="",transpose=FALSE,...) {
  .psm=if(is.list(PSM)) {
    PSM[[n]]
  } else {
    PSM
  }
  if(nm=="" && is.list(PSM))
    nm=names(PSM)[n]

  calls=get_calls(.psm)
  ann=get_ann_colors(calls$cl,verbose=TRUE)
  what=match.arg(what)
  cl = pheatmap:::cluster_mat(t(.psm),
                              method = "average",
                              distance = "euclidean" )
  #tree_col = clustering_callback(tree_col, t(psm))
  ## cl=hclust(as.dist(1-psm))
  obj=switch(what,
             psm=.psm,
             data=obsData)
  if(what=="data") {
    obj=obj/matrix(apply(abs(obj),2,max),nrow=nrow(obj),ncol=ncol(obj),byrow=TRUE)
    z=abs(obsData/sqrt(obsVars))
    obj[z<1.96]=0
    paletteLength <- 100
    myColor <- colorRampPalette((RColorBrewer::brewer.pal(n = 7, name = "RdYlBu")))(100)
    myColor[50]="#efefef"
                                        #colorRampPalette(c("yellow", "white", "blue"))(paletteLength)
                                        # length(breaks) == length(paletteLength) + 1
                                        # use floor and ceiling to deal with even/odd length pallettelengths
    myBreaks <- c(seq(min(obj), 0, length.out=ceiling(paletteLength/2) + 1),
                  seq(max(obj)/paletteLength, max(obj), length.out=floor(paletteLength/2)))
  } else {
    myColor=colorRampPalette((RColorBrewer::brewer.pal(n = 7,
                                                          name = "Blues")))(100)
    myBreaks=NA
  }
  if(transpose) {
  o=pheatmap((obj), #silent=TRUE,
             cluster_rows=cl,
             cluster_cols=what=="psm",
             clustering_method="average",
             color=myColor,
             breaks=myBreaks,
             annotation_row=as.data.frame(ann$ann)[,c("cluster"),drop=FALSE],
             annotation_colors=ann$colors,
             main=nm,
             ...)
  } else {
   o=pheatmap(t(obj), #silent=TRUE,
             cluster_cols=cl,
             cluster_rows=what=="psm",
             clustering_method="average",
             color=myColor,
             breaks=myBreaks,
             annotation_col=as.data.frame(ann$ann)[,c("cluster"),drop=FALSE],
             annotation_colors=ann$colors,
             main=nm,
             ...)
   }
           ## filename=if(show_plots) { NA } else { sub(".RData","_psm.png",savefile) })
}

readMeans=function(dirs) {
  means=lapply(dirs, function(d) {
    dim2=ifelse(grepl("all",d), 13, 3)
    tmp=scan(file.path(d,"clusterMeans.csv"),what="",sep="\n")
    v=lapply(tmp, function(x) matrix( scan(text = x, what = numeric(), quiet = TRUE, sep=","), ncol=dim2))
  })
}

readVars=function(dirs) {
  vars=lapply(dirs, function(d) {
    dim2=ifelse(grepl("all",d), 13, 3)
    tmp=scan(file.path(d,"clusterVars.csv"),what="",sep="\n")
    v=lapply(tmp, function(x) matrix( scan(text = x, what = numeric(), quiet = TRUE, sep=","), ncol=dim2))
  })
}
