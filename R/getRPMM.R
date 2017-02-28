getRPMM <- function(gset,mvps,
                    hm.main = NULL, row.main="EAC MVPs (N = 1033)",showmainlegend=TRUE, 
                    showannolegend=TRUE, mainlegendname="Beta-value",
                    annolegendname="rpmm cluster",mvp.orderlist=NULL){
  require("ComplexHeatmap")
  require("RPMM")
  bval <- as.matrix(getBeta(gset[rownames(gset) %in% mvps,])) 
  message("Running RPMM..")
  rpmm <- blcTree(t(bval), verbose = 0,
                  splitCriterion = blcSplitCriterionLevelWtdBIC)
  rpmmClass <- blcTreeLeafClasses(rpmm)
  a <- lapply(levels(rpmmClass), function(k) {colnames(bval)[which(rpmmClass == k)]}) # assign sample IDs from bval colnames
  names(a) <- levels(rpmmClass)
  m <- matrix(as.character(unlist(a)))
  # clu matrix dimensions contingent on cluster count (ie. 4x clusters shown below)
  x <- names(a)
  message("N = ",length(a)," clusters, continuing if 4 clusters present..")
  if(length(x)==4){
    message("4 clusters, returning cluster assignments...")
    clu <- matrix(c(rep(names(a)[[1]],length(a[[1]])),
                    rep(names(a)[[2]], length(a[[2]])),
                    rep(names(a)[[3]], length(a[[3]])),
                    rep(names(a)[[4]], length(a[[4]]))))
    sampleid <- unlist(a)
    returndf <- data.frame(sampleid,clu); rownames(returndf)<-NULL
  } else{
    return(message("ERROR: !n.RPMM.clusters==4"))
  }
  returndf$mlvl <- NA
  for(i in 1:4){
    whichgrpi <- which(colnames(bval) %in% as.character(returndf[returndf[,2]==x[i],1]))
    returndf[returndf[,2]==x[i],]$mlvl <- mean(rowMeans(bval[,whichgrpi]))
  }; returndf <- returndf[order(returndf$mlvl),] 
  mlvl <- as.factor(unique(returndf$mlvl))
  #mgrplvl <- levels(returndf$mlvl)
  col <- c(); 
  mcol <- c("MM","LM","IM","HM"); 
  #mcol<-mcol[rank(unique(returndf$mlvl))]
  for(i in 1:4){
    repi <- nrow(returndf[returndf$mlvl==mlvl[i],])
    col<-c(col,rep(mcol[i],repi))
  }; returndf$mcol <- col
  returndf[,1] <- as.character(returndf[,1]); 
  tbval <- as.data.frame(t(bval[,order(match(colnames(bval),returndf[,1]))]));
  
  if(!is.null(mvp.orderlist)){
    if(!length(intersect(mvp.orderlist,as.character(colnames(tbval))))==ncol(tbval)){
      return("ERROR: mvp.orderlist does not match bval probes.")
    }
    message("ordering probes on provided mvp list..")
    tbval <- tbval[,order(match(colnames(tbval),as.character(mvp.orderlist)))]
    if(!identical(colnames(tbval),mvp.orderlist)){
      return(message("ERROR: Matching order of probes to mvp list unsuccessful."))
    } else{
      # mvp IDs to be stored in return list
      return.mvplist = mvp.orderlist
    }
  } else{
    tbval <- tbval[,order(colMeans(tbval))]
    # mvplist to be stored in returned list
    message("ordering probes on mean values..")
    return.mvplist <- colnames(tbval)
  }
  if(identical(rownames(tbval),returndf[,1])){
    tbval$mcol <- returndf$mcol
  } else{
    return("ERROR: could not match returndf and tbval samples!")
  }
  # make heatmap, arguments: hm.main, row.main, mainlegendname, annoname, showmainlegend, showannolegend=FALSE
  hmval <- t(tbval[1:ncol(bval),1:nrow(bval)])
  message("ordering samples on mean methylation..")
  hmval <- hmval[,order(colMeans(hmval))]; returndf <- returndf[order(match(returndf[,1],colnames(hmval))),]
  if(!identical(returndf$sampleid,as.character(colnames(hmval)))){
    return(message("ERROR: reordering hmval map on sample means failed.."))
  }
  
  colside <- as.character(returndf$mcol)
  message("Relation of mean mgrp methylation to grp:\n",print(table(colside,returndf$mlvl)),"\n")
  ha = HeatmapAnnotation(df = data.frame(rpmmgrp=as.character(colside)), 
                         col = list(rpmmgrp = c("MM" =  "lightblue", "LM" = "gray","IM" = "coral","HM" = "yellow")),
                         show_legend = showannolegend,
                         name = annolegendname)
  breaks=seq(0,1,0.01)
  hmcol = colorRamp2(breaks,colorRampPalette(c("darkblue","yellow"))(n=length(breaks)))
  # heatmap to be stored in returned list
  hm1 <- Heatmap(hmval,
                 col=hmcol,
                 cluster_columns=FALSE,
                 column_dend_reorder = FALSE,
                 row_dend_reorder = FALSE,
                 cluster_rows = FALSE,
                 show_heatmap_legend = showmainlegend,
                 show_row_names = FALSE,
                 show_column_names = FALSE,
                 top_annotation = ha,
                 name=mainlegendname,
                 column_title = hm.main,  
                 row_title = row.main)
  # return as list
  return1 <- list(returndf,hm1,return.mvplist); 
  names(return1) <- c("rpmm_mgrpdf","heatmap","mvplist.hmorder")
  return(return1)
}
