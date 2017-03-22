# FUNCTION: exprmethyCor
# AUTHOR: SKM
# DESCRIPTION: correlation of gene expression and region methylation
# DEPENDANCIES: minfi
# ARGUMENTS
# cormethod = method to be used, options in helpfile of cor.test
# nprobecut = minimum probes mapping to gene of interest
# nsampcut = minimum sample count cutoff for available expr data (after filtering inf values) 
# annohm450 = HM450 annotation (all probes should also be available in gset)
# probelist = probes to use in investigation
# testgenes = list or vector of character gene identifiers to test for correlation
# totlen = end index on testgenes, must be less than length(testgenes)
# startlen = start index on testgenes, must be greater or equal to 1 
# exprset = counts or log2fc transformed expression data (colnames should match colnames of gset)
# gset = genomic ratio set containing beta values (colnames should match colnames of exprset)
# exprnames = optionally take subset 

exprmethyCor <- function(cormethod="spearman",nprobecut=3,nsampcut=15,
                         annohm450=anno450hm,
                         testgenes=testgenes.hm,totlen=length(testgenes.hm),startlen=1,
                         exprset,gset){
  # error messages
  if(totlen>length(testgenes)|totlen<startlen){
    return(message("ERROR: totlen is longer than length of testgenes or less than startlen."))
  }
  if(startlen>totlen){
    return(message("ERROR: startlen is greater than totlen."))
  }
  if(!length(intersect(colnames(exprset),colnames(gset)))>1){
    return(message("ERROR: not enough shared samples in methylation and expression sets to perform correlations!"))
  }
  
  testcor.hm <- as.data.frame(matrix(nrow=totlen,ncol=7)) 
  colnames(testcor.hm) <- c("gene","nsample.expr","n450probes","cor.estimate","cor.praw",
                            "meanexpr.log2fc","meanmethy")
  
  for(i in startlen:totlen){
    genei <- testgenes[i]; testcor.hm[i,1] <- genei # col1
    annoi <- annohm450[grep(genei,annohm450$UCSC_RefGene_Name),]
    
    expri <- exprset[rownames(exprset)==genei,]; expri <- expri[!is.infinite(expri)]
    exprnames.i <- names(expri)
    cond1 <- nrow(annoi)>=nprobecut; testcor.hm[i,3] <- nrow(annoi) # col3 
    cond2 <- length(expri)>=nsampcut; testcor.hm[i,2] <- length(expri) # col2
    if(cond1 & cond2){
      betai <- as.data.frame(getBeta(gset[rownames(gset) %in% rownames(annoi),
                                          colnames(gset) %in% exprnames.i]));
      betai <- betai[,order(match(colnames(betai),exprnames.i))]; 
      #message("sample data matched order: ", identical(colnames(betai),exprnames.i))
      
      methy <- colMeans(betai,na.rm=TRUE)
      cori <- cor.test(methy,expri,method=cormethod)
      testcor.hm[i,4] <- round(cori$estimate,5) # col4
      testcor.hm[i,5] <- round(cori$p.value,5) # col5
      testcor.hm[i,6] <- round(mean(expri),5) # col6
      testcor.hm[i,7] <- round(mean(methy),5) # col7
    } 
    if(!cond1){
      message("min hm450 probe count ",nprobecut," not met for gene ",i,", skipping cortest...")  
    }
    if(!cond2){
      message("min expression sample count ",nsampcut," not met for gene ",i,", skipping cortest...")  
    }
    message("completed ",i," of ",totlen," total (",round(100*(i/totlen),3),"%)")
    
  }
  return(testcor.hm)
  
}
