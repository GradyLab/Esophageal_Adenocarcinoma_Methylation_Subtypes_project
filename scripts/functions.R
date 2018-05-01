# script: a list of useful functions used in the manuscript code.

#======================================
#======================================
#======================================
# FUNCTION: mRgnSummary
# AUTHOR: SKM
# PURPOSE: Function that summarizes mean gene region methylation for pre-filtering.
# DEPENDANCIES: minfi
# ARGUMENTS: 
# testgenes = list of pre-filtered genes of interest
# annohm450 = annotation file, rownames contain all probes to search in gset
# testvar = character vector of sample groups, corresponding to colnames/sample order in gset
# gset = GenomicRatioSet or similar class, output from minfi
# nprobecutoff = minimum number of probes required to calculate a mean
# nacut = minimum fraction of all samples that must have available measurements for mean to be calculated
# roundnum = digits to round to (calls function 'round') 
# probesfilt = either "sd" or "var"; NULL defaults to taking a mean of all available probes in region (by sample group)
# varfraction = if probesfilt is "var", this sets fraction cutoff value for variance (by sample group)
# sdfraction = if probesfilt is "sd", this sets fraction cutoff for sd (by sample group)
# maxrange = index in testgenes to stop at, less than total length of testgenes

mRgnSummary <- function(testgenes,annohm450,testvar,gset,nprobecutoff=3,
                        nacut,roundnum=5,probesfilt=NULL,varfraction=1,
                        sdfraction=1,maxrange=length(testgenes)){
  
  # eval levels in samplegroups vector
  level1 <- as.character(levels(as.factor(testvar))[1]); level2 <- as.character(levels(as.factor(testvar))[2])
  
  # make return object
  return.methylvl <- matrix(nrow=maxrange,ncol=8); 
  colnames(return.methylvl)<-c("gene","nprobes",
                               paste0(level1,".meanmethy"),paste0(level2,".meanmethy"),
                               paste0(level1,".nprobes"),paste0(level2,".nprobes"),
                               paste0(level1,".meanmethy.sd"),paste0(level2,".meanmethy.sd"))
  if(!is.null(probesfilt)){
    if(probesfilt=="var"){
      colnames(return.methylvl)[3] <- paste0(colnames(return.methylvl)[3],".varfilt")
      colnames(return.methylvl)[4] <- paste0(colnames(return.methylvl)[4],".varfilt")
    }
    if(probesfilt=="sd"){
      colnames(return.methylvl)[3] <- paste0(colnames(return.methylvl)[3],".sdfilt")
      colnames(return.methylvl)[4] <- paste0(colnames(return.methylvl)[4],".sdfilt")
    } 
  }
  #
  for(i in 1:maxrange){
    genei <- paste0("(^|;)",testgenes[i],"(;.|$)"); return.methylvl[i,1] <- genei # genename col1
    probesi <- rownames(annohm450[grep(genei, annohm450$UCSC_RefGene_Name),]); return.methylvl[i,2]<-length(probesi) # nprobes col2
    betamethy.tn <- getBeta(gset[rownames(gset) %in% probesi,])
    has.na.probe <- mean(apply(betamethy.tn,2,function(x){length(x[is.na(x)])}))
    
    # process probe distribution filter options
    if(nrow(betamethy.tn)>=nprobecutoff){
      if(!is.null(probesfilt)){
        if(probesfilt=="var"){
          # lvl1
          level1.probemean <- rowMeans(betamethy.tn[,which(testvar==level1)],na.rm=TRUE); level1.var <-abs(var(as.numeric(level1.probemean)));
          level1.minmax <- c(mean(level1.probemean)-(varfraction*level1.var),mean(level1.probemean)+(varfraction*level1.var))
          lvl1retainprobe <- names(level1.probemean[which(level1.probemean>=level1.minmax[1] & level1.probemean<=level1.minmax[2])])
          if(length(lvl1retainprobe)>=nprobecutoff & has.na.probe<nacut){
            return.methylvl[i,3] <- round(mean(colMeans(betamethy.tn[rownames(betamethy.tn) %in% lvl1retainprobe,
                                                                     which(testvar==level1)])),roundnum)
          } else{
            return.methylvl[i,3] <- "NA.filt"
          }
          
          #lvl2
          level2.probemean <- rowMeans(betamethy.tn[,which(testvar==level2)],na.rm=TRUE); level2.var <-abs(var(as.numeric(level2.probemean)));
          level2.minmax <- c(mean(level2.probemean)-(varfraction*level2.var),mean(level2.probemean)+(varfraction*level2.var))
          lvl2retainprobe <- names(level2.probemean[which(level2.probemean>=level2.minmax[1] & level2.probemean<=level2.minmax[2])])
          
          if(length(lvl1retainprobe)>=nprobecutoff & has.na.probe<nacut){
            return.methylvl[i,4] <- round(mean(colMeans(betamethy.tn[rownames(betamethy.tn) %in% lvl2retainprobe,
                                                                     which(testvar==level2)])),roundnum)
          } else{
            return.methylvl[i,4] <- "NA.filt"
          }
          
        }
        if(probesfilt=="sd"){
          
          # level1
          level1.probemean <- rowMeans(betamethy.tn[,which(testvar==level1)],na.rm=TRUE); level1.sd <-abs(sd(as.numeric(level1.probemean)));
          return.methylvl[i,7] <- round(level1.sd,roundnum) # col7
          level1.minmax <- c(mean(level1.probemean)-(sdfraction*level1.sd),mean(level1.probemean)+(sdfraction*level1.sd))
          lvl1retainprobe <- names(level1.probemean[which(level1.probemean>=level1.minmax[1] & level1.probemean<=level1.minmax[2])])
          return.methylvl[i,5] <- length(lvl1retainprobe) # col5
          if(length(lvl1retainprobe)>=nprobecutoff & has.na.probe<nacut){
            return.methylvl[i,3] <- round(mean(colMeans(betamethy.tn[rownames(betamethy.tn) %in% lvl1retainprobe,
                                                                     which(testvar==level1)])),roundnum)
          } else{
            return.methylvl[i,3] <- "NA.filt"
          }
          
          # level2
          level2.probemean <- rowMeans(betamethy.tn[,which(testvar==level2)],na.rm=TRUE); level2.sd <-abs(sd(as.numeric(level2.probemean)));
          return.methylvl[i,8] <- round(level2.sd,roundnum) # col8
          level2.minmax <- c(mean(level2.probemean)-(sdfraction*level2.sd),mean(level2.probemean)+(sdfraction*level2.sd))
          lvl2retainprobe <- names(level2.probemean[which(level2.probemean>=level2.minmax[1] & level2.probemean<=level2.minmax[2])])
          return.methylvl[i,6] <- length(lvl2retainprobe) # col6
          if(length(lvl2retainprobe)>=nprobecutoff & has.na.probe<nacut){
            return.methylvl[i,4] <- round(mean(colMeans(betamethy.tn[rownames(betamethy.tn) %in% lvl2retainprobe,
                                                                     which(testvar==level2)])),roundnum)
          } else{
            return.methylvl[i,4] <- "NA.filt"
          }
          
        }
      } else {
        
        # level1
        return.methylvl[i,5] <- nrow(betamethy.tn) # col5
        level1.probemean <- rowMeans(betamethy.tn[,which(testvar==level1)]);
        lvl1retainprobe <- names(level1.probemean[is.na(level1.probemean)<nacut]) # apply na cutoff
        if(length(lvl1retainprobe)>=nprobecutoff){
          return.methylvl[i,3] <- round(mean(colMeans(betamethy.tn[rownames(betamethy.tn) %in% lvl1retainprobe,
                                                                   which(testvar==level1)])),roundnum)
        } else{
          return.methylvl[i,3] <- "NA.filt"
        }
        
        # level2
        return.methylvl[i,6] <- nrow(betamethy.tn) # col6
        level2.probemean <- rowMeans(betamethy.tn[,which(testvar==level2)]);
        lvl2retainprobe <- names(level2.probemean[is.na(level2.probemean)<nacut]) # apply na cutoff
        if(length(lvl2retainprobe)>=nprobecutoff){
          return.methylvl[i,4] <- round(mean(colMeans(betamethy.tn[rownames(betamethy.tn) %in% lvl2retainprobe,
                                                                   which(testvar==level2)])),roundnum)
        } else{
          return.methylvl[i,4] <- "NA.filt"
        }
        
        }
    }
    #
    message("gene ",i," of ",length(testgenes)," done (",round(100*(i/length(testgenes)),3),"%)")
    #
  }
  return(return.methylvl)
}


#======================================
#======================================
#======================================
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

#
