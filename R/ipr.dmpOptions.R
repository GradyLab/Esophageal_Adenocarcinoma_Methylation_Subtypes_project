ipr.dmpOptions <- function(df,getDMPs="ftest",varfilt=NULL,nsigprobes=1500,dmpVariablename,
                           dmpVariancefiltTissue){
  
  if(class(df) %in% c("MethylSet","GenomicMethylSet","GenomicRatioSet","RGChannelSet")){
    require(minfi)
    if(!getDMPs %in% c("wilcox","ftest")){
      return(message("ERROR: getDMPs test must be either 'wilcox' or 'ftest'"))
    }
    
    # dmp dataframe
    dfdmp <- data.frame(as.character(rownames(df)),
                        rep(0,nrow(df)))
    colnames(dfdmp) <- c("probeID","pval.unadj")
    tissuelvl <- levels(as.factor(eval(parse(text=paste0("df$",dmpVariablename)))))
    message("Detected tissue levels are:\n",paste(tissuelvl,sep=" ",collapse=";"))
    message("Expected hypervariant tissue group is: ",dmpVariancefiltTissue)
    whichhyper <- eval(parse(text=paste0("which(df$",dmpVariablename,"=='",dmpVariancefiltTissue,"')")))
    whichhypo <- eval(parse(text=paste0("which(df$",dmpVariablename,"!='",dmpVariancefiltTissue,"')")))
    
    ###
    if(getDMPs=="wilcox"){
      message("Performing DMP test: Wilcoxon Rank Sum...")
      dfdmp[,2] <- apply(as.matrix(getBeta(df)),1,function(x){
        return(wilcox.test(x[whichhyper],x[whichhypo])$p.value)
      })
      sigprobes <- dfdmp[order(dfdmp[,2]),1]
      sigprobes <- sigprobes[1:nsigprobes]
      message("Returning ",length(sigprobes)," significant DMPs.")
      dfdmp <- dfdmp[dfdmp[,1] %in% sigprobes,]
      dfdmp <- dfdmp[order(match(dfdmp[,1],sigprobes)),]
      
    }
    if(getDMPs=="ftest"){
      message("Performing DMP test: Ftest...")
      tissueindex <- eval(parse(text=paste0("df$",dmpVariablename)))
      dmpall <- dmpFinder(as.matrix(getBeta(df)),pheno=tissueindex,type=c("categorical"))
      dmpsig <- dmpall[1:nsigprobes,]
      sigprobes <- as.character(rownames(dmpsig))
      message("Returning ",length(sigprobes)," significant DMPs.")
      dfdmp <- dfdmp[dfdmp[,1] %in% as.character(sigprobes),]
      dfdmp <- dfdmp[order(match(dfdmp[,1],as.character(sigprobes))),]
      dfdmp[,2] <- dmpsig$qval
      colnames(dfdmp)[2] <- "q.value"                    
      
    }
    
    allvars <- apply(getBeta(df[rownames(df) %in% sigprobes,whichhyper]),1,var)
    message("DMP variances range from ",min(allvars)," to ",max(allvars))
    if(is.null(varfilt)){
      keepprobesvar <- allvars[allvars > quantile(allvars,seq(0,1,0.1))[4]]
      message("Excluding ",length(allvars)-length(keepprobesvar)," probes with variances <= ",quantile(allvars,seq(0,1,0.1))[4])
    } else{
      keepprobesvar <- allvars[allvars>varfilt]
      message("Excluding ",length(allvars)-length(keepprobesvar)," probes with variance <= ",varfilt)
    }
    
    dfdmp <- dfdmp[dfdmp[,1] %in% as.character(names(keepprobesvar)),]
    dfdmp <- dfdmp[order(match(dfdmp[,1],as.character(names(keepprobesvar)))),]
    dfdmp$VarfiltTissue.variance <- keepprobesvar
    
    # tissue means
    hypertissue.mean <- rowMeans(getBeta(df[rownames(df) %in% dfdmp[,1],whichhyper]))
    hypertissue.mean <- hypertissue.mean[order(match(names(hypertissue.mean),dfdmp[,1]))]
    hypotissue.mean <- rowMeans(getBeta(df[rownames(df) %in% dfdmp[,1],whichhypo]))
    hypotissue.mean <- hypotissue.mean[order(match(names(hypotissue.mean),dfdmp[,1]))]
    # df return, new variables
    dfdmp$Mean.varfiltTissue <- hypertissue.mean
    dfdmp$Mean.tissue2 <- hypotissue.mean
    dfdmp$Meandif.VFTminus2 <- dfdmp$Mean.varfiltTissue-dfdmp$Mean.tissue2
    dfdmp$probetype <- "NA"
    dfdmp[dfdmp$Meandif.VFTminus2>0,]$probetype <- "hyper"
    dfdmp[dfdmp$Meandif.VFTminus2<0,]$probetype <- "hypo"
    # get hyper and hypo lists
    hyper <- dfdmp[dfdmp$probetype=="hyper",1]
    hypo <- dfdmp[dfdmp$probetype=="hypo",1]
    return(dfdmp)
    
  } else{
    return(message("ERROR: Need a minfi set object to perform DMP analysis."))
  }  
  
}
