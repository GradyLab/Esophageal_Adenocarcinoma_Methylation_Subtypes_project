InfiniumPurifyR <- function(df,calcPurity=TRUE,
                            returnDMPonly=FALSE,plotunimod=FALSE,
                            hyper=NULL,hypo=NULL,
                            getDMPs=FALSE,dmpVariablename=NULL, # rename to dmpTissuecol
                            dmpVariancefiltTissue=NULL, # rename to dmpVarianceTissue
                            nsigprobes=NULL,varfilt=NULL,
                            plotunimod.main=NULL,plotunimod.ymax=NULL){
  
  if(!class(df) %in% c("data.frame","matrix")){
    if(!class(df) %in% c("MethylSet","GenomicMethylSet","GenomicRatioSet","RGChannelSet")){
      return(message("ERROR: DF variable not existant or not of valid class."))
    }
  }
  
  # 1. calculate DMPs?
  if(!getDMPs==FALSE){
    message("Setting up DMP test...")
    dfreturn1 <- ipr.dmpOptions(df=df,getDMPs=getDMPs,varfilt=varfilt,
                                nsigprobes=nsigprobes,
                                dmpVariablename=dmpVariablename,
                                dmpVariancefiltTissue=dmpVariancefiltTissue)
    hyper <- dfreturn1[dfreturn1$probetype=="hyper",1]
    hypo <- dfreturn1[dfreturn1$probetype=="hypo",1]
    if(returnDMPonly==TRUE){
      return(dfreturn1)
    }
  }
  
  # 2. calculate tumor purity?
  if(calcPurity==TRUE){
    dfreturn2 <- ipr.calcPurity(df=df,hyper=hyper,hypo=hypo)
    message("Success! Completed purity estimates for all samples.")  
  }
  
  # 3. plot unimodal distributions?
  if(plotunimod==TRUE){
    message("Plotting unimodal distributions...")
    ipr.densplotUnimod(df=df,hyper=hyper,hypo=hypo,main=plotunimod.main,plotunimod.ymax)
  }    
  
  # results summary/output
  message("Completed analysis. Mean estimated purity was: ",mean(dfreturn2[,2],na.rm=TRUE))
  if("dfreturn1" %in% ls()){
    dfreturn <- list(dfreturn1,dfreturn2)
    names(dfreturn) <- c("dmpresults","IP.estimates")
    return(dfreturn)
  } else{
    return(dfreturn2)
  }
  
}
