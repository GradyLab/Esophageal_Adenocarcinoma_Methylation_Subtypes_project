# Function Description: Estimate and depict tumor purity/enrichment using Zhang et al 2015's InfiniumPurify method.
# Author: SKM
# Original code ported to R from: https://bitbucket.org/zhengxiaoqi/infiniumpurify
# Variables:
#   df = methylation data (rows = markers, cols = samples)
#   calcPurity = whether to calculate purity
#   returnDMPonly = whether to only return differentially methylated probes (T vs. N)
#   plotunimod = whether to return plot of overlaid sample unimodal distributions for CpGs used in purity estimates
#   hyper = hypervariable (T vs. N) DMP list for purity estimations
#   hypo = hypomethylated (T vs. N) DMP list for purity estimations
#   getDMPs = whether to return list of DMPs (T vs. N, needs dmpVariablename and dmpVariancefiltTissue)
#   dmpVariablename = name of variable in pData("name-of-minfi-object")
#   dmpVariancefiltTissue = name of tissue in dmpVariablename variable that is hyperavariable/Tumor
#   nsigprobes= number of significant probes to return from DMP calculations
#   varfilt= filter for variance in DMP calculations
#   plotunimod.main = main title name for unimodal plot (needs plotunimod = TRUE) 
#   plotunimod.ymax= max for y-axis on unimodal plot (needs plotunimod = TRUE)

InfiniumPurifyR <- function(df,calcPurity=TRUE,
                            returnDMPonly=FALSE,plotunimod=FALSE,
                            hyper=NULL,hypo=NULL,
                            getDMPs=FALSE,dmpVariablename=NULL,
                            dmpVariancefiltTissue=NULL,
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
