ipr.calcPurity <- function(df,hyper,hypo){
  if(length(intersect(hyper,hypo))>0){
    return(message("ERROR: There is at least one probe in both hypo- and hyper-methylated CpG lists!"))
  }
  if(is.null(hyper)|is.null(hypo)){
    return(message("ERROR: No hyper or hypo probe IDs passed or inputted."))
  }
  # df: rownames need to be sample names/ID, colnames need to be probe IDs
  if(class(df) %in% c("MethylSet","GenomicMethylSet","GenomicRatioSet","RGChannelSet")){
    df <- t(getBeta(df))
  }
  dfreturn2 <- data.frame(rep("NA",nrow(df)),
                          rep(0,nrow(df)))
  colnames(dfreturn2) <- c("Sample_ID","IPR.purity")
  dfreturn2[,1] <- rownames(df)
  message("Computing kde modes of ",length(hyper)," hypermethylated and ",length(hypo)," hypomethylated DMPs for ",nrow(df)," samples...")
  message("Dim of df is: nrow=",nrow(df)," ncol=",ncol(df))
  
  for(samplei in 1:nrow(df)){
    # get sample ID and data
    dfi <- df[samplei,]
    isampleid <- as.character(rownames(df)[samplei])
    # get modified methylation levels
    imodmethyl <- ipr.dmpModmethy(betas=dfi,hyper=hyper,hypo=hypo)
    imodmethyl <- imodmethyl[!is.na(imodmethyl)]
    # get sample purity estimate
    idens.index <- unlist(density(imodmethyl)[1])
    idens.peaks <- unlist(density(imodmethyl)[2])
    ipurity <- idens.index[which(idens.peaks==max(idens.peaks))]
    dfreturn2[samplei,2] <- ipurity
    message("Sample ",isampleid," has estimated purity: ",round(ipurity,3),".")
    
    # progress bar
    ipr.progressBar(i=samplei,dur=nrow(df),unit="samples")

  }
  return(dfreturn2)
}
