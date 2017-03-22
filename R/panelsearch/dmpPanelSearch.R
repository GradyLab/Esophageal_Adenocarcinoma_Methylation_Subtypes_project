# variable descriptions:
  # gset is GenomicRatioSet of samples from which to perform DMP analysis and assemble/test panels
  # dmpprobelist is user-specified list of probes to draw from when assembling panels
  # dmplistrange is number of dmps from which to assemble random panels
  # testgroups is varname with identifiers for all cancer samples' mgroups
  # level1name is variable name/level to be assigned "1" (all others assigned 0) for DMP calc
  # comb.df.n, when !=NULL, is number of assembled panels to sample from all panel combinations
  # nprobespanel is total size of probe panels assembled at random from DMPs
  # nprobescutoff is minimum count of probes meeting betacutoff for sample to be "positive" in panel test
  # betacutoff is methylation value cutoff for probe to be considered "positive" in sample for panel test
  

dmpPanelSearch <- function(gset,dmplistrange=20,dmpprobelist=NULL,testgroups,level1name,comb.df.n=NULL,nprobespanel=5,nprobescutoff=3,betacutoff=0.5){
  require(minfi)
  # is gset a GRS?
  if(!class(gset)=="GenomicRatioSet"){
    return(print(paste0("Gset object needs to be of class 'GenomicRatioSet'"))
  }
  
  bvals.train <- getBeta(gset) 
  testgroups <- ifelse(testgroups %in% level1name,1,0)
  
  # Whether to run dmp test, and what probe list to use for panel assembly
  if(is.null(dmpprobelist)){
    print(paste0("Running DMP test..."))
    dmp.calc <- dmpFinder(bvals.train,testgroups,type="categorical")
    print(paste0("From DMP analysis, you are choosing the ",dmplistrange," top DMPs with qvalues ranging from ",
                 range(dmp.calc[1:dmplistrange,4])[1]," to ",range(dmp.calc[1:dmplistrange,4])[2]," for panel assembly"))
    # take top where samples are sorted by increasing q-val
    top.dmp.list <- rownames(dmp.calc[1:dmplistrange,])
    # assemble probe panels
    comb.df <- t(combn(c(top.dmp.list),nprobespanel))
    print(paste0("Out of ",dmplistrange," DMPs, ",nrow(comb.df)," panels of ",nprobespanel," probes were assembled"))
  }
  else{
    print(paste0("Using user-defined list of ",length(dmpprobelist)," probes for panel assembly."))
    top.dmp.list <- dmpprobelist
    # assemble probe panels
    comb.df <- t(combn(c(top.dmp.list),nprobespanel))
    print(paste0("Out of ",length(dmpprobelist)," provided probes, ",nrow(comb.df)," panels of ",nprobespanel," probes were assembled"))
  }
  
  # Choose random panel subset or whole panel (default)?
  if(!is.null(comb.df.n)){
    print(paste0("You are selecting ",comb.df.n," random panel(s) from ",nrow(comb.df)," total panel combinations"))
    comb.df <- comb.df[sample(nrow(comb.df),comb.df.n),] 
  }
  
  gset <- gset[rownames(gset) %in% top.dmp.list,]
  binarymatrix <- as.data.frame(t(getBeta(gset)))
  binarymatrix <- as.data.frame(ifelse(binarymatrix>=betacutoff,1,0))
  binarymatrix$ref <- testgroups
  data.df <- binarymatrix
  binarymatrix <- as.matrix(binarymatrix)
  ref <- as.numeric(testgroups) # assign numeric dichotomized reference variable (ie. HM=1; IM/LM/MM=0)
  tpfn <- sum(ref) # sum of positive state samples
  tnfp <- nrow(data.df)-tpfn # sum of negative state samples
  
  # return dataframe with scores and sens/spec (rows panels, col's patients and sens/spec)
  return.df <- as.data.frame(matrix(,nrow=nrow(comb.df),ncol=ncol(gset)+2))
  colnames(return.df) <- c(colnames(gset),"sens","spec")
  
  print(paste0("Calculating panel sensitivity and specificity..."))
  
  for (j in 1:nrow(comb.df)){
    panelmatrix <- binarymatrix[,colnames(binarymatrix) %in% comb.df[j,]]
    tp <- c()
    tn <- c()
    
    # determine scores according to panel test criteria (assemble character string then assess in console)
    matrixmodvar <- "paste0(ifelse(panelmatrix[,1]+"
    for(a in 3:ncol(comb.df)-1){
      matrixmodvar <- paste0(matrixmodvar,"panelmatrix[,",a,"]+")
    }
    matrixmodvar <- paste0(matrixmodvar,"panelmatrix[,",ncol(comb.df),"]>=",nprobescutoff,",1,0))")
    
    panel.calls <- as.numeric(eval(parse(text=matrixmodvar)))
    return.df[j,1:nrow(panelmatrix)] <- panel.calls
    for (k in 1:nrow(panelmatrix)){
      tp[k] <- ifelse(return.df[j,k]==1 & ref[k]==1, 1, 0)
      tn[k] <- ifelse(return.df[j,k]==0 & ref[k]==0, 1, 0)
    }
    return.df[j,nrow(binarymatrix)+1] <- sum(tp)/tpfn # calculate sensitivity
    return.df[j,nrow(binarymatrix)+2] <- sum(tn)/tnfp # calculate specificity
  }
  # assign rownames to results df (assemble character string then assess in console)
  rownamevar <- "paste0(comb.df[,1],';',"
  for(i in 2:ncol(comb.df)-1){
    rownamevar <- paste0(rownamevar,"comb.df[,",i,"],';',")
  }
  rownamevar <- paste0(rownamevar,"comb.df[,",ncol(comb.df),"])")
  
  rownames(return.df) <- eval(parse(text=rownamevar))
  
  print(paste0("Your search returned panels with maximum sensitivity of ",max(return.df$spec),
        " and maximum specificity of ",max(return.df$sens),
        ". Among optimal panels, there are a total of ",nrow(return.df[return.df$spec==max(return.df$spec),])," panels with maximum specificity, where the best panel shows specificity = ",max(return.df$spec)," and sensitivity = ",max(return.df[return.df$spec==max(return.df$spec),]$sens),
        ", and there are a total of ",nrow(return.df[return.df$sens==max(return.df$sens),])," panels with maximum sensitivity, where the best panel shows sensitivity = ",max(return.df$sens)," and specificity = ",max(return.df[return.df$sens==max(return.df$sens),]$spec),"."))
  
  print(paste0("...Done!"))
  return(return.df)
}
