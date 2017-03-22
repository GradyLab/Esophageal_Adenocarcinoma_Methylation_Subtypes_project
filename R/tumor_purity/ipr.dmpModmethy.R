ipr.dmpModmethy <- function(betas,hyper=hyper,hypo=hypo){
  
  dmpmod.cgretain <- c()
  for(i in 1:length(betas)){
    cpg <- as.character(names(betas)[i])
    beta.cpg <- as.numeric(betas[i])
    hypercg <- hyper[hyper %in% cpg]
    hypocg <- hypo[hypo %in% cpg]
    if(length(hypercg)>0 & length(hypocg)==0){
      dmpmod.cgretain[i] <- as.numeric(beta.cpg)
    }
    if(length(hypocg)>0 & length(hypercg)==0){
      dmpmod.cgretain[i] <- 1-as.numeric(beta.cpg)
    }
  }
  dmpmod.cgretain <- as.numeric(dmpmod.cgretain)
  message("Summary of modified methylation data: min=",round(min(dmpmod.cgretain),3),
          ", max=",round(max(dmpmod.cgretain),3),", mean=",round(mean(dmpmod.cgretain),3))

  return(dmpmod.cgretain)
  
}
