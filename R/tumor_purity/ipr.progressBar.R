ipr.progressBar <- function(i,dur,unit){
  
  stat.prognum <- round((i/dur)/2,2)*100
  stat.prog <- rep("=",stat.prognum)
  
  if(stat.prognum < dur){
    stat.left <- rep(".",50-stat.prognum)
  } else{
    stat.left <- ""
  }
  
  return(message("Progress: ",round(i/dur,2)*100,"%, ",i,"/",dur," ",unit,"\n||",
                 paste0(stat.prog),">>",paste0(stat.left),"|||"))
  
}
