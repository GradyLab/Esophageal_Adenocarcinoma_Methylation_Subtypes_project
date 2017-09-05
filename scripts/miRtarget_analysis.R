# miRNA target gene analysis
#library(miRNAtap)
library(org.Hs.eg.db)

lfct.na <- ifelse(is.infinite(lfct),NA,lfct)

#====================
# EXAMPLE: HM vs. IM
#====================
hmpat <- pdat[pdat$rpmm.mgroup=="HM",]$Participant
impat <- pdat[pdat$rpmm.mgroup=="IM",]$Participant

lfct.hmim <- lfct.na[,substr(colnames(lfct.na),9,12) %in% c(hmpat,impat)]
colnames(lfct.hmim) <- substr(colnames(lfct.hmim),9,12)

# get directionally accurate genes
mirtargets.list <- list()
for(i in 1:length(moi.hmim)){
  miri <- moi.hmim[i] # gsub("hsa-","",moi.hmlm[i])
  targetsi <- pci.df[pci.df$miR==miri,]$target.symbol
 
  targetexpr <- lfct.hmim[rownames(lfct.hmim) %in% targetsi,]
  
  dirmiri <- vpmir.hmim[vpmir.hmim$miR==moi.hmim[i],]$dif.hmim
  
  # filter on directionality of population expression 
  # (eg. if miR expr: group1>group2, retain targets where mRNA expr: group1<group2)
  if(dirmiri>0){
    targetfilter <- targetexpr[rowMeans(targetexpr[,hmpat])<rowMeans(targetexpr[,impat]),]
  }
  if(dirmiri<0){
    targetfilter <- targetexpr[rowMeans(targetexpr[,hmpat])>rowMeans(targetexpr[,impat]),]
  }
  
  # NA filter, <10 samples missing
  hmnact <- apply(targetfilter[,hmpat],1,function(x){length(x[is.na(x)])})
  imnact <- apply(targetfilter[,impat],1,function(x){length(x[is.na(x)])})
  targetfilter <- targetfilter[!(hmnact>10|imnact>10),]
  
  ttests.list <- list()
  for(j in 1:nrow(targetfilter)){
    ttests.list[[j]] <- t.test(targetfilter[j,hmpat],targetfilter[j,impat])$p.value
    message("finished ",j," or ",round(100*j/nrow(targetfilter)),"%")
  }
  
  targetdfi <- data.frame(HM.expr.target=rowMeans(targetfilter[,hmpat]),
                          IM.expr.target=rowMeans(targetfilter[,impat]),
                          miR.id=miri,
                          miR.dif=dirmiri,
                          targetgene.symbol=rownames(targetfilter),
                          targetgene.ttest.pval=unlist(ttests.list))
  
  # mirtargets.names <- targetsi
  namesi1 <- names(mirtargets.list)
  mirtargets.list[[i]] <- targetdfi
  names(mirtargets.list) <- c(namesi1,miri)
 
  message("finished ",i," miR or ",round(100*i/length(moi.hmim)),"%")
}


hmim.mirtargets.df <- as.data.frame(matrix(nrow=0,ncol=ncol(mirtargets.list[[1]])))
colnames(hmim.mirtargets.df) <- colnames(mirtargets.list[[1]])
for(k in 1:length(mirtargets.list)){
  hmim.mirtargets.df <- rbind(hmim.mirtargets.df,mirtargets.list[[k]])
}

#======================================
# Programmatically view volcano plot
#======================================
colplot <- as.character(hmim.mirtargets.df$miR.id)
miRuniquelist <- unique(hmim.mirtargets.df$miR.id)
col.lvl <- sample(colours(),length(miRuniquelist))

for(a in 1:length(miRuniquelist)){
  colplot[colplot==miRuniquelist[a]] <- col.lvl[a]
}; table(colplot)

xlimi <- c(-2,2); ylimi <- c(0,2.5)

jpeg("vplot-lfct-targetexpr-miR_HM-v-IM_eac-tcga.jpg",8,7,units="in",res=400)

par(oma = c(1, 1, 1, 12))
plot(hmim.mirtargets.df$HM.expr.target-hmim.mirtargets.df$IM.expr.target,
     -1*log10(hmim.mirtargets.df$targetgene.ttest.pval),
     xlim=xlimi,ylim=ylimi,col=colplot,pch=16,
     xlab="HM - IM",ylab = "-1*log(pval)")
abline(h=-1*log10(0.05),col="purple")
abline(v=-1,col="red"); abline(v=1,col="red")

# add labels to filtered genes
xfilt <- hmim.mirtargets.df[hmim.mirtargets.df$targetgene.ttest.pval<0.05 &
                              abs(hmim.mirtargets.df$HM.expr.target-hmim.mirtargets.df$IM.expr.target)>=1,]
text(xfilt$HM.expr.target-xfilt$IM.expr.target, 
     -1*log10(xfilt$targetgene.ttest.pval), 
     labels=xfilt$targetgene.symbol, cex= 0.7,pos=3)

# add legend from new blank plot window
par(fig = c(0, 1, 0, 1), oma = c(0, 25, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
legend("center",ncol=1,
       legend=c(as.character(miRuniquelist),"pval=0.05"),
       col=c(col.lvl,"purple"),
       pch=c(rep(16,length(miRuniquelist)),NA),
       lty=c(rep(NA,length(miRuniquelist)),1),
       cex=0.8)

dev.off()
