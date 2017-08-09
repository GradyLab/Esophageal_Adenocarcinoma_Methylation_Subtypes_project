# miRNA target gene analysis
# Author: Sean Maden
library(miRNAtap)
library(org.Hs.eg.db)

# lfct is log2FC converted data from Illumina HiSeq V2 RSEM-normalized counts
# Source: https://gdac.broadinstitute.org/

lfct.na <- ifelse(is.infinite(lfct),NA,lfct)

#=====================
# Example: HM vs. LM
#=====================

# pdat contains patient data linking IDs to methylation subtypes (HM/IM/LM/MM)
hmpat <- pdat[pdat$rpmm.mgroup=="HM",]$Participant
lmpat <- pdat[pdat$rpmm.mgroup=="LM",]$Participant

lfct.hmlm <- lfct.na[,substr(colnames(lfct.na),9,12) %in% c(hmpat,lmpat)]
colnames(lfct.hmlm) <- substr(colnames(lfct.hmlm),9,12)

# get directionally accurate genes
mirtargets.list <- list()
for(i in 1:length(moi.hmlm)){
  miri <- gsub("hsa-","",moi.hmlm[i])
  predi <- getPredictedTargets(miri,species="hsa")
  rankpred <- predi[,'rank_product']
  
  # targets given with entrez gene ids
  x <- org.Hs.egSYMBOL ; mapped_genes <- mappedkeys(x); xx <- as.list(x[mapped_genes]); #xx[1:5]
  xintpred <- xx[names(xx) %in% names(rankpred)]; xintpred <- xintpred[order(match(names(xintpred),names(rankpred)))]
  identical(names(xintpred),names(rankpred))
  
  names(rankpred) <- as.character(unlist(xintpred)) # finally map gene symbol to rank
  
  targetexpr <- lfct.hmlm[rownames(lfct.hmlm) %in% names(rankpred),]
  
  # filter on directionality (ie. for dirmiri>0, retain genes where HM expr < LM expr
  dirmiri <- vpmir.hmlm[vpmir.hmlm$miR==moi.hmlm[i],]$dif.hmlm
  if(dirmiri>0){
    targetfilter <- targetexpr[rowMeans(targetexpr[,hmpat])<rowMeans(targetexpr[,lmpat]),]
  }
  if(dirmiri<0){
    targetfilter <- targetexpr[rowMeans(targetexpr[,hmpat])>rowMeans(targetexpr[,lmpat]),]
  }
  
  # NA filter, <10 samples missing
  hmnact <- apply(targetfilter[,hmpat],1,function(x){length(x[is.na(x)])})
  lmnact <- apply(targetfilter[,lmpat],1,function(x){length(x[is.na(x)])})
  targetfilter <- targetfilter[!(hmnact>10|lmnact>10),]
  
  ttests.list <- list()
  for(j in 1:nrow(targetfilter)){
    ttests.list[[j]] <- t.test(targetfilter[j,hmpat],targetfilter[j,lmpat])$p.value
    message("finished ",j," or ",round(100*j/nrow(targetfilter)),"%")
  }
  
  targetdfi <- data.frame(HM.expr.target=rowMeans(targetfilter[,hmpat]),
                          lm.expr.target=rowMeans(targetfilter[,lmpat]),
                          miR.id=miri,
                          miR.dif=dirmiri,
                          targetgene.symbol=rownames(targetfilter),
                          targetgene.ttest.pval=unlist(ttests.list))
  
  mirtargets.names <- names(mirtargets.list)
  mirtargets.list[[i]] <- targetdfi
  names(mirtargets.list) <- c(mirtargets.names,miri)
  
  message("finished ",i," miR or ",round(100*i/length(moi.hmlm)),"%")
}


hmlm.mirtargets.df <- as.data.frame(matrix(nrow=0,ncol=ncol(mirtargets.list[[1]])))
colnames(hmlm.mirtargets.df) <- colnames(mirtargets.list[[1]])
for(k in 1:length(mirtargets.list)){
  hmlm.mirtargets.df <- rbind(hmlm.mirtargets.df,mirtargets.list[[k]])
}

#======================================
# Programmatically view volcano plot
#======================================

colplot <- as.character(hmlm.mirtargets.df$miR.id)
miRuniquelist <- unique(hmlm.mirtargets.df$miR.id)
col.lvl <- sample(colours(),length(miRuniquelist))
for(a in 1:length(miRuniquelist)){
  colplot[colplot==miRuniquelist[a]] <- col.lvl[a]
}; table(colplot)

jpeg("vplot-lfct-targetexpr-miR_hmlm_eac-tcga.jpg",8,7,units="in",res=400)

# define image margins, making room for legend at right
par(oma = c(1, 1, 1, 10))
plot(hmlm.mirtargets.df$HM.expr.target-hmlm.mirtargets.df$lm.expr.target,
     -1*log10(hmlm.mirtargets.df$targetgene.ttest.pval),
     xlim=c(-4,4),ylim=c(0,8),col=colplot,pch=16,
     xlab="HM - LM",ylab = "-1*log(pval)")
abline(h=-1*log10(0.05),col="purple")
abline(v=-1,col="red"); abline(v=1,col="red")
# add labels to filtered genes
xfilt <- hmlm.mirtargets.df[hmlm.mirtargets.df$targetgene.ttest.pval<0.05 &
                              abs(hmlm.mirtargets.df$HM.expr.target-hmlm.mirtargets.df$lm.expr.target)>=1,]
text(xfilt$HM.expr.target-xfilt$lm.expr.target, 
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

#
