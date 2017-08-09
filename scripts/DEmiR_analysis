# Integrated analysis of TCGA miR RNAseq data 
# Comparing high methylator samples to each other subtype
# Author: Sean Maden

# matmirna.log2fc converted from Illumina HiSeq V2 RNAseq RSEM-normalized expression file 
# Source: BROAD Firehose resource (https://gdac.broadinstitute.org/)

miR.lfc <- matmirna.log2fc[,substr(colnames(matmirna.log2fc),14,15)=="01"] # filter out normal samples
colnames(miR.lfc) <- substr(colnames(miR.lfc),9,12) # colnames as patient IDs
rownames(miR.lfc) <-gsub("\\|.*","",rownames(miR.lfc)) # rownames as miR IDs

# pdat is patient data for TCGA EACs
hmpat <-pdat[pdat$rpmm.mgroup=="HM",]$Participant # get patient ids for HM samples


#=============
# HM vs. IM
#=============

impat <- pdat[pdat$rpmm.mgroup=="IM",]$Participant # get IM patient IDs
mir.hmim <- miR.lfc[,colnames(miR.lfc) %in% c(hmpat,impat)]

# make a summary dataframe
mir.hmim.df <- data.frame(HM.mean=rowMeans(mir.hmim[,colnames(mir.hmim) %in% hmpat],na.rm=TRUE),
                          IM.mean=rowMeans(mir.hmim[,colnames(mir.hmim) %in% impat],na.rm=TRUE),
                          HM.nact=apply(mir.hmim[,colnames(mir.hmim) %in% hmpat],1,function(x){length(x[is.na(x)])}),
                          IM.nact=apply(mir.hmim[,colnames(mir.hmim) %in% impat],1,function(x){length(x[is.na(x)])}),
                          ttest.pval=rep(NA,nrow(mir.hmim)))
rownames(mir.hmim.df) <- rownames(mir.hmim)

# apply NA filter (<=10 missing samples)
mir.hmim.df <- mir.hmim.df[mir.hmim.df$HM.nact<=10 & mir.hmim.df$IM.nact<=10 ,]; dim(mir.hmim.df)

for(i in 1:nrow(mir.hmim.df)){
  
  # t-test of miR expression
  hmvali <- mir.hmim[rownames(mir.hmim.df)[i],hmpat]
  imvali <- mir.hmim[rownames(mir.hmim.df)[i],impat]
  hmvali <- hmvali[!is.na(hmvali)];imvali <- imvali[!is.na(imvali)]
  ti <- t.test(hmvali,imvali)
  
  mir.hmim.df[i,]$ttest.pval <- ti$p.value
}

# make volcano plot summary data frame
vpmir.hmim <- data.frame(miR=rownames(mir.hmim.df),
                         dif.hmim=mir.hmim.df$HM.mean-mir.hmim.df$IM.mean,
                         negl10p=-1*log10(mir.hmim.df$ttest.pval))

# HM vs. IM 
moi.hmim <- as.character(vpmir.hmim[abs(vpmir.hmim$dif.hmim)>=1 & vpmir.hmim$negl10p >= -1*log10(0.05),]$miR)

#===========
# HM vs. LM
#===========
lmpat <-pdat[pdat$rpmm.mgroup=="LM",]$Participant

mir.hmlm <- miR.lfc[,colnames(miR.lfc) %in% c(hmpat,lmpat)]
dim(mir.hmlm)
mir.hmlm.df <- data.frame(HM.mean=rowMeans(mir.hmlm[,colnames(mir.hmlm) %in% hmpat],na.rm=TRUE),
                          LM.mean=rowMeans(mir.hmlm[,colnames(mir.hmlm) %in% lmpat],na.rm=TRUE),
                          HM.nact=apply(mir.hmlm[,colnames(mir.hmlm) %in% hmpat],1,function(x){length(x[is.na(x)])}),
                          LM.nact=apply(mir.hmlm[,colnames(mir.hmlm) %in% lmpat],1,function(x){length(x[is.na(x)])}),
                          ttest.pval=rep(NA,nrow(mir.hmlm)))
rownames(mir.hmlm.df) <- rownames(mir.hmlm)

mir.hmlm.df <- mir.hmlm.df[mir.hmlm.df$HM.nact<=10 & mir.hmlm.df$LM.nact<=10 ,]; dim(mir.hmlm.df) # NA filter

for(i in 1:nrow(mir.hmlm.df)){
  hmvali <- mir.hmlm[rownames(mir.hmlm.df)[i],hmpat]
  lmvali <- mir.hmlm[rownames(mir.hmlm.df)[i],lmpat]
  hmvali <- hmvali[!is.na(hmvali)];lmvali <- lmvali[!is.na(lmvali)]
  ti <- t.test(hmvali,lmvali)
  mir.hmlm.df[i,]$ttest.pval <- ti$p.value
}

vpmir.hmlm <- data.frame(miR=rownames(mir.hmlm.df),
                         dif.hmlm=mir.hmlm.df$HM.mean-mir.hmlm.df$LM.mean,
                         negl10p=-1*log10(mir.hmlm.df$ttest.pval))

# HM vs. LM 
moi.hmlm <- as.character(vpmir.hmlm[abs(vpmir.hmlm$dif.hmlm)>=1 & vpmir.hmlm$negl10p >= -1*log10(0.05),]$miR)

#================
# HM vs. MM
#================
mmpat <-pdat[pdat$rpmm.mgroup=="MM",]$Participant

mir.hmmm <- miR.lfc[,colnames(miR.lfc) %in% c(hmpat,mmpat)]
dim(mir.hmmm)
mir.hmmm.df <- data.frame(HM.mean=rowMeans(mir.hmmm[,colnames(mir.hmmm) %in% hmpat],na.rm=TRUE),
                          MM.mean=rowMeans(mir.hmmm[,colnames(mir.hmmm) %in% mmpat],na.rm=TRUE),
                          HM.nact=apply(mir.hmmm[,colnames(mir.hmmm) %in% hmpat],1,function(x){length(x[is.na(x)])}),
                          LM.nact=apply(mir.hmmm[,colnames(mir.hmmm) %in% mmpat],1,function(x){length(x[is.na(x)])}),
                          ttest.pval=rep(NA,nrow(mir.hmmm)))
rownames(mir.hmmm.df) <- rownames(mir.hmmm)

mir.hmmm.df <- mir.hmmm.df[mir.hmmm.df$HM.nact<=10 & mir.hmmm.df$LM.nact<=10 ,]; dim(mir.hmmm.df) # NA filter

for(i in 1:nrow(mir.hmmm.df)){
  hmvali <- mir.hmmm[rownames(mir.hmmm.df)[i],hmpat]
  mmvali <- mir.hmmm[rownames(mir.hmmm.df)[i],mmpat]
  hmvali <- hmvali[!is.na(hmvali)];mmvali <- mmvali[!is.na(mmvali)]
  ti <- t.test(hmvali,mmvali)
  mir.hmmm.df[i,]$ttest.pval <- ti$p.value
}

vpmir.hmmm <- data.frame(miR=rownames(mir.hmmm.df),
                         dif.hmmm=mir.hmmm.df$HM.mean-mir.hmmm.df$MM.mean,
                         negl10p=-1*log10(mir.hmmm.df$ttest.pval))

# HM vs. MM 
moi.hmmm <- as.character(vpmir.hmmm[abs(vpmir.hmmm$dif.hmmm)>=1 & vpmir.hmmm$negl10p >= -1*log10(0.05),]$miR)

#=============================
# miR volcanos multiplot
#=============================
dev.off()

jpeg("volcano-miR-dif_hm-v-subtypes_eac-tcga.jpg",units="in",6,4,res=400)
par(mfrow=c(1,3),mar=c(6.5,4,3,1))

plot(vpmir.hmim$dif.hmim,
     vpmir.hmim$negl10p,
     xlab="HM - IM",
     ylab="-log10(p-val)",
     #main="HM vs. IM",
     pch=16,col=rgb(0.2,0.2,0.8,alpha=0.3),
     xlim=c(-2.5,2.5),
     ylim=c(0,4.5))
abline(v=-1,col="red");abline(v=1,col="red")
abline(h=-1*log10(0.05),col="purple")

plot(vpmir.hmlm$dif.hmlm,
     vpmir.hmlm$negl10p,
     xlab="HM - LM",
     #main="HM vs. LM",
     ylab="",
     pch=16,col=rgb(0.2,0.2,0.8,alpha=0.3),
     xlim=c(-2.5,2.5),
     ylim=c(0,4.5))
abline(v=-1,col="red");abline(v=1,col="red")
abline(h=-1*log10(0.05),col="purple")

plot(vpmir.hmmm$dif.hmmm,
     vpmir.hmmm$negl10p,
     xlab="HM - MM",
     ylab="",
     #main="HM vs. MM",
     pch=16,col=rgb(0.2,0.2,0.8,alpha=0.3),
     xlim=c(-2.5,2.5),
     ylim=c(0,4.5))
abline(v=-1,col="red");abline(v=1,col="red")
abline(h=-1*log10(0.05),col="purple")

# multiplot titles
mtext("Differential miR Expression",outer=TRUE,cex=1.2,line=-1.7)
mtext("Expression Difference (log2FC)",outer=TRUE,cex=0.9,line=-30)

dev.off()

#==========================
# venn diagrams of DEmiRs
#==========================

library(VennDiagram)

jpeg("venn-demiR_hm-v-subtypes_eac-tcga.jpg",4,4,units="in",res=400)
draw.triple.venn(area1=length(moi.hmim),
                 area2=length(moi.hmlm),
                 area3=length(moi.hmmm),
                 n12=length(intersect(moi.hmim,moi.hmlm)),
                 n23=length(intersect(moi.hmlm,moi.hmmm)),
                 n13=length(intersect(moi.hmim,moi.hmmm)),
                 n123=length(intersect(moi.hmim,intersect(moi.hmlm,moi.hmmm))),
                 category=c(paste0("HM-vs-IM\n(N=",length(moi.hmim),")"),
                            paste0("HM-vs-LM\n(N=",length(moi.hmlm),")"),
                            paste0("HM-vs-MM\n(N=",length(moi.hmmm),")")),
                 col=c("coral","gray","lightblue"),
                 alpha=c(0.2,0.2,0.2),
                 fill=c("coral","gray","lightblue"),
                 cat.fontfamily="arial",fontfamily="arial",
                 margin=0.01,cat.dist=0.1)
dev.off()

#
