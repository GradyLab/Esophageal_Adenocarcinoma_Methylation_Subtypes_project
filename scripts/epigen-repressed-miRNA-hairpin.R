#=======================================
# miRNA hairpin epigenetic repression
#=======================================
# read in and reformat miRNA hairpin locus expression
exp <- read.csv("tcga-eac-mirna-reads.csv",stringsAsFactors=FALSE)

rownames(exp) <- exp[,1];
exp <- exp[,grepl("reads_per_million_miRNA_mapped",as.character(exp[1,]))]; dim(exp)
exp <- exp[2:nrow(exp),]

# prefilter: at least 3 regional CpGs available
retainfilt <- c()
for(i in 1:nrow(xmirna)){
  if(length(unlist(strsplit(xmirna[i,]$cpg,";")))>=3){
    retainfilt <- c(retainfilt,i)
  }
  message(round(100*(i/nrow(xmirna)),3),"%")
}
mirna.retain <- intersect(rownames(exp),xmirna[retainfilt,]$mirna)
miranno <- xmirna[xmirna$mirna %in% mirna.retain,]
exp <- exp[rownames(exp) %in% mirna.retain,]
exp <- exp[order(match(rownames(exp),miranno$mirna)),]; identical(rownames(exp),miranno$mirna)

exp <- exp[,substr(colnames(exp),9,12) %in% pdat$Participant]

exp <- as.matrix(exp); class(exp) <- "numeric"
expn <- exp[,substr(colnames(exp),14,15)=="11"]; expt <- exp[,substr(colnames(exp),14,15)=="01"]
xn <- rowMeans(expn,na.rm=TRUE); summary(xn)

lfc.mhp <- expt
for(i in 1:nrow(lfc.mhp)){
  lfc.mhp[i,] <- log2((lfc.mhp[i,]+1)/xn[i])
  message(round(100*(i/nrow(lfc.mhp)),3),"%")
}

save(lfc.mhp,file="log2fcTN-miRNAhairpin_tcga-eac_firehose_RPMnorm.rda")

# set up correlation tests
gt <- gset.tcga.eac
gt <- gt[,order(match(gt$Participant,substr(colnames(lfc.mhp),9,12)))]
identical(gt$Participant,substr(colnames(lfc.mhp),9,12))

identical(gt[,gt$rpmm.mgroup=="HM"]$Participant,substr(colnames(lfc.mhp[,which(gt$rpmm.mgroup=="HM")]),9,12)) # subsetting doesn't change order matching

cordf.mhp <- as.data.frame(matrix(nrow=nrow(lfc.mhp),ncol=17))
colnames(cordf.mhp) <- c("mirna.hairpin.locus",
                         "HM.corp","HM.rho","HM.xmethy","HM.xexp",
                         "IM.corp","IM.rho","IM.xmethy","IM.xexp",
                         "LM.corp","LM.rho","LM.xmethy","LM.xexp",
                         "MM.corp","MM.rho","MM.xmethy","MM.xexp")

# run correlations for all, determine significant loci w. post-hoc filters
type <- c("HM","IM","LM","MM") # should be levels corresponding to variable 'rpmm.mgroup' in pData(gt)
ri <- list(c(2:5),c(6:9),c(10:13),c(14:17))
for(i in 1:nrow(exp)){
  cordf.mhp[i,1] <- rownames(exp)[i]
  cpgi <- intersect(unlist(strsplit(miranno[miranno$mirna==rownames(exp)[i],]$cpg,";")),rownames(gt))
  if(length(cpgi)>=3){
    for(j in 1:4){
      tj <- type[j]; whichtj <- gt$rpmm.mgroup==tj
      me.j <- as.numeric(colMeans(getBeta(gt[rownames(gt) %in% cpgi,whichtj])))
      ex.j <- as.numeric(lfc.mhp[i,substr(colnames(lfc.mhp),9,12) %in% gt[,whichtj]$Participant])
      
      ctj <- cor.test(ex.j,me.j,test="spearman")
      cordf.mhp[i,ri[[j]]] <- c(ctj$p.value,
                                ctj$estimate,
                                mean(me.j),
                                mean(ex.j,na.rm=TRUE))
    }
  }
  message(paste0(i," ",round(100*(i/nrow(exp)),3),"%"))
}

# clean up df
cdf <- cordf.mhp[!is.na(cordf.mhp[,2]),]; dim(cdf)
save(cdf,file="epigen-repgenes-mirna-hairpins_cor-exp-methy-summaries_tcga-eac.rda")

###

