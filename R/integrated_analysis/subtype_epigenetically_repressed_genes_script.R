#==================================================
# Epigenetically repressed genes in EAC subtypes
#==================================================

g <- gtcga[grep("TSS|5'UTR",getAnnotation(gtcga)$UCSC_RefGene_Group),] # gtcga a GenomeRatioSet including normal and tumor samples
mgenes <- unique(unlist(strsplit(getAnnotation(g)$UCSC_RefGene_Name,";"))) # get initial list of genes in methylation data

mdfall <- as.data.frame(matrix(nrow=length(mgenes),ncol=2)) # df to hold mean tumor and normal methylation across genes
rownames(mdfall) <- mgenes; colnames(mdfall) <- c("N","T")
bnorm <- getBeta(g[,g$Tissue_Type=="normal"])
btumor <- getBeta(g[,g$Tissue_Type=="tumor"])
aall <- as.data.frame(getAnnotation(g))

for(i in 1:nrow(mdfall)){
  genecgi <- rownames(aall[grep(paste0("(^|;)",rownames(mdfall)[i],"(;.|$)"),aall$UCSC_RefGene_Name),])
  if(length(genecgi)>=3){ # filter on availability of promoter CpGs in methylation data
    mdfall[i,1] <- mean(colMeans(bnorm[genecgi,])) # get mean promoter methylation at normal samples for gene i
    mdfall[i,2] <- mean(colMeans(btumor[genecgi,])) # get mean promoter methlation at tumor samples for gene i
    message("finished ",i," or ",round(i/nrow(mdfall),3),"%")
  }
}

m2 <- mdfall[!is.na(mdfall$N) & !is.na(mdfall$T),] # in case of NA methylation, filter that gene here
mgenesfilt <- rownames(m2[which(m2$T - m2$N >0 & m2$T>=0.2),]) # filter on properties of tumor mean methylation

#=========================
# Expression data filters
#=========================
# ct.t and ct.n are tables of tumor gene counts from RNAseq data (source: RSEM normalized genes, Illumina HiSeq v2 [http://firebrowse.org/?cohort=ESCA&download_dialog=true])

lfct <- ct.t 

for(i in 1:nrow(lfct)){ # convert tumor counts to log2FC
  ri <- as.numeric(ct.t[i,])
  lfct[i,] <- log2(as.numeric(ri)/mean(as.numeric(ri)))
}

eall <- lfct
nact <- apply(eall,1,function(x){length(x[is.na(x)])}) # count NAs (or samples with gene count = 0/log2FC = -Inf)

# final gene list to test
finalgenes <- intersect(rownames(eall[which(nact<= 0.4*87),]),mgenesfilt) # filter on less than 40% tumor samples missing

#================================
# Epigenetically repressed genes
#================================

# Using below helper function, iterate over subtypes
hm.repgenes.listobj <- get.repgenes("HM")
im.repgenes.listobj <- get.repgenes("IM")
lm.repgenes.listobj <- get.repgenes("LM")
mm.repgenes.listobj <- get.repgenes("MM")



# helper function
get.repgenes <- function(subtype,padj.cutoff=0.05,meanexp.cutoff= -2,rho.cutoff=0){
  return.list <- list()
  
  gs <- g[,g$rpmm.mgroup==subtype]
  
  medf <- as.data.frame(matrix(nrow=length(finalgenes),ncol=ncol(gs))) # methylation df to hold gene promoter methylation, by sample
  colnames(medf) <- gs$Participant
  rownames(medf) <- finalgenes
  
  aall <- as.data.frame(getAnnotation(gs))
          
  for(i in 1:length(finalgenes)){
    medf[i,] <- colMeans(getBeta(gs[rownames(aall[grep(paste0("(^|;)",rownames(medf)[i],"(;.|$)"),aall$UCSC_RefGene_Name),]),]))
  }
  
  exps <- lfct[rownames(lfct) %in% finalgenes,substr(colnames(lfct),9,12) %in% colnames(medf)]; 
  exps <- exps[order(match(rownames(exps),finalgenes)),]
  colnames(exps) <- substr(colnames(exps),9,12) # extract participant id from column names in expression data
  exps <- exps[,order(match(colnames(exps),colnames(medf)))]
  
  # ensure ordering is identical for methylation and expression
  identical(colnames(exps),colnames(medf));identical(rownames(medf),rownames(exps))
  
  cordf.s <- as.data.frame(matrix(nrow=nrow(medf),ncol=5))
  rownames(cordf.s) <- rownames(medf)
  colnames(cordf.s) <- c("mean.methy","mean.expr","nmissing.expr","spear.rho","spear.cor.pval")
  for(i in 1:nrow(cordf.s)){
    x <- as.numeric(exps[i,])
    rms <- !is.infinite(x) # create filter for missing samples (having -Inf log2FC)
    
    if(length(x[rms])>=10){
      m <- as.numeric(medf[i,]); m <- m[rms]
      
      cordf.s[i,1] <- mean(m)
      cordf.s[i,2] <- mean(as.numeric(x[rms],na.rm=TRUE)) 
      cordf.s[i,3] <- length(x[is.na(x)|is.infinite(x)])
      
      cori <- cor.test(m,x[rms],method="spearman")
      cordf.s[i,4] <- cori$estimate
      cordf.s[i,5] <- cori$p.value
    }
  }; cordf.s <- cordf.s[!is.na(cordf.s$spear.cor.pval),]
  
  cordf.s2 <- cordf.s[cordf.s$mean.expr <= meanexp.cutoff & cordf.s$spear.rho < rho.cutoff,]
  cordf.s2$padj <- p.adjust(cordf.s2$spear.cor.pval,method="BH")
  epigen.genes <- rownames(cordf.s2[cordf.s2$padj<=padj.cutoff,]) # get epigenetically repressed gene list
  
  return.list[[1]] <- cordf.s2; return.list[[2]] <- epigen.genes
  return.list[[3]] <- data.frame(padj.cutoff=padj.cutoff,meanexp.cutoff=meanexp.cutoff,correlation.cutoff=spear.cutoff)
  names(return.list) <- c("correlation_summary","epigenetically repressed genes","significance parameters")
  return(return.list)
}


#
