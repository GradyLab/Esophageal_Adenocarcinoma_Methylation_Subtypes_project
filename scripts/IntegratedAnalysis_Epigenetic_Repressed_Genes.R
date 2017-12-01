# Differential and Shared Epigenetically Repressed Genes Across EAC Methylation Subtypes
# Script author: Sean Maden

# Dependancies:
# github functions: https://github.com/metamaden/HM450_tools
# minfi, sva

load(paste0(getwd(),"/gtn.rda")) # gtn a GenomicRatioSet containing TCGA tumor and normal EACs, preprocessed as above. 
load(paste0(getwd(),"/log2fc_eac.rda") # log2fc_eac is tumor gene expression from Firehose-preprocessed dataset

### Prefilters: Expression data 
ttest.genesin <- matrix(nrow=nrow(log2fc.eac),ncol=6) 
colnames(ttest.genesin)<-c("geneID","inf.HM","inf.NHM","inf.IM","inf.LM","inf.MM")
colnames(log2fc.eac) <- substr(colnames(log2fc.eac),9,12) # all participant ids unique (no normal samples in expr dataset)
for(i in 1:nrow(log2fc.eac)){
  ttest.genesin[i,1] <- rownames(log2fc.eac)[i]; 
  genei.dat <- log2fc.eac[i,]
  ttest.genesin[i,2] <- length(genei.dat[is.infinite(genei.dat) & names(genei.dat) %in% hmid]) # HM
  ttest.genesin[i,3] <- length(genei.dat[is.infinite(genei.dat) & names(genei.dat) %in% nhmid]) # NHM
  ttest.genesin[i,4] <- length(genei.dat[is.infinite(genei.dat) & names(genei.dat) %in% imid]) # IM
  ttest.genesin[i,5] <- length(genei.dat[is.infinite(genei.dat) & names(genei.dat) %in% lmid]) # LM
  ttest.genesin[i,6] <- length(genei.dat[is.infinite(genei.dat) & names(genei.dat) %in% mmid]) # MM
} 
ttest.genesin <- as.data.frame(ttest.genesin) 
for(i in 2:6){ttest.genesin[,i]<-as.numeric(ttest.genesin[,i])} # makes df counting NAs

### Repeat below for each methylation subtype (Example shown for HM subtype samples)
exprhm <- log2fc.eac[,colnames(log2fc.eac) %in% hmid]; dim(exprhm) # 25 samples
testgenes.hm.nafilt <- as.character(ttest.genesin[ttest.genesin$inf.HM<16,1]) # NA count filter (<60% sample group size)
testgenes.hm.log2fcfilt <- rownames(exprhm[which(rowMeans(exprhm) <= -2),]) # -2 mean log2fc filter
testgenes.hm.filt <- intersect(testgenes.hm.nafilt,testgenes.hm.log2fcfilt) # genes that meet prefiltering criteria

### Prefilters: Methylation data
exprhm <- exprhm[,order(match(colnames(exprhm),hmid))]; 
gt.hm <- gtn[,gtn$Tissue_Code=="1" & gtn$Participant %in% hmid]; # gt.hm is tumor-only GenomicRatioSet
colnames(gt.hm) <- gt.hm$Participant 
gt.hm <- gt.hm[,order(match(colnames(gt.hm),hmid))]

gtn.hm <- gtn[,gtn$Participant %in% c(hmid)|gtn$Tissue_Code=="11"] # gtn.hm contains normal samples and subtype tumors
annohm450 <- as.data.frame(getAnnotation(gtn.hm))
annohm450.promoter=annohm450[grep("TSS|1stExon|5'UTR",annohm450$UCSC_RefGene_Group),] # isolate promoter region HM450 probes
testvar=gtn.hm$Tissue_Code; 
nprobecutoff=3; # min probe count needed to calculate region methylation, otherwise skipped
nacut.hm = 0.6*length(hmid) # above 60% n Tsubtype cutoff, region methylation not calculated
hm.msummary <- mRgnSummary(testgenes=testgenes.hm.filt,
                           annohm450=annohm450.promoter,
                           testvar=gtn.hm$Tissue_Code,
                           gset=gtn.hm,nprobecutoff,
                           nacut=nacut.hm,
                           roundnum=5,probesfilt=NULL,
                           maxrange=length(testgenes.hm.filt)) # calculate and filter region methylation (mean) by gene

hm.msum <- as.data.frame(hm.msummary); dim(hm.msum); for(i in 2:ncol(hm.msum)){hm.msum[,i]<-as.numeric(as.character(hm.msum[,i]))}
hm.msum <- hm.msum[!hm.msum[,2]<3,]; dim(hm.msum) # remove genes with less than 3 available HM450 promoter region probes

### Biological Filter: Methylation difference (T vs. NM)
methylationcutoff <- 0.2 # choose 0.2 betaval cut (where Tmethylation>Nmethylation, Beta-val scale)
testgenes.hm.mfilt <- as.character(hm.msum[abs(hm.msum[,3]-hm.msum[,4])>methylationcutoff & hm.msum[,3]>hm.msum[,4],1])

### Correlation testing
hm.cor <- exprmethyCor(cormethod="spearman",nprobecut=3,nsampcut=0.6*length(hmid),
                       annohm450=annohm450.promoter,
                       testgenes=testgenes.hm.mfilt,totlen=length(testgenes.hm.mfilt),startlen=1,
                       exprset=exprhm,gset=gt.hm) # Test of repression (correlation between promoter-region methylation and gene expression)

hm.cor$padj <- p.adjust(hm.cor[,5],method="BH") # Benjamini-Hotchberg adjustment for multiple testing
hm.corsig <- hm.cor[hm.cor$padj<0.05,]; 
hm.corsig.genes <- hm.corsig[,1]; hm.corsig.genes<-hm.corsig.genes[!is.na(hm.corsig.genes)] # List of IDs of significantly epigenetically repressed gene

# Repeat above for remaining EAC subtypes (IM, LM, MM)

#
