# Purpose: Get Epigenetically and Epitranscriptomically Repressed Genes
# Author: SKM

# Sample Methylation Subtypes
hmid <- c('A43E','A4OE','A4OI','A4OJ','A4ON','A4OX','A6DN','A6F8','A6FB','A6VZ',
          'A7BO','A893','A8NJ','A8NM','A8NR','A8NS','A8NT','A8NV','A8NW','A8WC',
          'A93C','A93E','A9GI','A9GL','A9GR')
imid <- c('A43I','A4OG','A4OH','A4OR','A4OU','A4OW','A5M8','A6DQ','A6FW','A6RE',
          'A6XG','A7RB','A88V','A88Y','A8EQ','A8NG','A8NH','A8NL','A8W8','A938',
          'A9GF','A9GH','A9GK','A9GM','AAAR','AASW')
lmid <- c('A43C','A4OF','A4OO','A4OP','A4OT','A4QS','A6BV','A6KZ','A6L6','A6XQ',
          'A6Y2','A7DP','A8NN','A8NU','A8WG','A939','A93D','A9CJ','A9GO','AA4D')
mmid <- c('A43M','A4OQ','A4OS','A6FH','A6L4','A6Y0','A88T','A891','A8NE','A8NF',
          'A8NI','A8W5','A9GG','A9GJ','A9GN','A9GQ')
nhmid <- c(mmid,lmid,imid)

#=====================
# Starting Datasets
# Methylation Data: gt/gtn (N = 87 EAC patients) from The Cancer Genome Atlas files.
# file location: https://portal.gdc.cancer.gov/

# Expression Data: log2fc from counts for N = 87 EAC patients, available on Firehose (step 1a).  
# file location: http://gdac.broadinstitute.org/runs/stddata__2016_01_28/data/ESCA/20160128
log2fc.eac <- mrna.eac.t.log2fc; colnames(log2fc.eac)<-substr(colnames(log2fc.eac),9,12) # expr for T samples

#================================================================
# Preprocessing and Filtering
# determine genes to include for ttest, count number of Inf vals
ttest.genesin <- matrix(nrow=nrow(log2fc.eac),ncol=6) 
colnames(ttest.genesin)<-c("geneID","inf.HM","inf.NHM","inf.IM","inf.LM","inf.MM")
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
for(i in 2:6){ttest.genesin[,i]<-as.numeric(ttest.genesin[,i])}

# Subtype-specific example show for HM samples:
# Filter 2a. n missing /inf <= 0.4*n_samples
testgenes.hm.nafilt <- as.character(ttest.genesin[ttest.genesin$inf.HM<(0.4*length(hmid)),1])
exprhm <- log2fc.eac[,colnames(log2fc.eac) %in% hmid]

#===================================================
# Biological relevance filters and correlation testing
# 3. Filter genes by mean subtype log2FC <= -2 
testgenes.hm.log2fcfilt <- rownames(exprhm[which(rowMeans(exprhm) <= -2),])
testgenes.hm.filt <- intersect(testgenes.hm.nafilt,testgenes.hm.log2fcfilt)

# 4. Methylation data filter on available CpGs and mean methylation
# gt and gtn are preprocessed/normalized GenomicRatioSets with Tumor and Tumor+Normal Matched samples
exprhm <- exprhm[,order(match(colnames(exprhm),hmid))]; identical(colnames(exprhm),hmid) 
gt.hm <- gt[,gt$Participant %in% hmid]
colnames(gt.hm) <- gt.hm$Participant 
gt.hm <- gt.hm[,order(match(colnames(gt.hm),hmid))] 
identical(colnames(gt.hm),colnames(exprhm)) # expr and methylation colnames should match
gtn.hm <- gtn[,gtn$Participant %in% c(hmid)|gtn$Tissue_Code=="11"] # gtn.hm contains normal samples and subtype tumors
annohm450 <- as.data.frame(getAnnotation(gtn.hm))
# isolate promoter region HM450 probes
annohm450.promoter=annohm450[grep("TSS|1stExon|5'UTR",annohm450$UCSC_RefGene_Group),] 
testvar=gtn.hm$Tissue_Code 
nprobecutoff=3 # available CpG count filter
nacut.hm = 0.6*length(hmid) # above 60% n Tsubtype cutoff, region methylation not calculated
# calculate and filter region methylation (mean) by gene
hm.msummary <- mRgnSummary(testgenes=testgenes.hm.filt,
                           annohm450=annohm450.promoter,
                           testvar=gtn.hm$Tissue_Code,
                           gset=gtn.hm,nprobecutoff,
                           nacut=nacut.hm,
                           roundnum=5,probesfilt=NULL,
                           maxrange=length(testgenes.hm.filt)) 

hm.msum <- as.data.frame(hm.msummary); dim(hm.msum); for(i in 2:ncol(hm.msum)){hm.msum[,i]<-as.numeric(as.character(hm.msum[,i]))}
# remove genes with less than 3 available HM450 promoter region probes
hm.msum <- hm.msum[!hm.msum[,2]<3,]; dim(hm.msum) 
methylationcutoff <- 0.2 # choose 0.2 betaval cut (where Tmethylation>Nmethylation, Beta-val scale)
testgenes.hm.mfilt <- as.character(hm.msum[abs(hm.msum[,3]-hm.msum[,4])>methylationcutoff & hm.msum[,3]>hm.msum[,4],1])

# 6. Run correlations expr:methylation at proximal promoters
hm.cor <- exprmethyCor(cormethod="spearman",nprobecut=3,nsampcut=0.6*length(hmid),
                       annohm450=annohm450.promoter,
                       testgenes=testgenes.hm.mfilt,totlen=length(testgenes.hm.mfilt),startlen=1,
                       exprset=exprhm,gset=gt.hm) # Test of repression (correlation between promoter-region methylation and gene expression)
hm.cor2 <- hm.cor[!is.na(hm.cor[,4]),]
hm.cor2$padj <- p.adjust(hm.cor2$cor.praw,method="BH")
hm.corgenes <- as.character(hm.cor2[hm.cor2$padj<0.05 & hm.cor2$meanexpr.log2fc< -2,]$gene)
dim(hm.cor2[hm.cor2$padj<0.05 & hm.cor2$cor.estimate<0,]) # epigenetically repressed genes

# 7. repeat above for all subtypes/sample groups
# 8. compare repressed genes

#************
