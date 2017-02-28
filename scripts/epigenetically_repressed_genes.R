# scripts for supplementary methods
# author: SKM
# initDate: 02Feb2017

# Dependancies:
# github functions: https://github.com/metamaden/HM450_tools
# minfi, sva


#=======================================================
# Supp Appendix #. TCGA Methylation Subtype Identifiers
#=======================================================
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

#======================================================================================
# Supp Fig #. Main Fig # - Differentially repressed genes across methylation subtypes
#======================================================================================
### Numbering below reflects step in methods flow chart (Supp. Fig. #)

### Steps 1b and 2b.
load(paste0(getwd(),"/gtn.rda")) # gtn a GenomicRatioSet containing TCGA tumor and normal EACs, preprocessed as above. 
# (optional): Batch correct subtype-filtered gtn ("gtn.hm" below) if clusters observed on batch, not tissue type

### 1a. Mapping and read depth optimization and RSEM filters as in (citation).
load(paste0(getwd(),"/log2fc_eac.rda") # log2fc_eac is tumor gene expression from Firehose-preprocessed dataset

### 2a. 
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

### Steps 3 and 4
exprhm <- exprhm[,order(match(colnames(exprhm),hmid))]; 
identical(colnames(exprhm),hmid) 
gt.hm <- gtn[,gtn$Tissue_Code=="1" & gtn$Participant %in% hmid]; # gt.hm is tumor-only GenomicRatioSet
colnames(gt.hm) <- gt.hm$Participant 
gt.hm <- gt.hm[,order(match(colnames(gt.hm),hmid))]
identical(colnames(gt.hm),colnames(exprhm)) # expr and methylation colnames should match

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

### Step 5
methylationcutoff <- 0.2 # choose 0.2 betaval cut (where Tmethylation>Nmethylation, Beta-val scale)
testgenes.hm.mfilt <- as.character(hm.msum[abs(hm.msum[,3]-hm.msum[,4])>methylationcutoff & hm.msum[,3]>hm.msum[,4],1])

### Step 6
hm.cor <- exprmethyCor(cormethod="spearman",nprobecut=3,nsampcut=0.6*length(hmid),
                       annohm450=annohm450.promoter,
                       testgenes=testgenes.hm.mfilt,totlen=length(testgenes.hm.mfilt),startlen=1,
                       exprset=exprhm,gset=gt.hm) # Test of repression (correlation between promoter-region methylation and gene expression)

hm.cor$padj <- p.adjust(hm.cor[,5],method="BH") # Benjamini-Hotchberg adjustment for multiple testing
hm.corsig <- hm.cor[hm.cor$padj<0.05,]; 
hm.corsig.genes <- hm.corsig[,1]; hm.corsig.genes<-hm.corsig.genes[!is.na(hm.corsig.genes)] # List of IDs of significantly epigenetically repressed gene

# Step 7. Repeat above for remaining groupsings and subtypes (IM, LM, MM, and NHM)
# Step 8. Comparing repressed genes (Main Fig. ##)
# Step 9. Pathway analyses of differentially repressed genes
