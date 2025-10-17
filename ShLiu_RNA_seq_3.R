#2024
#Libraries
###############################################
library(DESeq2)
library(ggplot2)
library(stringr)
library(vegan)
library(RColorBrewer)
??DESeq2

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("MatrixGenerics")
BiocManager::install("DESeq2")
#Functions
###############################################
#create PCOA plot 
plot_pcoa <-function(pcoa, conditions, shapeby, title) {
  eig<- eigenvals(pcoa)
  prop<-eig/sum(eig)
  PCOA1 <- paste("PCoA1 ",100*round(prop[1],3),"%")
  PCOA2 <- paste("PCoA2 ",100*round(prop[2],3),"%")
  
  pcoa.sum <- summary(pcoa)
  pcoa.sum.sites  <- data.frame(pcoa.sum$sites[,1:2])       # PC1 and PC2
  pcoa.sum.species  <- data.frame(pcoa.sum$species[,1:2])     # loadings for PC1 and PC2
  
  if(!is.null(conditions)) {
    pcoa.sum.sites["colorby"] <- conditions
  }
  else { pcoa.sum.sites["colorby"] = rep("1", nrow(pcoa.sum.sites)) }
  
  if(!is.null(shapeby)) {
    pcoa.sum.sites["shapeby"] <- shapeby
  }
  else { pcoa.sum.sites["shapeby"] = rep("1", nrow(pcoa.sum.sites)) }
  
  
  pcoa_plot <- ggplot(pcoa.sum.sites, aes(x=MDS1, y=MDS2)) + 
    geom_jitter(aes(color = colorby, shape = shapeby), size=4) +
    geom_hline(yintercept=0, linetype="dotted") +
    geom_vline(xintercept=0, linetype="dotted") +
    theme_bw() +
    theme( panel.grid.major.x = element_blank(),
           panel.grid.minor.x = element_blank(),
           panel.grid.major.y = element_blank(),
           panel.grid.minor.y = element_blank(),
           text = element_text(size=14), 
           legend.background = element_blank(), 
           legend.title = element_blank(), 
           legend.key = element_blank()) +
    ggtitle(title) + 
    xlab (PCOA1) + 
    ylab (PCOA2) + 
    stat_ellipse(level = 0.95, aes(group = colorby)) +
    coord_fixed()
  
  return(pcoa_plot)
}

CompileResults <- function(res, padj=0.05, l2fc=2, contrast) {
  #convert to dataframe
  res.df <- data.frame(res)
  #order the dataframe with most expressed genes first
  res.df <- res.df[order(res.df$log2FoldChange, decreasing = T), ]
  
  #Filter the significant DESeq2 Results
  res.df <- res.df[which(abs(res.df$log2FoldChange) > l2fc & res.df$padj < padj), ]
  res.df$Sequence <- factor(rownames(res.df), levels = rev(rownames(res.df)))
  
  #Change column order
  res.df <- res.df[c("Sequence", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")]
  
  #Add annotations to dataframe
  Annotations.contrast <- Annotations[match(rownames(res.df), Annotations$Transcript), ]
  res.df <- cbind(res.df, Annotations.contrast)
}


'%!in%' <- function(x,y)!('%in%'(x,y))


# Import Data
###############################################
setwd("/Users/Sherry/Desktop/RNA_bioinformatics")

#Import Metadata (prepared by Shuang Liu)
metadata <- read.delim("PS_gill_RNA_metafile3.txt", sep="\t")
tail(metadata)
head(metadata)
#import TPM data
TPM_quantmerge <- read.delim("TPM_quantmerge.filter.txt",
                             header=T)
rownames(TPM_quantmerge) <- TPM_quantmerge$Name; TPM_quantmerge$Name <- NULL

#Sort the dataset and keep the top 10,000 transcripts
TPM_quantmerge <- TPM_quantmerge[order(rowSums(TPM_quantmerge), decreasing = T), ]

#Remove any transcripts with < 1 mean TPM across samples
TPM_quantmerge.keep <- TPM_quantmerge[which(rowMeans(TPM_quantmerge) > 1.0), ]
nrow(TPM_quantmerge.keep)
#Import Annotations
Annotations <- read.delim("Annotation/est_transcripts95.annotation.tsv",
                          header = T)
Annotations$Transcript <- gsub("\\.p[1-9]+", "", Annotations$Peptide)
Annotations <- Annotations[!duplicated(Annotations$Transcript), ]
head(Annotations)

###############################################
#DeSeq2
###############################################

#import numReads data
quantmerge <- read.delim("est_numreads_quantmerge.txt",
                         header=T)
head(quantmerge)
rownames(quantmerge) <- quantmerge$Name; quantmerge$Name <- NULL

#Match numeric data with filtered TPM dataset
quantmerge <- quantmerge[which(rownames(quantmerge) %in% rownames(TPM_quantmerge.keep)), match(metadata$Sample, colnames(quantmerge))]
nrow(quantmerge)
head(quantmerge)
#Convert to integer
quantmerge <- as.matrix(quantmerge)
mode(quantmerge) <- "integer"
quantmerge[is.na(quantmerge)] <- 0
quantmerge.df <- data.frame(quantmerge)
write.table (quantmerge.df,
             "new2024/quantmergedata.txt",
             sep="\t", row.names=F, quote=F)
head(quantmerge)
#Run DESeq2 Analysis
#############################################
quantmerge.filter <- quantmerge[, which(colnames(quantmerge) %!in% c("CR4SW", "LCR1SW"))]
quantmerge.filter <- quantmerge.filter[which(rowSums(quantmerge.filter) != 0), ]
metadata.filter <- metadata[which(metadata$Sample %!in% c("CR4SW", "LCR1SW")), ]

dds <- DESeqDataSetFromMatrix(countData = quantmerge.filter, 
                              colData = metadata.filter, 
                              design = ~ Region + Treatment + Region:Treatment) #
dds <- DESeq(dds)

resultsNames(dds)

#Extract Normalized count data

dds.normcounts <- counts(dds, normalized=T)
normcounts.df <- data.frame(dds.normcounts)
normcounts.df$Sequence <- rownames(normcounts.df)
#nrow(dds.normcounts)
head(dds.normcounts)
write.table(normcounts.df,
            "new2024/normcounts.txt",
            sep="\t", row.names=F, quote=F)
sizeFactors(dds)

write.table(dds.normcounts,
            "new2024/dds.normcounts.txt",
            sep="\t", row.names=F, quote=F)
#Extract the DESeq2 Results and order (interaction term)
#############################################
#Extract Results
res.cRintercL <- results(dds, name="RegioncoastR.TreatmentSW")
res.iLintercL <- results(dds, name="RegioninlandL.TreatmentSW")
res.cRinteriL <- results(dds, contrast=list("RegioncoastR.TreatmentSW","RegioninlandL.TreatmentSW"))

#Compile Results

res.cRintercL.df <- CompileResults(res.cRintercL, contrast="cR.cL", padj=0.05, l2fc=2)
res.iLinercL.df <- CompileResults(res.iLintercL, contrast="iL.cL", padj=0.05, l2fc=2)
res.cRinteriL.df <- CompileResults(res.cRinteriL, contrast="cR.iL", padj=0.05, l2fc=2)
res.SW.AW.df <- CompileResults(res.SW.AW, contrast="SW.AW")

nrow(res.SW.AW.df[which(res.SW.AW.df$log2FoldChange > 0), ])
nrow(res.SW.AW.df[which(res.SW.AW.df$log2FoldChange < 0), ])

write.table(res.cRintercL.df,
            "new2024/results/cRinterCL.txt",
            sep="\t", row.names=F, quote=F)
write.table(res.cRintercL.df,
            "new2024/results/iLinterCL.txt",
            sep="\t", row.names=F, quote=F)
write.table(res.cRintercL.df,
            "new2024/results/cRinteriL.txt",
            sep="\t", row.names=F, quote=F)

#Isolate significant sequences
SigSequences <- c(rownames(res.cRintercL.df),
                  rownames(res.iLintercL.df),
                  rownames(res.cRinteriL.df))

SigSequences <- SigSequences[!duplicated(SigSequences)]
length(SigSequences)
write.table(SigSequences,
            "new2024/SigSequences.inter.txt",
            sep="\t", row.names=F, quote=F)
#Annotations for significant sequences (look for GO terms)
Sig.Annotations <- Annotations[which(Annotations$Transcript %in% SigSequences), ]
RefSeq <- gsub(" .*", "", Sig.Annotations$STICKLEBACK.UNIPROT)
UniProt <- gsub("tr\\|", "", Sig.Annotations$UNIPROT)
UniProt <- gsub("\\|.*", "", UniProt)
write.table(data.frame(Transcript = Sig.Annotations$Transcript,
                       RefSeq, 
                       UniProt, 
                       PFAM = Sig.Annotations$PFAMID, 
                       PGAMGO = Sig.Annotations$PFAMGO), 
            "new2024/Sig.Annotations.inter.tsv",
            sep="\t", row.names = F, quote = F)
##############################################################################
##Group the variables
dds$group <- factor(paste0(dds$Region, dds$Treatment))
design(dds) <- ~group
dds <- DESeq(dds)
resultsNames(dds)
#Extract Normalized count data for Group 
dds.normcounts.group <- counts(dds, normalized=T)
normcounts.group.df <- data.frame(dds.normcounts.group)
normcounts.group.df$Sequence <- rownames(normcounts.group.df)
#nrow(dds.normcounts)
head(dds.normcounts.group)
write.table(normcounts.df,
            "new2024/normcounts.group.txt",
            sep="\t", row.names=F, quote=F)
sizeFactors(dds)
##group the variables and compile the results######
res.CR.AW.SW <- results(dds, contrast=c("group", "coastRSW", "coastRAW"))
res.CR.AW.SW.df <- CompileResults(res.CR.AW.SW, contrast="SW.AW", padj=0.05, l2fc=2)
write.table(res.CR.AW.SW.df,
            "new2024/results/CR.AW.SW.txt",
            sep="\t", row.names=F, quote=F)

res.CL.AW.SW <- results(dds, contrast=c("group", "coastLSW", "coastLAW"))
res.CL.AW.SW.df <- CompileResults(res.CL.AW.SW, contrast="SW.AW", padj=0.05, l2fc=2)
write.table(res.CL.AW.SW.df,
            "new2024/results/CL.AW.SW.txt",
            sep="\t", row.names=F, quote=F)

res.IL.AW.SW <- results(dds, contrast=c("group", "inlandLSW", "inlandLAW"))
res.IL.AW.SW.df <- CompileResults(res.IL.AW.SW, contrast="SW.AW", padj=0.05, l2fc=2)
write.table(res.IL.AW.SW.df,
            "new2024/results/IL.AW.SW.txt",
            sep="\t", row.names=F, quote=F)


################################################################################
##likelihood ratio test for interaction term #########
ddslrt <- DESeq(dds, test = "LRT", reduced = ~ Region + Treatment)
reslrt <- results(ddslrt)
head(reslrt)
reslrt$Sequence <-factor(rownames(reslrt), levels = rev(rownames(reslrt)))
reslrt.df <- data.frame(reslrt)
reslrt.df <- reslrt.df[c("Sequence", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")]
head(reslrt.df)
reslrt.df$Sequence <- rownames(reslrt.df) 
write.table(reslrt.df,
            "new2024/results/interaction/LRTresults.inter.txt",
            sep="\t", row.names=F, quote=F)
Annotations.inter.lrt <- Annotations[match(rownames(reslrt.df), Annotations$Transcript), ]
head (Annotations.inter.lrt)
tail (Annotations.inter.lrt)
reslrt.inter.combine <- cbind(reslrt, Annotations.inter.lrt)
write.table(reslrt.inter.combine,
            "new2024/results/interaction/LRTresults.lrt.inter.txt",
            sep="\t", row.names=F, quote=F)
nrow(reslrt.inter.combine)

##group padj < 0.05
reslrt.keep <- reslrt.df[which(reslrt.df$padj < 0.05), ]
write.table(reslrt.keep,
            "new2024/LRTresultskeep.txt",
            sep="\t", row.names=F, quote=F)
### extract abs(log2fold change) >2
reslrt.keep3 <- reslrt.keep[which(abs(reslrt.keep$log2FoldChange) > 2), ]
head(reslrt.keep3)
nrow(reslrt.keep3)
reslrt.keep3.df <- data.frame(reslrt.keep3)
reslrt.keep3.df$Sequence <- factor(rownames(reslrt.keep3.df), levels = rev(rownames(reslrt.keep3.df)))
write.table(reslrt.keep3.df,
            "new2024/results/LRTresultskeep3.txt",
            sep="\t", row.names=F, quote=F)
##change the order of variables####
reslrt.keep <- reslrt.keep[c("Sequence", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")]
###Add annotations to dataframe
reslrt.keep.df <- data.frame(reslrt.keep)
reslrt.keep.df$Sequence <- factor(rownames(reslrt.keep.df), levels = rev(rownames(reslrt.keep.df)))
head(reslrt.keep.df)
tail(reslrt.keep.df)
Annotations.contrast <- Annotations[match(rownames(reslrt.keep.df), Annotations$Transcript), ]
head (Annotations.contrast)
tail (Annotations.contrast)
reslrt.keep.combine <- cbind(reslrt.keep, Annotations.contrast)
write.table(reslrt.keep.combine,
            "new2024/LRTresultskeepannot.txt",
            sep="\t", row.names=F, quote=F)
head(reslrt.keep.combine)
tail(reslrt.keep.combine)
nrow(reslrt.keep.df)
###annotation###########
sigseqlrt <- c(rownames(reslrt.keep.df))
length(sigseqlrt)
head(sigseqlrt)
sigseqlrt.sig <- sigseqlrt[!duplicated(sigseqlrt)]
length(sigseqlrt.sig)
sigseqlrt.sig.df <- data.frame(sigseqlrt.sig)
write.table(sigseqlrt.sig.df,
            "new2024/sigseqlrt.sig.txt",
            sep="\t", row.names=F, quote=F)
Sig.Annotations.inter <- Annotations[which(Annotations$Transcript %in% sigseqlrt.sig), ]
RefSeq <- gsub(" .*", "", Sig.Annotations.inter$STICKLEBACK.UNIPROT)
UniProt <- gsub("tr\\|", "", Sig.Annotations.inter$UNIPROT)
UniProt <- gsub("\\|.*", "", UniProt)

write.table(data.frame(Transcript = Sig.Annotations.inter$Transcript,
                       RefSeq, 
                       UniProt, 
                       PFAM = Sig.Annotations.inter$PFAMID, 
                       PGAMGO = Sig.Annotations.inter$PFAMGO), 
            "new2024/Sig.annotations.inter.tsv",
            sep="\t", row.names = F, quote = F)
### add annotation for the ones abs(log2fold chhange >2)
sigseqlrt2 <- c(rownames(reslrt.keep3.df))
length(sigseqlrt2)
head(sigseqlrt2)
sigseqlrt2.sig.df <- data.frame(sigseqlrt2)
write.table(sigseqlrt2.sig.df,
            "new2024/results/sigseqlrt.sig.txt",
            sep="\t", row.names=F, quote=F)
Sig.Annotations2.inter <- Annotations[which(Annotations$Transcript %in% sigseqlrt2), ]
RefSeq <- gsub(" .*", "", Sig.Annotations2.inter$STICKLEBACK.UNIPROT)
UniProt <- gsub("tr\\|", "", Sig.Annotations2.inter$UNIPROT)
UniProt <- gsub("\\|.*", "", UniProt)

write.table(data.frame(Transcript = Sig.Annotations2.inter$Transcript,
                       RefSeq, 
                       UniProt, 
                       PFAM = Sig.Annotations2.inter$PFAMID, 
                       PGAMGO = Sig.Annotations2.inter$PFAMGO), 
            "new2024/results/Sig.annotations2.inter.tsv",
            sep="\t", row.names = F, quote = F)

#####likelihood ratio test for region effect########################
quantmerge.filter <- quantmerge[, which(colnames(quantmerge) %!in% c("CR4SW", "LCR1SW"))]
quantmerge.filter <- quantmerge.filter[which(rowSums(quantmerge.filter) != 0), ]
metadata.filter <- metadata[which(metadata$Sample %!in% c("CR4SW", "LCR1SW")), ]

dds1 <- DESeqDataSetFromMatrix(countData = quantmerge.filter, 
                              colData = metadata.filter, 
                              design = ~ Region + Treatment) #
dds1 <- DESeq(dds1)
resultsNames(dds1)

ddslrtregion <- DESeq(dds1, test = "LRT", reduced = ~ Treatment)
reslrtregion <- results(ddslrtregion)
reslrtregion$Sequence <- factor(rownames(reslrtregion), levels = rev(rownames(reslrtregion)))
head(reslrtregion)
reslrtregion.df <- data.frame(reslrtregion)
reslrtregion.df$Sequence <- rownames(reslrtregion.df) 
write.table(reslrtregion.df,
            "new2024/LRTresultsregion.txt",
            sep="\t", row.names=F, quote=F)
reslrtregion.df<- reslrtregion.df[c("Sequence", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")]
Annotations.region.lrt <- Annotations[match(rownames(reslrtregion.df), Annotations$Transcript), ]
head (Annotations.region.lrt)
tail (Annotations.region.lrt)
reslrtregion.combine <- cbind(reslrtregion, Annotations.region.lrt)
tail(reslrtregion.combine)
write.table(reslrtregion.combine,
            "new2024/results/SplitByRegion/LRTresults.regionannot.txt",
            sep="\t", row.names=F, quote=F)
nrow(reslrtregion.combine)
##group padj < 0.05
reslrtregion.keep <- reslrtregion.df[which(reslrtregion.df$padj < 0.05), ]
write.table(reslrtregion.keep,
            "new2024/LRTresultskeepregion.txt",
            sep="\t", row.names=F, quote=F)
#### extracat the ones abs(log2foldchange)>2####
reslrtregion2 <- reslrtregion.keep.df[which(abs(reslrtregion.keep.df$log2FoldChange) >2), ]
nrow(reslrtregion2)
head (reslrtregion2)
reslrtregion2.keep.df <- data.frame (reslrtregion2)
reslrtregion2.keep.df$Sequence <- factor(rownames(reslrtregion2.keep.df), levels = rev(rownames(reslrtregion2.keep.df)))
write.table(reslrtregion2.keep.df,
            "new2024/results/LRTresultskeepregion2.txt",
            sep="\t", row.names=F, quote=F)
##change the order of variables####
reslrtregion.keep <- reslrtregion.keep[c("Sequence", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")]
###Add annotations to dataframe
reslrtregion.keep.df <- data.frame(reslrtregion.keep)
reslrtregion.keep.df$Sequence <- factor(rownames(reslrtregion.keep.df), levels = rev(rownames(reslrtregion.keep.df)))
nrow(reslrtregion.keep.df)
head(reslrtregion.keep.df)
tail(reslrtregion.keep.df)
Annotations.contrast <- Annotations[match(rownames(reslrtregion.keep.df), Annotations$Transcript), ]
head (Annotations.contrast)
tail (Annotations.contrast)
reslrtregion.keep.combine <- cbind(reslrtregion.keep, Annotations.contrast)
write.table(reslrtregion.keep.combine,
            "new2024/LRTresultskeepregionannot.txt",
            sep="\t", row.names=F, quote=F)

###annotation for region lrt###########
sigseqlrt.region <- c(rownames(reslrtregion.keep.df))
length(sigseqlrt.region)
head(sigseqlrt.region)
sigseqlrt.region.sig <- sigseqlrt.region[!duplicated(sigseqlrt.region)]
length(sigseqlrt.region.sig)
sigseqlrt.region.df <- data.frame(sigseqlrt.region)
write.table(sigseqlrt.region.df,
            "new2024/sigseqlrt.region.sig.txt",
            sep="\t", row.names=F, quote=F)
Sig.Annotations.region <- Annotations[which(Annotations$Transcript %in% sigseqlrt.region), ]
RefSeq <- gsub(" .*", "", Sig.Annotations.region$STICKLEBACK.UNIPROT)
UniProt <- gsub("tr\\|", "", Sig.Annotations.region$UNIPROT)
UniProt <- gsub("\\|.*", "", UniProt)

write.table(data.frame(Transcript = Sig.Annotations.region$Transcript,
                       RefSeq, 
                       UniProt, 
                       PFAM = Sig.Annotations.region$PFAMID, 
                       PGAMGO = Sig.Annotations.region$PFAMGO), 
            "new2024/Sig.annotations.region.tsv",
            sep="\t", row.names = F, quote = F)

###annotation for region lrt for the ones abs(log2foldchange) >2###########
sigseqlrt2.region <- c(rownames(reslrtregion2.keep.df))
length(sigseqlrt2.region)
head(sigseqlrt2.region)
sigseqlrt2.region.df <- data.frame(sigseqlrt2.region)
write.table(sigseqlrt2.region.df,
            "new2024/results/sigseqlrt2.region.sig.txt",
            sep="\t", row.names=F, quote=F)
Sig.Annotations2.region <- Annotations[which(Annotations$Transcript %in% sigseqlrt2.region), ]
RefSeq <- gsub(" .*", "", Sig.Annotations2.region$STICKLEBACK.UNIPROT)
UniProt <- gsub("tr\\|", "", Sig.Annotations2.region$UNIPROT)
UniProt <- gsub("\\|.*", "", UniProt)

write.table(data.frame(Transcript = Sig.Annotations2.region$Transcript,
                       RefSeq, 
                       UniProt, 
                       PFAM = Sig.Annotations2.region$PFAMID, 
                       PGAMGO = Sig.Annotations2.region$PFAMGO), 
            "new2024/results/Sig.annotations2.region.tsv",
            sep="\t", row.names = F, quote = F)
######likelihood ration test for salinity effect ##################
quantmerge.filter <- quantmerge[, which(colnames(quantmerge) %!in% c("CR4SW", "LCR1SW"))]
quantmerge.filter <- quantmerge.filter[which(rowSums(quantmerge.filter) != 0), ]
metadata.filter <- metadata[which(metadata$Sample %!in% c("CR4SW", "LCR1SW")), ]

dds2 <- DESeqDataSetFromMatrix(countData = quantmerge.filter, 
                               colData = metadata.filter, 
                               design = ~ Treatment + Region) #
dds2 <- DESeq(dds2)
resultsNames(dds2)

ddsltsalinity <- DESeq(dds2, test = "LRT", reduced = ~ Region)
resltsalinity <- results(ddsltsalinity)
resltsalinity.df <- data.frame(resltsalinity)
resltsalinity.df $Sequence <-rownames(resltsalinity.df)
ddsltsalinity.normcounts <- counts(ddsltsalinity, normalized=T)
head(resltsalinity.df)
ddsltsalinity.normcounts.df <- data.frame(ddsltsalinity.normcounts)
ddsltsalinity.normcounts.df$Sequence <- rownames(ddsltsalinity.normcounts)
#nrow(dds.normcounts)
write.table(normcounts.treat.df,
            "normcounts.treat.txt",
            sep="\t", row.names=F, quote=F)
resltsalinity.df <- data.frame(resltsalinity)
resltsalinity.df$Sequence <- rownames(resltsalinity.df) 
write.table(resltsalinity.df,
            "new2024/LRTresultssalinity.txt",
            sep="\t", row.names=F, quote=F)

write.table(ddsltsalinity.normcounts.df,
            "new2024/LRTresultssalinity.normcount.txt",
            sep="\t", row.names=F, quote=F)
## add annotations to LRT salinity######
resltsalinity.df <- resltsalinity.df[c("Sequence", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")]
Annotations.lrt.salinity <- Annotations[match(rownames(resltsalinity.df), Annotations$Transcript), ]
head (Annotations.lrt.salinity)
tail (Annotations.lrt.salinity)
resltsalinity.annotation.combine <- cbind(resltsalinity.df, Annotations.lrt.salinity)
write.table(resltsalinity.annotation.combine,
          "new2024/results/splitBySalinity/LRTresults.salinity.annot.txt",
           sep="\t", row.names=F, quote=F)
##group padj < 0.05
resltsalinity.keep <- resltsalinity.df[which(resltsalinity.df$padj < 0.05), ]
head(resltsalinity.keep)
write.table(resltsalinity.keep,
            "new2024/LRTresultskeepsalinity.txt",
            sep="\t", row.names=F, quote=F)

#### extract the ones abs(log2foldchange)>2
resltsalinity2 <- resltsalinity.keep[which(abs(resltsalinity.keep$log2FoldChange) >2 ), ]
nrow(resltsalinity2)
resltsalinity2.df <- data.frame (resltsalinity2)
resltsalinity2.df$Sequence <- factor(rownames(resltsalinity2.df), levels = rev(rownames(resltsalinity2.df)))
write.table(resltsalinity2.df,
            "new2024/results/LRTresultskeepsalinity2.txt",
            sep="\t", row.names=F, quote=F)

##change the order of variables####
resltsalinity.keep <- resltsalinity.keep[c("Sequence", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")]
###Add annotations to dataframe
#resltsalinity.keep.df <- data.frame(resltsalinity.keep)
#resltsalinity.keep.df$Sequence <- factor(rownames(resltsalinity.keep.df), levels = rev(rownames(resltsalinity.keep.df)))
nrow(resltsalinity.keep.df)
#head(resltsalinity.keep.df)
#tail(resltsalinity.keep.df)
#Annotations.contrast <- Annotations[match(rownames(resltsalinity.keep.df), Annotations$Transcript), ]
#head (Annotations.contrast)
#tail (Annotations.contrast)
#resltsalinity.keep.combine <- cbind(resltsalinity.keep, Annotations.contrast)
#write.table(resltsalinity.keep.combine,
#           "new2024/LRTresultskeepsalinityannot.txt",
#            sep="\t", row.names=F, quote=F)

###annotation for salinity lrt###########
sigseqlrt.salinity <- c(rownames(resltsalinity.keep.df))
length(sigseqlrt.salinity)
head(sigseqlrt.salinity)
sigseqlrt.salinity.sig <- sigseqlrt.salinity[!duplicated(sigseqlrt.salinity)]
length(sigseqlrt.salinity.sig)
sigseqlrt.salinity.df <- data.frame(sigseqlrt.salinity)
write.table(sigseqlrt.salinity.df,
            "new2024/sigseqlrt.salinity.sig.txt",
            sep="\t", row.names=F, quote=F)
Sig.Annotations.salinity <- Annotations[which(Annotations$Transcript %in% sigseqlrt.salinity), ]
RefSeq <- gsub(" .*", "", Sig.Annotations.salinity$STICKLEBACK.UNIPROT)
UniProt <- gsub("tr\\|", "", Sig.Annotations.salinity$UNIPROT)
UniProt <- gsub("\\|.*", "", UniProt)

write.table(data.frame(Transcript = Sig.Annotations.salinity$Transcript,
                       RefSeq, 
                       UniProt, 
                       PFAM = Sig.Annotations.salinity$PFAMID, 
                       PGAMGO = Sig.Annotations.salinity$PFAMGO), 
            "new2024/Sig.annotations.salinity.tsv",
            sep="\t", row.names = F, quote = F)

###annotation for salinity lrt for the ones abs(log2foldchange)>2###########
sigseqlrt2.salinity <- c(rownames(resltsalinity2.df))
length(sigseqlrt2.salinity)
head(resltsalinity2.df)
resltsalinity2.df$Sequence <- factor(rownames(resltsalinity2.df), levels = rev(rownames(resltsalinity2.df)))
sigseqlrt2.salinity.keep.df <- data.frame(sigseqlrt2.salinity)
sigseqlrt2.salinity.keep.df$Sequence <- factor(rownames(sigseqlrt2.salinity.keep.df), levels = rev(rownames(sigseqlrt2.salinity.keep.df)))
write.table(sigseqlrt2.salinity.keep.df,
            "new2024/results/sigseqlrt2.salinity.sig.txt",
            sep="\t", row.names=F, quote=F)
Sig.Annotations2.salinity <- Annotations[which(Annotations$Transcript %in% sigseqlrt2.salinity), ]
RefSeq <- gsub(" .*", "", Sig.Annotations2.salinity$STICKLEBACK.UNIPROT)
UniProt <- gsub("tr\\|", "", Sig.Annotations2.salinity$UNIPROT)
UniProt <- gsub("\\|.*", "", UniProt)
head(sigseqlrt2.salinity.keep.df)
write.table(data.frame(Transcript = Sig.Annotations2.salinity$Transcript,
                       RefSeq, 
                       UniProt, 
                       PFAM = Sig.Annotations2.salinity$PFAMID, 
                       PGAMGO = Sig.Annotations2.salinity$PFAMGO), 
            "new2024/results/Sig.annotations2.salinity.tsv",
            sep="\t", row.names = F, quote = F)

###try annother annotation########
Annotations.sal <- Annotations[match(rownames(resltsalinity2.df), Annotations$Transcript), ]
head (Annotations.sal)
tail (Annotations)
reslrtregion.keep.try.combine <- cbind(resltsalinity2.df, Annotations.sal)
write.table(reslrtregion.keep.try.combine,
            "new2024/LRTresultskeepregionannot.try.txt",
            sep="\t", row.names=F, quote=F)
### genes have higher expression level in FW compared in SW
resltsalinity3 <- resltsalinity.keep[which(resltsalinity.keep$log2FoldChange > 2), ]
nrow(resltsalinity3)
resltsalinity3.df <- data.frame (resltsalinity3)
resltsalinity3.df$Sequence <- factor(rownames(resltsalinity3.df), levels = rev(rownames(resltsalinity3.df)))
write.table(resltsalinity3.df,
            "new2024/results/splitBySalinity/LRTresultskeepsalinity3.txt",
            sep="\t", row.names=F, quote=F)
### annotation#########
resltsalinity3.df$Sequence <- rownames(resltsalinity3.df)
sigseqlrt3.salinity <- c(rownames(resltsalinity3.df))
length(sigseqlrt3.salinity)
head(sigseqlrt3.salinity)
resltsalinity3.df$Sequence <- factor(rownames(resltsalinity3.df), levels = rev(rownames(resltsalinity3.df)))
sigseqlrt3.salinity.keep.df <- data.frame(sigseqlrt3.salinity)
sigseqlrt3.salinity.keep.df$Sequence <- factor(rownames(sigseqlrt3.salinity.keep.df), levels = rev(rownames(sigseqlrt3.salinity.keep.df)))
write.table(sigseqlrt3.salinity.keep.df,
            "new2024/results/splitBySalinity/sigseqlrt3.salinity.sig.txt",
            sep="\t", row.names=F, quote=F)
sigseqlrt3.salinity <- read.delim("new2024/results/splitBySalinity/sigseqlrt3.salinity.sig.txt", header=T)
Sig.Annotations3.salinity <- Annotations[which(Annotations$Transcript %in% sigseqlrt3.salinity), ]
RefSeq <- gsub(" .*", "", Sig.Annotations3.salinity$STICKLEBACK.UNIPROT)
UniProt <- gsub("tr\\|", "", Sig.Annotations3.salinity$UNIPROT)
UniProt <- gsub("\\|.*", "", UniProt)
write.table(data.frame(Transcript = Sig.Annotations3.salinity$Transcript,
                       RefSeq, 
                       UniProt, 
                       PFAM = Sig.Annotations3.salinity$PFAMID, 
                       PGAMGO = Sig.Annotations3.salinity$PFAMGO), 
            "new2024/results/splitBySalinity/Sig.annotations.salinity.tsv",
            sep="\t", row.names = F, quote = F)
## another annotation#######
Annotations.salinity3 <- Annotations[match(rownames(resltsalinity3.df), Annotations$Transcript), ]
reslrtsalinity3.keep.try.combine <- cbind(resltsalinity3.df, Annotations.salinity3)
write.table(reslrtsalinity3.keep.try.combine,
            "new2024/results/splitBySalinity/LRTresultskeepsalinityannot3.txt",
            sep="\t", row.names=F, quote=F)
### Venn diagram########
library(ggplot2)
library(ggvenn)
a <- list('Treatment' = c(resltsalinity.df$Sequence),
          'Habitat' = c(reslrtregion.keep.df$Sequence),
          'Interaction' = c(reslrt.keep.df$Sequence))
ggvenn(a, c("Treatment", "Habitat", "Interaction"),
       fill_color = c("blue", "yellow", "red"),
       fill_alpha = 0.5,
       stroke_color = "black",
       stroke_alpha = 1,
       stroke_size = 1.5,
       set_name_size = 6,
       text_size = 4)
ggvenn(a)
#### Venn diagram for the ones abs(log2foldchange)>2##########
library(ggplot2)
library(ggvenn)
a <- list('Treatment' = c(resltsalinity2.df$Sequence),
          'Habitat' = c(reslrtregion2.keep.df$Sequence),
          'Interaction' = c(reslrt.keep3.df$Sequence))
ggvenn(a, c("Treatment", "Habitat", "Interaction"),
       fill_color = c("blue", "yellow", "red"),
       fill_alpha = 0.5,
       stroke_color = "black",
       stroke_alpha = 1,
       stroke_size = 1.5,
       set_name_size = 6,
       text_size = 4)
ggvenn(a)


##### in SW condition, region effect##########################
quantmerge.filter1 <- quantmerge[, which(colnames(quantmerge) %!in% c("CR4SW", "LCR1SW"))]
metadata.sw <- read.delim("PS_gill_RNA_metafile3_SW.txt", sep="\t")
quantmerge.sw <- quantmerge.filter1[, which(colnames(quantmerge.filter1) %in% metadata.sw$Sample)]
head(metadata.sw)
head(quantmerge.sw)
tail(metadata.sw)
ddssw <- DESeqDataSetFromMatrix(countData = quantmerge.sw, 
                               colData = metadata.sw, 
                               design = ~ Region) #
ddssw <- DESeq(ddssw)
resultsNames(ddssw)

#Extract Normalized count data
ddssw.normcounts <- counts(ddssw, normalized=T)
normcounts.sw.df <- data.frame(ddssw.normcounts)
normcounts.sw.df$Sequence <- rownames(normcounts.sw.df)
nrow(dds.normcounts)
head(dds.normcounts)
write.table(normcounts.sw.df,
            "new2024/normcounts.sw.txt",
            sep="\t", row.names=F, quote=F)
sizeFactors(ddsaw)
#Extract Results
res.cR.cL.sw <- results(ddssw, contrast = c("Region","coastR", "coastL"))
res.iL.cL.sw <- results(ddssw, contrast = c("Region","inlandL", "coastL"))
res.cR.iL.sw <- results(ddssw, contrast= c("Region","coastR", "inlandL"))

#Compile Results

res.cR.cL.sw.df <- CompileResults(res.cR.cL.sw, contrast="cR.cL", padj=0.05, l2fc=2)
res.iL.cL.sw.df <- CompileResults(res.iL.cL.sw, contrast="iL.cL", padj=0.05, l2fc=2)
res.cR.iL.sw.df <- CompileResults(res.cR.iL.sw, contrast="cR.iL", padj=0.05, l2fc=2)
#write tables######
write.table(res.cR.cL.sw.df,
            "new2024/results/cR.cL.sw.txt",
            sep="\t", row.names=F, quote=F)

write.table(res.iL.cL.sw.df,
            "new2024/results/iL.cL.sw.txt",
            sep="\t", row.names=F, quote=F)

write.table(res.cR.iL.sw.df,
            "new2024/results/cR.iL.sw.txt",
            sep="\t", row.names=F, quote=F)

###annotation in CR.CL.SW###########
sigseq.CR.CL.sw <- c(rownames(res.cR.cL.sw.df))
length(sigseq.CR.CL.sw)
head(sigseq.CR.CL.sw)
sigseq.CR.CL.sw.sig <- sigseq.CR.CL.sw[!duplicated(sigseq.CR.CL.sw)]
length(sigseq.CR.CL.sw.sig)
sigseq.CR.CL.sw.df <- data.frame(sigseq.CR.CL.sw)
write.table(sigseq.CR.CL.sw.df,
            "new2024/results/splitBySalinity/SW/sigseq.CR.CL.sw.sig.txt",
            sep="\t", row.names=F, quote=F)
Sig.Annotations.CR.CL.sw <- Annotations[which(Annotations$Transcript %in% sigseq.CR.CL.sw), ]
RefSeq <- gsub(" .*", "", Sig.Annotations.CR.CL.sw$STICKLEBACK.UNIPROT)
UniProt <- gsub("tr\\|", "", Sig.Annotations.CR.CL.sw$UNIPROT)
UniProt <- gsub("\\|.*", "", UniProt)

write.table(data.frame(Transcript = Sig.Annotations.CR.CL.sw$Transcript,
                       RefSeq, 
                       UniProt, 
                       PFAM = Sig.Annotations.CR.CL.sw$PFAMID, 
                       PGAMGO = Sig.Annotations.CR.CL.sw$PFAMGO), 
            "new2024/results/splitBySalinity/SW/Sig.annotations.CR.CL.sw.tsv",
            sep="\t", row.names = F, quote = F)
###annotation in IL.CL.SW######
sigseq.IL.CL.sw <- c(rownames(res.iL.cL.sw.df))
length(sigseq.IL.CL.sw)
head(sigseq.IL.CL.sw)
sigseq.IL.CL.sw.sig <- sigseq.IL.CL.sw[!duplicated(sigseq.IL.CL.sw)]
length(sigseq.IL.CL.sw.sig)
sigseq.IL.CL.sw.df <- data.frame(sigseq.IL.CL.sw)
write.table(sigseq.IL.CL.sw.df,
            "new2024/results/splitBySalinity/SW/sigseq.IL.CL.sw.sig.txt",
            sep="\t", row.names=F, quote=F)
Sig.Annotations.IL.CL.sw <- Annotations[which(Annotations$Transcript %in% sigseq.IL.CL.sw), ]
RefSeq <- gsub(" .*", "", Sig.Annotations.IL.CL.sw$STICKLEBACK.UNIPROT)
UniProt <- gsub("tr\\|", "", Sig.Annotations.IL.CL.sw$UNIPROT)
UniProt <- gsub("\\|.*", "", UniProt)

write.table(data.frame(Transcript = Sig.Annotations.IL.CL.sw$Transcript,
                       RefSeq, 
                       UniProt, 
                       PFAM = Sig.Annotations.IL.CL.sw$PFAMID, 
                       PGAMGO = Sig.Annotations.IL.CL.sw$PFAMGO), 
            "new2024/results/splitBySalinity/SW/Sig.annotations.IL.CL.sw.tsv",
            sep="\t", row.names = F, quote = F)
###annotation in CR.IL.SW######
sigseq.CR.IL.sw <- c(rownames(res.cR.iL.sw.df))
length(sigseq.CR.IL.sw)
head(sigseq.CR.IL.sw)
sigseq.CR.IL.sw.sig <- sigseq.CR.IL.sw[!duplicated(sigseq.CR.IL.sw)]
length(sigseq.CR.IL.sw.sig)
sigseq.CR.IL.sw.df <- data.frame(sigseq.CR.IL.sw)
write.table(sigseq.CR.IL.sw.df,
            "new2024/results/splitBySalinity/SW/sigseq.CR.IL.sw.sig.txt",
            sep="\t", row.names=F, quote=F)
Sig.Annotations.CR.IL.sw <- Annotations[which(Annotations$Transcript %in% sigseq.CR.IL.sw), ]
RefSeq <- gsub(" .*", "", Sig.Annotations.CR.IL.sw$STICKLEBACK.UNIPROT)
UniProt <- gsub("tr\\|", "", Sig.Annotations.CR.IL.sw$UNIPROT)
UniProt <- gsub("\\|.*", "", UniProt)

write.table(data.frame(Transcript = Sig.Annotations.CR.IL.sw$Transcript,
                       RefSeq, 
                       UniProt, 
                       PFAM = Sig.Annotations.CR.IL.sw$PFAMID, 
                       PGAMGO = Sig.Annotations.CR.IL.sw$PFAMGO), 
            "new2024/results/splitBySalinity/SW/Sig.annotations.CR.IL.sw.tsv",
            sep="\t", row.names = F, quote = F)

##one way analysis in SW###
ddslrt.sw <- DESeq(ddssw, test = "LRT", reduced = ~ 1)
reslrt.sw <- results(ddslrt.sw)
ddslrt.sw.normcounts <- counts(ddslrt.sw, normalized=T)
head(reslrt.sw)
reslrt.sw.df <- data.frame(reslrt.sw)
reslrt.sw.df$Sequence <-rownames(reslrt.sw.df)
reslrt.sw.df <- reslrt.sw.df[c("Sequence", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")]
head(reslrt.sw.df)
reslrt.sw.df <- data.frame(reslrt.sw)
reslrt.sw.df$Sequence <- rownames(reslrt.sw.df) 
write.table(reslrt.sw.df,
            "new2024/LRTresultsSW.txt",
            sep="\t", row.names=F, quote=F)
write.table(reslrt.sw.df,
            "new2024/results/splitBySalinity/SW/LRTresults.SW.txt",
            sep="\t", row.names=F, quote=F)
### add annotation########
Annotations.lrt.sw <- Annotations[match(rownames(reslrt.sw.df), Annotations$Transcript), ]
head (Annotations.lrt.sw)
tail (Annotations.lrt.sw)
reslrt.sw.combine <- cbind(reslrt.sw.df, Annotations.lrt.sw)
write.table(reslrt.sw.combine,
            "new2024/results/splitBySalinity/SW/LRTresults.annot.sw.txt",
            sep="\t", row.names=F, quote=F)
nrow(reslrt.sw.combine)

##group padj < 0.05
reslrt.sw.keep <- reslrt.sw.df[which(reslrt.sw.df$padj < 0.05), ]
write.table(reslrt.sw.keep,
            "new2024/LRTresultskeep.sw.txt",
            sep="\t", row.names=F, quote=F)
#### extract for the ones abs(log2foldchange)>2####
reslrt2.sw.keep <- reslrt.sw.keep.df[which(abs(reslrt.sw.keep.df$log2FoldChange) >2), ]
nrow(reslrt2.sw.keep)
reslrt2.sw.keep.df <- data.frame(reslrt2.sw.keep)
reslrt2.sw.keep.df$Sequence <- factor(rownames(reslrt2.sw.keep.df), levels = rev(rownames(reslrt2.sw.keep.df)))
write.table(reslrt2.sw.keep,
            "new2024/results/LRTresultskeep2.sw.txt",
            sep="\t", row.names=F, quote=F)

##change the order of variables####
reslrt.sw.keep <- reslrt.sw.keep[c("Sequence", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")]
###Add annotations to dataframe
reslrt.sw.keep.df <- data.frame(reslrt.sw.keep)
reslrt.sw.keep.df$Sequence <- factor(rownames(reslrt.sw.keep.df), levels = rev(rownames(reslrt.sw.keep.df)))
nrow(reslrt.sw.keep.df)
head(reslrt.sw.keep.df)
tail(reslrt.sw.keep.df)
Annotations.contrast.sw <- Annotations[match(rownames(reslrt.sw.keep.df), Annotations$Transcript), ]
head (Annotations.contrast.sw)
tail (Annotations.contrast.sw)
reslrt.sw.keep.combine <- cbind(reslrt.sw.keep, Annotations.contrast.sw)
write.table(reslrt.sw.keep.combine,
            "new2024/LRTresultskeepannot.sw.txt",
            sep="\t", row.names=F, quote=F)

###annotation for lrt in SW better than ablove methods for annotation###############
sigseq.reslrt.sw <- c(rownames(reslrt.sw.keep.df))
length(sigseq.reslrt.sw)
head(sigseq.reslrt.sw)
sigseq.reslrt.sw.sig <- sigseq.reslrt.sw[!duplicated(sigseq.reslrt.sw)]
length(sigseq.reslrt.sw.sig)
sigseq.reslrt.sw.df <- data.frame(sigseq.reslrt.sw)
write.table(sigseq.reslrt.sw.df,
            "new2024/results/splitBySalinity/SW/sigseq.reslrt.sw.sig.txt",
            sep="\t", row.names=F, quote=F)
Sig.Annotations.lrtregion.sw <- Annotations[which(Annotations$Transcript %in% sigseq.reslrt.sw), ]
RefSeq <- gsub(" .*", "", Sig.Annotations.lrtregion.sw$STICKLEBACK.UNIPROT)
UniProt <- gsub("tr\\|", "", Sig.Annotations.lrtregion.sw$UNIPROT)
UniProt <- gsub("\\|.*", "", UniProt)

write.table(data.frame(Transcript = Sig.Annotations.lrtregion.sw$Transcript,
                       RefSeq, 
                       UniProt, 
                       PFAM = Sig.Annotations.lrtregion.sw$PFAMID, 
                       PGAMGO = Sig.Annotations.lrtregion.sw$PFAMGO), 
            "new2024/results/splitBySalinity/SW/Sig.annotations.lrtregion.sw.tsv",
            sep="\t", row.names = F, quote = F)

###annotation for lrt in SW (for the ones abs(log2foldchange)>2 better than ablove methods for annotation###############
sigseq.reslrt2.sw <- c(rownames(reslrt2.sw.keep.df))
length(sigseq.reslrt2.sw)
head(sigseq.reslrt2.sw)
sigseq.reslrt2.sw.df <- data.frame(sigseq.reslrt2.sw)
write.table(sigseq.reslrt2.sw.df,
            "new2024/results/splitBySalinity/SW/sigseq.reslrt2.sw.sig.txt",
            sep="\t", row.names=F, quote=F)
Sig.Annotations.lrtregion2.sw <- Annotations[which(Annotations$Transcript %in% sigseq.reslrt2.sw), ]
RefSeq <- gsub(" .*", "", Sig.Annotations.lrtregion2.sw$STICKLEBACK.UNIPROT)
UniProt <- gsub("tr\\|", "", Sig.Annotations.lrtregion2.sw$UNIPROT)
UniProt <- gsub("\\|.*", "", UniProt)

write.table(data.frame(Transcript = Sig.Annotations.lrtregion2.sw$Transcript,
                       RefSeq, 
                       UniProt, 
                       PFAM = Sig.Annotations.lrtregion2.sw$PFAMID, 
                       PGAMGO = Sig.Annotations.lrtregion2.sw$PFAMGO), 
            "new2024/results/splitBySalinity/SW/Sig.annotations.lrtregion2.sw.tsv",
            sep="\t", row.names = F, quote = F)
###in FW condition, region effect#############################
quantmerge.filter2 <- quantmerge[, which(colnames(quantmerge) %!in% c("CR4SW", "LCR1SW"))]
metadata.aw <- read.delim("PS_gill_RNA_metafile3_AW.txt", sep="\t")
quantmerge.aw <- quantmerge.filter2[, which(colnames(quantmerge.filter2) %in% metadata.aw$Sample)]
head(metadata.aw)
head(quantmerge.aw)
tail(metadata.aw)
ddsaw <- DESeqDataSetFromMatrix(countData = quantmerge.aw, 
                                colData = metadata.aw, 
                                design = ~ Region) #
ddsaw <- DESeq(ddsaw)
resultsNames(ddsaw)

#Extract Normalized count data
ddsaw.normcounts <- counts(ddsaw, normalized=T)
normcounts.aw.df <- data.frame(ddsaw.normcounts)
normcounts.aw.df$Sequence <- rownames(normcounts.aw.df)
nrow(ddsaw.normcounts)
head(dds.normcounts)
write.table(normcounts.aw.df,
            "new2024/normcounts.aw.txt",
            sep="\t", row.names=F, quote=F)
sizeFactors(ddsaw)
## add annotations###########


#Extract Results
res.cR.cL.aw <- results(ddsaw, contrast = c("Region","coastR", "coastL"))
res.iL.cL.aw <- results(ddsaw, contrast = c("Region","inlandL", "coastL"))
res.cR.iL.aw <- results(ddsaw, contrast= c("Region","coastR", "inlandL"))

#Compile Results

res.cR.cL.aw.df <- CompileResults(res.cR.cL.aw, contrast="cR.cL", padj=0.05, l2fc=2)
res.iL.cL.aw.df <- CompileResults(res.iL.cL.aw, contrast="iL.cL", padj=0.05, l2fc=2)
res.cR.iL.aw.df <- CompileResults(res.cR.iL.aw, contrast="cR.iL", padj=0.05, l2fc=2)
#write tables######
write.table(res.cR.cL.aw.df,
            "new2024/results/cR.cL.aw.txt",
            sep="\t", row.names=F, quote=F)

write.table(res.iL.cL.aw.df,
            "new2024/results/iL.cL.aw.txt",
            sep="\t", row.names=F, quote=F)

write.table(res.cR.iL.aw.df,
            "new2024/results/cR.iL.aw.txt",
            sep="\t", row.names=F, quote=F)
###annotation in CR.CL.AW###########
sigseq.CR.CL.aw <- c(rownames(res.cR.cL.aw.df))
length(sigseq.CR.CL.aw)
head(sigseq.CR.CL.aw)
sigseq.CR.CL.aw.sig <- sigseq.CR.CL.aw[!duplicated(sigseq.CR.CL.aw)]
length(sigseq.CR.CL.aw.sig)
sigseq.CR.CL.aw.df <- data.frame(sigseq.CR.CL.aw)
write.table(sigseq.CR.CL.aw.df,
            "new2024/results/splitBySalinity/SW/sigseq.CR.CL.aw.sig.txt",
            sep="\t", row.names=F, quote=F)
Sig.Annotations.CR.CL.aw <- Annotations[which(Annotations$Transcript %in% sigseq.CR.CL.aw), ]
RefSeq <- gsub(" .*", "", Sig.Annotations.CR.CL.aw$STICKLEBACK.UNIPROT)
UniProt <- gsub("tr\\|", "", Sig.Annotations.CR.CL.aw$UNIPROT)
UniProt <- gsub("\\|.*", "", UniProt)

write.table(data.frame(Transcript = Sig.Annotations.CR.CL.aw$Transcript,
                       RefSeq, 
                       UniProt, 
                       PFAM = Sig.Annotations.CR.CL.aw$PFAMID, 
                       PGAMGO = Sig.Annotations.CR.CL.aw$PFAMGO), 
            "new2024/results/splitBySalinity/FW/Sig.annotations.CR.CL.aw.tsv",
            sep="\t", row.names = F, quote = F)
###annotation in IL.CL.SW######
sigseq.IL.CL.aw <- c(rownames(res.iL.cL.aw.df))
length(sigseq.IL.CL.aw)
head(sigseq.IL.CL.aw)
sigseq.IL.CL.aw.sig <- sigseq.IL.CL.aw[!duplicated(sigseq.IL.CL.aw)]
length(sigseq.IL.CL.aw.sig)
sigseq.IL.CL.aw.df <- data.frame(sigseq.IL.CL.aw)
write.table(sigseq.IL.CL.aw.df,
            "new2024/results/splitBySalinity/FW/sigseq.IL.CL.aw.sig.txt",
            sep="\t", row.names=F, quote=F)
Sig.Annotations.IL.CL.aw <- Annotations[which(Annotations$Transcript %in% sigseq.IL.CL.aw), ]
RefSeq <- gsub(" .*", "", Sig.Annotations.IL.CL.aw$STICKLEBACK.UNIPROT)
UniProt <- gsub("tr\\|", "", Sig.Annotations.IL.CL.aw$UNIPROT)
UniProt <- gsub("\\|.*", "", UniProt)

write.table(data.frame(Transcript = Sig.Annotations.IL.CL.aw$Transcript,
                       RefSeq, 
                       UniProt, 
                       PFAM = Sig.Annotations.IL.CL.aw$PFAMID, 
                       PGAMGO = Sig.Annotations.IL.CL.aw$PFAMGO), 
            "new2024/results/splitBySalinity/FW/Sig.annotations.IL.CL.aw.tsv",
            sep="\t", row.names = F, quote = F)
###annotation in CR.IL.SW######
sigseq.CR.IL.aw <- c(rownames(res.cR.iL.aw.df))
length(sigseq.CR.IL.aw)
head(sigseq.CR.IL.aw)
sigseq.CR.IL.aw.sig <- sigseq.CR.IL.aw[!duplicated(sigseq.CR.IL.aw)]
length(sigseq.CR.IL.aw.sig)
sigseq.CR.IL.aw.df <- data.frame(sigseq.CR.IL.aw)
write.table(sigseq.CR.IL.aw.df,
            "new2024/results/splitBySalinity/FW/sigseq.CR.IL.aw.sig.txt",
            sep="\t", row.names=F, quote=F)
Sig.Annotations.CR.IL.aw <- Annotations[which(Annotations$Transcript %in% sigseq.CR.IL.aw), ]
RefSeq <- gsub(" .*", "", Sig.Annotations.CR.IL.aw$STICKLEBACK.UNIPROT)
UniProt <- gsub("tr\\|", "", Sig.Annotations.CR.IL.aw$UNIPROT)
UniProt <- gsub("\\|.*", "", UniProt)

write.table(data.frame(Transcript = Sig.Annotations.CR.IL.aw$Transcript,
                       RefSeq, 
                       UniProt, 
                       PFAM = Sig.Annotations.CR.IL.aw$PFAMID, 
                       PGAMGO = Sig.Annotations.CR.IL.aw$PFAMGO), 
            "new2024/results/splitBySalinity/FW/Sig.annotations.CR.IL.aw.tsv",
            sep="\t", row.names = F, quote = F)

##one way analysis in FW###
ddslrt.aw <- DESeq(ddsaw, test = "LRT", reduced = ~ 1)
reslrt.aw <- results(ddslrt.aw)
head(reslrt.aw)
reslrt.aw.df <- data.frame(reslrt.aw)
reslrt.aw.df$Sequence <- rownames(reslrt.aw.df)
reslrt.aw.df <- reslrt.aw.df[c("Sequence", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")]
write.table(reslrt.aw.df,
            "new2024/results/splitBySalinity/FW/LRTresults.AW.txt",
            sep="\t", row.names=F, quote=F)

### add annotation####
Annotations.lrt.aw <- Annotations[match(rownames(reslrt.aw.df), Annotations$Transcript), ]
head (Annotations.lrt.aw)
tail (reslrt.aw.df)
reslrt.aw.combine <- cbind(reslrt.aw.df, Annotations.lrt.aw)
write.table(reslrt.aw.combine,
            "new2024/results/splitBySalinity/FW/LRTresults.annot.aw.txt",
            sep="\t", row.names=F, quote=F)

nrow(reslrt.aw.combine)
##group padj < 0.05
reslrt.aw.keep <- reslrt.aw.df[which(reslrt.aw.df$padj < 0.05), ]
write.table(reslrt.aw.keep,
            "new2024/LRTresultskeep.aw.txt",
            sep="\t", row.names=F, quote=F)

#### extract the ones abs(log2foldchange)>2 ########
reslrt2.aw.keep <- reslrt.aw.keep.df[which(abs(reslrt.aw.keep.df$log2FoldChange) >2 ), ]
nrow(reslrt2.aw.keep)
reslrt2.aw.keep.df <- data.frame(reslrt2.aw.keep)
write.table(reslrt2.aw.keep,
            "new2024/results/LRTresultskeep2.aw.txt",
            sep="\t", row.names=F, quote=F)
##change the order of variables####
reslrt.aw.keep <- reslrt.aw.keep[c("Sequence", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")]
###Add annotations to dataframe
reslrt.aw.keep.df <- data.frame(reslrt.aw.keep)
reslrt.aw.keep.df$Sequence <- factor(rownames(reslrt.aw.keep.df), levels = rev(rownames(reslrt.aw.keep.df)))
head(reslrt.aw.keep.df)
tail(reslrt.aw.keep.df)
Annotations.contrast.aw <- Annotations[match(rownames(reslrt.aw.keep.df), Annotations$Transcript), ]
head (Annotations.contrast.aw)
tail (Annotations.contrast.aw)
reslrt.aw.keep.combine <- cbind(reslrt.aw.keep, Annotations.contrast.aw)
write.table(reslrt.aw.keep.combine,
            "new2024/LRTresultskeepannot.aw.txt",
            sep="\t", row.names=F, quote=F)
head(reslrt.aw.keep.combine)
tail(reslrt.aw.keep.combine)
###annotation for lrt in FW better than ablove methods for annotation###############
sigseq.reslrt.aw <- c(rownames(reslrt.aw.keep.df))
length(sigseq.reslrt.aw)
head(sigseq.reslrt.aw)
sigseq.reslrt.aw.sig <- sigseq.reslrt.aw[!duplicated(sigseq.reslrt.aw)]
length(sigseq.reslrt.aw.sig)
sigseq.reslrt.aw.df <- data.frame(sigseq.reslrt.aw)
write.table(sigseq.reslrt.aw.df,
            "new2024/results/splitBySalinity/FW/sigseq.reslrt.aw.sig.txt",
            sep="\t", row.names=F, quote=F)
Sig.Annotations.lrtregion.aw <- Annotations[which(Annotations$Transcript %in% sigseq.reslrt.aw), ]
RefSeq <- gsub(" .*", "", Sig.Annotations.lrtregion.aw$STICKLEBACK.UNIPROT)
UniProt <- gsub("tr\\|", "", Sig.Annotations.lrtregion.aw$UNIPROT)
UniProt <- gsub("\\|.*", "", UniProt)

write.table(data.frame(Transcript = Sig.Annotations.lrtregion.aw$Transcript,
                       RefSeq, 
                       UniProt, 
                       PFAM = Sig.Annotations.lrtregion.aw$PFAMID, 
                       PGAMGO = Sig.Annotations.lrtregion.aw$PFAMGO), 
            "new2024/results/splitBySalinity/FW/Sig.annotations.lrtregion.aw.tsv",
            sep="\t", row.names = F, quote = F)

###annotation for lrt in FW better than above methods for annotation###############
sigseq.reslrt2.aw <- c(rownames(reslrt2.aw.keep.df))
length(sigseq.reslrt2.aw)
head(sigseq.reslrt2.aw)
sigseq.reslrt2.aw.df <- data.frame(sigseq.reslrt2.aw)
write.table(sigseq.reslrt2.aw.df,
            "new2024/results/splitBySalinity/FW/sigseq.reslrt2.aw.sig.txt",
            sep="\t", row.names=F, quote=F)
Sig.Annotations.lrtregion2.aw <- Annotations[which(Annotations$Transcript %in% sigseq.reslrt2.aw), ]
RefSeq <- gsub(" .*", "", Sig.Annotations.lrtregion2.aw$STICKLEBACK.UNIPROT)
UniProt <- gsub("tr\\|", "", Sig.Annotations.lrtregion2.aw$UNIPROT)
UniProt <- gsub("\\|.*", "", UniProt)

write.table(data.frame(Transcript = Sig.Annotations.lrtregion2.aw$Transcript,
                       RefSeq, 
                       UniProt, 
                       PFAM = Sig.Annotations.lrtregion2.aw$PFAMID, 
                       PGAMGO = Sig.Annotations.lrtregion2.aw$PFAMGO), 
            "new2024/results/splitBySalinity/FW/Sig.annotations.lrtregion2.aw.tsv",
            sep="\t", row.names = F, quote = F)

#### pair-wise comparison of DE ################################################
##DE analysis for CR vs CL ##########
quantmerge.filterCR.CL <- quantmerge[, which(colnames(quantmerge) %!in% c("CR4SW", "LCR1SW"))]
metadata.coastal <- read.delim("PS_gill_RNA_metafile_coastal.txt", sep="\t")
quantmerge.coastal2 <- quantmerge.filterCR.CL[, which(colnames(quantmerge.filterCR.CL) %in% metadata.coastal$Sample)]
head (quantmerge.filterCR.CL)
head(metadata.coastal)
head(quantmerge.coastal2)
tail(metadata.coastal.filter)
ddscoastal <- DESeqDataSetFromMatrix(countData = quantmerge.coastal2, 
                                colData = metadata.coastal, 
                                design = ~ Region + Treatment + Region:Treatment) #
ddscoastal <- DESeq(ddscoastal)
resultsNames(ddscoastal)

#Extract Normalized count data
ddscoastal.normcounts <- counts(ddscoastal, normalized=T)
normcounts.coastal.df <- data.frame(ddscoastal.normcounts)
normcounts.coastal.df$Sequence <- rownames(normcounts.coastal.df)
write.table(normcounts.coastal.df,
            "new2024/results/coastal/normcounts.coastal.txt",
            sep="\t", row.names=F, quote=F)
sizeFactors(ddscoastal)
#Extract Results
res.cR.cL.inter <- results(ddscoastal, name = c("RegioncoastR.TreatmentSW"))
#Compile Results
res.cR.cL.inter.df <- CompileResults(res.cR.cL.inter, contrast="cR.cL", padj=0.05, l2fc=2)
write.table(res.cR.cL.inter.df,
            "new2024/results/coastal/res.cR.cL.inter.txt",
            sep="\t", row.names=F, quote=F)
res.cR.cL.inter.df$Sequence <- factor(rownames(res.cR.cL.inter.df), levels = rev(rownames(res.cR.cL.inter.df)))

###annotation ###############
sigseq.res.cR.cL.inter <- c(rownames(res.cR.cL.inter.df))
length(sigseq.res.cR.cL.inter)
head(sigseq.res.cR.cL.inter)
sigseq.res.cR.cL.inter.sig <- sigseq.res.cR.cL.inter[!duplicated(sigseq.res.cR.cL.inter)]
length(sigseq.res.cR.cL.inter.sig)
sigseq.res.cR.cL.inter.sig.df <- data.frame(sigseq.res.cR.cL.inter.sig)
write.table(sigseq.res.cR.cL.inter.sig.df,
            "new2024/results/coastal/sigseq.res.cR.cL.inter.sig.txt",
            sep="\t", row.names=F, quote=F)
Sig.Annotations.coastal.inter <- Annotations[which(Annotations$Transcript %in% sigseq.res.cR.cL.inter.sig), ]
RefSeq <- gsub(" .*", "", Sig.Annotations.coastal.inter$STICKLEBACK.UNIPROT)
UniProt <- gsub("tr\\|", "", Sig.Annotations.coastal.inter$UNIPROT)
UniProt <- gsub("\\|.*", "", UniProt)

write.table(data.frame(Transcript = Sig.Annotations.coastal.inter$Transcript,
                       RefSeq, 
                       UniProt, 
                       PFAM = Sig.Annotations.coastal.inter$PFAMID, 
                       PGAMGO = Sig.Annotations.coastal.inter$PFAMGO), 
            "new2024/results/coastal/Sig.Annotations.coastal.inter.tsv",
            sep="\t", row.names = F, quote = F)
######## regional effect in coastal (CR and CL)#############
ddscoastalregion <- DESeqDataSetFromMatrix(countData = quantmerge.coastal2, 
                                     colData = metadata.coastal, 
                                     design = ~ Region) #
ddscoastalregion <- DESeq(ddscoastalregion)
resultsNames(ddscoastalregion)
#Extract Normalized count data
ddscoastalregion.normcounts <- counts(ddscoastalregion, normalized=T)
normcounts.coastalregion.df <- data.frame(ddscoastalregion.normcounts)
normcounts.coastalregion.df$Sequence <- rownames(normcounts.coastalregion.df)
write.table(normcounts.coastalregion.df,
            "new2024/results/coastal/normcounts.coastalregion.txt",
            sep="\t", row.names=F, quote=F)
sizeFactors(ddscoastal)
#Extract Results
res.cR.cL.region <- results(ddscoastalregion, name = c("Region_coastR_vs_coastL"))
#Compile Results
res.cR.cL.region.df <- CompileResults(res.cR.cL.region, contrast="cR.cL", padj=0.05, l2fc=2)
write.table(res.cR.cL.region.df,
            "new2024/results/coastal/res.cR.cL.region.txt",
            sep="\t", row.names=F, quote=F)
nrow(res.cR.cL.region.df)
res.cR.cL.region.df$Sequence <- factor(rownames(res.cR.cL.region.df), levels = rev(rownames(res.cR.cL.region.df)))
###annotation ###############
sigseq.res.cR.cL.region <- c(rownames(res.cR.cL.region.df))
length(sigseq.res.cR.cL.region)
head(sigseq.res.cR.cL.region)
sigseq.res.cR.cL.region.sig.df <- data.frame(sigseq.res.cR.cL.region)
write.table(sigseq.res.cR.cL.region.sig.df,
            "new2024/results/coastal/sigseq.res.cR.cL.region.sig.txt",
            sep="\t", row.names=F, quote=F)
Sig.Annotations.coastal.region <- Annotations[which(Annotations$Transcript %in% sigseq.res.cR.cL.region), ]
RefSeq <- gsub(" .*", "", Sig.Annotations.coastal.region$STICKLEBACK.UNIPROT)
UniProt <- gsub("tr\\|", "", Sig.Annotations.coastal.region$UNIPROT)
UniProt <- gsub("\\|.*", "", UniProt)

write.table(data.frame(Transcript = Sig.Annotations.coastal.region$Transcript,
                       RefSeq, 
                       UniProt, 
                       PFAM = Sig.Annotations.coastal.region$PFAMID, 
                       PGAMGO = Sig.Annotations.coastal.region$PFAMGO), 
            "new2024/results/coastal/Sig.Annotations.coastal.region.tsv",
            sep="\t", row.names = F, quote = F)
### salinity effect on coastal region########
ddscoastalsalinity <- DESeqDataSetFromMatrix(countData = quantmerge.coastal2, 
                                           colData = metadata.coastal, 
                                           design = ~ Treatment) #
ddscoastalsalinity <- DESeq(ddscoastalsalinity)
resultsNames(ddscoastalsalinity)
#Extract Normalized count data
ddscoastalsalinity.normcounts <- counts(ddscoastalsalinity, normalized=T)
normcounts.coastalsalinity.df <- data.frame(ddscoastalsalinity.normcounts)
normcounts.coastalsalinity.df$Sequence <- rownames(normcounts.coastalsalinity.df)
write.table(normcounts.coastalsalinity.df,
            "new2024/results/coastal/normcounts.coastalsalinity.txt",
            sep="\t", row.names=F, quote=F)
sizeFactors(ddscoastalsalinity)
#Extract Results
res.cR.cL.salinity <- results(ddscoastalsalinity, name = c("Treatment_SW_vs_AW"))
#Compile Results
res.cR.cL.salinity.df <- CompileResults(res.cR.cL.salinity, contrast="cR.cL", padj=0.05, l2fc=2)
write.table(res.cR.cL.salinity.df,
            "new2024/results/coastal/res.cR.cL.salinity.txt",
            sep="\t", row.names=F, quote=F)
nrow(res.cR.cL.salinity.df)
res.cR.cL.salinity.df$Sequence <- factor(rownames(res.cR.cL.salinity.df), levels = rev(rownames(res.cR.cL.salinity.df)))
###annotation ###############
sigseq.res.cR.cL.salinity <- c(rownames(res.cR.cL.salinity.df))
length(sigseq.res.cR.cL.salinity)
head(sigseq.res.cR.cL.salinity)
sigseq.res.cR.cL.salinity.sig.df <- data.frame(sigseq.res.cR.cL.salinity)
write.table(sigseq.res.cR.cL.salinity.sig.df,
            "new2024/results/coastal/sigseq.res.cR.cL.salinity.sig.txt",
            sep="\t", row.names=F, quote=F)
Sig.Annotations.coastal.salinity <- Annotations[which(Annotations$Transcript %in% sigseq.res.cR.cL.salinity), ]
RefSeq <- gsub(" .*", "", Sig.Annotations.coastal.salinity$STICKLEBACK.UNIPROT)
UniProt <- gsub("tr\\|", "", Sig.Annotations.coastal.salinity$UNIPROT)
UniProt <- gsub("\\|.*", "", UniProt)

write.table(data.frame(Transcript = Sig.Annotations.coastal.salinity$Transcript,
                       RefSeq, 
                       UniProt, 
                       PFAM = Sig.Annotations.coastal.salinity$PFAMID, 
                       PGAMGO = Sig.Annotations.coastal.salinity$PFAMGO), 
            "new2024/results/coastal/Sig.Annotations.coastal.salinity.tsv",
            sep="\t", row.names = F, quote = F)
### Venn diagram for CR and CL########
library(ggplot2)
library(ggvenn)
a <- list('Treatment' = c(res.cR.cL.salinity.df$Sequence),
          'Habitat' = c(res.cR.cL.region.df$Sequence),
          'Interaction' = c(res.cR.cL.inter.df$Sequence))
ggvenn(a, c("Treatment", "Habitat", "Interaction"),
       fill_color = c("green","purple", "red"),
       fill_alpha = 0.5,
       stroke_color = "black",
       stroke_alpha = 1,
       stroke_size = 1.5,
       set_name_size = 6,
       text_size = 4)
ggvenn(a)

dev.off()
#### DE analysis for lake populations (CL and IL)#########################
metadata.lake <- read.delim("PS_gill_RNA_metafile_lake.txt", sep="\t")
quantmerge.lake <- quantmerge[, which(colnames(quantmerge) %in% metadata.lake$Sample)]
head (quantmerge.lake)
head(metadata.lake)
ddslake <- DESeqDataSetFromMatrix(countData = quantmerge.lake, 
                                     colData = metadata.lake, 
                                     design = ~ Region + Treatment + Region:Treatment) #
ddslake <- DESeq(ddslake)
resultsNames(ddslake)

#Extract Normalized count data
ddslake.normcounts <- counts(ddslake, normalized=T)
normcounts.lake.df <- data.frame(ddslake.normcounts)
normcounts.lake.df$Sequence <- rownames(normcounts.lake.df)
write.table(normcounts.lake.df,
            "new2024/results/lake/normcounts.lake.txt",
            sep="\t", row.names=F, quote=F)
sizeFactors(ddslake)
#Extract Results
res.cL.iL.inter <- results(ddslake, name = c("RegioninlandL.TreatmentSW"))
#Compile Results
res.cL.iL.inter.df <- CompileResults(res.cL.iL.inter, contrast="cR.iL", padj=0.05, l2fc=2)
write.table(res.cL.iL.inter.df,
            "new2024/results/lake/res.cL.iL.inter.txt",
            sep="\t", row.names=F, quote=F)
res.cL.iL.inter.df$Sequence <- factor(rownames(res.cL.iL.inter.df), levels = rev(rownames(res.cL.iL.inter.df)))
###annotation ###############
sigseq.res.cL.iL.inter <- c(rownames(res.cL.iL.inter.df))
length(sigseq.res.cL.iL.inter)
head(sigseq.res.cL.iL.inter)
sigseq.res.cL.iL.inter.sig.df <- data.frame(sigseq.res.cL.iL.inter)
write.table(sigseq.res.cL.iL.inter.sig.df,
            "new2024/results/lake/sigseq.res.cL.iL.inter.sig.txt",
            sep="\t", row.names=F, quote=F)
Sig.Annotations.lake.inter <- Annotations[which(Annotations$Transcript %in% sigseq.res.cL.iL.inter), ]
RefSeq <- gsub(" .*", "", Sig.Annotations.lake.inter$STICKLEBACK.UNIPROT)
UniProt <- gsub("tr\\|", "", Sig.Annotations.lake.inter$UNIPROT)
UniProt <- gsub("\\|.*", "", UniProt)

write.table(data.frame(Transcript = Sig.Annotations.lake.inter$Transcript,
                       RefSeq, 
                       UniProt, 
                       PFAM = Sig.Annotations.lake.inter$PFAMID, 
                       PGAMGO = Sig.Annotations.lake.inter$PFAMGO), 
            "new2024/results/lake/Sig.Annotations.lake.inter.tsv",
            sep="\t", row.names = F, quote = F)
######## regional effect in lake#############
ddslakeregion <- DESeqDataSetFromMatrix(countData = quantmerge.lake, 
                                           colData = metadata.lake, 
                                           design = ~ Region) #
ddslakeregion <- DESeq(ddslakeregion)
resultsNames(ddslakeregion)
#Extract Normalized count data
ddslakeregion.normcounts <- counts(ddslakeregion, normalized=T)
normcounts.lakeregion.df <- data.frame(ddslakeregion.normcounts)
normcounts.lakeregion.df$Sequence <- rownames(normcounts.lakeregion.df)
write.table(normcounts.lakeregion.df,
            "new2024/results/lake/normcounts.lakeregion.txt",
            sep="\t", row.names=F, quote=F)
sizeFactors(ddslake)
#Extract Results
res.cL.iL.region <- results(ddslakeregion, name = c("Region_inlandL_vs_coastL"))
#Compile Results
res.cL.iL.lake.df <- CompileResults(res.cL.iL.region, contrast="cL.iL", padj=0.05, l2fc=2)
write.table(res.cL.iL.lake.df,
            "new2024/results/lake/res.cL.iL.region.txt",
            sep="\t", row.names=F, quote=F)
nrow(res.cL.iL.lake.df)
res.cL.iL.lake.df$Sequence <- factor(rownames(res.cL.iL.lake.df), levels = rev(rownames(res.cL.iL.lake.df)))
###annotation ###############
sigseq.res.cL.iL.region <- c(rownames(res.cL.iL.lake.df))
length(sigseq.res.cL.iL.region)
head(sigseq.res.cL.iL.region)
sigseq.res.cL.iL.region.sig.df <- data.frame(sigseq.res.cL.iL.region)
write.table(sigseq.res.cL.iL.region.sig.df,
            "new2024/results/lake/sigseq.res.cL.iL.region.sig.txt",
            sep="\t", row.names=F, quote=F)
Sig.Annotations.lake.region <- Annotations[which(Annotations$Transcript %in% sigseq.res.cL.iL.region), ]
RefSeq <- gsub(" .*", "", Sig.Annotations.lake.region$STICKLEBACK.UNIPROT)
UniProt <- gsub("tr\\|", "", Sig.Annotations.lake.region$UNIPROT)
UniProt <- gsub("\\|.*", "", UniProt)

write.table(data.frame(Transcript = Sig.Annotations.lake.region$Transcript,
                       RefSeq, 
                       UniProt, 
                       PFAM = Sig.Annotations.lake.region$PFAMID, 
                       PGAMGO = Sig.Annotations.lake.region$PFAMGO), 
            "new2024/results/lake/Sig.Annotations.lake.region.tsv",
            sep="\t", row.names = F, quote = F)
### salinity effect on lake region########
ddslakesalinity <- DESeqDataSetFromMatrix(countData = quantmerge.lake, 
                                             colData = metadata.lake, 
                                             design = ~ Treatment) #
ddslakesalinity <- DESeq(ddslakesalinity)
resultsNames(ddslakesalinity)
#Extract Normalized count data
ddslakesalinity.normcounts <- counts(ddslakesalinity, normalized=T)
normcounts.lakesalinity.df <- data.frame(ddslakesalinity.normcounts)
normcounts.lakesalinity.df$Sequence <- rownames(normcounts.lakesalinity.df)
write.table(normcounts.lakesalinity.df,
            "new2024/results/lake/normcounts.lakesalinity.txt",
            sep="\t", row.names=F, quote=F)
sizeFactors(ddscoastalsalinity)
#Extract Results
res.cL.iL.salinity <- results(ddslakesalinity, name = c("Treatment_SW_vs_AW"))
#Compile Results
res.cL.iL.salinity.df <- CompileResults(res.cL.iL.salinity, contrast="cL.iL", padj=0.05, l2fc=2)
write.table(res.cL.iL.salinity.df,
            "new2024/results/lake/res.cL.iL.salinity.txt",
            sep="\t", row.names=F, quote=F)
nrow(res.cL.iL.salinity.df)
res.cL.iL.salinity.df$Sequence <- factor(rownames(res.cL.iL.salinity.df), levels = rev(rownames(res.cL.iL.salinity.df)))
###annotation ###############
sigseq.res.cL.iL.salinity <- c(rownames(res.cL.iL.salinity.df))
length(sigseq.res.cL.iL.salinity)
head(sigseq.res.cL.iL.salinity)
sigseq.res.cL.iL.salinity.sig.df <- data.frame(sigseq.res.cL.iL.salinity)
write.table(sigseq.res.cL.iL.salinity.sig.df,
            "new2024/results/lake/sigseq.res.cL.iL.salinity.sig.txt",
            sep="\t", row.names=F, quote=F)
Sig.Annotations.lake.salinity <- Annotations[which(Annotations$Transcript %in% sigseq.res.cL.iL.salinity), ]
RefSeq <- gsub(" .*", "", Sig.Annotations.lake.salinity$STICKLEBACK.UNIPROT)
UniProt <- gsub("tr\\|", "", Sig.Annotations.lake.salinity$UNIPROT)
UniProt <- gsub("\\|.*", "", UniProt)

write.table(data.frame(Transcript = Sig.Annotations.lake.salinity$Transcript,
                       RefSeq, 
                       UniProt, 
                       PFAM = Sig.Annotations.lake.salinity$PFAMID, 
                       PGAMGO = Sig.Annotations.lake.salinity$PFAMGO), 
            "new2024/results/lake/Sig.Annotations.lake.salinity.tsv",
            sep="\t", row.names = F, quote = F)
#######Venn diagram for CL and IL#############
library(ggplot2)
library(ggvenn)

a <- list('Treatment' = c(res.cL.iL.salinity.df$Sequence),
          'Habitat' = c(res.cL.iL.lake.df$Sequence),
          'Interaction' = c(res.cL.iL.inter.df$Sequence))
ggvenn(a, c("Treatment", "Habitat", "Interaction"),
       fill_color = c("green","yellow", "blue"),
       fill_alpha = 0.5,
       stroke_color = "black",
       stroke_alpha = 1,
       stroke_size = 1.5,
       set_name_size = 6,
       text_size = 4)
ggvenn(a)

dev.off()
#### DE analysis for CR VS IL populations#########################
quantmerge.filterCR.IL <- quantmerge[, which(colnames(quantmerge) %!in% c("CR4SW", "LCR1SW"))]
metadata.CRIL <- read.delim("PS_gill_RNA_metafile_CR_IL.txt", sep="\t")
metadata.CRIL.filter <- metadata.CRIL[which(metadata.CRIL$Sample %!in% c("CR4SW", "LCR1SW")), ]
quantmerge.CRIL <- quantmerge.filterCR.IL[, which(colnames(quantmerge.filterCR.IL) %in% metadata.CRIL$Sample)]
head (quantmerge.CRIL)
head(metadata.CRIL)
ddsCRIL <- DESeqDataSetFromMatrix(countData = quantmerge.CRIL, 
                                  colData = metadata.CRIL.filter, 
                                  design = ~ Region + Treatment + Region:Treatment) #
ddsCRIL <- DESeq(ddsCRIL)
resultsNames(ddsCRIL)

#Extract Normalized count data
ddsCRIL.normcounts <- counts(ddsCRIL, normalized=T)
normcounts.CRIL.df <- data.frame(ddsCRIL.normcounts)
normcounts.CRIL.df$Sequence <- rownames(normcounts.CRIL.df)
write.table(normcounts.CRIL.df,
            "new2024/results/CRVSIL/normcounts.CRIL.txt",
            sep="\t", row.names=F, quote=F)
sizeFactors(ddsCRIL)
#Extract Results
res.cR.iL.inter <- results(ddsCRIL, name = c("RegioninlandL.TreatmentSW"))
#Compile Results
res.cR.iL.inter.df <- CompileResults(res.cR.iL.inter, contrast="cR.iL", padj=0.05, l2fc=2)
write.table(res.cR.iL.inter.df,
            "new2024/results/CRVSIL/res.cR.iL.inter.txt",
            sep="\t", row.names=F, quote=F)

res.cR.iL.inter.df$Sequence <- factor(rownames(res.cR.iL.inter.df), levels = rev(rownames(res.cR.iL.inter.df)))
###annotation ###############
sigseq.res.cR.iL.inter <- c(rownames(res.cR.iL.inter.df))
length(sigseq.res.cR.iL.inter)
head(sigseq.res.cR.iL.inter)
sigseq.res.cR.iL.inter.sig.df <- data.frame(sigseq.res.cR.iL.inter)
write.table(sigseq.res.cR.iL.inter.sig.df,
            "new2024/results/CRVSIL/sigseq.res.cR.iL.inter.sig.txt",
            sep="\t", row.names=F, quote=F)
Sig.Annotations.CRIL.inter <- Annotations[which(Annotations$Transcript %in% sigseq.res.cR.iL.inter), ]
RefSeq <- gsub(" .*", "", Sig.Annotations.CRIL.inter$STICKLEBACK.UNIPROT)
UniProt <- gsub("tr\\|", "", Sig.Annotations.CRIL.inter$UNIPROT)
UniProt <- gsub("\\|.*", "", UniProt)

write.table(data.frame(Transcript = Sig.Annotations.CRIL.inter$Transcript,
                       RefSeq, 
                       UniProt, 
                       PFAM = Sig.Annotations.CRIL.inter$PFAMID, 
                       PGAMGO = Sig.Annotations.CRIL.inter$PFAMGO), 
            "new2024/results/CRVSIL/Sig.Annotations.CRIL.inter.tsv",
            sep="\t", row.names = F, quote = F)
######## regional effect in CR VS IL#############
ddsCRILregion <- DESeqDataSetFromMatrix(countData = quantmerge.CRIL, 
                                        colData = metadata.CRIL.filter, 
                                        design = ~ Region) #
ddsCRILregion <- DESeq(ddsCRILregion)
resultsNames(ddsCRILregion)
#Extract Normalized count data
ddsCRILregion.normcounts <- counts(ddsCRILregion, normalized=T)
normcounts.CRILregion.df <- data.frame(ddsCRILregion.normcounts)
normcounts.CRILregion.df$Sequence <- rownames(normcounts.CRILregion.df)
write.table(normcounts.CRILregion.df,
            "new2024/results/CRVSIL/normcounts.CRILregion.txt",
            sep="\t", row.names=F, quote=F)
sizeFactors(ddsCRILregion)
#Extract Results
res.cR.iL.region <- results(ddsCRILregion, name = c("Region_inlandL_vs_coastR"))
#Compile Results
res.cR.iL.region.df <- CompileResults(res.cR.iL.region, contrast="cR.iL", padj=0.05, l2fc=2)
write.table(res.cR.iL.region.df,
            "new2024/results/CRVSIL/res.cR.iL.region.txt",
            sep="\t", row.names=F, quote=F)
nrow(res.cR.iL.region.df)
res.cR.iL.region.df$Sequence <- factor(rownames(res.cR.iL.region.df), levels = rev(rownames(res.cR.iL.region.df)))
###annotation ###############
sigseq.res.cR.iL.region <- c(rownames(res.cR.iL.region.df))
length(sigseq.res.cR.iL.region)
head(sigseq.res.cR.iL.region)
sigseq.res.cR.iL.region.sig.df <- data.frame(sigseq.res.cR.iL.region)
write.table(sigseq.res.cR.iL.region.sig.df,
            "new2024/results/CRVSIL/sigseq.res.cR.iL.region.sig.txt",
            sep="\t", row.names=F, quote=F)
Sig.Annotations.CRIL.region <- Annotations[which(Annotations$Transcript %in% sigseq.res.cR.iL.region), ]
RefSeq <- gsub(" .*", "", Sig.Annotations.CRIL.region$STICKLEBACK.UNIPROT)
UniProt <- gsub("tr\\|", "", Sig.Annotations.CRIL.region$UNIPROT)
UniProt <- gsub("\\|.*", "", UniProt)

write.table(data.frame(Transcript = Sig.Annotations.CRIL.region$Transcript,
                       RefSeq, 
                       UniProt, 
                       PFAM = Sig.Annotations.CRIL.region$PFAMID, 
                       PGAMGO = Sig.Annotations.CRIL.region$PFAMGO), 
            "new2024/results/CRVSIL/Sig.Annotations.CRIL.region.tsv",
            sep="\t", row.names = F, quote = F)
### salinity effect on CR and IL########
ddsCRILsalinity <- DESeqDataSetFromMatrix(countData = quantmerge.CRIL, 
                                          colData = metadata.CRIL.filter, 
                                          design = ~ Treatment) #
ddsCRILsalinity <- DESeq(ddsCRILsalinity)
resultsNames(ddsCRILsalinity)
#Extract Normalized count data
ddsCRILsalinity.normcounts <- counts(ddsCRILsalinity, normalized=T)
normcounts.CRILsalinity.df <- data.frame(ddsCRILsalinity.normcounts)
normcounts.CRILsalinity.df$Sequence <- rownames(normcounts.CRILsalinity.df)
write.table(normcounts.CRILsalinity.df,
            "new2024/results/CRVSIL/normcounts.CRILsalinity.txt",
            sep="\t", row.names=F, quote=F)
sizeFactors(ddsCRILalsalinity)
#Extract Results
res.cR.iL.salinity <- results(ddsCRILsalinity, name = c("Treatment_SW_vs_AW"))
#Compile Results
res.cR.iL.salinity.df <- CompileResults(res.cR.iL.salinity, contrast="cR.iL", padj=0.05, l2fc=2)
write.table(res.cR.iL.salinity.df,
            "new2024/results/CRVSIL/res.cR.iL.salinity.txt",
            sep="\t", row.names=F, quote=F)
nrow(res.cR.iL.salinity.df)
res.cR.iL.salinity.df$Sequence <- factor(rownames(res.cR.iL.salinity.df), levels = rev(rownames(res.cR.iL.salinity.df)))
###annotation ###############
sigseq.res.cR.iL.salinity <- c(rownames(res.cR.iL.salinity.df))
length(sigseq.res.cR.iL.salinity)
head(sigseq.res.cR.iL.salinity)
sigseq.res.cR.iL.salinity.sig.df <- data.frame(sigseq.res.cR.iL.salinity)
write.table(sigseq.res.cR.iL.salinity.sig.df,
            "new2024/results/CRVSIL/sigseq.res.cR.iL.salinity.sig.txt",
            sep="\t", row.names=F, quote=F)
Sig.Annotations.CRIL.salinity <- Annotations[which(Annotations$Transcript %in% sigseq.res.cR.iL.salinity), ]
RefSeq <- gsub(" .*", "", Sig.Annotations.CRIL.salinity$STICKLEBACK.UNIPROT)
UniProt <- gsub("tr\\|", "", Sig.Annotations.CRIL.salinity$UNIPROT)
UniProt <- gsub("\\|.*", "", UniProt)

write.table(data.frame(Transcript = Sig.Annotations.CRIL.salinity$Transcript,
                       RefSeq, 
                       UniProt, 
                       PFAM = Sig.Annotations.CRIL.salinity$PFAMID, 
                       PGAMGO = Sig.Annotations.CRIL.salinity$PFAMGO), 
            "new2024/results/CRVSIL/Sig.Annotations.CRIL.salinity.tsv",
            sep="\t", row.names = F, quote = F)
#######Venn diagram for CR and IL#############
library(ggplot2)
library(ggvenn)

a <- list('Treatment' = c(res.cR.iL.salinity.df$Sequence),
          'Habitat' = c(res.cR.iL.region.df$Sequence),
          'Interaction' = c(res.cR.iL.inter.df$Sequence))
ggvenn(a, c("Treatment", "Habitat", "Interaction"),
       fill_color = c("green","red", "blue"),
       fill_alpha = 0.5,
       stroke_color = "black",
       stroke_alpha = 1,
       stroke_size = 1.5,
       set_name_size = 6,
       text_size = 4)
ggvenn(a)

dev.off()

#######################################################################################
#Heatmap of Significant Results for CR and IL interaction effect ######
#######################################################################################
ddsCRIL.normcounts <- read.delim("new2024/results/CRVSIL/normcounts.CRIL.txt", sep="\t")
res.cR.iL.inter.df <- read.delim("new2024/results/CRVSIL/res.cR.iL.inter.txt", sep="\t")
metadata.CRIL <- read.delim("PS_gill_RNA_metafile_CR_IL.txt", sep="\t")
metadata.CRIL.filter <- metadata.CRIL[which(metadata.CRIL$Sample %!in% c("CR4SW", "LCR1SW")), ]
head(normcounts.CRIL.df)
head(res.cR.iL.inter.df)
head(metadata.CRIL.filter)
ddsCRIL.normcounts <- ddsCRIL.normcounts [which(rownames(ddsCRIL.normcounts) %in% res.cR.iL.inter.df$Sequence), ]

colnames(ddsCRIL.normcounts) == metadata.CRIL.filter$Sample
dim(ddsCRIL.normcounts)

#Remove genes with more than half zero
ddsCRIL.normcounts <- ddsCRIL.normcounts[apply(ddsCRIL.normcounts == 0, 1, sum) <= 39, ]
#dim(dds.normcounts.sig)

ddsCRIL.normcounts2 <- log2(ddsCRIL.normcounts)
ddsCRIL.normcounts2[is.infinite(ddsCRIL.normcounts2)] <- 0
ddsCRIL.normcounts2[is.na(ddsCRIL.normcounts2)] <- 0

write.table(cbind(Sequence = rownames(ddsCRIL.normcounts2), ddsCRIL.normcounts2),
            "new2024/results/CRVSIL/ddsCRIL.normcounts2.txt",
            row.names = F,
            quote = F,
            sep="\t")

ddsCRIL.normcounts2 <- read.table("new2024/results/CRVSIL/ddsCRIL.normcounts2.txt", sep="\t", header = T)
rownames(ddsCRIL.normcounts2) <- ddsCRIL.normcounts2$Sequence; ddsCRIL.normcounts2$Sequence <- NULL

head(ddsCRIL.normcounts2)
my_palette <- colorRampPalette(c("#18043d","#f7f7f7", "#e66101"))(n = 299)
png("/Users/Sherry/Desktop/heatmapt.CRIL.inter.png",
    width = 5*500,
    height = 5*500,
    res = 500,
    pointsize = 8)
#png("/Users/Sherry/Desktop/heatmapt.common.region.png")
metadata.CRIL.filter$RegionTreatment <- factor(paste(metadata.CRIL.filter$Region, metadata.CRIL.filter$Treatment))
#head(metadata.filter)
n <- nlevels(metadata.CRIL.filter$RegionTreatment)
dat.col <- data.frame(RegionTreatment=unique(metadata.CRIL.filter$RegionTreatment),
                      RegionColors=c("#8DD3C7", "#30b09a",
                                     "#FFFFB3", "#b5b538"))
                                   
#SampleCol <- as.character(merge(metadata.filter,dat.col)$RegionColors)
SampleCol <- as.character(dat.col[match(metadata.CRIL.filter$RegionTreatment, dat.col$RegionTreatment), ]$RegionColors)
gplots::heatmap.2(as.matrix(ddsCRIL.normcounts2),
                  main = "Abundace of DE Transcripts", # heat map title
                  notecol="black",      # change font color of cell labels to black
                  density.info="none",  # turns off density plot inside color legend
                  trace="none",  # turns off trace lines inside the heat map
                  scale="column",
                  dendrogram="none",
                  Colv=FALSE, 
                  col=my_palette,
                  cexCol = 0.6, 
                  labRow = F,
                  lhei =c(1,5),
                  ColSideColors = SampleCol,
                  distfun = function(x) dist(x,method = 'euclidian'),
                  margins =c(4,5),
)


par(lend = 2)           # square line ends for the color legend

legend(x=-0.05, y=0.98,      # location of the legend on the heatmap plot
       legend = as.character(dat.col$RegionTreatment), # category labels
       col = as.character(dat.col$RegionColors),  # color key
       lty= 1,             # line style
       lwd = 8,            # line width
       cex=0.8,
       bty="n"
)
dev.off()

#######################################################################################
#Heatmap of Significant Results for CL and IL interaction effect ######
#######################################################################################
ddslake.normcounts <- ddslake.normcounts [which(rownames(ddslake.normcounts) %in% res.cL.iL.inter.df$Sequence), ]
colnames(ddslake.normcounts) == metadata.lake$Sample
dim(ddslake.normcounts)

#Remove genes with more than half zero
ddslake.normcounts <- ddslake.normcounts[apply(ddslake.normcounts == 0, 1, sum) <= 25, ]
#dim(dds.normcounts.sig)

ddslake.normcounts2 <- log2(ddslake.normcounts)
ddslake.normcounts2[is.infinite(ddslake.normcounts2)] <- 0
ddslake.normcounts2[is.na(ddslake.normcounts2)] <- 0
write.table(cbind(Sequence = rownames(ddslake.normcounts2), ddslake.normcounts2),
            "new2024/results/lake/ddslake.normcounts2.txt",
            row.names = F,
            quote = F,
            sep="\t")

ddslake.normcounts2 <- read.table("new2024/results/lake/ddslake.normcounts2.txt", sep="\t", header = T)
rownames(ddslake.normcounts2) <- ddslake.normcounts2$Sequence; ddslake.normcounts2$Sequence <- NULL
head(ddslake.normcounts2)
nrow(ddslake.normcounts2)
my_palette <- colorRampPalette(c("#18043d","#f7f7f7", "#e66101"))(n = 299)
png("/Users/Sherry/Desktop/heatmapt.lake.inter.png",
    width = 5*500,
    height = 5*500,
    res = 500,
    pointsize = 8)
#png("/Users/Sherry/Desktop/heatmapt.common.region.png")
metadata.lake$RegionTreatment <- factor(paste(metadata.lake$Region, metadata.lake$Treatment))
#head(metadata.filter)
n <- nlevels(metadata.lake$RegionTreatment)
dat.col <- data.frame(RegionTreatment=unique(metadata.lake$RegionTreatment),
                      RegionColors=c("#8DD3C7", "#30b09a",
                                     "#FFFFB3", "#b5b538"))

#SampleCol <- as.character(merge(metadata.filter,dat.col)$RegionColors)
SampleCol <- as.character(dat.col[match(metadata.lake$RegionTreatment, dat.col$RegionTreatment), ]$RegionColors)
gplots::heatmap.2(as.matrix(ddslake.normcounts2),
                  main = "Abundace of DE Transcripts", # heat map title
                  notecol="black",      # change font color of cell labels to black
                  density.info="none",  # turns off density plot inside color legend
                  trace="none",  # turns off trace lines inside the heat map
                  scale="column",
                  dendrogram="none",
                  Colv=FALSE, 
                  col=my_palette,
                  cexCol = 0.6, 
                  labRow = F,
                  lhei =c(1,5),
                  ColSideColors = SampleCol,
                  distfun = function(x) dist(x,method = 'euclidian'),
                  margins =c(4,5),
)


par(lend = 2)           # square line ends for the color legend

legend(x=-0.05, y=0.98,      # location of the legend on the heatmap plot
       legend = as.character(dat.col$RegionTreatment), # category labels
       col = as.character(dat.col$RegionColors),  # color key
       lty= 1,             # line style
       lwd = 8,            # line width
       cex=0.8,
       bty="n"
)

dev.off()

#######################################################################################
#Heatmap of Significant Results for CR and CL interaction effect ######
#######################################################################################
ddscoastal.normcounts <- ddscoastal.normcounts [which(rownames(ddscoastal.normcounts) %in% res.cR.cL.inter.df$Sequence), ]
colnames(ddscoastal.normcounts) == metadata.coastal$Sample
dim(ddscoastal.normcounts)

#Remove genes with more than half zero
ddscoastal.normcounts <- ddscoastal.normcounts[apply(ddscoastal.normcounts == 0, 1, sum) <= 27, ]
#dim(dds.normcounts.sig)

ddscoastal.normcounts2 <- log2(ddscoastal.normcounts)
ddscoastal.normcounts2[is.infinite(ddscoastal.normcounts2)] <- 0
ddscoastal.normcounts2[is.na(ddscoastal.normcounts2)] <- 0
write.table(cbind(Sequence = rownames(ddscoastal.normcounts2), ddscoastal.normcounts2),
            "new2024/results/coastal/ddscoastal.normcounts2.txt",
            row.names = F,
            quote = F,
            sep="\t")

ddscoastal.normcounts2 <- read.table("new2024/results/coastal/ddscoastal.normcounts2.txt", sep="\t", header = T)
rownames(ddscoastal.normcounts2) <- ddscoastal.normcounts2$Sequence; ddscoastal.normcounts2$Sequence <- NULL
head(ddscoastal.normcounts2)
nrow(ddscoastal.normcounts2)
my_palette <- colorRampPalette(c("#18043d","#f7f7f7", "#e66101"))(n = 299)
png("/Users/Sherry/Desktop/heatmapt.coastal.inter.png",
    width = 5*500,
    height = 5*500,
    res = 500,
    pointsize = 8)
#png("/Users/Sherry/Desktop/heatmapt.common.region.png")
metadata.coastal$RegionTreatment <- factor(paste(metadata.coastal$Region, metadata.coastal$Treatment))
#head(metadata.filter)
n <- nlevels(metadata.coastal$RegionTreatment)
dat.col <- data.frame(RegionTreatment=unique(metadata.coastal$RegionTreatment),
                      RegionColors=c("#8DD3C7", "#30b09a",
                                     "#FFFFB3", "#b5b538"))

#SampleCol <- as.character(merge(metadata.filter,dat.col)$RegionColors)
SampleCol <- as.character(dat.col[match(metadata.coastal$RegionTreatment, dat.col$RegionTreatment), ]$RegionColors)
gplots::heatmap.2(as.matrix(ddscoastal.normcounts2),
                  main = "Abundace of DE Transcripts", # heat map title
                  notecol="black",      # change font color of cell labels to black
                  density.info="none",  # turns off density plot inside color legend
                  trace="none",  # turns off trace lines inside the heat map
                  scale="column",
                  dendrogram="none",
                  Colv=FALSE, 
                  col=my_palette,
                  cexCol = 0.6, 
                  labRow = F,
                  lhei =c(1,5),
                  ColSideColors = SampleCol,
                  distfun = function(x) dist(x,method = 'euclidian'),
                  margins =c(4,5),
)


par(lend = 2)           # square line ends for the color legend

legend(x=-0.05, y=0.98,      # location of the legend on the heatmap plot
       legend = as.character(dat.col$RegionTreatment), # category labels
       col = as.character(dat.col$RegionColors),  # color key
       lty= 1,             # line style
       lwd = 8,            # line width
       cex=0.8,
       bty="n"
)


dev.off()

#######################################################################################
#Heatmap of Significant Results for interaction effect across all habitats######
#######################################################################################
dds.normcounts <- dds.normcounts[which(rownames(dds.normcounts) %in% reslrt.keep3.df$Sequence), ]
colnames(dds.normcounts) == metadata.filter$Sample
dim(dds.normcounts)

#Remove genes with more than half zero
dds.normcounts <- dds.normcounts[apply(dds.normcounts == 0, 1, sum) <= 39, ]
#dim(dds.normcounts.sig)

dds.normcounts2 <- log2(dds.normcounts)
dds.normcounts2[is.infinite(dds.normcounts2)] <- 0
dds.normcounts2[is.na(dds.normcounts2)] <- 0
write.table(cbind(Sequence = rownames(dds.normcounts2), dds.normcounts2),
            "new2024/results/interaction/dds.normcounts2.txt",
            row.names = F,
            quote = F,
            sep="\t")

dds.normcounts2 <- read.table("new2024/results/interaction/dds.normcounts2.txt", sep="\t", header = T)
rownames(dds.normcounts2) <- dds.normcounts2$Sequence; dds.normcounts2$Sequence <- NULL
head(dds.normcounts2)
nrow(dds.normcounts2)
my_palette <- colorRampPalette(c("#18043d","#f7f7f7", "#e66101"))(n = 299)
png("/Users/Sherry/Desktop/heatmapt.inter.png",
    width = 5*500,
    height = 5*500,
    res = 500,
    pointsize = 8)
#png("/Users/Sherry/Desktop/heatmapt.common.region.png")
metadata.filter$RegionTreatment <- factor(paste(metadata.filter$Region, metadata.filter$Treatment))
#head(metadata.filter)
n <- nlevels(metadata.filter$RegionTreatment)
dat.col <- data.frame(RegionTreatment=unique(metadata.filter$RegionTreatment),
                      RegionColors=c("#8DD3C7", "#30b09a",
                                     "#FFFFB3", "#b5b538",
                                     "#BEBADA", "#564d91"))

#SampleCol <- as.character(merge(metadata.filter,dat.col)$RegionColors)
SampleCol <- as.character(dat.col[match(metadata.filter$RegionTreatment, dat.col$RegionTreatment), ]$RegionColors)
gplots::heatmap.2(as.matrix(dds.normcounts2),
                  main = "Abundace of DE Transcripts", # heat map title
                  notecol="black",      # change font color of cell labels to black
                  density.info="none",  # turns off density plot inside color legend
                  trace="none",  # turns off trace lines inside the heat map
                  scale="column",
                  dendrogram="none",
                  Colv=FALSE, 
                  col=my_palette,
                  cexCol = 0.6, 
                  labRow = F,
                  lhei =c(1,5),
                  ColSideColors = SampleCol,
                  distfun = function(x) dist(x,method = 'euclidian'),
                  margins =c(4,5),
)


par(lend = 2)           # square line ends for the color legend

legend(x=-0.05, y=0.98,      # location of the legend on the heatmap plot
       legend = as.character(dat.col$RegionTreatment), # category labels
       col = as.character(dat.col$RegionColors),  # color key
       lty= 1,             # line style
       lwd = 8,            # line width
       cex=0.8,
       bty="n"
)


dev.off()

#### salinity comparison of DE within a region/habitat################################################
## Salinity DE analysis for CR  ##########
metadata.CR <- read.delim("PS_gill_RNA_metafile_CR.txt", sep="\t")

quantmerge.filterCR <- quantmerge[, which(colnames(quantmerge) %!in% c("CR4SW", "LCR1SW"))]
metadata.CR <- read.delim("PS_gill_RNA_metafile_CR.txt", sep="\t")
quantmerge.CR2 <- quantmerge.filterCR[, which(colnames(quantmerge.filterCR) %in% metadata.CR$Sample)]
head(quantmerge.CR2)
tail(quantmerge.CR2)
ddsCR <- DESeqDataSetFromMatrix(countData = quantmerge.CR2, 
                                     colData = metadata.CR, 
                                     design = ~ Treatment) #
ddsCR <- DESeq(ddsCR)
resultsNames(ddsCR)

#Extract Normalized count data
ddsCR.normcounts <- counts(ddsCR, normalized=T)
normcounts.CR.df <- data.frame(ddsCR.normcounts)
normcounts.CR.df$Sequence <- rownames(normcounts.CR.df)
write.table(normcounts.CR.df,
            "new2024/results/coastal/normcounts.CR.txt",
            sep="\t", row.names=F, quote=F)
sizeFactors(ddsCR)
#Extract Results
res.cR.salinity <- results(ddsCR, name = c("Treatment_SW_vs_AW"))
#Compile Results
res.cR.salinity.df <- CompileResults(res.cR.salinity, contrast="SW.FW", padj=0.05, l2fc=2)
write.table(res.cR.salinity.df,
            "new2024/results/coastal/res.cR.salinity.txt",
            sep="\t", row.names=F, quote=F)
res.cR.salinity.df$Sequence <- factor(rownames(res.cR.salinity.df), levels = rev(rownames(res.cR.salinity.df)))

## Salinity DE analysis for CL  ##########
metadata.CL <- read.delim("PS_gill_RNA_metafile_CL.txt", sep="\t")
quantmerge.filterCL <- quantmerge[, which(colnames(quantmerge) %!in% c("CR4SW", "LCR1SW"))]
quantmerge.CL2 <- quantmerge.filterCL[, which(colnames(quantmerge.filterCL) %in% metadata.CL$Sample)]
head(quantmerge.CL2)
tail(quantmerge.CL2)
ddsCL <- DESeqDataSetFromMatrix(countData = quantmerge.CL2, 
                                colData = metadata.CL, 
                                design = ~ Treatment) #
ddsCL <- DESeq(ddsCL)
resultsNames(ddsCL)

#Extract Normalized count data
ddsCL.normcounts <- counts(ddsCL, normalized=T)
normcounts.CL.df <- data.frame(ddsCL.normcounts)
normcounts.CL.df$Sequence <- rownames(normcounts.CL.df)
write.table(normcounts.CL.df,
            "new2024/results/coastal/normcounts.CL.txt",
            sep="\t", row.names=F, quote=F)
sizeFactors(ddsCL)
#Extract Results
res.cL.salinity <- results(ddsCL, name = c("Treatment_SW_vs_AW"))
#Compile Results
res.cL.salinity.df <- CompileResults(res.cL.salinity, contrast="SW.FW", padj=0.05, l2fc=2)
write.table(res.cL.salinity.df,
            "new2024/results/coastal/res.cL.salinity.txt",
            sep="\t", row.names=F, quote=F)
res.cL.salinity.df$Sequence <- factor(rownames(res.cL.salinity.df), levels = rev(rownames(res.cL.salinity.df)))
nrow(res.cR.salinity.df)

## Salinity DE analysis for IL  ##########
metadata.IL <- read.delim("PS_gill_RNA_metafile_IL.txt", sep="\t")
quantmerge.filterIL <- quantmerge[, which(colnames(quantmerge) %!in% c("CR4SW", "LCR1SW"))]
quantmerge.IL2 <- quantmerge.filterIL[, which(colnames(quantmerge.filterIL) %in% metadata.IL$Sample)]
head(quantmerge.IL2)
tail(quantmerge.IL2)
ddsIL <- DESeqDataSetFromMatrix(countData = quantmerge.IL2, 
                                colData = metadata.IL, 
                                design = ~ Treatment) #
ddsIL <- DESeq(ddsIL)
resultsNames(ddsIL)

#Extract Normalized count data
ddsIL.normcounts <- counts(ddsIL, normalized=T)
normcounts.IL.df <- data.frame(ddsIL.normcounts)
normcounts.IL.df$Sequence <- rownames(normcounts.IL.df)
write.table(normcounts.IL.df,
            "new2024/results/coastal/normcounts.IL.txt",
            sep="\t", row.names=F, quote=F)
sizeFactors(ddsIL)
#Extract Results
res.iL.salinity <- results(ddsIL, name = c("Treatment_SW_vs_AW"))
#Compile Results
res.iL.salinity.df <- CompileResults(res.iL.salinity, contrast="SW.FW", padj=0.05, l2fc=2)
write.table(res.iL.salinity.df,
            "new2024/results/coastal/res.iL.salinity.txt",
            sep="\t", row.names=F, quote=F)
res.iL.salinity.df$Sequence <- factor(rownames(res.iL.salinity.df), levels = rev(rownames(res.iL.salinity.df)))
nrow(res.iL.salinity.df)

#######################################################################################
#Heatmap of Significant Results for regional effect of all ######
#######################################################################################
dds.normcounts <- read.delim("new2024/normcounts.txt", sep="\t")
sig.region <- read.table("new2024/results/Region/sigseqlrt2.region.sig.txt", sep="\t", header = F)$V1
metadata <- read.delim("PS_gill_RNA_metafile6.txt", sep="\t")
metadata.filter <- metadata[which(metadata$Sample %!in% c("CR4SW", "LCR1SW")), ]
head(dds.normcounts)
head(sig.region)
head(metadata.filter)
tail(sig.region)
ddsregion.normcounts <- dds.normcounts[which(rownames(dds.normcounts) %in% sig.region), ]
head(ddsregion.normcounts)
colnames(ddsregion.normcounts) == metadata.filter$Sample
dim(ddsregion.normcounts)

#Remove genes with more than half zero
#ddsregion.normcounts <- ddsregion.normcounts[apply(ddsregion.normcounts == 0, 1, sum) <= 39, ]
#dim(ddsregion.normcounts)

ddsregion.normcounts2 <- log2(ddsregion.normcounts)
ddsregion.normcounts2[is.infinite(ddsregion.normcounts2)] <- 0
ddsregion.normcounts2[is.na(ddsregion.normcounts2)] <- 0
head (ddsregion.normcounts2)
write.table(cbind(Sequence = rownames(ddsregion.normcounts2), ddsregion.normcounts2),
            "new2024/results/Region/ddsregion.normcounts2.txt",
            row.names = F,
            quote = F,
            sep="\t")

ddsregion.normcounts2 <- read.table("new2024/results/Region/ddsregion.normcounts2.txt", sep="\t", header = T)
rownames(ddsregion.normcounts2) <- ddsregion.normcounts2$Sequence; ddsregion.normcounts2$Sequence <- NULL

head(ddsregion.normcounts2)
my_palette <- colorRampPalette(c("#18043d","#f7f7f7", "#e66101"))(n = 299)
png("/Users/Sherry/Desktop/heatmapt.region.png",
    width = 5*500,
    height = 5*500,
    res = 500,
    pointsize = 8)
#png("/Users/Sherry/Desktop/heatmapt.region.png")
metadata.filter$RegionTreatment <- factor(paste(metadata.filter$Region, metadata.filter$Treatment))
#head(metadata.filter)
n <- nlevels(metadata.filter$RegionTreatment)
#n<- nlevels(metadata.filter$Region)
dat.col <- data.frame(RegionTreatment=unique(metadata.filter$RegionTreatment),
                      RegionColors=c("#8DD3C7", "#30b09a",
                                     "#FFFFB3", "#b5b538",
                                     "#BEBADA", "#564d91"))
#dat.col <- data.frame(Region=unique(metadata.filter$Region),
                      #RegionColors=c("#8DD3C7", "#30b09a",
                          #           "#FFFFB3", "#b5b538",
                          #           "#BEBADA", "#564d91"))
#SampleCol <- as.character(merge(metadata.filter,dat.col)$RegionColors)
#SampleCol <- as.character(merge(metadata.filter,dat.col)$RegionColors)
SampleCol <- as.character(dat.col[match(metadata.filter$RegionTreatment, dat.col$RegionTreatment), ]$RegionColors)
gplots::heatmap.2(as.matrix(ddsregion.normcounts2),
                  main = "Abundace of DE Transcripts", # heat map title
                  notecol="black",      # change font color of cell labels to black
                  density.info="none",  # turns off density plot inside color legend
                  trace="none",  # turns off trace lines inside the heat map
                  scale="column",
                  dendrogram="none",
                  Colv=FALSE, 
                  col=my_palette,
                  cexCol = 0.9, 
                  labRow = F,
                  lhei =c(1,5),
                  ColSideColors = SampleCol,
                  distfun = function(x) dist(x,method = 'euclidian'),
                  margins =c(4,5),
)


par(lend = 2)           # square line ends for the color legend

legend(x=-0.05, y=0.98,      # location of the legend on the heatmap plot
       legend = as.character(dat.col$Region), # category labels
       col = as.character(dat.col$RegionColors),  # color key
       lty= 1,             # line style
       lwd = 8,            # line width
       cex=1.5,
       bty="n"
)
dev.off()

#######################################################################################
#Heatmap of Significant Results for salinity effect of all ######
#######################################################################################
#dds.normcounts <- read.delim("new2024/normcounts.txt", sep="\t")
sig.treat <- read.table("new2024/results/splitBySalinity/sigseqlrt3.salinity.sig.txt", sep="\t")$V1
metadata <- read.delim("PS_gill_RNA_metafile5.txt", sep="\t")
metadata.filter <- metadata[which(metadata$Sample %!in% c("CR4SW", "LCR1SW")), ]
head(dds.normcounts)
head(sig.region)
head(metadata.filter)
tail(sig.region)
ddstreat.normcounts <- dds.normcounts[which(rownames(dds.normcounts) %in% sig.treat), ]
head(ddstreat.normcounts)
colnames(ddstreat.normcounts) == metadata.filter$Sample
dim(ddstreat.normcounts)

#Remove genes with more than half zero
#ddsregion.normcounts <- ddsregion.normcounts[apply(ddsregion.normcounts == 0, 1, sum) <= 39, ]
#dim(ddsregion.normcounts)

ddstreat.normcounts2 <- log2(ddstreat.normcounts)
ddstreat.normcounts2[is.infinite(ddstreat.normcounts2)] <- 0
ddstreat.normcounts2[is.na(ddstreat.normcounts2)] <- 0
head (ddstreat.normcounts2)
write.table(cbind(Sequence = rownames(ddstreat.normcounts2), ddstreat.normcounts2),
            "new2024/results/splitBySalinity/ddstreat.normcounts2.txt",
            row.names = F,
            quote = F,
            sep="\t")


ddstreat.normcounts2 <- read.table("new2024/results/splitBySalinity/ddstreat.normcounts2.txt", sep="\t", header = T)

rownames(ddstreat.normcounts2) <- ddstreat.normcounts2$Sequence; ddstreat.normcounts2$Sequence <- NULL

head(ddstreat.normcounts2)
my_palette <- colorRampPalette(c("#18043d","#f7f7f7", "#e66101"))(n = 299)
png("/Users/Sherry/Desktop/heatmapt.treat.png",
    width = 5*500,
    height = 5*500,
    res = 500,
    pointsize = 8)
#png("/Users/Sherry/Desktop/heatmapt.treat.png")
metadata.filter$RegionTreatment <- factor(paste(metadata.filter$Region, metadata.filter$Treatment))
#head(metadata.filter)
n <- nlevels(metadata.filter$RegionTreatment)
#n<- nlevels(metadata.filter$Region)
dat.col <- data.frame(RegionTreatment=unique(metadata.filter$RegionTreatment),
                      RegionColors=c("#8DD3C7", "#30b09a",
                                     "#FFFFB3", "#b5b538",
                                     "#BEBADA", "#564d91"))
#dat.col <- data.frame(Region=unique(metadata.filter$Region),
#RegionColors=c("#8DD3C7", "#30b09a",
#           "#FFFFB3", "#b5b538",
#           "#BEBADA", "#564d91"))
#SampleCol <- as.character(merge(metadata.filter,dat.col)$RegionColors)
#SampleCol <- as.character(merge(metadata.filter,dat.col)$RegionColors)
SampleCol <- as.character(dat.col[match(metadata.filter$RegionTreatment, dat.col$RegionTreatment), ]$RegionColors)
gplots::heatmap.2(as.matrix(ddstreat.normcounts2),
                  main = "Abundace of DE Transcripts", # heat map title
                  notecol="black",      # change font color of cell labels to black
                  density.info="none",  # turns off density plot inside color legend
                  trace="none",  # turns off trace lines inside the heat map
                  scale="column",
                  dendrogram="none",
                  Colv=FALSE, 
                  col=my_palette,
                  cexCol = 0.9, 
                  labRow = F,
                  lhei =c(1,5),
                  ColSideColors = SampleCol,
                  distfun = function(x) dist(x,method = 'euclidian'),
                  margins =c(4,5),
)


par(lend = 2)           # square line ends for the color legend

legend(x=-0.05, y=0.98,      # location of the legend on the heatmap plot
       legend = as.character(dat.col$Region), # category labels
       col = as.character(dat.col$RegionColors),  # color key
       lty= 1,             # line style
       lwd = 8,            # line width
       cex=1.5,
       bty="n"
)
dev.off()

#### Venn diagram for salinity effect in each habita##########
library(ggplot2)
library(ggvenn)
salinity2 <- list('coastal river' = c(res.cR.salinity.df$Sequence),
          'coastal lake' = c(res.cL.salinity.df$Sequence),
          'interior lake' = c(res.iL.salinity.df$Sequence))
ggvenn(salinity2, c("coastal river", "coastal lake", "interior lake"),
       fill_color = c("purple", "yellow", "green"),
       fill_alpha = 0.5,
       stroke_color = "black",
       stroke_alpha = 1,
       stroke_size = 1.5,
       set_name_size = 6,
       text_size = 4)
ggvenn(salinity2)
res.cR.salinity <- read.delim("new2024/results/coastal/res.cR.salinity.txt", sep="\t")
res.cR.salinity.df <- data.frame(res.cR.salinity)
res.cR.salinity.df$Sequence <- rownames(res.cR.salinity.df)
res.cL.salinity <- read.delim("new2024/results/coastal/res.cL.iL.region.txt", sep="\t")
res.cL.salinity.df <- data.frame(res.cL.salinity)
res.cL.salinity.df$Sequence <- rownames(res.cL.salinity.df)
res.iL.salinity <- read.delim("new2024/results/coastal/res.iL.salinity.txt", sep="\t")
res.iL.salinity.df <- data.frame(res.iL.salinity)
res.iL.salinity.df$Sequence <- rownames(res.iL.salinity.df)

#### Venn diagram for habiat effect in each pair##########
library(ggplot2)
library(ggvenn)
habitat <- list('coastal river_ coastal lake' = c(res.cR.cL.region.df$Sequence),
                  'interior lake_coastal lake' = c(res.cL.iL.region.df$Sequence),
                  'coastal river_interior lake' = c(res.cR.iL.region.df$Sequence))
ggvenn(habitat, c("coastal river_ coastal lake", "interior lake_coastal lake", "coastal river_interior lake"),
       fill_color = c("purple", "green", "grey"),
       fill_alpha = 0.5,
       stroke_color = "black",
       stroke_alpha = 1,
       stroke_size = 1.5,
       set_name_size = 6,
       text_size = 4)
ggvenn(habitat)
res.cR.cL.region <- read.delim("new2024/results/coastal/sigseq.res.cR.cL.region.sig.txt", sep="\t")
res.cR.cL.region.df <- data.frame(res.cR.cL.region)
res.cR.cL.region.df$Sequence <- rownames(res.cR.cL.region.df)
res.cL.iL.region <- read.delim("new2024/results/lake/res.cL.iL.region.txt", sep="\t")
res.cL.iL.region.df <- data.frame(res.cL.iL.region)
res.cL.iL.region.df$Sequence <- rownames(res.cL.iL.region.df)
res.cR.iL.region <- read.delim("new2024/results/CRVSIL/res.cR.iL.region.txt", sep="\t")
res.cR.iL.region.df <- data.frame(res.cR.iL.region)
res.cR.iL.region.df$Sequence <- rownames(res.cR.iL.region.df)

#### Venn diagram for interaction effect in each pair##########
library(ggplot2)
library(ggvenn)
inter <- list('coastal river_ coastal lake' = c(res.cR.cL.inter.df$Sequence),
                'interior lake_coastal lake' = c(res.cL.iL.inter.df$Sequence),
                'coastal river_interior lake' = c(res.cR.iL.inter.df$Sequence))
ggvenn(inter, c("coastal river_ coastal lake", "interior lake_coastal lake", "coastal river_interior lake"),
       fill_color = c("purple", "yellow", "red"),
       fill_alpha = 0.5,
       stroke_color = "black",
       stroke_alpha = 1,
       stroke_size = 1.5,
       set_name_size = 6,
       text_size = 4)
ggvenn(inter)
res.cR.cL.inter <- read.delim("new2024/results/coastal/res.cR.cL.inter.txt", sep="\t")
res.cR.cL.inter.df <- data.frame(res.cR.cL.inter)
res.cR.cL.inter.df$Sequence <- rownames(res.cR.cL.inter.df)
res.cL.iL.inter <- read.delim("new2024/results/lake/res.cL.iL.inter.txt", sep="\t")
res.cL.iL.inter.df <- data.frame(res.cL.iL.inter)
res.cL.iL.inter.df$Sequence <- rownames(res.cL.iL.inter.df)
res.cR.iL.inter <- read.delim("new2024/results/CRVSIL/res.cR.iL.inter.txt", sep="\t")
res.cR.iL.inter.df <- data.frame(res.cR.iL.inter)
res.cR.iL.inter.df$Sequence <- rownames(res.cR.iL.inter.df)
