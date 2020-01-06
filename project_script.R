library(dplyr)
library(DESeq2)
library(pheatmap)
library(ggplot2)

# importing data
counts <- read.table("Input_data/unprocessed_counts.tab", row.names = 1, check.names = F)
pheno <- read.table("Input_data/pheno.tab",row.names = 1, check.names = F)
stat <- read.table("Input_data/stat.tab", row.names = 1, check.names = F)

# clean sample names
colnames(counts) <- gsub(".sorted.bam$", "", colnames(counts)) # removes files extension
colnames(counts) <- gsub("^121317.", "", colnames(counts)) # removes tag from lib samples

## removal of samples
counts <- counts %>% dplyr::select(-contains("unmatched"))
cor.table <- cor(counts, method = "spearman") # spearman cor table

# cor heat map: pre-removal
pheatmap(cor.table, fontsize_row = 7, fontsize_col = 7)
dev.copy(png,'Output/Figures/Fig1-pre-removal.png') # saving options for figure
dev.off()

# sample removals
pheno.dis <- c("L6.TTAGGC", "L3.TTAGGC") # pheno discrepancy: gender and no pheno
low.cor <- names(which(apply(cor.table, 2, mean) < 0.88))
remove <- c(low.cor, pheno.dis) # remove samples with cor > .88
counts <- counts %>% dplyr::select(-one_of(remove)) 

# cor heat map: post-removal
cor.table <- cor(counts, method = "spearman") # spearman cor table
pheatmap(cor.table, fontsize_row = 7, fontsize_col = 7)
dev.copy(png,'Output/Figures/Fig2-post-removal.png') # saving options for figure
dev.off()

#read stat for removed samples
stat <- stat %>% pivot_longer(col= 2:ncol(stat), 
                              names_to = "sample" , 
                              values_to = "count") %>% filter(count > 0)


remove.stat <- stat[which(stat$sample %in% low.cor),]
ggplot(remove.stat, aes(sample, count, fill = Status)) + geom_bar(stat = "identity", position = "fill") 
dev.copy(png,'Output/Figures/Fig3-stat-removal-ratio.png') # saving options for figure
dev.off()
options(scipen=10000) # gets rid of scietific notations
ggplot(remove.stat, aes(sample, count, fill = Status)) + geom_bar(stat = "identity")
dev.copy(png,'Output/Figures/Fig4-stat-removal.png') # saving options for figure
options(scipen=1)

# saving options
write.table(counts, "Output/Data/counts.tab", sep = "\t")

## DESeq2 analysis
counts <- counts[,match(rownames(pheno), colnames(counts))] # ordering counts with pheno
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = pheno,
                              design = ~  pheno) # experiment declaration 

# PCA plot with VST transform 
vsd <- vst(dds)
pcaData <- DESeq2::plotPCA(vsd, intgroup = c("pheno","batch_number"))
# indiv labels
ggplot(pcaData$data)  + 
  geom_text(aes(PC1, PC2, color = batch_number), label = pheno$indiv) + # with replicate labels
  ylab(pcaData$labels$y) +
  xlab(pcaData$labels$x)
# batch_number labels
ggplot(pcaData$data)  + 
  geom_point(aes(PC1, PC2, color = pheno, shape = batch_number)) + # with replicate labels
  ylab(pcaData$labels$y) +
  xlab(pcaData$labels$x)
dev.copy(png,'Output/Figures/Fig5-PCA.png') # saving options for figure
dev.off()

# differential expression 
dds <- DESeq(dds)
plotDispEsts(dds) # dispersion plot
dev.copy(png,'Output/Figures/dispersion.png') # saving options for figure
dev.off()

res <- results(dds, contrast = c("pheno", "C", "S")) #specify group comparisons here
summary(res)

# top results
top <- rownames(res)[which(res$padj < 0.05)] # sig genes
plotCounts(dds, "RPL11", "pheno") # plot normalized counts of specific genes by group

# MA plot
plotMA(res, ylim=c(-4,4)) 

# skrinkage estimator
resultsNames(dds) # select coef for srinkage
resLFC <- lfcShrink(dds, coef = 4, type = "apeglm")
plotMA(resLFC, ylim=c(-1,1))
