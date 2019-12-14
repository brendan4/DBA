library(dplyr)
library(DESeq2)

# importing data
counts <- read.table("Input_data/unprocessed_counts.tab", row.names = 1, check.names = F)
pheno <- read.table("Input_data/pheno.tab",row.names = 1, check.names = F)

# clean sample names
colnames(counts) <- gsub(".sorted.bam$", "", colnames(counts)) # removes files extension
colnames(counts) <- gsub("^121317.", "", colnames(counts)) # removes tag from lib samples

# removal of samples
counts <- counts %>% dplyr::select(-contains("unmatched"))
cor.table <- cor(counts, method = "spearman") # spearman cor table
pheno.dis <- c("L6.TTAGGC", "L3.TTAGGC") # pheno discrepancy
remove <- c(names(which(apply(cor.table, 2, sum) < 
                          ncol(cor.table)*.88)), pheno.dis) # remove samples with cor > .88
counts <- counts %>% dplyr::select(-one_of(remove)) 

# saving options
write.table(counts, "Output/Data/counts.tab", sep = "\t")

## DESeq2 analysis
dds <- DESeqDataSetFromMatrix(countData = counts.t,
                              colData = pheno,
                              design = ~  pheno) # experiment declaration 

# PCA plot with VST transform 
vsd <- vst(dds)
pcaData <- DESeq2::plotPCA(vsd, intgroup = c("pheno"))
ggplot(pcaData$data)  + # with replicate labels
  geom_text(aes(PC1, PC2, color = pheno), label = pheno$Replicates)+
  ylab(pcaData$labels$y)+
  xlab(pcaData$labels$x)

# differential expression 
dds <- DESeq(dds)
plotDispEsts(dds)

res <- results(dds, contrast = c("pheno", "C", "W")) #specify group comparisons here
summary(res)

# top results
top <- rownames(res)[which(res$padj < 0.00005)]
plotCounts(dds, "RPL11", "pheno") # plot normalized counts of specific genes by group

# MA plot
plotMA(res, ylim=c(-4,4)) 

# skrinkage estimator
resultsNames(dds) # select coef for srinkage
resLFC <- lfcShrink(dds, coef = 4, type = "apeglm")
plotMA(resLFC, ylim=c(-1,1))
