library(tidyverse)
library(DESeq2)
library(pheatmap)
library(ggplot2)
library(readr)
library(propr)
source ('./mapCode-colors.r')

# importing data
counts <- read.table("Input_data/unprocessed_counts.tab", row.names = 1, check.names = F)
pheno <- read.table("Input_data/pheno.tab",row.names = 1, check.names = F)
stat <- read.table("Input_data/stat.tab", row.names = 1, check.names = F)
pheno.data <- read_csv("Output/data/pheno_out.txt")

# clean sample names
colnames(counts) <- gsub(".sorted.bam$", "", colnames(counts)) # removes files extension
colnames(counts) <- gsub("^121317.", "", colnames(counts)) # removes tag from lib samples

# matching pheno.data names to count data
idx <- which(gsub(".\\d*$", "",colnames(counts)) %in%
               gsub("^121317.", "", pheno.data$library_ID)) # match earlier hiseq rows
pheno.data[which(pheno.data$lane=="earlier_HiseqRun"), 
           "index_sequence"] <- colnames(counts)[idx] # using counts naming convention for earlier hiseq rows

## Step: removal of samples
counts <- counts %>% dplyr::select(-contains("unmatched")) # remove unmatched cols
cor.table <- cor(counts, method = "spearman") # spearman cor table


# cor heat map: pre-removal; see if we can remove the dendrogram from one side
breaksList = seq(0, 1, by = .01)
h1 = pheatmap(cor.table, fontsize_row = 7, fontsize_col = 7, breaks = breaksList,
              color = color.palette2(length(breaksList)) )
dev.copy(pdf,'Output/Figures/Fig1-pre-removal.pdf') # saving options for figure
dev.off()

# sample removals
pheno.dis <- c("L6.TTAGGC", "L3.TTAGGC") # pheno discrepancy: gender and no pheno
low.cor <- names(which(apply(cor.table, 2, mean) < 0.9)) # samples names with cor > 0.9
remove <- c(low.cor, pheno.dis) 
counts <- counts %>% dplyr::select(-one_of(remove)) # removing samples from counts
pheno.data <- pheno.data[-which(!pheno.data$index_sequence %in% colnames(counts)),] # removing samples from pheno.data

# cor heat map: post-removal
cor.table <- cor(counts, method = "spearman") # spearman cor table
pheatmap(cor.table, fontsize_row = 7, fontsize_col = 7)
dev.copy(pdf,'Output/Figures/post-removal.pdf') # saving options for figure
dev.off()

# read stat for removed samples
stat <- stat %>% pivot_longer(col= 2:ncol(stat), 
                              names_to = "sample" , 
                              values_to = "count") %>% # piviting data for plotting
  filter(count > 0) # filtering out stats with no counts 

# ratio counting stats
remove.stat <- stat[which(stat$sample %in% low.cor),] # subsetting removed samples 
ggplot(remove.stat, aes(sample, count, fill = Status)) + 
  geom_bar(stat = "identity", position = "fill") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) # plotting by ratio
dev.copy(pdf,'Output/Figures/Fig2-stat-removal-ratio.pdf') # saving options for figure
dev.off()

# counting stats
options(scipen=10000) # gets rid of scientific notations
ggplot(remove.stat, aes(sample, count, fill = Status)) + 
  geom_bar(stat = "identity") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) # plotting out of ratio
dev.copy(pdf,'Output/Figures/stat-removal.pdf') # saving options for figure
options(scipen=1)
dev.off()

# saving options
write.table(counts, "Output/Data/counts.tab", sep = "\t")

# adding phenotype to pheno.data
pheno.data <- pheno.data[match(rownames(pheno), pheno.data$index_sequence),] # matching with pheno
pheno.data$pheno <- pheno$pheno # adding pheno data
pheno.data$batch_number <- factor(pheno.data$batch_number) # factoring batch_number

## Step: propd analysis
# Extract counts for ribosomal proteins to analyze differential proportionality
rp = grep ( "^RP[L|S]", row.names(counts), perl =T ) 
rps6k = grep ( "RPS6K", row.names(counts)) 
dash = grep ("-", row.names(counts))
rp = setdiff(rp, rps6k) 
rp = setdiff(rp, dash)
# Will remove RPS4Y1, Y2 and X for proportionality
to_remove = c(82, 83, 87, 88)
# Will keep all the RP-like proteins for now. 
rp = rp[-to_remove]
row.names(counts)[rp] 
rp_counts = t( counts[rp, ] ) 

#droping s from pheno.data
idx <- which(pheno.data$pheno == "S")
pheno.data[idx] <- "C"
pheno.data$pheno <- factor(pheno.data$pheno)

# There is an option to add voom weights using weighted = T
# Unclear to me whether this would be appropriate if using rp counts alone
# alpha is a positive and its value is used. 
pd = propd(rp_counts, as.character(pheno.data$pheno), alpha= NA, p = 1000 ) 
# Figure out the difference between disjointed and emergent
# We seem to like emergent based on the description
theta_e <- setEmergent(pd)
theta_e = updateCutoffs(theta_e, cutoff = seq(0.01, 0.91, 0.3))
tab <- getResults(theta_e)
# Plot RPL4 and RPS3
plot(theta_e@counts[, grep ("^RPL4$", colnames(theta_e@counts)  ) ], 
     theta_e@counts[, grep ("^RPS3$", colnames(theta_e@counts)  ) ], 
     col = ifelse(theta_e@group == "C", "red", "blue"))

grp1 <- theta_e@group == "C"
grp2 <- theta_e@group == "W"
abline(a = 0, b = theta_e@counts[grp1, grep ("^RPS3$", colnames(theta_e@counts)  ) ] /
         theta_e@counts[grp1, grep ("^RPL4$", colnames(theta_e@counts)  ) ], col = "red")
abline(a = 0, b = theta_e@counts[grp2, grep ("^RPS3$", colnames(theta_e@counts)  ) ] / 
         theta_e@counts[grp2, grep ("^RPL4$", colnames(theta_e@counts)  ) ], col = "blue")

plot(theta_e@counts[, grep ("^RPL4$", colnames(theta_e@counts)  )] / 
       theta_e@counts[, grep ("^RPS3$", colnames(theta_e@counts)  )],
     col = ifelse(theta_e@group == "C", "red", "blue"))

# # RPL11 is not in the top list
# parallel(theta_e, include = "RPL4")


## Step: DESeq2 analysis - collasping replicates - no batch effects
counts <- counts[,match(pheno.data$index_sequence, colnames(counts))] # ordering counts with pheno
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = pheno.data,
                              design = ~  pheno) # experiment declaration 
counts(dds)

# PCA plot with VST transform 
vsd <- vst(dds)
# assay(vsd) gives the transformed counts data
pcaData <- DESeq2::plotPCA(vsd, intgroup = c("pheno","batch_number"))

# indiv labels
ggplot(pcaData$data)  + 
  geom_text(aes(PC1, PC2, color = pheno), label = pheno.data$individual_number) + # with replicate labels
  ylab(pcaData$labels$y) +
  xlab(pcaData$labels$x)
dev.copy(pdf,'Output/Figures/Fig3-PCA-individual.pdf') # saving options for figure
dev.off()

# batch_number labels
ggplot(pcaData$data)  + 
  geom_point(aes(PC1, PC2, color = batch_number)) + # with batch shapes
  ylab(pcaData$labels$y) +
  xlab(pcaData$labels$x) + 
  scale_color_manual(values = c("#CC79A7", "#0072B2", "#F0E442", "#D55E00"))
dev.copy(pdf,'Output/Figures/Fig4-PCA-batch-shape.pdf') # saving options for figure
dev.off()

# differential expression 
dds <- collapseReplicates(dds, groupby = dds$individual_number)
collapsed.counts <- "counts"(dds)
dds <- DESeq(dds)
plotDispEsts(dds) # dispersion plot
dev.copy(pdf,'Output/Figures/dispersion.pdf') # saving options for figure
dev.off()

res <- results(dds, contrast = c("pheno", "C", "W")) #specify group comparisons here
summary(res) # summary of results for group comparison

# top results
top <- rownames(res)[which(res$padj < 0.05)] # view names of sig genes 
plotCounts(dds, "RPL11", "pheno") # plot normalized counts of specific genes by group

# MA plot
plotMA(res, ylim=c(-4,4)) 

# skrinkage estimator
resultsNames(dds) # select coef for srinkage
resLFC <- lfcShrink(dds, coef = 3, type = "apeglm")
plotMA(resLFC, ylim=c(-1,1))


data(iris)
keep <- iris$Species %in% c("setosa", "versicolor")
counts <- iris[keep, 1:4] * 10
group <- ifelse(iris[keep, "Species"] == "setosa", "A", "B")
group
keep
