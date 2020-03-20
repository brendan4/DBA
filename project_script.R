#!diagnostics off
library(tidyr)
library(dplyr)
library(DESeq2)
library(pheatmap)
library(ggplot2)
library(readr)
library(propr)
library(variancePartition)
library(limma)
library(edgeR)
source ('./mapCode-colors.r') # loading colors

# importing data
counts <- read.table("Input_data/unprocessed_counts.tab", row.names = 1, check.names = F)
pheno <- read.table("Input_data/pheno.tab",row.names = 1, check.names = F)
stat <- read.table("Input_data/stat.tab", row.names = 1, check.names = F)
pheno.data <- read_csv("Output/data/pheno_out.txt")
pheno.data <- as.data.frame(pheno.data)

# cleaning sample names
colnames(counts) <- gsub(".sorted.bam$", "", colnames(counts)) # removes files extension
colnames(counts) <- gsub("^121317.", "", colnames(counts)) # removes nonunique tag from lib samples
counts <- counts %>% dplyr::select(-contains("unmatched")) # remove unmatched cols


# matching pheno.data lib names to count data convention
idx <- which(gsub(".\\d*$", "",colnames(counts)) %in%
               gsub("^121317.", "", pheno.data$library_ID)) # match earlier hiseq rows
pheno.data[which(pheno.data$lane=="earlier_HiseqRun"), 
           "index_sequence"] <- colnames(counts)[idx] # using counts naming convention for earlier hiseq rows

# forming tech rep information 
pheno.data <- pheno.data %>% mutate (tech = gsub("^\\d+_", "", pheno.data$library_ID) ) 
counts <- counts[,-which(!colnames(counts) %in% pheno.data$index_sequence)] # not in pheno table
rownames(stat) <- stat[,1]
stat <-  stat[, -which(!colnames(stat) %in% pheno.data$index_sequence)] #protecting stat -> shift to rownames
stat["Status"] <- rownames(stat) # move stat out of rownames
stat <- stat[c(ncol(stat),seq(ncol(stat)-1))] # reordering stat


# creating better names for samples
indiv <- unique(pheno.data$individual_number)
for( i in 1:length(indiv) ) { # cycle through each individual
  idx <- which(pheno.data$individual_number %in% indiv[i]) # match all entries that belong to a individual
  samples <- unique(pheno.data$tech[idx]) # unique samples 
  for( j in 1:length(samples) ) { # cycle through each sample for that individual
    idx1 <- which(pheno.data$tech %in% samples[j]) # obtain idx for all entries with belong to a sample
    name <- paste(indiv[i], "_s", j, sep = "") # forming new name 
    if (length(idx1) > 1) {
      techs <- rep(paste(name, "_t", seq(length(idx1)), sep = ""), length(idx1) )
      pheno.data$individual_number[idx1] <- techs
    } else {
      pheno.data$individual_number[idx1] <- name
    }
    
  }
}


# unifying names
counts <- counts[,match(pheno.data$index_sequence, colnames(counts))] # ordering counts with pheno
colnames(counts) <- pheno.data$individual_number # replacing colnames of counts 
stat[2:ncol(stat)] <- stat[,match(pheno.data$index_sequence, colnames(stat) )] # ordering stat with pheno
colnames(stat)[2:ncol(stat)] <- pheno.data$individual_number

## Step: removal of samples
cor.table <- cor(counts, method = "spearman") # spearman cor table
mat_cluster_cols <- hclust(dist(t(counts)))


# cor heat map: pre-removal
breaksList = seq(0, 1, by = .01)
pheatmap(cor.table, 
         fontsize_row = 7, 
         fontsize_col = 7, 
         breaks = breaksList,
         treeheight_row = 0,
         cluster_cols = mat_cluster_cols,
         cluster_rows = mat_cluster_cols,
         color = color.palette2(length(breaksList)) )


dev.copy(pdf,'Output/Figures/Fig1-pre-removal.pdf') # saving options for figure
dev.off()

# sample removals
pheno.sex <- pheno.data$individual_number[which(pheno.data$index_sequence %in% "L6.TTAGGC")]# pheno discrepancy: gender and no pheno
low.cor <- names(which(apply(cor.table, 2, mean) < 0.9)) # samples names with cor > 0.9
remove <- c(low.cor, pheno.sex) 
counts <- counts %>% dplyr::select(-one_of(remove)) # removing samples from counts



#boxplot of counts
colsum <- as.data.frame( colSums(counts))
colsum$name <- rownames(colsum)
colnames(colsum) <- c("count", "names")
boxp1 <- ggplot(colsum) + 
  geom_boxplot(outlier.colour = "red", outlier.shape = 1, aes(y = count)) +
  ylab("Assigned Counts") +
  theme(axis.text.x = element_blank()) 
boxp1

dev.copy(pdf,'Output/Figures/assigned-boxplot.pdf') # saving options for figure
dev.off()

# cor heat map: post-removal
cor.table <- cor(counts, method = "spearman") # spearman cor table
breaksList = seq(0, 1, by = .01)
pheatmap(cor.table, 
         fontsize_row = 7, 
         fontsize_col = 7, 
         breaks = breaksList,
         treeheight_row = 0, # remove row dendrograms
         na_col = "white",
         color = color.palette2(length(breaksList)) )

dev.copy(pdf,'Output/Figures/post-removal.pdf') # saving options for figure
dev.off()

# read stat for removed samples
stat <- stat %>% pivot_longer(col= 2:ncol(stat), 
                              names_to = "sample" , 
                              values_to = "count") %>% # piviting data for plotting
  filter(count > 0) %>% # filtering out stats with no counts 
  filter(Status != "Unassigned_Ambiguity") # filter out low rep stat 

# mean statistics of kept samples
kept.stat <- stat[which(!stat$sample %in% low.cor), ] %>% 
  group_by(Status) %>% 
  summarize(count = mean(count)) %>% 
  mutate(sample = "Mean")
remove.stat <- stat[which(stat$sample %in% low.cor), ] # statistics for remove samples only

# ratio of counting stats
ggplot() + 
  geom_bar(data = remove.stat, 
           aes(sample, count, fill = Status), 
           stat = "identity", position = "fill") + 
  geom_bar(data = kept.stat, 
           aes(sample, count, fill = Status), 
           stat = "identity", position = "fill") +
  geom_text(data = remove.stat %>%
              group_by(sample) %>% 
              summarize(sums = sum(count)), 
            aes(x = sample, label = sums, y = 1.025), 
            size = 4, color = "gray25") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual("legend", values = c("Assigned" = "#7F0000", 
                                         "Unassigned_NoFeatures" = "#EF6548",
                                         "Unassigned_Unmapped" = "#FDD49E")) # colors from brewer.heat

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

# removal of samples from pheno
pheno.data <- pheno.data[-which(!pheno.data$individual_number %in% colnames(counts)),] # removing samples from pheno.data
# adding phenotype to pheno.data
pheno.data <- pheno.data[match(rownames(pheno), pheno.data$index_sequence),] # matching with pheno
pheno.data$pheno <- pheno$pheno # adding pheno data
pheno.data$batch_number <- factor(pheno.data$batch_number) # factoring batch_number

# collapsing replicates without deseq
counts <- counts[,match(pheno.data$individual_number, colnames(counts))] # ordering counts with pheno
sp <- split(seq(along = pheno.data$tech), pheno.data$tech) # splitting and obtaining index of tech for each sample 
countdata <- sapply(sp, function(i) rowSums(counts[ , i, drop = FALSE])) # summing the rows of count data based on tech status 
idx <- sapply(sp, function(i) i[1]) # obtaining idx for each sample
colnames(countdata) <- pheno.data$individual_number[idx] # replacing names
pheno.coll <- pheno.data[idx,] # collapsing pheno data
pheno.coll <- pheno.coll [order(pheno.coll$pheno), ]
countdata<- countdata[,match(pheno.coll$individual_number, colnames(countdata))] # ordering counts with pheno

## Step: propd analysis
# Extract counts for ribosomal proteins to analyze differential proportionality
rp = grep ( "^RP[L|S]", row.names(countdata), perl =T ) 
rps6k = grep ( "RPS6K", row.names(countdata)) 
dash = grep ("-", row.names(countdata))
rp = setdiff(rp, rps6k) 
rp = setdiff(rp, dash)
row.names(countdata)[rp] 
grep ( "L1$", row.names(countdata)[rp] , perl = T )
# Will remove RPS4Y1, Y2 and X for proportionality
to_remove = c(16, 17, 25, 50, 57, 82, 83, 84, 87, 88)
# Will keep all the RP-like proteins for now. 
rp = rp[-to_remove]
row.names(countdata)[rp] 
rp_counts = t( countdata[rp, ] ) 
rowMeans(countdata[rp, ])


#OPTIONAL
# droping s from pheno.data
idx <- which(pheno.coll$pheno == "S")
pheno.coll$pheno[idx] <- "C"
pheno.coll$pheno <- factor(pheno.coll$pheno)

#OPTIONAL
# changing C to mutation
idx <- which(pheno.coll$pheno == "C")
pheno.coll$pheno <- as.character(pheno.coll$pheno)
pheno.coll$pheno[idx] <- "c.396+3A>G"

#OPTIONAL
# wiltype to non-C 
idx <- which(pheno.coll$pheno == "W")
pheno.coll$pheno <- as.character(pheno.coll$pheno)
pheno.coll$pheno[idx] <- "Non-Carr"
pheno.coll$pheno <- factor(pheno.coll$pheno)

# There is an option to add voom weights using weighted = T
# Unclear to me whether this would be appropriate if using rp counts alone
# alpha is a positive and its value is used. 
pd = propd(rp_counts, as.character(pheno.coll$pheno), alpha= NA, p = 1000 )  

# Figure out the difference between disjointed and emergent
# We seem to like emergent based on the description
theta_e <- setEmergent(pd)
theta_e = updateCutoffs(theta_e, cutoff = seq(0.05, 0.95, 0.3))
tab <- getResults(theta_e)

# Plot RPL27 and RPS3A
plot(theta_e@counts[, grep ("^RPL27$", colnames(theta_e@counts)  ) ], 
     theta_e@counts[, grep ("^RPS3A$", colnames(theta_e@counts)  ) ], 
     col = ifelse(theta_e@group == "C", "red", "blue"))

grp1 <- theta_e@group == "C"
grp2 <- theta_e@group == "W"
abline(a = 0, b = theta_e@counts[grp1, grep ("^RPS3A$", colnames(theta_e@counts)  ) ] /
         theta_e@counts[grp1, grep ("^RPL27$", colnames(theta_e@counts)  ) ], col = "red")
abline(a = 0, b = theta_e@counts[grp2, grep ("^RPS3A$", colnames(theta_e@counts)  ) ] / 
         theta_e@counts[grp2, grep ("^RPL27$", colnames(theta_e@counts)  ) ], col = "blue")


# gene pair selection
gene.pair.list <- data.frame(partner = 1:10, pair = 1:10)
for (i in 1:10){
 gene.pair.list[i, ] <- tab[i, c("Partner", "Pair")]
}

# plots gene relationships
plot.ratio <- function (gene.list){
  for( i in 1:nrow(gene.pair.list)){
    var.pair <- data.frame(ratio = theta_e@counts[, grep (paste("^",gene.pair.list[i,1],"$", sep= ""), colnames(theta_e@counts)  )] / 
                             theta_e@counts[, grep (paste("^",gene.pair.list[i,2],"$", sep= ""), colnames(theta_e@counts)  )], 
                           pheno = pheno.coll$pheno)
    print(ggplot(var.pair) + 
      geom_point(aes(x = seq(1:nrow(var.pair)),
                     y = ratio ,
                     color = pheno),
                 size = 2) +
      xlab("Sample Index") +
      ylab("Ratio") +
      ggtitle(paste(gene.pair.list[i,1], "/",  gene.pair.list[i,2])))
  }
  
}

plot.ratio(gene.list = gene.pair.list)
#preparing data for ggplot
var.pair <- data.frame(ratio = theta_e@counts[, grep ("^RPL27$", colnames(theta_e@counts)  )] / 
                         theta_e@counts[, grep ("^RPS3A$", colnames(theta_e@counts)  )], 
                       pheno = pheno.coll$pheno)
ggplot(var.pair) + 
  geom_point(aes(x = seq(1:nrow(var.pair)),
                 y = ratio ,
                 color = pheno)) +
  xlab("Sample Index") +
  ylab("Ratio") +
  ggtitle("RPL27 / RPS3A")

# base r plot
plot(theta_e@counts[, grep ("^RPL27$", colnames(theta_e@counts)  )] / 
       theta_e@counts[, grep ("^RPS3A$", colnames(theta_e@counts)  )],
     col = ifelse(theta_e@group == "C", "red", "blue"))

# basic plots 
plot(colMeans(rp_counts[which(pheno.coll$pheno == "C"), ]), colMeans(rp_counts[which(pheno.coll$pheno == "W"), ]))
boxplot(rp_counts[which(pheno.coll$pheno == "C"), "RPL11"], rp_counts[which(pheno.coll$pheno == "W"), "RPL11"])
ggplot(tab, aes(x = lrv1, y = lrv2)) + geom_point() + geom_abline(intercept = 0, color = "red")
plot(tab[,"lrv1"], tab[,"lrv2"],col = ifelse( grepl("L", tab$Partner) & 
                                                grepl("L", tab$Pair)  == F, "Red", "black"))

sum(colMeans(rp_counts[which(pheno.coll$pheno == "C"), ])/ colMeans(rp_counts[which(pheno.coll$pheno == "W") , ]) < 1)

# making data long format for ggplot 
rp_long <- as.data.frame(rp_counts)
rp_long$indiv <- rownames(rp_long)
rp_long <- rp_long %>% pivot_longer( cols = 1:(ncol(rp_long) -1), names_to = "genes", values_to = c("counts"))
rp_long <- rp_long %>% mutate(pheno = case_when(  rp_long$indiv %in% pheno.coll$individual_number[which(pheno.coll$pheno == "C")] ~ "C",
                                                  rp_long$indiv %in% pheno.coll$individual_number[which(pheno.coll$pheno == "W")] ~ "W",
                                                  rp_long$indiv %in% pheno.coll$individual_number[which(pheno.coll$pheno == "S")] ~ "S"))
# box plot of genes
rp_long %>% filter(genes == 'RPL11') %>%
  ggplot(aes(x = pheno, y = counts)) + 
  labs(title = "RPL11 counts", y = "Counts", x = "Phenotype") + 
  geom_boxplot(aes(fill = pheno), outlier.colour="black", outlier.shape=1, outlier.size=3)

dev.copy(pdf,'Output/Figures/RPL11-counts.pdf') # saving options for figure
dev.off()


parallel(theta_e, include = "RPL27")

## To generate the graph of variance, we will need to calculate variance of log-ratios for the two groups
## VLR (X, Y) = var (log (X/Y))
## We will do a nested for-loop if there is a more clever algorithm feel free to implement
gene.num <- ncol(rp_counts)
carrier_wt_vlr = matrix(nrow = (gene.num*(gene.num-1) )/ 2, ncol = 2, dimnames = list(c(rep("NA", gene.num*(gene.num-1)/2)) ,c("C", "W")))

row_id = 1
vlr = function (x,y) { 
  return ( var ( log(x/y) ) ) 
}
for ( rp_ind in 1:dim(rp_counts)[2] ) { 
  for (rp2_ind in seq( (rp_ind+1), dim(rp_counts)[2]) )  {
    pair_name = paste( colnames(rp_counts)[rp_ind], colnames(rp_counts)[rp2_ind], sep = "_" )
    row.names(carrier_wt_vlr)[row_id] = pair_name
    carrier_wt_vlr[row_id, "C"] = vlr ( rp_counts[pheno.coll$pheno == "C",rp_ind], rp_counts[pheno.coll$pheno == "C",rp2_ind])
    carrier_wt_vlr[row_id, "W"] =  vlr ( rp_counts[pheno.coll$pheno == "W",rp_ind], rp_counts[pheno.coll$pheno == "W",rp2_ind]) 
    row_id = row_id + 1
  }
}


# OPTIONAL: droping s from pheno.data
idx <- which(pheno.data$pheno == "S")
pheno.data$pheno[idx] <- "C"
pheno.data$pheno <- factor(pheno.data$pheno)

# changing C to mutation
idx <- which(pheno.data$pheno == "C")
pheno.data$pheno <- as.character(pheno.data$pheno)
pheno.data$pheno[idx] <- "c.396+3A>G"


# wiltype to non-C 
idx <- which(pheno.data$pheno == "W")
pheno.data$pheno <- as.character(pheno.data$pheno)
pheno.data$pheno[idx] <- "Non-Carr"
pheno.data$pheno <- factor(pheno.data$pheno)


## Step: DESeq2 analysis - collasping replicates - no batch effects
counts <- counts[,match(pheno.data$individual_number, colnames(counts))] # ordering counts with pheno
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = pheno.data,
                              design = ~pheno) # experiment declaration 
counts(dds)

# PCA plot with VST transform 
vsd <- vst(dds)
# assay(vsd) gives the transformed counts data
pcaData <- DESeq2::plotPCA(vsd, intgroup = c("pheno","batch_number"))

# indiv labels
ggplot(pcaData$data)  + 
  geom_text(aes(PC1, PC2, color = pheno), label = pheno.data$individual_number) + # with replicate labels
  ylab(pcaData$labels$y) +
  xlab(pcaData$labels$x) + 
  scale_color_manual(values = c( "#d11f12", "#547294", "#F0E442"))
dev.copy(pdf,'Output/Figures/Fig3-PCA-individual.pdf') # saving options for figure
dev.off()

# batch_number labels
ggplot(pcaData$data)  + 
  geom_point(aes(PC1, PC2, color = batch_number), size = 3) + # with batch shapes
  ylab(pcaData$labels$y) +
  xlab(pcaData$labels$x) + 
  scale_color_manual(values = c("#CC79A7", "#0072B2", "#F0E442", "#D55E00"))
dev.copy(pdf,'Output/Figures/Fig4-PCA-batch-shape.pdf') # saving options for figure
dev.off()

# differential expression 
dds <- collapseReplicates(dds, groupby = dds$tech)
collapsed.counts <- "counts"(dds)
dds <- DESeq(dds)
plotDispEsts(dds) # dispersion plot
dev.copy(pdf,'Output/Figures/dispersion.pdf') # saving options for figure
dev.off()

res <- results(dds, contrast = c("pheno", "W", "C")) #specify group comparisons here
summary(res) # summary of results for group comparison
resultsNames(dds)

# filtering 
plot(metadata(res)$filterNumRej, type="b", xlab="quantiles of 'baseMean'",
     ylab="number of rejections")
dev.copy(pdf,'Output/Figures/DESeq-flitering.pdf') # saving options for figure
dev.off()

# top results
top <- rownames(res)[which(res$padj < 0.05)] # view names of sig genes 
plotCounts(dds, "RPL11", "pheno") # plot normalized counts of specific genes by group
dev.copy(pdf,'Output/Figures/RPL11-pheno.pdf') # saving options for figure
dev.off()

write.table(x = data.frame(gene = rownames(res)[which(res$padj < 0.05)], 
                           logfold = res$log2FoldChange[which(res$padj < 0.05)],  
                           padj = res$padj[which(res$padj < 0.05)], basemean = res$baseMean[which(res$padj < 0.05)]), 
            file = "gene_names.csv", sep = ",")
deseq.data = data.frame(gene = rownames(res)[which(res$padj < 0.05)], 
               logfold = res$log2FoldChange[which(res$padj < 0.05)],  
               padj = res$padj[which(res$padj < 0.05)], basemean = res$baseMean[which(res$padj < 0.05)])

# ploting specific genes
norm <- counts(dds, normalized=TRUE)
gene.plot <- t(norm[ which(rownames(norm) %in% c("RPL11", "CDK11A", "CRIM1", "GATA2", "MCC")), ])
carr.pheno <- pheno.data[which(pheno.data$pheno == "C"), ]
carr.genes <- gene.plot[which( rownames(gene.plot) %in% carr.pheno$tech), ]
carr.genes <- as.data.frame(carr.genes)
carr.means <- as.data.frame(t(colMeans(carr.genes)))

wild.pheno <- pheno.data[which(pheno.data$pheno == "W"), ]
wild.genes <- gene.plot[which( rownames(gene.plot) %in% wild.pheno$tech), ]
wild.genes <- as.data.frame(wild.genes)
wild.means <- as.data.frame(t(colMeans(wild.genes)))

ggplot() + geom_point(data = wild.genes, aes(y = RPL11, x = "W"), color = "#547294") + 
  geom_point(data = carr.genes, aes(y= RPL11, x = "c.396+3A>G"), color = "#d11f12") + 
  geom_bar(data = wild.means, aes(y = RPL11, "W"), stat = "identity", color = "#547294", alpha = .2)+
  geom_bar(data = carr.means, aes(y = RPL11, "c.396+3A>G"), stat = "identity", color = "#d11f12", alpha = .2) +
  labs(x = "Genotype", y = "RPL11 Normalized")+
  theme_classic()
dev.copy(pdf,'Output/Figures/RPL11.pdf') # saving options for figure
dev.off()

ggplot() + geom_point(data = wild.genes, aes(y = CDK11A, x = "Non-carriers"), color = "#547294") + 
  geom_point(data = carr.genes, aes(y= CDK11A, x = "c.396+3A>G"), color = "#d11f12") + 
  geom_bar(data = wild.means, aes(y = CDK11A, "Non-carriers"), stat = "identity", color = "#547294", alpha = .2)+
  geom_bar(data = carr.means, aes(y = CDK11A, "c.396+3A>G"), stat = "identity", color = "#d11f12", alpha = .2) +
  labs(x = "Genotype", y = "CDK11A Normalized")+
  theme_classic()
dev.copy(pdf,'Output/Figures/CDK11A.pdf') # saving options for figure
dev.off()

ggplot() + geom_point(data = wild.genes, aes(y = CRIM1, x = "Non-carriers"), color = "#547294") + 
  geom_point(data = carr.genes, aes(y= CRIM1, x = "c.396+3A>G"), color = "#d11f12") + 
  geom_bar(data = wild.means, aes(y = CRIM1, "Non-carriers"), stat = "identity", color = "#547294", alpha = .2)+
  geom_bar(data = carr.means, aes(y = CRIM1, "c.396+3A>G"), stat = "identity", color = "#d11f12", alpha = .2) +
  labs(x = "Genotype", y = "CRIM1 Normalized")+
  theme_classic()
dev.copy(pdf,'Output/Figures/CRIM1.pdf') # saving options for figure
dev.off()

ggplot() + geom_point(data = wild.genes, aes(y = GATA2, x = "Non-carriers"), color = "#547294") + 
  geom_point(data = carr.genes, aes(y= GATA2, x = "c.396+3A>G"), color = "#d11f12") + 
  geom_bar(data = wild.means, aes(y = GATA2, "Non-carriers"), stat = "identity", color = "#547294", alpha = .2)+
  geom_bar(data = carr.means, aes(y = GATA2, "c.396+3A>G"), stat = "identity", color = "#d11f12", alpha = .2) +
  labs(x = "Genotype", y = "GATA2 Normalized")+
  theme_classic()
dev.copy(pdf,'Output/Figures/GATA2.pdf') # saving options for figure
dev.off()

ggplot() + geom_point(data = wild.genes, aes(y = MCC, x = "Non-carriers"), color = "#547294") + 
  geom_point(data = carr.genes, aes(y= MCC, x = "c.396+3A>G"), color = "#d11f12") + 
  geom_bar(data = wild.means, aes(y = MCC, "Non-carriers"), stat = "identity", color = "#547294", alpha = .2)+
  geom_bar(data = carr.means, aes(y = MCC, "c.396+3A>G"), stat = "identity", color = "#d11f12", alpha = .2) +
  labs(x = "Genotype", y = "MCC Normalized")+
  theme_classic()
dev.copy(pdf,'Output/Figures/MCC.pdf') # saving options for figure
dev.off()


# MA plot
DESeq2::plotMA(res, ylim=c(-4,4)) 
text(600,3, label = "Non-carriers")
text(600,-3, label = "c.396+3A>G" , col = "dodgerblue4")
dev.copy(pdf,'Output/Figures/MA-plot.pdf') # saving options for figure
dev.off()

# skrinkage estimator
DESeq2::resultsNames(dds) # select coef for srinkage
resLFC <- DESeq2::lfcShrink(dds, coef = "pheno_W_vs_C", type = "ashr")
DESeq2::plotMA(resLFC, ylim=c(-3,3))


## Adding dream/voom/limma for mixed effects model of differential expression 
# Will use countdata -> technical replicates collapsed 
# Matching phenotypes are in pheno.coll

# 1/3 of the samples is a good starting point
isexpr = rowSums(cpm (countdata) > 1 ) >= 1
geneExpr = DGEList( countdata[isexpr, ])
geneExpr = calcNormFactors( geneExpr, method = "TMM")

design_DE = as.data.frame ( pheno.coll[,c(3, 9)] ) 
design_DE$Individual = sapply ( strsplit(design_DE$individual_number , split = "_"), "[[" , 1 ) 
design_DE$StatusSubType = as.character( design_DE$pheno ) 
design_DE$StatusSubType[11:14 ] = rep ( "S" , 4 ) 
design_DE = design_DE [, -1]
row.names(design_DE) = colnames(countdata)

## We will apply dream to fit the mixed effects model 
form <- ~ pheno + (1|Individual)
vobjDream = voomWithDreamWeights(geneExpr, form, design_DE, plot = T)
fitmm = dream(vobjDream, form, design_DE)
fitmm$design

topTable ( fitmm, coef = "phenoW", number = 20)

# We can alternatively fit with three levels
form <- ~ StatusSubType + (1|Individual)
vobjDream = voomWithDreamWeights(geneExpr, form, design_DE, plot = T)
fitmm = dream(vobjDream, form, design_DE)
head ( fitmm$design ) 
limma <- as.data.frame (topTable ( fitmm, coef = 2, number = 10000) )
limma$gene_names <- rownames(limma)
colnames(limma) <- paste("dream", colnames(limma), sep = "_")
limma.both <- limma[which(rownames(limma) %in% deseq.data$gene), ]

#exporting limma data
merged.tables <- merge(deseq.data, limma.both, by.x = "gene", by.y = "dream_gene_names")
write.table(merged.tables, file = "dseq-dream.csv", sep = ",")



## STEP: ribosome plots
# collapsing replicates without deseq
counts <- counts[,match(pheno.data$individual_number, colnames(counts))] # ordering counts with pheno
sp <- split(seq(along = pheno.data$tech), pheno.data$tech) # splitting and obtaining index of tech for each sample 
countdata <- sapply(sp, function(i) rowSums(counts[ , i, drop = FALSE])) # summing the rows of count data based on tech status 
idx <- sapply(sp, function(i) i[1]) # obtaining idx for each sample
colnames(countdata) <- pheno.data$individual_number[idx] # replacing names
pheno.coll <- pheno.data[idx,] # collapsing pheno data
pheno.coll <- pheno.coll [order(pheno.coll$pheno), ]

# plotting ribosomal genes
norm <- counts(dds, normalized=TRUE)
norm <- norm[,match(pheno.coll$tech, colnames(norm))] # ordering counts with pheno

rp = grep ( "^RP[L|S]", row.names(norm), perl =T ) 
rps6k = grep ( "RPS6K", row.names(norm)) 
dash = grep ("-", row.names(norm))

rp = setdiff(rp, rps6k) 
rp = setdiff(rp, dash)

row.names(norm)[rp] 
grep ( "L1$", row.names(norm)[rp] , perl = T )
# Will remove RPS4Y1, Y2 and X for proportionality
to_remove = c(16, 17, 25, 50, 57, 82, 83, 84, 87, 88)
# Will keep all the RP-like proteins for now. 
rp = rp[-to_remove]
rps = grep ( "^RPS", row.names(norm)[rp], perl =T ) 
rpl = grep ( "^RPL", row.names(norm)[rp], perl =T ) 

row.names(norm)[rp[rpl]] 
row.names(norm)[rp[rps]]

rpl_counts <- norm[rp[rpl], ]
rps_counts <- norm[rp[rps], ]


wild.idx <- pheno.coll$tech[which(pheno.coll$pheno == "W")]
car.idx <- pheno.coll$tech[which(pheno.coll$pheno == "C")]

# large subunit 
par(mfrow=c(1,2))
boxplot(colSums(rpl_counts[, which(colnames(rpl_counts) %in% wild.idx)]), ylim = c(60000, 150000))
title("Noncarrier - large subunit")
boxplot(colSums(rpl_counts[, which(colnames(rpl_counts) %in% car.idx)]), ylim = c(60000, 150000))
title("c.396+3A>G - large subunit")
dev.copy(pdf,'Output/Figures/large-ribo-boxplot.pdf') # saving options for figure
dev.off()

# small subunit 
par(mfrow=c(1,2))
boxplot(colSums(rps_counts[,which(colnames(rps_counts) %in% wild.idx)]), yim = c(40000, 80000))
title("Noncarrier - small subunit")
boxplot(colSums(rps_counts[,which(colnames(rps_counts) %in% car.idx)]), yim = c(40000, 80000))
title("c.396+3A>G- small subunit")
dev.copy(pdf,'Output/Figures/small-ribo-boxplot.pdf') # saving options for figure
dev.off()

rpl_counts <- rpl_counts[,match(pheno.coll$tech, colnames(rpl_counts))] # ordering counts with pheno
rpl_long <- as.data.frame(rpl_counts)
colnames(rpl_long) <- pheno.coll$individual_number
rpl_long$indiv <- rownames(rpl_long)
rpl_long <- rpl_long %>% pivot_longer(cols = 1:(ncol(rpl_long)-1), names_to = "sample", values_to = c("counts"))
rpl_long<- rpl_long %>% mutate(pheno = case_when(rpl_long$sample %in% pheno.coll$individual_number[which(pheno.coll$pheno == "C")] ~ "C",
                                                  rpl_long$sample %in% pheno.coll$individual_number[which(pheno.coll$pheno == "W")] ~ "W",
                                                  rpl_long$sample %in% pheno.coll$individual_number[which(pheno.coll$pheno == "S")] ~ "S"))

rpl_long %>% ggplot(aes(x = indiv, y = counts)) + 
  geom_point() + 
  facet_grid(~pheno) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
