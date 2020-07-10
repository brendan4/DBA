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
library(forcats)
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
pheno.data$indiv <- pheno.data$individual_number #preserve indiv data


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

stat <- stat %>% mutate(Status = 
                          case_when(Status == "Assigned" ~ "Mapped", 
                                    Status == "Unassigned_NoFeatures" ~ "Mapped", 
                                    Status == "Unassigned_Unmapped" ~ "Unmapped"))

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
  scale_fill_manual("legend", values = c("Mapped" = "#7F0000",
                                         "Unmapped" = "#FDD49E")) # colors from brewer.heat
dev.copy(pdf,'Output/Figures/Fig2-stat-removal-ratio.pdf') # saving options for figure
dev.off()

# counting stats: non ratio
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

#boxplot of counts
colsum <- as.data.frame(colSums(counts))
colsum$name <- rownames(colsum)
colnames(colsum) <- c("count", "names")
colsum <- colsum[match(pheno.data$individual_number, colsum$names ),] 
colsum$pheno <- pheno.data$pheno

# keep s
colsum <- colsum %>% mutate(pheno = case_when(pheno == "W" ~ "Non-Carr", pheno == "C" ~ "Carr", pheno == "S" ~ "DBA"))
colsum <- colsum %>% mutate(temp = case_when(pheno == "Non-Carr" ~ "W", pheno == "Carr" | pheno == 'DBA' ~ "M"))

# OPTIONAL: drop S
colsum <- colsum[-which(colsum$pheno %in% "S"), ] 
colsum <- colsum %>% mutate(pheno = case_when(pheno == "W" ~ "Non-Carr", pheno == "C" ~ "c.396+3A>G"))

boxp1 <- ggplot(colsum) + 
  geom_boxplot(outlier.colour = "red", outlier.shape = 1, aes(y = count)) +
  ylab("Assigned Counts") +
  theme(axis.text.x = element_blank()) + facet_grid(~pheno)
boxp1
dev.copy(pdf,'Output/Figures/assigned-boxplot.pdf') # saving options for figure
dev.off()

dotplot <- ggplot(colsum) +
  geom_point(aes(x = 1,y = count, fill = pheno, color = temp, group = pheno), 
             size = 2, shape = 21, position = position_dodge(width = .5))+
  theme(axis.text.x = element_blank(), axis.title.x = element_blank())+
  scale_fill_manual(values = c( "salmon1", "red3","gray55")) +
  scale_color_manual(values = c('red3','black'))

dotplot
dev.copy(pdf,'Output/Figures/assigned-dotplot.pdf') # saving options for figure
dev.off()

# collapsing replicates without deseq
collapse.rep <- function(counts, pheno.data) {
  counts <- counts[,match(pheno.data$individual_number, colnames(counts))] # ordering counts with pheno
  sp <- split(seq(along = pheno.data$tech), pheno.data$tech) # splitting and obtaining index of tech for each sample 
  countdata <- sapply(sp, function(i) rowSums(counts[ , i, drop = FALSE])) # summing the rows of count data based on tech status 
  idx <- sapply(sp, function(i) i[1]) # obtaining idx for each sample
  colnames(countdata) <- pheno.data$individual_number[idx] # replacing names
  pheno.coll <- pheno.data[idx,] # collapsing pheno data
  pheno.coll <- pheno.coll [order(pheno.coll$pheno), ]
  countdata<- countdata[,match(pheno.coll$individual_number, colnames(countdata))] # ordering counts with pheno
  data <- list(countdata, pheno.coll)
  return(data)
}


## Step: propd analysis
# collapsing replicates without deseq
data <- collapse.rep(counts, pheno.data)
countdata <- data[[1]]
pheno.coll <- data[[2]]

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
# converting S to C from pheno.data

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


#Droping S from rp_counts and pheno.coll
idx <- which(pheno.coll$pheno == "S")
rp_counts <- rp_counts[-idx,]
pheno.coll <- pheno.coll[-idx,]

# There is an option to add voom weights using weighted = T
# Unclear to me whether this would be appropriate if using rp counts alone
# alpha is a positive and its value is used. 
pd = propd(rp_counts, as.character(pheno.coll$pheno), alpha= NA, p = 1000 , weighted = F)  


# Figure out the difference between disjointed and emergent
# We seem to like emergent based on the description
theta_e <- setEmergent(pd)
theta_e = updateCutoffs(theta_e, cutoff = seq(0.05, 0.95, 0.3))
tab <- getResults(theta_e)
theta_e@fdr
updateF(theta_e)

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
carrier_wt_vlr = matrix(nrow = (gene.num*(gene.num-1) )/ 2, 
                        ncol = 2, 
                        dimnames = list(c(rep("NA", 
                                              gene.num*(gene.num-1)/2)),
                                        c("C", "W")))

row_id = 1
vlr = function (x,y) { 
  return ( var ( log(x/y) ) ) 
}
#rp2_ind = 2
for ( rp_ind in 1:dim(rp_counts)[2] ) { 
  for (rp2_ind in seq( (rp_ind+1), dim(rp_counts)[2]) )  {
    pair_name = paste( colnames(rp_counts)[rp_ind], colnames(rp_counts)[rp2_ind], sep = "_" )
    row.names(carrier_wt_vlr)[row_id] = pair_name
    carrier_wt_vlr[row_id, "C"] = vlr ( rp_counts[pheno.coll$pheno == "C",rp_ind], rp_counts[pheno.coll$pheno == "C",rp2_ind])
    carrier_wt_vlr[row_id, "W"] =  vlr ( rp_counts[pheno.coll$pheno == "W",rp_ind], rp_counts[pheno.coll$pheno == "W",rp2_ind]) 
    row_id = row_id + 1
  }
}
vlr = function (x,y) { 
  return ( var ( log(x/y) ) ) 
}
mlr = function (x,y) { 
  return ( mean ( log(x/y) ) ) 
}
wtm <- function (x, w){
  return ( sum (x * w) / sum(w) )
}
wtv <- function (x,w){
  xbar = wtm(x, w)
  return ( sum(w * (x - xbar) ** 2) * (sum(w) / ((sum(w) ** 2) - sum((w ** 2)))) )
}

#OLD METHOD - SEE propd_sampling: propE without propr package
ct <- t(countdata)
ct <- ct[-which(pheno.coll.all$pheno == "S"),]
pheno.coll.all <- pheno.coll # saving for later use if needed
pheno.coll <- pheno.coll[-which(pheno.coll$pheno == "S"),]
group <- as.character(pheno.coll$pheno)
per <- 100

# obtaining weights : calculate with zeros
design <- matrix(0, nrow = nrow(ct), ncol = 2)
design[group == unique(group)[1], 1] <- 1
design[group == unique(group)[2], 2] <- 1
v <- limma::voom(t(ct), design = design)
weights <- t(v$weights)
W <- weights 

# replacing the zeros
if (any(as.matrix(countdata) == 0) & is.na(alpha)) {
  message("Alert: Replacing 0s with next smallest value.")
  zeros <- ct == 0
  ct[zeros] <- min(ct[!zeros])}


# groups
group1 <- group == unique(group)[1]
group2 <- group == unique(group)[2]
n1 <- sum(group1)
n2 <- sum(group2)

# NOT FOR WEIGHTED propd 
p1 <- n1 - 1
p2 <- n2 - 1
p <- n1 + n2 - 1

# permutation building
message("Alert: Fixing permutations to active random seed.")
permutes <- as.data.frame(matrix(0, nrow = nrow(ct), 
                                 ncol = per))
for (col in 1:ncol(permutes)) permutes[, col] <- sample(1:nrow(ct))

rpl11.idx <- which(colnames(ct) %in% "RPL11")
c = c()
w = c()
gene = c()
genes = c()
lrv = c()
c.lrm = c()
w.lrm = c()

# non weighted
for (i in 1:ncol(ct)){
  gene[i] <- colnames(ct)[i]
  genes[i] <- paste(colnames(ct)[i],'_',colnames(ct)[rpl11.idx], sep ="")
  c[i] = vlr(ct[pheno.coll$pheno == "C",i], 
      ct[pheno.coll$pheno == "C", rpl11.idx])
  w[i] = vlr(ct[pheno.coll$pheno == "W",i], 
      ct[pheno.coll$pheno == "W", rpl11.idx])
  lrv[i] = vlr(ct[which(pheno.coll$pheno %in% c("C","W")),i], 
               ct[which(pheno.coll$pheno %in% c("C","W")),rpl11.idx])
  c.lrm[i] = mlr(ct[pheno.coll$pheno == "C",i], 
                 ct[pheno.coll$pheno == "C", rpl11.idx])
  w.lrm[i] = mlr(ct[pheno.coll$pheno == "W",i], 
                 ct[pheno.coll$pheno == "W", rpl11.idx])
}


# run for weighted
p1 <- c()
p2 <- c()
p <- c()

# weighed values 
for (i in 1:ncol(ct)){
  gene[i] <- colnames(ct)[i]
  genes[i] <- paste(colnames(ct)[i],'_',colnames(ct)[rpl11.idx], sep ="")
  c.wij = W[pheno.coll$pheno == "C",i] * W[pheno.coll$pheno == "C",rpl11.idx]
  c[i] = wtv(log(ct[pheno.coll$pheno == "C", i] / ct[pheno.coll$pheno == "C", rpl11.idx]), c.wij)
  
  w.wij = W[pheno.coll$pheno == "W",i] * W[pheno.coll$pheno == "W",rpl11.idx]
  w[i] = wtv(log(ct[pheno.coll$pheno == "W", i] / ct[pheno.coll$pheno == "W", rpl11.idx]), w.wij)
  
  wij = W[which(pheno.coll$pheno %in% c("C","W")),i] * W[which(pheno.coll$pheno %in% c("C","W")),rpl11.idx]
  lrv[i] = wtv(log(ct[which(pheno.coll$pheno %in% c("C","W")), i] / ct[which(pheno.coll$pheno %in% c("C","W")), rpl11.idx]), wij)
  
  c.lrm[i] = mlr(ct[pheno.coll$pheno == "C",i], 
                 ct[pheno.coll$pheno == "C", rpl11.idx])
  w.lrm[i] = mlr(ct[pheno.coll$pheno == "W",i], 
                 ct[pheno.coll$pheno == "W", rpl11.idx])
 
  # calculating omega for weights 
  c.n = sum(c.wij)
  p1[i] = c.n - sum(c.wij**2)/ c.n 
  
  w.n = sum(w.wij)
  p2[i] = w.n - sum(w.wij**2)/ w.n
  
  n = sum(wij)
  p[i] = n - sum(wij**2)/ n
  
}


theta <- (p1 * c + p2 * w)/(p * lrv)
theta_e <- 1 - pmax(p1 * c, p2 * w)/(p * lrv)
theta_f <- pmax(p1 * c, p2 * w)/(p * lrv)
theta_g <- pmin(p1 * c, p2 * w)/(p * lrv)
ct.rpl11 <- data.frame(gene = gene, genes = genes, c= c, w = w, lrv = lrv,
                       c.lrm = c.lrm, w.lrm = w.lrm,
                       theta = theta, theta_e = theta_e, theta_f = theta_f, 
                       theta_g = theta_g )


# Permutations
cutoff = seq(0.05, 0.95, 0.3)
FDR <- as.data.frame(matrix(0, nrow = length(cutoff), ncol = 4))
colnames(FDR) <- c("cutoff", "randcounts", "truecounts", 
                   "FDR")
FDR$cutoff <- cutoff
per <- ncol(permutes)

# weighted permutations: make sure LRV = weighed version
for (k in 1:per) {
  shuffle <- permutes[, k]
  
  c = c()
  w = c()
  p1 <- c()
  p2 <- c()
  p <- c()
  lrv <- c()
  count.perm <- ct[shuffle, ]
  design <- matrix(0, nrow = nrow(count.perm), ncol = 2)
  design[group == unique(group)[1], 1] <- 1
  design[group == unique(group)[2], 2] <- 1
  v <- limma::voom(t(count.perm), design = design)
  W <- t(v$weights)
  # weighed values 
  for (i in 1:ncol(count.perm)){
  
    c.wij = W[pheno.coll$pheno == "C",i] * W[pheno.coll$pheno == "C",rpl11.idx]
    c[i] = wtv(log(count.perm[pheno.coll$pheno == "C", i] / count.perm[pheno.coll$pheno == "C", rpl11.idx]), c.wij)
    
    w.wij = W[pheno.coll$pheno == "W",i] * W[pheno.coll$pheno == "W",rpl11.idx]
    w[i] = wtv(log(count.perm[pheno.coll$pheno == "W", i] / count.perm[pheno.coll$pheno == "W", rpl11.idx]), w.wij)
    
    wij = W[which(pheno.coll$pheno %in% c("C","W")),i] * W[which(pheno.coll$pheno %in% c("C","W")),rpl11.idx]
    lrv[i] = wtv(log(count.perm[which(pheno.coll$pheno %in% c("C","W")), i] / count.perm[which(pheno.coll$pheno %in% c("C","W")), rpl11.idx]), wij)
    
    # calculating omega for weights 
    c.n = sum(c.wij)
    p1[i] = c.n - sum(c.wij**2)/ c.n 
    
    w.n = sum(w.wij)
    p2[i] = w.n - sum(w.wij**2)/ w.n
    
    n = sum(wij)
    p[i] = n - sum(wij**2)/ n
    

    
  }
  theta_e <- 1 - pmax(p1 * c, p2 * w)/(p * lrv)
  pkt <- theta_e
  for (cut in 1:nrow(FDR)) {
    FDR[cut, "randcounts"] <- FDR[cut, "randcounts"] + 
      sum(pkt < FDR[cut, "cutoff"], na.rm = T)
  }
}

FDR$randcounts <- FDR$randcounts/per
for (cut in 1:nrow(FDR)) {
  FDR[cut, "truecounts"] <- sum(ct.rpl11$theta_e < 
                                  FDR[cut, "cutoff"], na.rm = T)
  FDR[cut, "FDR"] <- FDR[cut, "randcounts"]/FDR[cut, "truecounts"]
}

# OLD METHOD - SEE propd_sampling.r: nonweighted permutations
cutoff = seq(0.05, 0.95, 0.3)
FDR <- as.data.frame(matrix(0, nrow = length(cutoff), ncol = 4))
colnames(FDR) <- c("cutoff", "randcounts", "truecounts", 
                   "FDR")
FDR$cutoff <- cutoff
per <- ncol(permutes)
p1 <- n1 - 1
p2 <- n2 - 1
p <- n1 + n2 - 1

for (k in 1:per) {
  shuffle <- permutes[, k]
  
  c = c()
  w = c()
  
  count.perm <- ct[shuffle, ]

  # weighed values 
  for (i in 1:ncol(count.perm)){
   
    c[i] = vlr(count.perm[pheno.coll$pheno == "C",i], 
               count.perm[pheno.coll$pheno == "C", rpl11.idx])
    w[i] = vlr(count.perm[pheno.coll$pheno == "W",i], 
               count.perm[pheno.coll$pheno == "W", rpl11.idx])
    
    
  }
  theta_e <- 1 - pmax(p1 * c, p2 * w)/(p * lrv.old)
  pkt <- theta_e
  for (cut in 1:nrow(FDR)) {
    FDR[cut, "randcounts"] <- FDR[cut, "randcounts"] + 
      sum(pkt < FDR[cut, "cutoff"], na.rm = T)
  }
}

FDR$randcounts <- FDR$randcounts/per
for (cut in 1:nrow(FDR)) {
  FDR[cut, "truecounts"] <- sum(ct.rpl11$theta_e < 
                                  FDR[cut, "cutoff"], na.rm = T)
  FDR[cut, "FDR"] <- FDR[cut, "randcounts"]/FDR[cut, "truecounts"]
}

# non-weight vs weighted
non.weight <- ct.rpl11.old[order(ct.rpl11.old$theta_e), "gene"]
weight <- ct.rpl11[order(ct.rpl11$theta_e), "gene"]

plot(rowMeans(countdata[which(rownames(countdata) %in% non.weight[1:100]),]),
rowMeans(countdata[which(rownames(countdata) %in% weight[1:100]),]))

median(countdata[which(rownames(countdata) %in% non.weight[1:100]),])
median(countdata[which(rownames(countdata) %in% weight[1:100]),])

  
plot(log(ct[, grep ("^RPL11$", colnames(ct)  )] / 
       ct[, grep ("^WRAP53$", colnames(ct)  )]),
     col = ifelse(group == "C", "red", "blue"))

# gene pair selection: FOR MULTIPLE PLOT GENERATION
gene.pair.list <- data.frame(gene = ct.rpl11.old[ct.rpl11.old$theta_e < .1 & !is.na(ct.rpl11.old$theta_e )
                                             , "gene"], rpl11 = "RPL11")

# plots gene relationships
plot.ratio <- function (gene.list){
  for( i in 1:nrow(gene.pair.list)){
    var.pair <- data.frame(ratio = ct[, grep (paste("^",gene.pair.list[i,1],"$", sep= ""), colnames(ct)  )] / 
                             ct[, grep (paste("^",gene.pair.list[i,2],"$", sep= ""), colnames(ct)  )], 
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

# PROPD pair plots
var.pair <- data.frame(ratio = ct[, grep (paste("^","RPL4","$", sep= ""), colnames(ct)  )] / 
                         ct[, grep (paste("^","RPL11","$", sep= ""), colnames(ct)  )], 
                       Genotype = pheno.coll$pheno)
var.pair <- var.pair %>% mutate(Genotype = 
                      case_when(Genotype == "C" ~ "c.396+3A", 
                                Genotype == "W" ~ "Non-Carr"))
ggplot(var.pair) + 
        geom_point(aes(x = seq(1:nrow(var.pair)),
                       y = ratio ,
                       color = Genotype,
                       fill = Genotype), size = 2,  shape = 21) +
        xlab("Sample Index") +
        ylab("Ratio") +
        ggtitle(paste("RPL4", "/",  "RPL11")) +
  theme_classic() + 
  scale_color_manual(values = c("red3", "black")) + 
  scale_fill_manual(values = c("lightsalmon", "gray55"))
   
# direction of proportions 
ct.rpl11.old %>% 
  mutate(direction = 
           case_when(ct.rpl11$c > ct.rpl11$w ~ "C", ct.rpl11$c < ct.rpl11$w ~ "W")) %>% 
  group_by(direction) %>% summarise(count = n())

# filter for top results
ct.rpl11.old %>% filter(theta_e < 0.15) %>%
  mutate(direction = 
           case_when(c > w ~ "C", c < w ~ "W")) %>% 
  group_by(direction) %>% summarise(count = n())



## DESeq data prep
# THIS MUST BE RUN FOR DESeq
pheno.data.alt <- pheno.data # make changes to pheno.data.alt

# OPTION 1: droping s from pheno.data
idx <- which(pheno.data.alt$pheno == "S")
pheno.data.alt$pheno[idx] <- "C"
pheno.data.alt$pheno <- factor(pheno.data.alt$pheno)

# OPTION 2: Sick individuals as carriers
# THIS METHOD IS WHAT IS USED FOR THE PAPER
pheno.data.alt$pheno[which(pheno.data.alt$pheno == "S")] <- "C"

# OPTIONAL: changing C to mutation
idx <- which(pheno.data.alt$pheno == "C")
pheno.data.alt$pheno <- as.character(pheno.data.alt$pheno)
pheno.data.alt$pheno[idx] <- "c.396+3A" #OPTION 1
pheno.data.alt$pheno[idx] <- "Carr" #OPTION 2

# OPTIONAL: wildtype to non-C 
idx <- which(pheno.data.alt$pheno == "W")
pheno.data.alt$pheno <- as.character(pheno.data.alt$pheno)
pheno.data.alt$pheno[idx] <- "Non-Carr"
pheno.data.alt$pheno <- factor(pheno.data.alt$pheno)

# OPTIONAL: S to DBA 
idx <- which(pheno.data.alt$pheno == "S")
pheno.data.alt$pheno <- as.character(pheno.data.alt$pheno)
pheno.data.alt$pheno[idx] <- "DBA"

## Step: DESeq2 analysis - collasping replicates - no batch effects
counts <- counts[,match(pheno.data$individual_number, colnames(counts))] # ordering counts with pheno
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = pheno.data.alt,
                              design = ~pheno) # experiment declaration 
counts(dds)

# PCA PLOTS: with VST transform 
vsd <- vst(dds)
# assay(vsd) gives the transformed counts data
pcaData <- DESeq2::plotPCA(vsd, intgroup = c("pheno","batch_number"))

# PCA: indiv labels
ggplot(pcaData$data)  + 
  geom_text(aes(PC1, PC2, color = pheno), label = pheno.data$individual_number) + # with replicate labels
  ylab(pcaData$labels$y) +
  xlab(pcaData$labels$x) + 
  scale_color_manual(name = "Pheno",values = c( "red3", "salmon1", "gray55"))
dev.copy(pdf,'Output/Figures/Fig3-PCA-individual.pdf') # saving options for figure
dev.off()

# PCA: indiv point version
ggplot(pcaData$data)  + 
  geom_point(aes(PC1, PC2, color = pheno, fill = pheno), shape = 21) + # with replicate labels
  ylab(pcaData$labels$y) +
  xlab(pcaData$labels$x) + 
  scale_color_manual(name = "delete_me",values = c( "red3", "red3", "Black")) + 
  scale_fill_manual(name = "Pheno",values = c( "salmon1", "red3", "gray55"))
dev.copy(pdf,'Output/Figures/PCA-individual-dot.pdf') # saving options for figure
dev.off()

# PCA: batch_number labels
ggplot(pcaData$data)  + 
  geom_point(aes(PC1, PC2, color = batch_number), size = 3) + # with batch shapes
  ylab(pcaData$labels$y) +
  xlab(pcaData$labels$x) + 
  scale_color_manual(values = c("#CC79A7", "#0072B2", "#F0E442", "#D55E00")) + theme_classic()
dev.copy(pdf,'Output/Figures/Fig4-PCA-batch-shape.pdf') # saving options for figure
dev.off()

# ploting PC3 and PC4 must be done manually
vv <- rowVars(assay(vsd))
ntop <- 500 
select <- order(vv, decreasing = T)[seq_len(min(ntop,length(vv)))] # 500 most var genes

pca <- prcomp(t(assay(vsd)[select,]))
PC1 <- pca$x[,1]
PC2 <- pca$x[,2]
PC3 <- pca$x[,3]
pcaData<- data.frame(pc1= PC2,pc2 =PC3)

#PCA variation 
pca.var <- pca$sdev^2
pca.var <- round(pca.var/sum(pca.var)*100, 1)

# PCA: plot PC3 and PC4 by pheno 
ggplot(pcaData) + geom_text(aes(pc1,pc2, color= pheno.data$pheno), 
                            label = rownames(pcaData))+
  ylab(paste("PC4:", pca.var[4], '% variance')) +
  xlab(paste("PC3", pca.var[3], '% variance')) + 
  scale_color_manual(values = c( "#d11f12", "#03AB11","#FF7700"))
dev.copy(pdf,'Output/Figures/PCA2-3-individual.pdf') # saving options for figure
dev.off()

# PCA: plot PC3 and PC4 by bacth_number
ggplot(pcaData) + geom_point(aes(pc1,pc2, color = pheno.data$batch_number), size = 3)+
  ylab(paste("PC4:", pca.var[4], '% variance')) +
  xlab(paste("PC3", pca.var[3], '% variance')) +
  scale_color_manual(values = c("#CC79A7", "#0072B2", "#F0E442", "#D55E00"))+ theme_classic()
dev.copy(pdf,'Output/Figures/PCA2-3-batch-shape.pdf') # saving options for figure
dev.off()

# differential expression: collapsing technical: HIGHLY RECOMMENDED
dds <- collapseReplicates(dds, groupby = dds$tech) # collapsing by technical 
collapsed.counts <- "counts"(dds)
dds <- DESeq(dds)

# dispersion plot
plotDispEsts(dds) 
dev.copy(pdf,'Output/Figures/dispersion.pdf') # saving options for figure
dev.off()

# results: group comparisons 
res <- results(dds, contrast = c("pheno", "W", "C")) #specify group comparisons here
summary(res) # summary of results for group comparison
resultsNames(dds)

# top results
top <- rownames(res)[which(res$padj < 0.01)] # view names of sig genes 
plotCounts(dds, "RPL11", "pheno") # plot normalized counts of specific genes by group
dev.copy(pdf,'Output/Figures/RPL11-pheno.pdf') # saving options for figure
dev.off()

# INDIVIDUAL GEN PLOTS: extracting data from results
deseq.data = data.frame(gene = rownames(res)[which(res$padj < 0.01)], 
                        logfold = res$log2FoldChange[which(res$padj < 0.01)],  
                        padj = res$padj[which(res$padj < 0.01)], 
                        basemean = res$baseMean[which(res$padj < 0.01)])


# MA PLOT: plots results of differential expression analysis - assumes DESeq has been ran
DESeq2::plotMA(res, ylim=c(-4,4))  # DESeq ma plot
text(600,3, label = "Non-carriers") # add text labels not figure
text(600,-3, label = "c.396+3A>G" , col = "dodgerblue4")
dev.copy(pdf,'Output/Figures/MA-plot.pdf') # saving options for figure
dev.off()

# skrinkage estimator
DESeq2::resultsNames(dds) # select coef for srinkage
# attempts to reduce low count noise in plots
resLFC <- DESeq2::lfcShrink(dds, 
                            coef = "pheno_W_vs_C", 
                            type = "ashr") 
DESeq2::plotMA(resLFC, ylim=c(-3,3))

# allocating srinkaged values for diff genes
SV.data <- data.frame(gene = rownames(resLFC)[which(resLFC$padj < .05)], 
       logfold = resLFC$log2FoldChange[which(resLFC$padj < .05)],  
       padj = resLFC$padj[which(resLFC$padj < .05)])

## Adding dream/voom/limma for mixed effects model of differential expression 
# Will use countdata -> technical replicates collapsed 
# Matching phenotype are in pheno.coll

# 1/3 of the samples is a good starting point
isexpr = rowSums(cpm (countdata) > 1 ) >= 1
geneExpr = DGEList( countdata[isexpr, ])
geneExpr = calcNormFactors( geneExpr, method = "TMM")

design_DE = as.data.frame ( pheno.coll[,c(3, 10)] ) #keep pheno and long individual name
design_DE$Individual = sapply ( strsplit(design_DE$individual_number , split = "_"), "[[" , 1 ) 
design_DE$StatusSubType = as.character( design_DE$pheno ) 
design_DE$pheno[11:14 ] = rep ( "C" , 4 ) # this is dependent on data used: change pheno to "C" if no sick or statussubtype to "S" if no sick
design_DE = design_DE [, -1]
row.names(design_DE) = colnames(countdata)

## We will apply dream to fit the mixed effects model 
form <- ~ pheno + (1|Individual)
vobjDream = voomWithDreamWeights(geneExpr, form, design_DE, plot = T)
fitmm = dream(vobjDream, form, design_DE)
fitmm$design
topTable ( fitmm, coef = "phenoW", number = 20)
limma.SCvN <- as.data.frame (topTable ( fitmm,coef = "phenoW", number = 10000) )


# We can alternatively fit with three levels
form <- ~ StatusSubType + (1|Individual)
vobjDream = voomWithDreamWeights(geneExpr, form, design_DE, plot = T)
fitmm = dream(vobjDream, form, design_DE)
fitmm$design  
limma <- as.data.frame (topTable ( fitmm, coef = "StatusSubTypeW", number = 10000) )
limma$gene_names <- rownames(limma)
colnames(limma) <- paste("dream", colnames(limma), sep = "_")
limma.both <- limma[which(rownames(limma) %in% deseq.data$gene), ]

# exporting limma data
merged.tables <- merge(deseq.data, limma.both, by.x = "gene", by.y = "dream_gene_names")
merged.tables <- merged.tables[-grep("^ENSG.+",merged.tables$gene), ] # removing non-annotated gene
write.table(merged.tables, file = "Output/data/dseq-dream.csv", sep = ",")

cdiff.data <- read_csv("Output/data/dseq-dream.csv") # reading in combined diff
plot.diff <- data.frame(gene = cdiff.data$gene[which(cdiff.data$padj < 0.01 & cdiff.data$dream_P.Value < 0.05)],
                        logfold = cdiff.data$logfold[which(cdiff.data$padj < 0.01 & cdiff.data$dream_P.Value < 0.05)],
                        basemean = cdiff.data$basemean[which(cdiff.data$padj < 0.01 & cdiff.data$dream_P.Value < 0.05)]) # enforcing cutoffs
diff.results <- data.frame(gene = rownames(res)[which(res$padj < 1)], 
                           logfold = res$log2FoldChange[which(res$padj < 1)],  
                           padj = res$padj[which(res$padj < 1)], 
                           basemean = res$baseMean[which(res$padj < 1)])

# simple MA plot 
ggplot(diff.results) +
  geom_point(aes(x=log(basemean+0.1), y=logfold), alpha =.75, color = "gray45")+
  geom_label(data = plot.diff, aes(x=log(basemean+0.1), y=logfold, label = gene), 
             color = "red", size = 1.5)+
  theme_classic()

ggplot(diff.results) +
  geom_point(aes(x=log(basemean+0.1), y=logfold), alpha =.75, color = "gray45")+
  geom_point(data = plot.diff, aes(x=log(basemean+0.1), y=logfold), 
             color = "red", size = 1.5)+
  theme_classic()

# GENE BOXPLOTS: expression of genes by phenotype

# data is prepped and separated by phenotype
norm<- counts(dds, normalized=TRUE) # PLEASE MAKE SURE DDS IS COLLASPED
gene.plot <- t(norm[ which(rownames(norm) %in% plot.diff$gene), ]) # must provide gene names here
carr.pheno <- pheno.data[which(pheno.data.alt$pheno == "C"), ]
carr.genes <- gene.plot[which( rownames(gene.plot) %in% carr.pheno$tech), ]
carr.genes <- as.data.frame(carr.genes)
carr.means <- as.data.frame(t(colMeans(carr.genes)))

wild.pheno <- pheno.data[which(pheno.data.alt$pheno == "W"), ]
wild.genes <- gene.plot[which( rownames(gene.plot) %in% wild.pheno$tech), ]
wild.genes <- as.data.frame(wild.genes)
wild.means <- as.data.frame(t(colMeans(wild.genes)))

# GENE BOXPLOTS: must prep data above
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


ggplot() + geom_point(data = wild.genes, aes(y = GATA2, x = "Non-carriers"), color = "#547294") + 
  geom_point(data = carr.genes, aes(y= GATA2, x = "c.396+3A>G"), color = "#d11f12") + 
  geom_bar(data = wild.means, aes(y = GATA2, "Non-carriers"), stat = "identity", color = "#547294", alpha = .2)+
  geom_bar(data = carr.means, aes(y = GATA2, "c.396+3A>G"), stat = "identity", color = "#d11f12", alpha = .2) +
  labs(x = "Genotype", y = "GATA2 Normalized")+
  theme_classic()
dev.copy(pdf,'Output/Figures/GATA2.pdf') # saving options for figure
dev.off()

ggplot() + geom_point(data = wild.genes, aes(y = ENSG00000285976  , x = "Non-carriers"), color = "#547294") + 
  geom_point(data = carr.genes, aes(y= ENSG00000285976  , x = "c.396+3A>G"), color = "#d11f12") + 
  geom_bar(data = wild.means, aes(y = ENSG00000285976  , "Non-carriers"), stat = "identity", color = "#547294", alpha = .2)+
  geom_bar(data = carr.means, aes(y = ENSG00000285976  , "c.396+3A>G"), stat = "identity", color = "#d11f12", alpha = .2) +
  labs(x = "Genotype", y = "ENSG00000285976 Normalized")+
  theme_classic()
dev.copy(pdf,'Output/Figures/ENSG00000285976.pdf') # saving options for figure
dev.off()



# FRY?ROAST collapse - gene set analysis
# preping the data
data <- collapse.rep(counts, pheno.data)
countdata <- data[[1]]
pheno.coll <- data[[2]]

dds.coll <- collapseReplicates(dds, groupby = dds$indiv)
collapsed.counts <- "counts"(dds.coll)

# subsetting RP data for fry
rp = grep ( "^RP[L|S]", row.names(countdata), perl =T ) 
rps6k = grep ( "RPS6K", row.names(countdata)) 
dash = grep ("-", row.names(countdata))
rp = setdiff(rp, rps6k) 
rp = setdiff(rp, dash)
row.names(countdata)[rp] 
grep ( "L1$", row.names(countdata)[rp] , perl = T )
# Will remove RPS4Y1, Y2 and X for proportionality
to_remove = c(16, 17, 25, 29, 50, 57, 82, 83, 84, 87, 88)
# Will keep all the RP-like proteins for now. 
rp = rp[-to_remove]
row.names(countdata)[rp]
isexpr = rowSums(cpm (countdata) > 1 ) >= 2
geneExpr = DGEList( countdata[isexpr, ])
geneExpr = calcNormFactors( geneExpr, method = "TMM")

#design for fry
design_DE = as.data.frame ( pheno.coll[,c(9, 10)] ) 
design_DE$Individual = design_DE$indiv 
design_DE$StatusSubType = as.character( design_DE$pheno ) 
design_DE = design_DE [, -1]
countdata <- countdata[,match(pheno.coll$indiv, colnames(countdata))]
design_DE$Individual == colnames(countdata)
row.names(design_DE) = colnames(countdata)

#FRY/ROAST 
#use rp_counts for propr analysis
design <- model.matrix(~0+design_DE$pheno)
colnames(design) <- levels(factor(design_DE$pheno))
y <- estimateDisp(geneExpr, robust = T, design = design)

cVw <- makeContrasts(C-W, levels = design)
cVw.fry <- fry(y, index= rownames(rp_counts), contrast = cVw, design=design)

cVs <- makeContrasts(C-S, levels = design)
cVs.fry <- fry(y, index= rownames(rp_counts), contrast = cVs, design=design)

sVw <- makeContrasts(S-W, levels = design)
sVw.fry <- fry(y, index= rownames(rp_counts), contrast = sVw, design=design)


## RIBOSOME PLOTS
# collapsing replicates without deseq
data <- collapse.rep(counts, pheno.data)
countdata <- data[[1]]
pheno.coll <- data[[2]]


# norm counts by pheno plot : all counts 
norm<- counts(dds, normalized=TRUE)
wild<- colSums(norm[,which(colnames(norm) %in% pheno.coll$tech[pheno.coll$pheno == "W"])])
carr<- colSums(norm[,which(colnames(norm) %in% pheno.coll$tech[pheno.coll$pheno == "C"])])
sick<- colSums(norm[,which(colnames(norm) %in% pheno.coll$tech[pheno.coll$pheno == "S"])])
sums <- data.frame(count =colSums(norm))
sums$pheno[which(rownames(sums) %in% pheno.coll$tech[pheno.coll$pheno == "W"])] <- "W"
sums$pheno[which(rownames(sums) %in% pheno.coll$tech[pheno.coll$pheno == "C"])] <- "C"
sums$pheno[which(rownames(sums) %in% pheno.coll$tech[pheno.coll$pheno == "S"])] <- "S"
sums <- sums %>% mutate(pheno = case_when(pheno == "W" ~ "Non-Carr", pheno == "C" ~ "Carr", pheno == "S" ~ "DBA"))
# keep s in distrubutions
sums <- sums %>% mutate(temp = case_when(pheno == "Non-Carr" ~ "W", pheno == "Carr" | pheno == 'DBA' ~ "M"))

# boxplot of distrubution of normalized counts 
ggplot(sums) + 
  geom_boxplot(outlier.colour = "red", outlier.shape = 1, aes(y = count)) +
  ylab("Normalized Counts") +
  theme(axis.text.x = element_blank()) + facet_grid(~pheno)
dev.copy(pdf,'Output/Figures/allcounts-boxplot.pdf') # saving options for figure
dev.off()

# dotplot of distrubution of normalized counts
dotplot <- ggplot(sums) +
  geom_point(aes(x = 1,y = count, fill = pheno, color = temp), size = 2, shape = 21, 
             position = position_dodge(width = .5))+
  ylab("Normalized Counts") +
  theme(axis.text.x = element_blank(), axis.title.x = element_blank())+
  scale_fill_manual(values = c( "salmon1", "red3" ,"gray55")) +
  scale_color_manual(values = c('red3','black'))

dotplot
dev.copy(pdf,'Output/Figures/allcounts-dotplot.pdf') # saving options for figure
dev.off()

# plotting ribosomal genes
norm <- counts(dds, normalized=TRUE)
norm <- norm[,match(pheno.coll$tech, colnames(norm))] # ordering counts with pheno

# subsetting Ribo genes
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

# PRIMATIVE: large subunit boxplots
par(mfrow=c(1,2))
boxplot(colSums(rpl_counts[, which(colnames(rpl_counts) %in% wild.idx)]), ylim = c(60000, 150000))
title("Noncarrier - large subunit")
boxplot(colSums(rpl_counts[, which(colnames(rpl_counts) %in% car.idx)]), ylim = c(60000, 150000))
title("c.396+3A>G - large subunit")
dev.copy(pdf,'Output/Figures/large-ribo-boxplot.pdf') # saving options for figure
dev.off()

# PRIMATIVE: RPL11 boxplots
par(mfrow=c(1,2))
boxplot(rpl_counts[which(rownames(rpl_counts) %in% c("RPL11")),
                   which(colnames(rpl_counts) %in% wild.idx)], ylim = c(900, 2750))
boxplot(rpl_counts[which(rownames(rpl_counts) %in% c("RPL11")),
                   which(colnames(rpl_counts) %in% car.idx)], ylim = c(900, 2750))


# PRIMATIVE: small subunit boxplots
par(mfrow=c(1,2))
boxplot(colSums(rps_counts[,which(colnames(rps_counts) %in% wild.idx)]), yim = c(40000, 80000))
title("Noncarrier - small subunit")
boxplot(colSums(rps_counts[,which(colnames(rps_counts) %in% car.idx)]), yim = c(40000, 80000))
title("c.396+3A>G- small subunit")
dev.copy(pdf,'Output/Figures/small-ribo-boxplot.pdf') # saving options for figure
dev.off()


# DATAPREP: long ribo gene datat
rpl_counts <- rpl_counts[,match(pheno.coll$tech, colnames(rpl_counts))] # ordering counts with pheno
rpl_long <- as.data.frame(rpl_counts)
colnames(rpl_long) <- pheno.coll$individual_number
rpl_long$means <- rowMeans(rpl_long)
rpl_long$gene <- rownames(rpl_long)
rpl_long <- rpl_long %>% pivot_longer(cols = 1:(ncol(rpl_long)-2), names_to = "sample", values_to = c("counts"))
# adding phenotype data to each row
rpl_long<- rpl_long %>% mutate(pheno = case_when(rpl_long$sample %in% pheno.coll$individual_number[which(pheno.coll$pheno == "C")] ~ "C",
                                                  rpl_long$sample %in% pheno.coll$individual_number[which(pheno.coll$pheno == "W")] ~ "W",
                                                  rpl_long$sample %in% pheno.coll$individual_number[which(pheno.coll$pheno == "S")] ~ "S"))
# WARNING: NORMALIZES DATA
rpl_long <- rpl_long %>% mutate(counts = counts/means)

# NOT OPTIONAL: name change
rpl_long <- rpl_long %>% mutate(pheno = case_when(rpl_long$pheno == "C" ~ "Carr",
                                     rpl_long$pheno == "S" ~ "S",
                                     rpl_long$pheno == "W" ~ "Non-Carr"))

# Adding Genotype data to each row
rpl_long <- rpl_long %>% mutate(geno = case_when(rpl_long$pheno == "Carr" | rpl_long$pheno == "S" ~ "M", 
                                                 rpl_long$pheno == "Non-Carr" ~ "W") )

# Ordering data for clearity
ordering.long <- rpl_long %>% group_by(gene) %>% 
  summarise(counts = median(counts)) %>% 
  mutate(gene = fct_reorder(gene, counts))
long.levels <- levels(ordering.long$gene)

rpl_long$gene <- factor(rpl_long$gene, levels=long.levels)

# large subunit ribosomal gene dotplot
rpl_long_plot <- rpl_long %>% ggplot(aes(x = gene, y = counts, fill = pheno, group = pheno, color = geno)) +
  #geom_bar(stat= "identity", position=position_dodge(), alpha = .1) + 
  geom_point(position = position_dodge(width = .5), shape = 21)+ 
  scale_fill_manual(values = c( "lightsalmon", "gray55","red3")) +
  scale_color_manual(values = c("M"="red3","W"='black',
                                "Carr"= "lightsalmon","Non-Carr"= "gray55","S"="red3")) +
  stat_summary(aes(color = pheno), fun.y = median, fun.ymin = median, fun.ymax = median,
               geom = "crossbar", width = 0.5)+
  theme_classic()
rpl_long_plot + theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.text=element_text(size=8))
dev.copy(pdf,'Output/Figures/large-ribo-dotplot.pdf') # saving options for figure
dev.off()


# DATAPREP: small ribo gene data
rps_counts <- rps_counts[,match(pheno.coll$tech, colnames(rps_counts))] # ordering counts with pheno
rps_long <- as.data.frame(rps_counts)
colnames(rps_long) <- pheno.coll$individual_number
rps_long$means <- rowMeans(rps_long)
rps_long$gene <- rownames(rps_long)
rps_long <- rps_long %>% pivot_longer(cols = 1:(ncol(rps_long)-2), names_to = "sample", values_to = c("counts"))
rps_long<- rps_long %>% mutate(pheno = case_when(rps_long$sample %in% pheno.coll$individual_number[which(pheno.coll$pheno == "C")] ~ "C",
                                                 rps_long$sample %in% pheno.coll$individual_number[which(pheno.coll$pheno == "W")] ~ "W",
                                                 rps_long$sample %in% pheno.coll$individual_number[which(pheno.coll$pheno == "S")] ~ "S"))
# WARNING: NORMALIZES DATA
rps_long <- rps_long %>% mutate(counts = counts/means)

#NOT OPTIONAL WITH GENOTYPE DATA
rps_long <- rps_long %>% mutate(pheno = case_when(rps_long$pheno == "C" ~ "Carr",
                                      rps_long$pheno == "S" ~ "S",
                                      rps_long$pheno == "W" ~ "Non-Carr"))

rps_long <- rps_long %>% mutate(geno = case_when(rps_long$pheno == "Carr" | rps_long$pheno == "S" ~ "M", 
                                     rps_long$pheno == "Non-Carr" ~ "W") ) 

#factors for sorting
ordering.short <- rps_long %>% group_by(gene) %>% 
  summarise(counts = median(counts)) %>% 
  mutate(gene = fct_reorder(gene, counts))
short.levels <- levels(ordering.short$gene)

rps_long$gene <- factor(rps_long$gene, levels=short.levels)

# small subunit ribosomal gene dotplot
rps_long %>% ggplot(aes(x = gene, y = counts, fill = pheno, group = pheno, color = geno)) +
  #geom_bar(stat= "identity", position=position_dodge(), alpha = .1) + 
  geom_point(shape = 21, position = position_dodge(width = .5))+ 
  theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.text=element_text(size=8))+
  ylab("counts")+
  scale_fill_manual(values = c( "lightsalmon", "gray55","red3")) +
  scale_color_manual(values = c("M"="red3","W"='black',
                                "Carr"= "lightsalmon","Non-Carr"= "gray55","S"="red3")) +
  stat_summary(aes(color = pheno), fun.y = median, fun.ymin = median, fun.ymax = median,
               geom = "crossbar", width = 0.5) +
  theme_classic()
dev.copy(pdf,'Output/Figures/small-ribo-dotplot.pdf') # saving options for figure
dev.off()

# MAKING sub and RPL11 boxplots
# drop sick 
rpl_long_plot <- rpl_long %>% 
  filter(gene != "RPL11") %>%
  group_by(sample, pheno) %>% 
  summarize(counts = sum(counts)) %>% 
  filter(pheno != "S") %>% mutate(type = "Large")
RPL11_plot <- rpl_long %>% 
  filter(gene == "RPL11") %>%
  filter(pheno != "S") %>% 
  group_by(sample, pheno) %>% 
  summarise(counts = sum(counts)) %>% mutate(type = "RPL11")
rps_long_plot <- rps_long %>% 
  group_by(sample, pheno) %>% 
  summarize(counts = sum(counts)) %>% 
  filter(pheno != "S") %>% mutate(type = "Small")

# ribo boxplots
ribo.plot <- rbind(rpl_long_plot, rps_long_plot, RPL11_plot)
ribo.plot <- ribo.plot %>% mutate(pheno = case_when(pheno == "Carr" ~ "c.396+3A>G", pheno == "Non-Carr" ~ "Noncarrier"))

#ploting ribo boxplots
ggplot(ribo.plot) + 
  geom_boxplot( aes(x = pheno, y = counts)) + 
  facet_wrap(~type, scales = "free") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme_classic()
dev.copy(pdf,'Output/Figures/ribo-boxplot.pdf') # saving options for figure
dev.off()

# individual gene boxplots
RPS19_plot <- rpl_long %>% 
  filter(gene == "RPL26") %>%
  filter(pheno != "S") %>% 
  group_by(sample, pheno) %>% 
  summarise(counts = sum(counts)) %>% mutate(type = "RPL26")
RPS19_plot <- RPS19_plot %>% mutate(pheno = case_when(pheno == "Carr" ~ "c.396+3A>G", pheno == "Non-Carr" ~ "Noncarrier"))

ggplot(RPS19_plot) + 
  geom_boxplot( aes(x = pheno, y = counts)) + 
  theme_classic() + ggtitle("RPL26")

# RPL11 dot plot
RPL11_dot <- rpl_long %>% 
  filter(gene == "RPL11") %>% 
  group_by(sample, pheno, geno) %>% 
  summarise(counts = sum(counts)) %>% mutate(type = "RPL11")

dotplot <- ggplot(RPL11_dot) +
  geom_point(aes(x = 1, 
                 y = counts, 
                 fill = pheno, 
                 color = geno, 
                 group = pheno), 
             size = 2, 
             shape = 21, 
             position = position_dodge(width = .5))+
  theme(axis.text.x = element_blank(), 
        axis.title.x = element_blank())+
  scale_fill_manual(values = 
                      c( "salmon1", "gray55", "red3")) +
  scale_color_manual(values = c("M"="red3","W"='black',
                                "Carr"= "lightsalmon","Non-Carr"= "gray55","S"="red3"))
dotplot
dev.copy(pdf,'Output/Figures/rpl11-dotplot.pdf') # saving options for figure
dev.off()

# RPL11 dot plot
RPS29_dot <- rps_long %>% 
  filter(gene == "RPS29") %>% 
  group_by(sample, pheno) %>% 
  summarise(counts = sum(counts)) %>% mutate(type = "RPL11")

dotplot <- ggplot(RPS29_dot) +
  geom_point(aes(x = 1, 
                 y = counts, 
                 fill = pheno,
                 group = pheno), 
             size = 2, 
             shape = 21, 
             position = position_dodge(width = .5))+
  theme(axis.text.x = element_blank(), 
        axis.title.x = element_blank())+
  scale_fill_manual(values = 
                      c( "salmon1", "gray55", "red3")) 
dotplot


#investigating outliner
ggplot(RPL11_dot) +
  geom_text(aes(x = 1, 
                 y = counts, 
                 color = geno, 
                 group = pheno,
                label = sample), 
             size = 2, 
             position = position_dodge(width = .5))+
  theme(axis.text.x = element_blank(), 
        axis.title.x = element_blank())+
  scale_color_manual(values = c("M"="red3","W"='black',
                                "Carr"= "lightsalmon",
                                "Non-Carr"= "gray55","S"="red3"))


rps_long <- as.data.frame(rps_counts)
colnames(rps_long) <- pheno.coll$individual_number
rps_long$means <- rowMeans(rps_long)
rps_long$gene <- rownames(rps_long)
rps_long <- rps_long %>% pivot_longer(cols = 1:(ncol(rps_long)-2), names_to = "sample", values_to = c("counts"))
rps_long<- rps_long %>% mutate(pheno = case_when(rps_long$sample %in% pheno.coll$individual_number[which(pheno.coll$pheno == "C")] ~ "C",
                                                 rps_long$sample %in% pheno.coll$individual_number[which(pheno.coll$pheno == "W")] ~ "W",
                                                 rps_long$sample %in% pheno.coll$individual_number[which(pheno.coll$pheno == "S")] ~ "S"))
rps_long <- rps_long %>% mutate(counts = counts/means)
rps_long <- rps_long %>% mutate(pheno = case_when(rps_long$pheno == "C" ~ "Carr",
                                                  rps_long$pheno == "S" ~ "S",
                                                  rps_long$pheno == "W" ~ "Non-Carr"))

rps_long <- rps_long %>% mutate(geno = case_when(rps_long$pheno == "Carr" | rps_long$pheno == "S" ~ "M", 
                                                 rps_long$pheno == "Non-Carr" ~ "W") )
# cdk11a plot 
cdk11a <- data.frame(count = norm[which(rownames(norm) %in% "CDK11A"),])
cdk11a$tech <- rownames(cdk11a)
cdk11a <- cdk11a[match(pheno.coll$tech, rownames(cdk11a)),] # ordering counts with pheno
cdk11a$indiv <- pheno.coll$individual_number
cdk11a$pheno <- pheno.coll$pheno
cdk11a <- cdk11a %>% 
  mutate(geno = case_when(cdk11a$pheno == "S" | cdk11a$pheno == "C" ~ "M", 
                     cdk11a$pheno == "W" ~ "W"))

ggplot(data = cdk11a) +
geom_text(aes(x = 1, 
              y = count, 
              group = pheno,
              label = indiv,
              color = geno), 
          size = 3, 
          position = position_dodge(width = .5))+
  theme(axis.text.x = element_blank(), 
        axis.title.x = element_blank())+
  scale_color_manual(values = c("M"="red3","W"='black',
                                "Carr"= "lightsalmon",
                                "Non-Carr"= "gray55","S"="red3")) +
  ggtitle("CDK11a")

# RPL11 plot 
rpl11 <- data.frame(count = norm[which(rownames(norm) %in% "RPL11"),])
rpl11$tech <- rownames(rpl11)
rpl11 <- rpl11[match(pheno.coll$tech, rownames(rpl11)),] # ordering counts with pheno
rpl11$indiv <- pheno.coll$individual_number
rpl11$pheno <- pheno.coll$pheno
rpl11 <- rpl11 %>% 
  mutate(geno = case_when(rpl11$pheno == "S" | rpl11$pheno == "C" ~ "M", 
                          rpl11$pheno == "W" ~ "W"))

ggplot(data = rpl11) +
  geom_text(aes(x = 1, 
                y = count, 
                group = pheno,
                label = indiv,
                color = geno), 
            size = 3, 
            position = position_dodge(width = .5))+
  theme(axis.text.x = element_blank(), 
        axis.title.x = element_blank())+
  scale_color_manual(values = c("M"="red3","W"='black',
                                "Carr"= "lightsalmon",
                                "Non-Carr"= "gray55","S"="red3")) +
  ggtitle("RPL11")

# indiv plots 
norm.plot <- norm[,match(pheno.coll$tech, colnames(norm))]
colnames(norm.plot) <- pheno.coll$individual_number
norm.plot <- data.frame(total = colSums(norm.plot), 
                        pheno = pheno.coll$pheno, 
                        indiv = pheno.coll$individual_number)
norm.plot <- norm.plot %>% mutate(geno = 
                       case_when(norm.plot$pheno == "C" | norm.plot$pheno == "S"~ "M",
                                 norm.plot$pheno == "W" ~ "W"))
ggplot(data = norm.plot)+
  geom_text(aes(x = 1, 
                           y = total, 
                           group = pheno,
                           label = indiv,
                           color = geno), 
                       size = 3, 
                       position = position_dodge(width = .5))+
  theme(axis.text.x = element_blank(), 
        axis.title.x = element_blank())+
  scale_color_manual(values = c("M"="red3","W"='black',
                                "Carr"= "lightsalmon",
                                "Non-Carr"= "gray55","S"="red3"))


# bfactors export
rps_counts <- rps_counts[,match(pheno.coll$tech, colnames(rps_counts))] # ordering counts with pheno
rps <- as.data.frame(rps_counts)
colnames(rps) <- pheno.coll$individual_number

rpl_counts <- rpl_counts[,match(pheno.coll$tech, colnames(rpl_counts))] # ordering counts with pheno
rpl <- as.data.frame(rpl_counts)
colnames(rpl) <- pheno.coll$individual_number
full <- rbind(rpl, rps)

full <- full[,match(pheno.coll$individual_number, colnames(full))]
wild.idx <- which(pheno.coll$pheno %in% "W") 
carr.idx <- which(pheno.coll$pheno %in% "C")

bfactors <- as.data.frame(rowMeans(full[,wild.idx,])/rowMeans(full[,carr.idx]))
write.table(bfactors, sep = ",", file = "bfactors.txt")

# nonskrined values for bfactors
deseq.data = data.frame( logfold = res$log2FoldChange[which(rownames(res) %in% rownames(full))], 
                         row.names = rownames(res)[which(rownames(res) %in% rownames(full))],
                         res$padj[which(rownames(res) %in% rownames(full))])

# here is the Skinage values for bfactor
SV.data <- data.frame(row.names = rownames(resLFC)[which(rownames(resLFC) %in% rownames(full))], 
                      logfold = resLFC$log2FoldChange[which(rownames(resLFC) %in% rownames(full))])
write.table(SV.data, sep = ",", file = "bfactors.txt")

library(plotrix)
f <- colorRamp(c("white", "blue"))
range01 <- function(x){(x-min(x))/(max(x)-min(x))}
colors <- f(range01(bfactors$`rowMeans(full[, carr.idx])/rowMeans(full[, wild.idx, ])`))
color.scale(range01(bfactors$`rowMeans(full[, carr.idx])/rowMeans(full[, wild.idx, ])`))

# pdf reporting
BiocManager::install("ReportingTools")
library("ReportingTools")
des2Report <- HTMLReport(shortName = 'RNAseq_analysis_with_DESeq2',
                         title = 'RNA-seq analysis of differential expression using DESeq2',
                         reportDirectory = "./report")
## This might take a while
publish(dds,des2Report, pvalueCutoff=0.05,
        annotation.db="biomaRt", factor = colData(dds)$pheno,
        reportDir="./reports")
finish(des2Report)


write.csv(SV.data,file = "ribo-values.csv")

#ploting CvS
norm <- counts(dds, normalized=TRUE)
norm <- norm[,match(pheno.coll$tech, colnames(norm))] # ordering counts with pheno
sick.idx <- which(pheno.coll$pheno == "S")
c.idx <- which(pheno.coll$pheno == "C")
w.idx <- which(pheno.coll$pheno == "W")
II.8.idx <- which(pheno.coll$indiv == "II.8")
III.3_s2.idx <- which(pheno.coll$individual_number == "III.3_s2")
III.3_s2 <- norm[,III.3_s2.idx]
II.8 <- norm[,II.8.idx]
sick <- rowMeans(norm[,sick.idx ])
carrier <- rowMeans(norm[,c.idx ])
wild <-rowMeans(norm[,w.idx ])
svc.plot <- log(data.frame(sick = sick, carrier = carrier, 
                           wildtype = wild, II.8 = II.8, III.3_s2 = III.3_s2) + 1)

# sig gene names
plot.diff <- cdiff.data$gene[which(cdiff.data$padj < 0.01 & cdiff.data$dream_P.Value < 0.05)]
sick.diff <- rowMeans(norm[which(rownames(norm) %in% plot.diff), sick.idx])
c.diff <- rowMeans(norm[which(rownames(norm) %in% plot.diff), c.idx])
w.diff <- rowMeans(norm[which(rownames(norm) %in% plot.diff), w.idx])
eight.diff <- norm[which(rownames(norm) %in% plot.diff), II.8.idx]
three.diff <- norm[which(rownames(norm) %in% plot.diff), III.3_s2.idx]
svc.diff <- log(data.frame(sick = sick.diff, carrier = c.diff, 
                           wildtype =w.diff, eight= eight.diff, three = three.diff) +1)

svc.plot %>% ggplot() + 
  geom_point(aes(x = sick, y = carrier), alpha = 0.5) + 
  geom_point(data = svc.diff, aes(x = sick, 
                                 carrier), 
            color = 'red3')+
  theme_classic()

svc.plot %>% ggplot() + 
  geom_point(aes(x = wildtype, y = carrier), alpha = 0.5) + 
  geom_point(data = svc.diff, aes(x = wildtype, 
                                 carrier), 
            color = 'red3')+
  theme_classic()

svc.plot %>% ggplot() + 
  geom_point(aes(x = wildtype, y = sick), alpha = 0.5) + 
  geom_point(data = svc.diff, aes(x = wildtype, 
                                  sick), 
             color = 'red3')+
  theme_classic()

svc.plot %>% ggplot() + 
  geom_point(aes(x = II.8, y = carrier), alpha = 0.5) + 
  geom_point(data = svc.diff, aes(x = eight, 
                                  carrier), 
             color = 'red3')+
  theme_classic()

svc.plot %>% ggplot() + 
  geom_point(aes(x =III.3_s2 , y = II.8), alpha = 0.5) + 
  geom_point(data = svc.diff, aes(x = three, 
                                  eight), 
             color = 'red3')+
  theme_classic()
