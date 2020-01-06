library(biomaRt)
library(Rsubread)

# listing files for first histat run 
main <- "/media/brendan/Elements/dba/DBA_121317" # root folder that contains bam files for early histat runs 
setwd(main)
all.folders <- list.files() # lists all folders in that are bams
bam.files <- c() # holder of bam files 
for(i in 1:length(all.folders)){ # cycles through folder and adds bam file 
  setwd(paste(main, all.folders[i],sep = "/"))
  bam.files <- c(bam.files, paste(main, all.folders[i],list.files(pattern = ".bam$"), sep = "/"))
  setwd(main)
}

# listing file for all other runs 
setwd("..") # up one level 
main2 <-"/media/brendan/Elements/dba/main_batch" # root folder that contains bam files for other runs 
setwd(main2)
all.folders <- list.files()
for(i in 1:length(all.folders)){
  setwd(paste(main2, all.folders[i],sep = "/"))
  bam.files <- c(bam.files, paste(main2, all.folders[i],list.files(pattern = ".bam$"), sep = "/"))
  setwd(main2)
}

setwd("..") # up one level 
setwd("annotation") # opens annotation folder
anno <- list.files(pattern = ".gz$", full.names = TRUE) # sets annotation location 
fc <- featureCounts(bam.files, annot.ext= anno, isGTFAnnotationFile = T,isPairedEnd=TRUE) # conducts counting
counts <- fc$counts
stat <- fc$stat
write.table(stat, "Input_data/stat.tab", sep = "\t")

#gene name conversions
ids <- rownames(counts) # gene names from counts
ids <- gsub("\\.\\d*$", "", genes) # remove version tag
ensembl_database <- useEnsembl(biomart="ensembl", 
                      dataset="hsapiens_gene_ensembl") # open ensembl 
results <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"),
                 filter = "ensembl_gene_id", 
                 values = ids, mart = ensembl_database) # use genes to obtain hgnc_symbol

#for each entry in results find corresponding tag
for(i in 1:nrow(results)){
  idx <- grep(paste(results$ensembl_gene_id[i], "\\.\\d*$", sep= ""), rownames(counts)) # find corresponding tag
  if (nchar(results$hgnc_symbol[i]) == 0){ # if gene result is bank replace with id
    rownames(counts)[idx] <- results$ensembl_gene_id[i]}
  else if (results$hgnc_symbol[i] %in% rownames(counts)){ # if symbol already used replaced with id
    rownames(counts)[idx] <- results$ensembl_gene_id[i] 
  }else{ # else replace with gene symbol
    rownames(counts)[idx] <- results$hgnc_symbol[i]
  }
}
write.table(counts, "unprocessed_counts.tab", sep = "\t")

