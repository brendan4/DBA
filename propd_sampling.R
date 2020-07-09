library("dplyr")
library("ggpubr") 
library("ggplot2")

# collapsing technical replicates and ordering pheno with count
data <- collapse.rep(counts, pheno.data)
countdata <- data[[1]]
pheno.coll <- data[[2]]

# collapsed stats
length(unique(pheno.coll$indiv[pheno.coll$pheno == "C"])) # unique C individuals
length(unique(pheno.coll$indiv[pheno.coll$pheno == "W"])) # unique WT individuals
pheno.coll$indiv[which(duplicated(pheno.coll$indiv))] # duplicated individauls

# perparing data
ct <- t(countdata) # we transpose for calculations
pheno.coll.all <- pheno.coll # saving for later use if needed

# OPTION 1: removing sick individuals
ct <- ct[-which(pheno.coll.all$pheno == "S"),] # optional drop sick
pheno.coll <- pheno.coll[-which(pheno.coll$pheno == "S"),] #optional drop sick

# OPTION 2: Sick indivduals as carriers
pheno.coll$pheno[which(pheno.coll$pheno == "S")] <- "C"

# checking individauls facts after converting S to C
length(unique(pheno.coll$indiv[pheno.coll$pheno == "C"]))
length(unique(pheno.coll$indiv[pheno.coll$pheno == "W"]))
pheno.coll$indiv[which(duplicated(pheno.coll$indiv))]
length(which(pheno.coll$pheno == "C"))
length(which(pheno.coll$pheno == "W"))

### TESTING OPTION 1
per <- 10 # set number of permutations
test <- permute.samples(ct, pheno.coll, per)

per <- 10 # set number of permutations
test <- permute.samples(ct, pheno.coll, per, weighted = T)

# direction of proportions 
dir <-test %>% # directions of all results
  mutate(direction = 
           case_when(c > w ~ "C", c < w ~ "W")) %>% 
  group_by(direction) %>% summarise(count = n())

test %>% filter(theta_e < 0.15) %>% # direction of top results
  mutate(direction = 
           case_when(c > w ~ "C", c < w ~ "W")) %>% 
  group_by(direction) %>% summarise(count = n())

# ploting 
# gene pair selection: FOR MULTIPLE PLOT GENERATION
gene.pair.list <- data.frame(gene = test[test$theta_e < .05 & # set treshhold here
                                                   !is.na(test$theta_e )
                                                 , "gene"], rpl11 = "RPL11")
plot.ratio(gene.list = gene.pair.list)

# testing direction functions
test2 <- gen.cutoff(ct,pheno.coll,
                    per = 50, 
                    directions = T)

# testing cutoff functions 
test3 <- gen.cutoff(ct,pheno.coll, 
                    per = 50, 
                    directions = F,
                    do.cutoff = T, 
                    prop.values = test,
                    cutoff = c(0.025, 0.05, 0.10, 0.15, 0.20, 0.30, 0.50, 0.60, 0.80, 0.95))

# generating gene list
gene.list <- test$gene[which(test$theta_e < 0.10)]
write(gene.list, file = "genedata.txt", sep = ",")

### TESTING OPTION 2
per <- 50 # set number of permutations
test <- permute.samples(ct, pheno.coll, per, sample.size = 14)

# generate cutoff for propd data
gen.cutoff <- function(ct, pheno.coll, per, 
                       directions = F, do.cutoff = T, 
                       cutoff = seq(0.05, 0.95, 0.3),
                       prop.values = NA, sample.size = 10){
  
  # randomly mixes samples and calculates directions
  if (directions){
    # building permutation table
    permutes.all <- as.data.frame(matrix(0, nrow = nrow(ct), 
                                     ncol = per))
    for (col in 1:ncol(permutes.all)) permutes.all[, col] <- sample(nrow(ct), replace = T)
    
    #directions table
    direction.data <- as.data.frame(matrix(0, nrow = 2, 
                                           ncol = per))
    rownames(direction.data) <- c("C", "W")
    
    for(k in 1:per){
      shuffle <- permutes.all[,k]
      c.idx <- which(pheno.coll$pheno == "C")
      w.idx <- which(pheno.coll$pheno == "W")
      
      # building permutation table
      shuffle2 <- sample(shuffle[w.idx], replace = T, size = sample.size) # random selection 
      
      per.count <- ct[c(shuffle2,shuffle[c.idx]),]
      per.pheno <- pheno.coll[c(shuffle2,shuffle[c.idx]),]
      per.pheno[1:sample.size, "pheno"] <- "W"
      per.pheno[sample.size:sample.size*2, "pheno"] <- "C"
      propd.data <- propd.calc(per.count, per.pheno)
        
      direction <- propd.data %>% 
          mutate(direction = 
                   case_when(c > w ~ "C", c < w ~ "W")) %>% 
          group_by(direction) %>% summarise(count = n()) %>% na.omit()
        
      direction.data["C",k] <- direction[which(direction$direction == "C"), "count"]
      print(direction.data["C",k])
      direction.data["W",k] <- direction[which(direction$direction == "W"), "count"]
        
     
    }
    return(direction.data)
  }
  # generates the cut off for propd via FDR
  if (do.cutoff){
    if(is.na(prop.values)){
      stop("Must provide percalculated propd values as prop.values")
    }
    # permutation building
    message("Alert: Fixing permutations to active random seed.")
    permutes.cut <- as.data.frame(matrix(0, nrow = nrow(ct), 
                                         ncol = per))
    for (col in 1:ncol(permutes.cut)) permutes.cut[, col] <- sample(1:nrow(ct))
    # Permutations
  
    FDR <- as.data.frame(matrix(0, nrow = length(cutoff), ncol = 4))
    colnames(FDR) <- c("cutoff", "randcounts", "truecounts", 
                       "FDR")
    FDR$cutoff <- cutoff

    for (k in 1:per) {
      
      shuffle <- permutes.cut[, k] 
      wt <- sample(w.idx, replace = T, size = sample.size) # random selection of 10 wildtype
      pheno.coll2 <- pheno.coll[c(1:sample.size, wt),] # collapse pheno data
      shuffle <- shuffle[c(1:sample.size, wt)] # rebuild shuffle
      
      # groups
      group <- as.character(pheno.coll2$pheno)
      group1 <- group == unique(group)[1]
      group2 <- group == unique(group)[2]
      n1 <- sum(group1)
      n2 <- sum(group2)
      
      # NOT FOR WEIGHTED propd 
      p1 <- n1 - 1
      p2 <- n2 - 1
      p <- n1 + n2 - 1
      
      c = c()
      w = c()
      lrv = c()
      
      count.perm <- ct[shuffle, ]
      
      # replacing the zeros
      if (any(as.matrix(count.perm) == 0)) {
        message("Alert: Replacing 0s with next smallest value.")
        zeros <- count.perm == 0
        count.perm[zeros] <- min(count.perm[!zeros])}
      
      # weighed values 
      for (i in 1:ncol(count.perm)){
        
        c[i] = vlr(count.perm[pheno.coll2$pheno == "C",i], 
                   count.perm[pheno.coll2$pheno == "C", rpl11.idx])
        w[i] = vlr(count.perm[pheno.coll2$pheno == "W",i], 
                   count.perm[pheno.coll2$pheno == "W", rpl11.idx])
        lrv[i] = vlr(count.perm[which(pheno.coll2$pheno %in% c("C","W")),i], 
                     count.perm[which(pheno.coll2$pheno %in% c("C","W")),rpl11.idx])
        
        
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
      FDR[cut, "truecounts"] <- sum(prop.values$theta_e < 
                                      FDR[cut, "cutoff"], na.rm = T)
      FDR[cut, "FDR"] <- FDR[cut, "randcounts"]/FDR[cut, "truecounts"]
    }
    return(FDR)
  }

}

# permutes WT data 
permute.samples <- function(ct, pheno.coll, per, weighted = F, sample.size= 10){
  c.idx <- which(pheno.coll$pheno == "C")
  w.idx <- which(pheno.coll$pheno == "W")
  if (weighted){
    # groups
    group <- as.character(pheno.coll$pheno)
    
    # obtaining weights : calculate with zeros
    design <- matrix(0, nrow = nrow(ct), ncol = 2)
    design[group == unique(group)[1], 1] <- 1
    design[group == unique(group)[2], 2] <- 1
    v <- limma::voom(t(ct), design = design)
    weights <- t(v$weights)
    W <- weights 
  }
  
  # building permutation table
  permutes <- as.data.frame(matrix(0, nrow = sum(pheno.coll$pheno == "C"), 
                                   ncol = per))
  for (col in 1:ncol(permutes)) permutes[, col] <- sample(w.idx, replace = T, size = sample.size) # random selection 
  
  holder <- rep(0, time = ncol(ct))
  avg.data <- data.frame(c = holder, w = holder, lrv = holder,
                         c.lrm = holder, w.lrm = holder,
                         theta = holder, theta_e = holder, theta_f = holder, 
                         theta_g = holder)
  # calculate permutation
  for(k in 1:per){
    shuffle <- permutes[,k]
    per.count <- ct[c(shuffle,c.idx),]
    per.pheno <- pheno.coll[c(shuffle,c.idx),]
    
    if (weighted){
      per.W <- W[c(shuffle,c.idx),]
      propd.data <- propd.calc(per.count, per.pheno, weighted, W = per.W)
    } else {
      propd.data <- propd.calc(per.count, per.pheno)
    }
    
    # suming data with pervious
    for( item in 1:ncol(avg.data)){
      name <- colnames(avg.data)[item]
      avg.data[,name] <- avg.data[,name] + propd.data[,name]
      
    }
  }
  avg.data <- avg.data/per
  avg.data$gene  <- propd.data$gene
  return(avg.data)
  
}


propd.calc <- function(ct, pheno.coll, weighted = F, W = NA, gene.fix = "RPL11"){
  if (weighted){
    
    # replacing the zeros
    if (any(as.matrix(ct) == 0)) {
      message("Alert: Replacing 0s with next smallest value.")
      zeros <- ct == 0
      ct[zeros] <- min(ct[!zeros])}
    
    # run for weighted
    p1 <- c()
    p2 <- c()
    p <- c()
    
    # forming data holders
    rpl11.idx <- which(colnames(ct) %in% gene.fix)
    c = c()
    w = c()
    gene = c()
    genes = c()
    lrv = c()
    c.lrm = c()
    w.lrm = c()
    
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
    
    # calculates theta values
    theta <- (p1 * c + p2 * w)/(p * lrv)
    theta_e <- 1 - pmax(p1 * c, p2 * w)/(p * lrv)
    theta_f <- pmax(p1 * c, p2 * w)/(p * lrv)
    theta_g <- pmin(p1 * c, p2 * w)/(p * lrv)
    ct.rpl11 <- data.frame(gene = gene, genes = genes, c= c, w = w, lrv = lrv,
                           c.lrm = c.lrm, w.lrm = w.lrm,
                           theta = theta, theta_e = theta_e, theta_f = theta_f, 
                           theta_g = theta_g )
    return(ct.rpl11)
    
    
  } else {
  # replacing the zeros
  if (any(as.matrix(ct) == 0)) {
    message("Alert: Replacing 0s with next smallest value.")
    zeros <- ct == 0
    ct[zeros] <- min(ct[!zeros])}
  
  # groups
  group <- as.character(pheno.coll$pheno)
  group1 <- group == unique(group)[1]
  group2 <- group == unique(group)[2]
  n1 <- sum(group1)
  n2 <- sum(group2)
  
  # NOT FOR WEIGHTED propd 
  p1 <- n1 - 1
  p2 <- n2 - 1
  p <- n1 + n2 - 1
  
  # forming data holders
  rpl11.idx <- which(colnames(ct) %in% gene.fix)
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
  
  theta <- (p1 * c + p2 * w)/(p * lrv)
  theta_e <- 1 - pmax(p1 * c, p2 * w)/(p * lrv)
  theta_f <- pmax(p1 * c, p2 * w)/(p * lrv)
  theta_g <- pmin(p1 * c, p2 * w)/(p * lrv)
  ct.rpl11 <- data.frame(gene = gene, genes = genes, c= c, w = w, lrv = lrv,
                         c.lrm = c.lrm, w.lrm = w.lrm,
                         theta = theta, theta_e = theta_e, theta_f = theta_f, 
                         theta_g = theta_g )
  return(ct.rpl11)
  }
}


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

# Var log ratio
vlr = function (x,y) { 
  return ( var ( log(x/y) ) ) 
}


# plots gene relationships
plot.ratio <- function (gene.list, plot.indiv = F){
  if (plot.indiv){
    ct <- t(ct)
    ct <- ct[,match(pheno.coll$individual_number, colnames(ct))] # ordering counts with pheno
    sp <- split(seq(along = pheno.coll$indiv), pheno.coll$indiv) # splitting and obtaining index of tech for each sample 
    ct <- sapply(sp, function(i) rowMeans(ct[ , i, drop = FALSE])) # summing the rows of count data based on tech status 
    idx <- sapply(sp, function(i) i[1]) # obtaining idx for each sample
    colnames(ct) <- pheno.coll$individual_number[idx] # replacing names
    pheno.coll <- pheno.coll[idx,] # collapsing pheno data
    pheno.coll <- pheno.coll [order(pheno.coll$pheno), ]
    ct<- ct[,match(pheno.coll$individual_number, colnames(ct))] # ordering counts with pheno
    ct <- t(ct)
  } 

  for( i in 1:nrow(gene.list)){
    var.pair <- data.frame(ratio = ct[, grep (paste("^",gene.pair.list[i,1],"$", sep= ""), colnames(ct)  )] / 
                             ct[, grep (paste("^",gene.pair.list[i,2],"$", sep= ""), colnames(ct)  )], 
                           Genotype = pheno.coll$pheno)
    var.pair <- var.pair %>% mutate(Genotype = 
                                      case_when(Genotype == "C" ~ "c.396+3A", 
                                                Genotype == "W" ~ "Non-Carr"))
    print(ggplot(var.pair) + 
            geom_point(aes(x = seq(1:nrow(var.pair)),
                           y = ratio,
                           color = Genotype,
                           fill = Genotype), 
                       size = 2,  
                       shape = 21) +
            xlab("Sample Index") +
            ylab("Ratio") +
            ggtitle(paste(gene.pair.list[i,1], "/",  gene.pair.list[i,2]))+
            theme_classic() + 
            scale_color_manual(values = c("red3", "black")) + 
            scale_fill_manual(values = c("lightsalmon", "gray55")))
  }
  
}

# sample permutations

individual.samples <- function(ct, pheno.coll, per= 10, fix.gene = "RPL11"){
  # creating data holder
  holder <- rep(0, time = ncol(ct))
  avg.data <- data.frame(c = holder, w = holder, lrv = holder,
                         c.lrm = holder, w.lrm = holder,
                         theta = holder, theta_e = holder, theta_f = holder, 
                         theta_g = holder)
  
  # setting idx by phenotype
  c.idx <- which(pheno.coll$pheno == "C")
  w.idx <- which(pheno.coll$pheno == "W")
  
  # finding unique samples
  c.indiv <- unique(pheno.coll$indiv[c.idx]) 
  w.indiv <- unique(pheno.coll$indiv[w.idx])
  
  # maximal possible sample size
  sample.size <- min(length(c.indiv), length(w.indiv))
  
  # forming permations table
  permutes <- as.data.frame(matrix(0, nrow = sample.size * 2, 
                                   ncol = per))
  
  for (col in 1:ncol(permutes)){
    per.w <- sample(w.indiv, size = sample.size, replace = F) # per select wildtype
    sampled.c <- c()
    sampled.w <- c()
    
    # selecting one sample per individual
    for( i in 1:sample.size){
      curr.c <- which(pheno.coll$indiv %in% c.indiv[i])
      curr.w <- which(pheno.coll$indiv %in% per.w[i])
      
      if (length(curr.c) == 1){
        sampled.c[i] <- curr.c
      } else {
        sampled.c[i] <- sample(curr.c, size = 1)
      }
      
      if (length(curr.w) == 1){
        sampled.w[i] <- curr.w
      } else {
        sampled.w[i] <- sample(curr.w, size = 1)
      }
      
    }
    
    # post-condition c sampled must = unqiue c
    if (sum(pheno.coll$indiv[sampled.c] != c.indiv) != 0){
      stop("unique groups don't match")
    }
    permutes[, col] <- c(sampled.c, sampled.w)
  }
  
  # calculate values for each permutation
  for( i in 1:per){
    shuffle <- permutes[,i]
    per.pheno <- pheno.coll[shuffle,]
    per.ct <- ct[shuffle,]
    if (sum(rownames(per.ct) != per.pheno$individual_number) != 0){
      stop("samples do not match between pheno and ct. Make sure ct and pheno ordered properly")
    }
    
    propd.data <- propd.calc(per.ct, per.pheno, gene.fix = fix.gene)
    
    # suming data with pervious
    for( item in 1:ncol(avg.data)){
      name <- colnames(avg.data)[item]
      avg.data[,name] <- avg.data[,name] + propd.data[,name]
    }
  }
  avg.data <- avg.data/per
  avg.data$gene  <- propd.data$gene
  return(avg.data)
}

indiv.results <- individual.samples(ct, pheno.coll, per = 50)


# direction of proportions 
dir <- indiv.results %>% # directions of all results
  mutate(direction = 
           case_when(c > w ~ "C", c < w ~ "W")) %>% 
  group_by(direction) %>% summarise(count = n())

top <- indiv.results %>% filter(theta_e < 0.15) %>% # direction of top results
  mutate(direction = 
           case_when(c > w ~ "C", c < w ~ "W")) %>% 
  group_by(direction) %>% summarise(count = n())

# forming bar plots for directionality
# binomial test for directions
dir.holder <- c()
dir.holder[1] <- dir$count[1]
dir.holder[2] <- dir$count[2]
sum(dir.holder)

binom.test(dir.holder[1], sum(dir.holder))
binom.test(dir.holder) # same as above 
binom.test(top$count, sum(top$count)) # for top results


stat.test <- data.frame(p = " p < 2.2e-16",
                        y.position = max(dir$count) + max(dir$count)/10, 
                        group1 = 1, group2 = 2)

dir <- na.omit(dir) %>% 
  mutate(Genotype = 
           case_when(direction == "C" ~ "c.396+3A", direction == "W" ~ "Noncarr"))

dir.p <- dir %>% ggplot() + geom_bar(aes(x = Genotype, y = count, color = Genotype), 
                                     stat = "identity", fill= "white")  +
  stat_pvalue_manual(
    data = stat.test, label = "p",
    xmin = "group1", xmax = "group2",
    y.position = "y.position") + 
  scale_color_manual(values = c("Red3", "Black")) + 
  xlab("Direction of Variance") + ylab("Number of Pairs")
dir.p + theme_classic()


# ploting top results 
if (nrow(top)==1){
  # creating w row
  top <- rbind(top, c("W", 0, "Noncarr"))
  top$count <- as.numeric(top$count)
  
}
stat.test <- data.frame(p = " p < 2.2e-16",
                        y.position = max(top$count) + max(top$count)/10, 
                        group1 = 1, group2 = 2)

top <- top %>% mutate(Genotype = 
           case_when(direction == "C" ~ "c.396+3A", direction == "W" ~ "Non-Carr"))

dir.p <- top %>% ggplot() + geom_bar(aes(x = Genotype, y = count, color = Genotype), 
                                     stat = "identity", fill= "white")  +
  stat_pvalue_manual(
    data = stat.test, label = "p",
    xmin = "group1", xmax = "group2",
    y.position = "y.position") + 
  scale_color_manual(values = c("Red3", "Black")) + 
  xlab("Direction of Variance") + ylab("Number of Pairs")
dir.p + theme_classic()


# update cuttoffs for individuals
# we must first form equal groups size for between C and W
# setting idx by phenotype
gen.cutoff.indiv <- function(ct, pheno.coll, per, cutoff, prop.values){
  c.idx <- which(pheno.coll$pheno == "C")
  w.idx <- which(pheno.coll$pheno == "W")
  
  # finding unique samples
  c.indiv <- unique(pheno.coll$indiv[c.idx]) 
  w.indiv <- unique(pheno.coll$indiv[w.idx])
  
  # maximal possible sample size
  sample.size <- min(length(c.indiv), length(w.indiv))
  
  # forming FDR table
  FDR <- as.data.frame(matrix(0, nrow = length(cutoff), ncol = 4))
  colnames(FDR) <- c("cutoff", "randcounts", "truecounts", 
                     "FDR")
  FDR$cutoff <- cutoff
  
  rpl11.idx <- which(colnames(ct) %in% "RPL11") # RPL11 index 
  
  for( i in 1:per){
    #must be repeat with every cycle
    per.w <- sample(w.indiv, size = sample.size, replace = F) # per select wildtype
    
    # once we have equal group size we can mix C and W, then calculate
    per.sample <- sample(c(c.indiv,per.w), size = sample.size*2)
    # c = 1:sample.size, w = sample.size: 2* sample.size
    c.idx.per <- 1:sample.size
    w.idx.per <- (1+sample.size):(2*sample.size)
    
    c.indiv.per <- per.sample[c.idx.per]
    w.indiv.per <- per.sample[w.idx.per]
    sampled.c <- c()
    sampled.w <- c()
    
    # selecting individuals
    for (i in 1:sample.size){
      curr.c <- which(pheno.coll$indiv %in% c.indiv.per[i])
      curr.w <- which(pheno.coll$indiv %in% w.indiv.per[i])
      
      if (length(curr.c) == 1){ # if 1 individual don't sample
        sampled.c[i] <- curr.c
      } else {
        sampled.c[i] <- sample(curr.c, size = 1)
      }
      
      if (length(curr.w) == 1){
        sampled.w[i] <- curr.w
      } else {
        sampled.w[i] <- sample(curr.w, size = 1)
      }
      
    }
    
    per.ct <- ct[c(sampled.c, sampled.w),]
    per.pheno <- pheno.coll[c(sampled.c, sampled.w),]
    per.pheno$pheno[1:sample.size] <- "C"
    per.pheno$pheno[(sample.size + 1):(sample.size*2)] <- "W"
    
    # groups
    group <- as.character(per.pheno$pheno)
    group1 <- group == unique(group)[1]
    group2 <- group == unique(group)[2]
    n1 <- sum(group1)
    n2 <- sum(group2)
    
    # NOT FOR WEIGHTED propd 
    p1 <- n1 - 1
    p2 <- n2 - 1
    p <- n1 + n2 - 1
    
    c = c()
    w = c()
    lrv = c()
    
    count.perm <- per.ct
    
    # replacing the zeros
    if (any(as.matrix(count.perm) == 0)) {
      zeros <- count.perm == 0
      count.perm[zeros] <- min(count.perm[!zeros])}
    
    # weighed values 
    for (i in 1:ncol(count.perm)){
      
      c[i] = vlr(count.perm[per.pheno$pheno == "C",i], 
                 count.perm[per.pheno$pheno == "C", rpl11.idx])
      w[i] = vlr(count.perm[per.pheno$pheno == "W",i], 
                 count.perm[per.pheno$pheno == "W", rpl11.idx])
      lrv[i] = vlr(count.perm[which(per.pheno$pheno %in% c("C","W")),i], 
                   count.perm[which(per.pheno$pheno %in% c("C","W")),rpl11.idx])
      
      
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
    FDR[cut, "truecounts"] <- sum(prop.values$theta_e < 
                                    FDR[cut, "cutoff"], na.rm = T)
    FDR[cut, "FDR"] <- FDR[cut, "randcounts"]/FDR[cut, "truecounts"]
  }
  return(FDR)
}

cutoff = c(0.025, 0.05, 0.10, 0.15, 0.20, 0.30, 0.50, 0.60, 0.80, 0.95) # gen cutoffs
prop.values <- indiv.results # must provide compare results
FDR.indiv <- gen.cutoff.indiv(ct, pheno.coll, 
                              per = 10, 
                              cutoff = cutoff, 
                              prop.values = prop.values)

# individuals: CvsW
# assumes sick individuals are removed: see OPTION 1

indiv.results.nos <- individual.samples(ct, pheno.coll, per = 100)

# ordering results for indiv.results
order.results <- indiv.results.nos[order(indiv.results.nos$theta_e),]

# direction of proportions 
dir <- indiv.results.nos %>% # directions of all results
  mutate(direction = 
           case_when(c > w ~ "C", c < w ~ "W")) %>% 
  group_by(direction) %>% summarise(count = n())

indiv.results.nos %>% filter(theta_e < 0.15) %>% # direction of top results
  mutate(direction = 
           case_when(c > w ~ "C", c < w ~ "W")) %>% 
  group_by(direction) %>% summarise(count = n())

# binomial test for directions
dir.holder <- c()
dir.holder[1] <- dir$count[1]
dir.holder[2] <- dir$count[2]
sum(dir.holder)

binom.test(dir.holder[1], sum(dir.holder))
binom.test(dir.holder) # same as above 

# preping data for plotting
stat.test <- data.frame(p = " p < 2.2e-16",
                        y.position = max(dir$count) + max(dir$count)/10, 
                        group1 = 1, group2 = 2)

dir <- na.omit(dir) %>% 
  mutate(phenotype = 
           case_when(direction == "C" ~ "Carr", direction == "W" ~ "Non-Carr"))

# plot for directionality
dir.p <- dir %>% ggplot() + geom_bar(aes(x = phenotype, y = count, color = phenotype), 
                                     stat = "identity", fill= "white")  +
  stat_pvalue_manual(
    data = stat.test, label = "p",
    xmin = "group1", xmax = "group2",
    y.position = "y.position") + 
  scale_color_manual(values = c("lightsalmon", "gray55")) + 
  xlab("Direction of Variance") + ylab("Number of Pairs")
dir.p + theme_classic()


# generating cutoff 
cutoff = c(0.025, 0.05, 0.10, 0.15, 0.20, 0.30, 0.50, 0.60, 0.80, 0.95) # gen cutoffs
FDR.indiv.nos <- gen.cutoff.indiv(ct, pheno.coll, 
                              per = 100, 
                              cutoff = cutoff, 
                              prop.values = indiv.results.nos)


# UNDER WORKS : new.method <- indiv.results
individual.samples<- function(ct, pheno.coll, per= 10){
  # creating data holder
  holder <- rep(0, time = ncol(ct))
  avg.data <- data.frame(c = holder, w = holder, lrv = holder,
                         c.lrm = holder, w.lrm = holder,
                         theta = holder, theta_e = holder, theta_f = holder, 
                         theta_g = holder)
  
  # setting idx by phenotype
  c.idx <- which(pheno.coll$pheno == "C")
  w.idx <- which(pheno.coll$pheno == "W")
  
  # finding unique samples
  c.indiv <- unique(pheno.coll$indiv[c.idx]) 
  w.indiv <- unique(pheno.coll$indiv[w.idx])
  
  # maximal possible sample size
  sample.size <- min(length(c.indiv), length(w.indiv))
  
  # forming permations table
  permutes <- as.data.frame(matrix(0, nrow = sample.size * 2, 
                                   ncol = per))
  
  for (col in 1:ncol(permutes)){
    per.w <- sample(w.indiv, size = sample.size, replace = F) # per select wildtype
    per.w <- sample(per.w, size = sample.size, replace = T)
    per.c <- sample(c.indiv, size = sample.size, replace = T)
    sampled.c <- c()
    sampled.w <- c()
    
    
    # selecting one sample per individual
    for( i in 1:sample.size){
      curr.c <- which(pheno.coll$indiv %in% per.c[i])
      curr.w <- which(pheno.coll$indiv %in% per.w[i])
      
      if (length(curr.c) == 1){
        sampled.c[i] <- curr.c
      } else {
        sampled.c[i] <- sample(curr.c, size = 1)
      }
      
      if (length(curr.w) == 1){
        sampled.w[i] <- curr.w
      } else {
        sampled.w[i] <- sample(curr.w, size = 1)
      }
      
    }
    
    # post-condition c sampled must = unqiue c
    if (sum(pheno.coll$indiv[sampled.c] %in% c.indiv) != sample.size){
      stop("unique groups don't match")
    }
    permutes[, col] <- c(sampled.c, sampled.w)
  }
  
  # calculate values for each permutation
  for( i in 1:per){
    shuffle <- permutes[,i]
    per.pheno <- pheno.coll[shuffle,]
    per.ct <- ct[shuffle,]
    if (sum(rownames(per.ct) != per.pheno$individual_number) != 0){
      stop("samples do not match between pheno and ct. Make sure ct and pheno ordered properly")
    }
    
    propd.data <- propd.calc(per.ct, per.pheno)
    
    # suming data with pervious
    for( item in 1:ncol(avg.data)){
      name <- colnames(avg.data)[item]
      avg.data[,name] <- avg.data[,name] + propd.data[,name]
    }
  }
  avg.data <- avg.data/per
  avg.data$gene  <- propd.data$gene
  return(avg.data)
}


# exporting data for funcassociate
clean.results <- indiv.results[-grep("^ENSG.+" ,indiv.results$gene ), ]
ordered <- clean.results[order(clean.results$theta_e), ]

write.table(clean.results$gene, file='entire.tsv', 
            quote=FALSE, sep='\t', col.names = F,
            row.names = F)
write.table(ordered$gene [which(ordered$theta_e < 0.15)], 
            file='filtered-ordered.tsv', quote=FALSE,
            sep='\t', col.names = F, row.names = F)
write.table(ordered$gene, file='ordered.tsv',
            quote=FALSE, sep='\t', col.names = F, 
            row.names = F)

# filtering results for direction of C
c.dir <- ordered %>% 
  mutate(direction = 
           case_when(c > w ~ "C", c < w ~ "W")) %>% 
  filter(direction == "C") 

write.table(c.dir$gene, file='c-dir.tsv',
            quote=FALSE, sep='\t', col.names = F, 
            row.names = F)

# exporting data for individual no sick 
clean.results <- indiv.results.nos[-grep("^ENSG.+" ,indiv.results.nos$gene ), ]
ordered <- clean.results[order(clean.results$theta_e), ]

write.table(ordered$gene [which(ordered$theta_e < 0.15)], 
            file='filtered-ordered-nos.tsv', quote=FALSE,
            sep='\t', col.names = F, row.names = F) # ordered and Filtered 
write.table(ordered$gene, file='ordered-nos.tsv',
            quote=FALSE, sep='\t', col.names = F, 
            row.names = F) #  just ordered data 

# filtering results for direction of C
c.dir.nos <- ordered %>% 
  mutate(direction = 
           case_when(c > w ~ "C", c < w ~ "W")) %>% 
  filter(direction == "C") 

write.table(c.dir.nos$gene, file='c-dir.nos.tsv',
            quote=FALSE, sep='\t', col.names = F, 
            row.names = F) # writing c dir for funcA

# direction of wildtype
w.dir.nos <- ordered %>% 
  mutate(direction = 
           case_when(c > w ~ "C", c < w ~ "W")) %>% 
  filter(direction == "W") 

write.table(w.dir.nos$gene, file='w-dir.nos.tsv',
            quote=FALSE, sep='\t', col.names = F, 
            row.names = F) # writing for funcA

# loading funcA data
func.anno <- read.table("func-anno.tsv")
rownames(func.anno)[which(func.anno$GO.0044429 > 0)]
gene.pair.list <- data.frame(gene =  rownames(func.anno)[which(func.anno$GO.0044429 > 0)], 
                             rpl11 = "RPL11")
gene.pair.list <- data.frame (gene = rownames(func.anno)[which(func.anno$GO.0044429 > 0 & 
                                              rownames(func.anno) %in% 
                                              rownames(func.anno[1:57,]))], rpl11 = "RPL11")
plot.ratio(gene.list = gene.pair.list) # ploting mito 

func.results <- read.table("funcassociate_results.tsv",col.names = 1)
colnames(func.results) <- unlist(func.results[1, ,drop= T])
func.results <- func.results[-1,]

#exporting mito theta e values
mito.prop <- indiv.results[which(indiv.results$gene %in% rownames(func.anno)[which(func.anno$GO.0044429 > 0 & rownames(func.anno) %in% rownames(func.anno[1:57,]))]),]
mito.prop <- mito.prop[,colnames(mito.prop) %in% c("gene", "theta_e")]
colnames(mito.prop)
write.csv(mito.prop, file = "mito-theta_e.csv")



#FRY?ROAST collapse
counts <- counts[,match(pheno.data$individual_number, colnames(counts))] # ordering counts with pheno
sp <- split(seq(along = pheno.data$indiv), pheno.data$indiv) # splitting and obtaining index of tech for each sample 
countdata <- sapply(sp, function(i) rowSums(counts[ , i, drop = FALSE])) # summing the rows of count data based on tech status 
idx <- sapply(sp, function(i) i[1]) # obtaining idx for each sample
colnames(countdata) <- pheno.data$indiv[idx] # replacing names
pheno.coll <- pheno.data[idx,] # collapsing pheno data
pheno.coll <- pheno.coll [order(pheno.coll$pheno), ]

dds.coll <- collapseReplicates(dds, groupby = dds$indiv)
collapsed.counts <- "counts"(dds.coll)
#mito data for fry
sub.idx <- which(rownames(countdata) %in% rownames(func.anno)[which(func.anno$GO.0044429 > 0 & rownames(func.anno) %in% rownames(func.anno[1:57,]))])



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

# change gen.set here
gene.set <- rownames(func.anno)[which(func.anno$GO.0044429 > 0 & rownames(func.anno) %in% rownames(func.anno[1:57,]))]

cVw <- makeContrasts(C-W, levels = design)
cVw.fry <- fry(y, index= gene.set, contrast = cVw, design=design)

cVs <- makeContrasts(C-S, levels = design)
cVs.fry <- fry(y, index= gene.set, contrast = cVs, design=design)

sVw <- makeContrasts(S-W, levels = design)
sVw.fry <- fry(y, index= gene.set, contrast = sVw, design=design)

#investigaing the number of zeros
plot(log(countdata[,colnames(countdata) %in% "III.6_s1"] + 1),log(countdata[,colnames(countdata) %in% "II.8_s1"] +1) )
sum(countdata[,colnames(countdata) %in% "III.6_s1"] == 0)
sum(countdata[,colnames(countdata) %in% "II.8_s1"] == 0)
sum.zero <- function(x){
  sum(x == 0)
}
apply(countdata, 2, sum.zero)

# droping II.8_s1
ct<- countdata[,-which(colnames(countdata) %in% "II.8_s1")]
ct <- t(ct)
pheno.coll <- pheno.coll[-which(pheno.coll$individual_number %in% "II.8_s1"),]
indiv.results.drop <- individual.samples(ct, pheno.coll, per = 100)

# plotting some results
gene.pair.list <- data.frame(gene =  indiv.results.drop$gene[order(indiv.results.drop$theta_e)][1:10], 
                             rpl11 = "RPL11")
plot.ratio(gene.list = gene.pair.list)


# CDK11A and propd
CDK11a <- individual.samples(ct, pheno.coll, per = 50, fix.gene = "CDK11A")

# direction of proportions 
dir <- CDK11a %>% # directions of all results
  mutate(direction = 
           case_when(c > w ~ "C", c < w ~ "W")) %>% 
  group_by(direction) %>% summarise(count = n())

top <- CDK11a %>% filter(theta_e < 0.25) %>% # direction of top results
  mutate(direction = 
           case_when(c > w ~ "C", c < w ~ "W")) %>% 
  group_by(direction) %>% summarise(count = n())

# forming bar plots for directionality
# binomial test for directions
dir.holder <- c()
dir.holder[1] <- dir$count[1]
dir.holder[2] <- dir$count[2]
sum(dir.holder)

binom.test(dir.holder[1], sum(dir.holder))
binom.data <- binom.test(dir.holder) # same as above 

stat.test <- data.frame(p = round(binom.data$p.value, 2),
                        y.position = max(dir$count) + max(dir$count)/10, 
                        group1 = 1, group2 = 2)

dir <- na.omit(dir) %>% 
  mutate(Genotype = 
           case_when(direction == "C" ~ "c.396+3A", direction == "W" ~ "Noncarr"))

dir.p <- dir %>% ggplot() + geom_bar(aes(x = Genotype, y = count, color = Genotype), 
                                     stat = "identity", fill= "white")  +
  stat_pvalue_manual(
    data = stat.test, label = "p",
    xmin = "group1", xmax = "group2",
    y.position = "y.position") + 
  scale_color_manual(values = c("Red3", "Black")) + 
  xlab("Direction of Variance") + ylab("Number of Pairs")
dir.p + theme_classic()


# ploting top results 
top.binom <- binom.test(top$count, sum(top$count)) # for top results
if (nrow(top)==1){
  # creating w row
  top <- rbind(top, c("W", 0, "Noncarr"))
  top$count <- as.numeric(top$count)
  
}
stat.test <- data.frame(p = round(top.binom$p.value, 2),
                        y.position = max(top$count) + max(top$count)/10, 
                        group1 = 1, group2 = 2)

top <- top %>% mutate(Genotype = 
                        case_when(direction == "C" ~ "c.396+3A", direction == "W" ~ "Non-Carr"))

dir.p <- top %>% ggplot() + geom_bar(aes(x = Genotype, y = count, color = Genotype), 
                                     stat = "identity", fill= "white")  +
  stat_pvalue_manual(
    data = stat.test, label = "p",
    xmin = "group1", xmax = "group2",
    y.position = "y.position") + 
  scale_color_manual(values = c("Red3", "Black")) + 
  xlab("Direction of Variance") + ylab("Number of Pairs")
dir.p + theme_classic()

# exporting data for funcassociate
clean.results <- CDK11a[-grep("^ENSG.+" ,CDK11a$gene ), ]
ordered <- clean.results[order(clean.results$theta_e), ]


write.table(ordered$gene [which(ordered$theta_e < 0.25)], 
            file='filtered-ordered-CDK11A.tsv', quote=FALSE,
            sep='\t', col.names = F, row.names = F)


length(unique (pheno.coll$indiv[which(pheno.coll$pheno=="W")]))


# diff prop visual aid
gene <- "RPL4"
w.var <- var(log(ct[pheno.coll$pheno == "W",which(colnames(ct) %in% gene)]/ct[pheno.coll$pheno == "W",which(colnames(ct) %in% "RPL11")]))
c.var <- var(log(ct[pheno.coll$pheno == "C",which(colnames(ct) %in% gene)]/ct[pheno.coll$pheno == "C",which(colnames(ct) %in% "RPL11")]))
t.var <- var(log(ct[,which(colnames(ct) %in% gene)]/ct[,which(colnames(ct) %in% "RPL11")]))
line.plot <- data.frame(variation = c(w.var,c.var,t.var),type = c("Noncarr", "Carr", "Total"), min.line = c(0,0,0))
ggplot(line.plot, aes(x=type, y = variation)) + 
  geom_linerange(aes(x=type, y=NULL, ymin=min.line, ymax=variation, color = type)) + 
  theme_classic() + 
  ggtitle(gene) + 
  scale_color_manual(values = c("lightsalmon", "gray55", "black"))


line.plot


# PROPD pair plots with variation lines
var.pair <- data.frame(ratio = ct[, grep (paste("^","FIBP","$", sep= ""), colnames(ct)  )] / 
                         ct[, grep (paste("^","RPL11","$", sep= ""), colnames(ct)  )], 
                       Genotype = pheno.coll$pheno)
var.pair <- var.pair %>% mutate(Genotype = 
                                  case_when(Genotype == "C" ~ "c.396+3A", 
                                            Genotype == "W" ~ "Non-Carr"))
min.line <- c(min(var.pair$ratio[which(var.pair$Genotype == "c.396+3A")]), min(var.pair$ratio[which(var.pair$Genotype == "Non-Carr")]) )
max.line <- c(max(var.pair$ratio[which(var.pair$Genotype == "c.396+3A")]), max(var.pair$ratio[which(var.pair$Genotype == "Non-Carr")]))
x <- c(length(which(var.pair$Genotype == "c.396+3A"))/2, length(var.pair$Genotype) - length(which(var.pair$Genotype == "Non-Carr"))/2)

lines <- data.frame(min.line, max.line, x, Genotype = c("c.396+3A","Non-Carr"))

ggplot(var.pair) + 
  geom_point(aes(x = seq(1:nrow(var.pair)),
                 y = ratio ,
                 color = Genotype,
                 fill = Genotype), size = 2,  shape = 21) +
  xlab("Sample Index") +
  ylab("Ratio") +
  ggtitle(paste("FIBP", "/",  "RPL11")) +
  theme_classic() + 
  geom_linerange(aes(x=x, y=NULL, ymin=min.line, ymax=max.line, color = Genotype), data=lines)+
  scale_color_manual(values = c("red3", "black")) + 
  scale_fill_manual(values = c("lightsalmon", "gray55")) 
