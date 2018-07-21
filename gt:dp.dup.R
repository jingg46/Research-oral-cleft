library(VariantAnnotation)

##target sequencing data
target <- readVcf("/users/jli/target.dup.vcf", "hg19")
target.gt <- geno(target)$GT
target.gt <- as.data.frame(target.gt, stringsAsFactors = F)
target.gt[target.gt == "."] <- NA


gmkf <- readVcf("/users/jli/gmkf.dup.vcf", "hg19")
gmkf.gt <- geno(gmkf)$GT
gmkf.gt <- as.data.frame(gmkf.gt, stringsAsFactors = F)

gmkf.gt[gmkf.gt == "0/."] <- NA
gmkf.gt[gmkf.gt == "1/."] <- NA
gmkf.gt[gmkf.gt == "."] <- NA

temp <- rowRanges(gmkf)
temp1 <- temp@ranges
gmkf.position <- temp1@start

##map rs name and position together for gmkf data
rs.name <- names(temp)
length(unique(rs.name))
gmkf.rs.position <- data.frame(rs.name, gmkf.position)

temp <- rowRanges(target)
temp1 <- temp@ranges
target.position <- temp1@start

###check for duplicate
length(unique(target.position))
length(unique(gmkf.position))


###get rid of ALL duplicated positions
gmkf.dup.inx <- c(which(duplicated(gmkf.position)), which(duplicated(gmkf.position))-1)
gmkf.dup.inx <- sort(gmkf.dup.inx)
gmkf.gt <- gmkf.gt[-gmkf.dup.inx,]
gmkf.position <- gmkf.position[-gmkf.dup.inx]

target.inx <- which(target.position %in% gmkf.position)
gmkf.inx <- which(gmkf.position %in% target.position)

gmkf.position <- gmkf.position[gmkf.inx]

gmkf.rs.position <- gmkf.rs.position[gmkf.rs.position$gmkf.position %in% gmkf.position,]
write.table(gmkf.rs.position$rs.name, "gmkf.dup.rs.name.txt",sep = " ", quote = FALSE)

gmkf.gt <- gmkf.gt[gmkf.inx,]
target.gt <- target.gt[target.inx,]

##get duplicate positions for snp plot
##gmkf.dup.position <- gmkf.position[which(gmkf.position %in% target.position)]
##gmkf.dup.position <- data.frame(cbind(rep("chr8", length(gmkf.dup.position)), gmkf.dup.position))
##write.table(gmkf.dup.position, "gmkf.dup.position.txt", quote = F, row.names = F, col.names = F, sep = " ")

combine.id <- read.csv("/users/jli/combine.id.csv", stringsAsFactors = F)

###change gmkf id into target id
gmkf.gt.id <- colnames(gmkf.gt)
inx = rep(NA, length(gmkf.gt.id))
for (i in 1:length(gmkf.gt.id)){
  inx[i] <- which(combine.id$gmkf.id == gmkf.gt.id[i])
}
colnames(gmkf.gt) <- combine.id$target.id[inx]

###order colnames in both data sets
gmkf.gt <- gmkf.gt[ ,order(names(gmkf.gt))]
target.gt <- target.gt[ ,order(names(target.gt))]

###check 
###by position

position.count <- rowSums(gmkf.gt!= target.gt, na.rm=T)

###check abnormal positions
list.pos <- sort(position.count, decreasing = T)[1:10]
##other than 213 (rs1119880, gmkf: 8:129320002_C/T, target: 8:129320002_C/A), all (error rate larget than 10%) are from the same positions


###check maf for abnormal positions

maf <- function(x){
  alt <- (sum(x == "0/1", na.rm = T) + sum(x == "1/1", na.rm = T) *2)/((402 - sum(is.na(x)))*2)
  ref <- (sum(x == "0/1", na.rm = T) + sum(x == "0/0", na.rm = T) *2)/((402 - sum(is.na(x)))*2)
  if(alt > ref){
   paste(c("ref", ref)) 
  }
  else{
    paste(c("alt", alt))
  }
}


for (i in 1:length(list.pos)){
  print(c(maf(target.gt[which(position.count == list.pos[i]),]), maf(gmkf.gt[which(position.count == list.pos[i]),])))
}###only 1 position has different allels as minor allel (8:129950945_G/A)

###by sample
sample.count <- colSums(gmkf.gt != target.gt, na.rm=T)

pdf("hist.dup.gt.by.position.pdf")
hist(position.count, main = "Histogram of mismatch in GT by position", breaks = 30, xlab = "Number of pairs with mismatch calls")
dev.off()
pdf("hist.dup.gt.by.sample.pdf")
hist(sample.count, main = "Histogram of mismatch in GT by sample", breaks = 30, xlab = "Number of mismatch calls per pair")
dev.off()


###check DP
target.dp <- geno(target)$DP
target.dp <- as.data.frame(target.dp, stringsAsFactors = F)
target.dp[target.dp == "."] <- NA

gmkf.dp <- geno(gmkf)$DP
gmkf.dp <- as.data.frame(gmkf.dp, stringsAsFactors = F)
gmkf.dp[gmkf.dp == "."] <- NA

temp <- rowRanges(gmkf)
temp1 <- temp@ranges
gmkf.position <- temp1@start

temp <- rowRanges(target)
temp1 <- temp@ranges
target.position <- temp1@start

gmkf.dup.inx <- c(which(duplicated(gmkf.position)), which(duplicated(gmkf.position))-1)
gmkf.dup.inx <- sort(gmkf.dup.inx)
gmkf.dp <- gmkf.dp[-gmkf.dup.inx,]
gmkf.position <- gmkf.position[-gmkf.dup.inx]

target.inx <- which(target.position %in% gmkf.position)
gmkf.inx <- which(gmkf.position %in% target.position)

gmkf.dp <- gmkf.dp[gmkf.inx,]
target.dp <- target.dp[target.inx,]

gmkf.dp.id <- colnames(gmkf.dp)
inx = rep(NA, length(gmkf.dp.id))
for (i in 1:length(gmkf.dp.id)){
  inx[i] <- which(combine.id$gmkf.id == gmkf.dp.id[i])
}
colnames(gmkf.dp) <- combine.id$target.id[inx]

###order colnames in both data sets
gmkf.dp <- gmkf.dp[ ,order(names(gmkf.dp))]
target.dp <- target.dp[ ,order(names(target.dp))]

###check 
###by position
position.count <- rowSums(gmkf.dp - target.dp, na.rm=T)

###by sample
sample.count <- rep(NA, ncol(gmkf.dp))
for (i in 1:ncol(gmkf.dp)){
  sample.count[i] <- mean((gmkf.dp[,i]-mean(gmkf.dp[,i], na.rm = T))/sd(gmkf.dp[,i], na.rm = T) 
                          - (target.dp[,i] - mean(target.dp[,i], na.rm = T))/sd(target.dp[,i], na.rm = T), na.rm=T)
}

pdf("hist.dup.dp.by.position.pdf")
hist(position.count, main = "Histogram of mismatch in DP by position")
dev.off()
pdf("hist.dup.dp.by.sample.pdf")
hist(sample.count, main = "Histogram of mismatch in DP by sample")
dev.off()


###check for gq
target.gq <- geno(target)$GQ
target.gq <- as.data.frame(target.gq, stringsAsFactors = F)
target.gq[target.gq == "."] <- NA

gmkf.gq <- geno(gmkf)$GQ
gmkf.gq <- as.data.frame(gmkf.gq, stringsAsFactors = F)
gmkf.gq[gmkf.gq == "."] <- NA

temp <- rowRanges(gmkf)
temp1 <- temp@ranges
gmkf.position <- temp1@start

temp <- rowRanges(target)
temp1 <- temp@ranges
target.position <- temp1@start

gmkf.dup.inx <- c(which(duplicated(gmkf.position)), which(duplicated(gmkf.position))-1)
gmkf.dup.inx <- sort(gmkf.dup.inx)
gmkf.gq <- gmkf.gq[-gmkf.dup.inx,]
gmkf.position <- gmkf.position[-gmkf.dup.inx]

target.inx <- which(target.position %in% gmkf.position)
gmkf.inx <- which(gmkf.position %in% target.position)

gmkf.gq <- gmkf.gq[gmkf.inx,]
target.gq <- target.gq[target.inx,]

gmkf.gq.id <- colnames(gmkf.gq)
inx = rep(NA, length(gmkf.gq.id))
for (i in 1:length(gmkf.gq.id)){
  inx[i] <- which(combine.id$gmkf.id == gmkf.gq.id[i])
}
colnames(gmkf.gq) <- combine.id$target.id[inx]

###order colnames in both data sets
gmkf.gq <- gmkf.gq[ ,order(names(gmkf.gq))]
target.gq <- target.gq[ ,order(names(target.gq))]

###check overall gq distribution for WGS sample-wise
gmkf.avg.gq <- colMeans(gmkf.gq, na.rm = T)
gmkf.avg.gq <- as.numeric(gmkf.avg.gq)
pdf("hist.gmkf.dup.gq.by.sample.pdf")
hist(gmkf.avg.gq, main = "Histogram of GQ distribution for WGS by sample", breaks = 30)
dev.off()

###check 
###by position
position.count <- rowSums(gmkf.gq - target.gq, na.rm=T)

###by sample
sample.count <- rep(NA, ncol(gmkf.gq))
for (i in 1:ncol(gmkf.gq)){
  sample.count[i] <- mean((gmkf.dp[,i]-mean(gmkf.dp[,i], na.rm = T))/sd(gmkf.dp[,i], na.rm = T) 
       - (target.dp[,i] - mean(target.dp[,i], na.rm = T))/sd(target.dp[,i], na.rm = T), na.rm=T)
}

##abnormal sample (< -0.002, > 0.002): "H_ME-DS10829_1-DS10829_1" "H_ME-DS10829_2-DS10829_2"
##3"H_ME-DS10895_1-DS10895_1" "H_ME-DS11323_3-DS11323_3" "H_ME-DS11337_2-DS11337_2"

pdf("hist.dup.gq.by.position.pdf")
hist(position.count, main = "Histogram of mismatch in GQ by position")
dev.off()
pdf("hist.dup.gq.by.sample.pdf")
hist(sample.count, main = "Histogram of mismatch in GQ by sample")
dev.off()


########gmkf data############
###check normalized dp by mismatch in GT across individual

inx <- mapply(function(x,y) which(x != y), gmkf.gt, target.gt)
avg.mismatch.dp <- vector(length = ncol(gmkf.dp))
avg.match.dp <- vector(length = ncol(gmkf.dp))

for (i in 1: ncol(gmkf.dp)){
  if (length(unlist(inx[i]) != 0)){
    a <- unlist(inx[i])
    avg.mismatch.dp[i] <- mean(gmkf.dp[a,i], na.rm = T)
    avg.match.dp[i] <- mean(gmkf.dp[-a,i], na.rm = T)
  }
  else{
    avg.mismatch.dp[i] <- NA
    avg.match.dp[i] <- mean(gmkf.dp[,i], na.rm = T)
  }
}

avg.dp <- data.frame(cbind(c(avg.match.dp, avg.mismatch.dp),  c(rep("match", length(avg.match.dp)), rep("mismatch", length(avg.mismatch.dp)))),stringsAsFactors=FALSE)
avg.dp[,1] <- as.numeric(avg.dp[,1])

##create summary statistics for avg.dp based on gt match
summary(avg.dp[avg.dp$X2 == "match",1])
sd(avg.dp[avg.dp$X2 == "match",1])
summary(avg.dp[avg.dp$X2 == "mismatch",1])
sd(avg.dp[avg.dp$X2 == "mismatch",1], na.rm = T)

library(ggplot2)
pdf("hist.avg.mismatch.dp.pdf")
ggplot(avg.dp, aes(x=X1, fill=X2)) + 
  geom_histogram(alpha=0.4, position="identity") +
  ggtitle("Histogram of average DP based on mismatch GT by sample") +
  xlab("Avg DP per individual") +
  scale_fill_discrete(name = "GT match type")
dev.off()

###check normalized gq by mismatch in GT across individual

avg.mismatch.gq <- vector(length = ncol(gmkf.gq))
avg.match.gq <- vector(length = ncol(gmkf.gq))

for (i in 1: ncol(gmkf.gq)){
  if (length(unlist(inx[i]) != 0)){
    a <- unlist(inx[i])
    avg.mismatch.gq[i] <- mean(gmkf.gq[a,i], na.rm = T)
    avg.match.gq[i] <- mean(gmkf.gq[-a,i], na.rm = T)
  }
  else{
    avg.mismatch.gq[i] <- NA
    avg.match.gq[i] <- mean(gmkf.gq[,i], na.rm = T)
  }
}

avg.gq <- data.frame(cbind(c(avg.match.gq, avg.mismatch.gq),  c(rep("match", length(avg.match.gq)), rep("mismatch", length(avg.mismatch.gq)))),stringsAsFactors=FALSE)
avg.gq[,1] <- as.numeric(avg.gq[,1])

##create summary statistics for avg.gq based on gt match
summary(avg.gq[avg.gq$X2 == "match",1])
sd(avg.gq[avg.gq$X2 == "match",1])
summary(avg.gq[avg.gq$X2 == "mismatch",1])
sd(avg.gq[avg.gq$X2 == "mismatch",1], na.rm = T)

library(ggplot2)
pdf("hist.avg.mismatch.gq.pdf")
ggplot(avg.gq, aes(x=X1, fill=X2)) + 
  geom_histogram(alpha=0.4, position="identity") +
  ggtitle("Histogram of average GQ based on mismatch GT by sample") +
  xlab("Avg GQ per individual") +
  scale_fill_discrete(name = "GT match type")
dev.off()

###########targeted data############
###check normalized dp by mismatch in GT across individual

inx <- mapply(function(x,y) which(x != y), target.gt, gmkf.gt)
avg.mismatch.dp <- vector(length = ncol(target.dp))
avg.match.dp <- vector(length = ncol(target.dp))

for (i in 1: ncol(target.dp)){
  if (length(unlist(inx[i]) != 0)){
    a <- unlist(inx[i])
    avg.mismatch.dp[i] <- mean(target.dp[a,i], na.rm = T)
    avg.match.dp[i] <- mean(target.dp[-a,i], na.rm = T)
  }
  else{
    avg.mismatch.dp[i] <- NA
    avg.match.dp[i] <- mean(target.dp[,i], na.rm = T)
  }
}

avg.dp <- data.frame(cbind(c(avg.match.dp, avg.mismatch.dp),  c(rep("match", length(avg.match.dp)), rep("mismatch", length(avg.mismatch.dp)))),stringsAsFactors=FALSE)
avg.dp[,1] <- as.numeric(avg.dp[,1])

##create summary statistics for avg.dp based on gt match
summary(avg.dp[avg.dp$X2 == "match",1])
sd(avg.dp[avg.dp$X2 == "match",1])
summary(avg.dp[avg.dp$X2 == "mismatch",1])
sd(avg.dp[avg.dp$X2 == "mismatch",1], na.rm = T)

library(ggplot2)
pdf("hist.avg.mismatch.dp.target.pdf")
ggplot(avg.dp, aes(x=X1, fill=X2)) + 
  geom_histogram(alpha=0.4, position="identity") +
  ggtitle("Histogram of average DP based on mismatch GT by sample (Targeted Data)") +
  xlab("Avg DP per individual") +
  scale_fill_discrete(name = "GT match type")
dev.off()

###check normalized gq by mismatch in GT across individual

avg.mismatch.gq <- vector(length = ncol(target.gq))
avg.match.gq <- vector(length = ncol(target.gq))

for (i in 1: ncol(target.gq)){
  if (length(unlist(inx[i]) != 0)){
    a <- unlist(inx[i])
    avg.mismatch.gq[i] <- mean(target.gq[a,i], na.rm = T)
    avg.match.gq[i] <- mean(target.gq[-a,i], na.rm = T)
  }
  else{
    avg.mismatch.gq[i] <- NA
    avg.match.gq[i] <- mean(target.gq[,i], na.rm = T)
  }
}

avg.gq <- data.frame(cbind(c(avg.match.gq, avg.mismatch.gq),  c(rep("match", length(avg.match.gq)), rep("mismatch", length(avg.mismatch.gq)))),stringsAsFactors=FALSE)
avg.gq[,1] <- as.numeric(avg.gq[,1])

##create summary statistics for avg.gq based on gt match
summary(avg.gq[avg.gq$X2 == "match",1])
sd(avg.gq[avg.gq$X2 == "match",1])
summary(avg.gq[avg.gq$X2 == "mismatch",1])
sd(avg.gq[avg.gq$X2 == "mismatch",1], na.rm = T)

library(ggplot2)
pdf("hist.avg.mismatch.gq.target.pdf")
ggplot(avg.gq, aes(x=X1, fill=X2)) + 
  geom_histogram(alpha=0.4, position="identity") +
  ggtitle("Histogram of average GQ based on mismatch GT by sample (Targeted Data)") +
  xlab("Avg GQ per individual") +
  scale_fill_discrete(name = "GT match type")
dev.off()



###look at individual level, pick some random samples, compare DP based on match in GT
sample.inx <- sample(ncol(gmkf.gq), 5, replace = F)

i=1
a <- unlist(inx[sample.inx[i]])
individual.gq <- data.frame(gmkf.gq[,sample.inx[i]])
individual.gq$mark <- "match"
individual.gq$mark[a] <- "mismatch"
colnames(individual.gq)[1] <- "dp"

pdf("sample.gq.comparison.pdf")
ggplot(individual.gq, aes(x=dp, fill=mark)) + 
  geom_histogram(alpha=0.4, position="identity") +
  ggtitle("Histogram of sample-level GQ based on mismatch GT by sample 1") +
  xlab("GQ for this sample across all positions") +
  scale_fill_discrete(name = "GT match type")
dev.off()

i=2
a <- unlist(inx[sample.inx[i]])
individual.gq <- data.frame(gmkf.gq[,sample.inx[i]])
individual.gq$mark <- "match"
individual.gq$mark[a] <- "mismatch"
colnames(individual.gq)[1] <- "dp"

pdf("sample.gq.comparison2.pdf")
ggplot(individual.gq, aes(x=dp, fill=mark)) + 
  geom_histogram(alpha=0.4, position="identity") +
  ggtitle("Histogram of sample-level GQ based on mismatch GT by sample 2") +
  xlab("GQ for this sample across all positions") +
  scale_fill_discrete(name = "GT match type")
dev.off()

i=3
a <- unlist(inx[sample.inx[i]])
individual.gq <- data.frame(gmkf.gq[,sample.inx[i]])
individual.gq$mark <- "match"
individual.gq$mark[a] <- "mismatch"
colnames(individual.gq)[1] <- "dp"

pdf("sample.gq.comparison3.pdf")
ggplot(individual.gq, aes(x=dp, fill=mark)) + 
  geom_histogram(alpha=0.4, position="identity") +
  ggtitle("Histogram of sample-level GQ based on mismatch GT by sample 3") +
  xlab("GQ for this sample across all positions") +
  scale_fill_discrete(name = "GT match type")
dev.off()

i=4
a <- unlist(inx[sample.inx[i]])
individual.gq <- data.frame(gmkf.gq[,sample.inx[i]])
individual.gq$mark <- "match"
individual.gq$mark[a] <- "mismatch"
colnames(individual.gq)[1] <- "dp"

pdf("sample.gq.comparison4.pdf")
ggplot(individual.gq, aes(x=dp, fill=mark)) + 
  geom_histogram(alpha=0.4, position="identity") +
  ggtitle("Histogram of sample-level GQ based on mismatch GT by sample 4") +
  xlab("GQ for this sample across all positions") +
  scale_fill_discrete(name = "GT match type")
dev.off()

i=5
a <- unlist(inx[sample.inx[i]])
individual.gq <- data.frame(gmkf.gq[,sample.inx[i]])
individual.gq$mark <- "match"
individual.gq$mark[a] <- "mismatch"
colnames(individual.gq)[1] <- "dp"

pdf("sample.gq.comparison5.pdf")
ggplot(individual.gq, aes(x=dp, fill=mark)) + 
  geom_histogram(alpha=0.4, position="identity") +
  ggtitle("Histogram of sample-level GQ based on mismatch GT by sample 5") +
  xlab("GQ for this sample across all positions") +
  scale_fill_discrete(name = "GT match type")
dev.off()







