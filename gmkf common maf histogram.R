library(VariantAnnotation)

##target sequencing data
target <- readVcf("/users/jli/biallelic.target.vcf", "hg19")
#Warning message:
#In .bcfHeaderAsSimpleList(header) :
#duplicate keys in header will be forced to unique rownames
target.gt <- geno(target)$GT
target.gt <- as.data.frame(target.gt, stringsAsFactors = F)
target.gt[target.gt == "."] <- NA
##all gt: 0/0, 0/1, 1/1, .
maf.t <- function(x){
  alt <- (sum(x == "0/1", na.rm = T) + sum(x == "1/1", na.rm = T) *2)/((722 - sum(is.na(x)))*2)
  ref <- (sum(x == "0/1", na.rm = T) + sum(x == "0/0", na.rm = T) *2)/((722 - sum(is.na(x)))*2)
  return(min(c(alt, ref)))
}

target.maf <- apply(target.gt, 1, maf.t)

##target rare vs. common by maf
temp <- rowRanges(target)
temp1 <- temp@ranges
target.position <- temp1@start
target.position <- target.position[which(!is.na(target.maf))]
target.maf <- target.maf[!is.na(target.maf)]
target.rare.position <- target.position[which(target.maf < 0.01)]
target.common.position <- target.position[which(target.maf >= 0.01)]


##gmkf sequencing data
gmkf <- readVcf("/dcl01/beaty/data/gmkf/euro/vcfs/filtered/8q24.recode.vcf", "hg19")
gmkf.gt <- geno(gmkf)$GT
gmkf.gt <- as.data.frame(gmkf.gt, stringsAsFactors = F)
#sum(gmkf.gt == "0/.")  612395
#sum(gmkf.gt == "1/.")  2394

gmkf.gt[gmkf.gt == "0/."] <- NA
gmkf.gt[gmkf.gt == "1/."] <- NA
gmkf.gt[gmkf.gt == "."] <- NA

maf.g <- function(x){
  alt <- (sum(x == "0/1", na.rm = T) + sum(x == "1/1", na.rm = T) *2)/((996 - sum(is.na(x)))*2)
  ref <- (sum(x == "0/1", na.rm = T) + sum(x == "0/0", na.rm = T) *2)/((996 - sum(is.na(x)))*2)
  return(min(c(alt, ref)))
}

gmkf.maf <- apply(gmkf.gt, 1, maf.g)

##target rare vs. common by maf
temp <- rowRanges(gmkf)
temp1 <- temp@ranges
gmkf.position <- temp1@start
gmkf.position <- gmkf.position[which(!is.na(gmkf.maf))]
gmkf.maf <- gmkf.maf[!is.na(gmkf.maf)]
gmkf.rare.position <- gmkf.position[which(gmkf.maf < 0.01)]
gmkf.common.position <- gmkf.position[which(gmkf.maf >= 0.01)]

##non-overlapping maf of gmkf common
gmkf.common.maf <- gmkf.maf[gmkf.maf >= 0.01]
nonoverlap.gmkf.common.maf <- gmkf.common.maf[which(!gmkf.common.position %in% target.position)]
pdf("non-overlapping.gmkf.common.maf.pdf")
hist(nonoverlap.gmkf.common.maf, main = "Histogram of maf of gmkf common variants (non-overlappkng regions)", breaks = 30)










