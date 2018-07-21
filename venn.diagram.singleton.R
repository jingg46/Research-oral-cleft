library(VariantAnnotation)

##target sequencing data
target <- readVcf("/users/jli/target.parent.vcf", "hg19")
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
gmkf <- readVcf("/users/jli/gmkf.parent.vcf", "hg19")
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

##venn diagram
library("VennDiagram")
pdf("target.rare.venn_diagram.pdf")
venn.plot <- venn.diagram(list(target.rare.position, gmkf.position), NULL, fill=c("red", "green"), 
                          alpha=c(0.5,0.5), cex = 2, cat.fontface=4, category.names=c("TS (rare)", "Gmkf"))
grid.draw(venn.plot)
dev.off()

pdf("target.common.venn_diagram.pdf")
venn.plot <- venn.diagram(list(target.common.position, gmkf.position), NULL, fill=c("red", "green"), 
                          alpha=c(0.5,0.5), cex = 2, cat.fontface=4, category.names=c("TS (common)", "WGS"), 
                          cat.just = list(c(1,0), c(0,0)))
grid.draw(venn.plot)
dev.off()

pdf("gmkf.rare.venn_diagram.pdf")
venn.plot <- venn.diagram(list(gmkf.rare.position, target.position), NULL, fill=c("green", "red"), 
                          alpha=c(0.5,0.5), cex = 2, cat.fontface=4, category.names=c("WGS (rare)", "TS"),
                          cat.just = list(c(1,0), c(0,0)))
grid.draw(venn.plot)
dev.off()

pdf("gmkf.common.venn_diagram.pdf")
venn.plot <- venn.diagram(list(gmkf.common.position, target.position), NULL, fill=c("green", "red"), 
                          alpha=c(0.5,0.5), cex = 2, cat.fontface=4, category.names=c("WGS (common)", "TS"),
                          cat.just = list(c(1,0), c(0,0)))
grid.draw(venn.plot)
dev.off()

##singleton
target.nonoverlap.gt <- target.gt[which(!target.rare.position %in% gmkf.position),]
count.variant <- function(x){
  (sum(x == "1/1", na.rm = T) + sum(x == "0/1", na.rm = T))
}

target.count.variant <- apply(target.nonoverlap.gt, 1, count.variant)
sum(target.count.variant == 1)  

gmkf.nonoverlap.gt <- gmkf.gt[which(!gmkf.rare.position %in% target.position),]

gmkf.count.variant <- apply(gmkf.nonoverlap.gt, 1, count.variant)
sum(gmkf.count.variant == 1)   











