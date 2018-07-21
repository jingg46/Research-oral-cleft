library(VariantAnnotation)
target.common <- readVcf("/users/jli/target.8q24.common.recode.vcf.recode.vcf", "hg19")
temp <- rowRanges(target.common)
temp1 <- temp@ranges
target.common.position <- temp1@start

gmkf.common <- readVcf("/users/jli/gmkf.8q24.common.recode.vcf.recode.vcf", "hg19")
temp <- rowRanges(gmkf.common)
temp1 <- temp@ranges
gmkf.common.position <- temp1@start

library("VennDiagram")
pdf("common.venn_diagram.pdf")
venn.plot <- venn.diagram(list(target.common.position, gmkf.common.position), NULL, fill=c("red", "green"), 
                          alpha=c(0.5,0.5), cex = 2, cat.fontface=4, category.names=c("Target common variants", "Gmkf common variants"))
grid.draw(venn.plot)
dev.off()