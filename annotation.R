#$ -l mem_free=10G,h_vmem=11G
#$ -l h_fsize=1000G

library(data.table)
c8.anno <- fread("/dcl01/beaty/data/gmkf/euro/anno/c8.annotation.txt", header = T, sep = "\t")
c8.anno <- c8.anno[c8.anno$StartPosition >= 129295896 & c8.anno$StartPosition <= 130354946, ]
write.table(c8.anno, "/dcl01/beaty/data/gmkf/euro/anno/c8.annotation.8124.txt")
