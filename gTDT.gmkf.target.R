

################################################################################

#Sign into cluster
ssh jli@jhpce01.jhsph.edu

#remove dup id from targeted vcf
#module load vcftools
#vcftools --vcf /users/jli/biallelic.target.vcf --remove target.dup.id.txt --out target.no.dup.vcf --recode
#mv target.no.dup.vcf.recode.vcf target.no.dup.vcf

#Allocate memory
qrsh -l mem_free=30G,h_vmem=31G
R

################load geno file for gmkf data ##################
load("/dcl01/beaty/data/gmkf/euro/vcfs/filtered/chr8_geno_filt_12_06_2017.rda")

###############targeted data############
library(VariantAnnotation)
library(trio)


target <- readVcf("/users/jli/target.no.dup.vcf","hg19")
target.ped <- read.csv("/users/jli/target.ped.no.dup.csv")
target.ped <- target.ped[,-1]    ##remove first row (contain row number )

###drop IDs that are not in complete trios
toDrop<-c("H_ME-DS10781_1-DS10781_1", "H_ME-DS11205_1-DS11205_1")
target<-target[,setdiff(colnames(target), toDrop)]

###get target vcf rs name mapped with position number
temp <- rowRanges(target)
temp1 <- temp@ranges
target.position <- temp1@start
target.rs.name <- names(temp)

target.rs.position <- data.frame(cbind(target.position, target.rs.name), stringsAsFactors = F)
write.csv(target.rs.position, "target.rs.position.csv", row.names = F)

target.ped<-subset(target.ped, !pid %in% toDrop)
target.geno <- vcf2geno(target, target.ped)
#Warning message:
#In vcf2geno(target, target.ped) :
#  Since NA of the SNVs were monomorphic, these SNVs were removed.

###get gmkf vcf rs name mapped with position number
gmkf <- readVcf("/users/jli/gmkf.dup.vcf", "hg19")
temp <- rowRanges(gmkf)
temp1 <- temp@ranges
gmkf.position <- temp1@start
gmkf.rs.name <- names(temp)

gmkf.rs.position <- data.frame(cbind(gmkf.position, gmkf.rs.name), stringsAsFactors = F)
write.csv(gmkf.rs.position, "gmkf.rs.position.csv", row.names = F)

################filtering for target geno file########################
## which rows of the genotype matrix contain the parents?
parRows<-c(seq(1, nrow(target.geno), by=3), seq(2, nrow(target.geno), by=3))

## check mafs and fraction of missing calls
fracMissing<-colSums(is.na(target.geno[parRows,]))/length(parRows)
summary(fracMissing)
sum(fracMissing > 0.05)
target.geno.filter <- target.geno[,-which(fracMissing > 0.05)]   ###remove snp with missing calls > 0.05

mafs<-colSums(target.geno.filter[parRows, ], na.rm=TRUE)/(2*colSums(!is.na(target.geno.filter[parRows,])))
sum(mafs < 0.01, na.rm = T)
target.geno.filter <- target.geno.filter[, -which(mafs < 0.01)]


## Calculate Hardy-Weinberg equlilbrium p-value
## allele frequencies for non-ref (p) and ref (q) alleles
parRows<-c(seq(1, nrow(target.geno.filter), by=3), seq(2, nrow(target.geno.filter), by=3))   ###re-generate index for parents
mafs<-colSums(target.geno.filter[parRows, ], na.rm=TRUE)/(2*colSums(!is.na(target.geno.filter[parRows,])))
p<-mafs
q<-1-mafs
n<-colSums(!is.na(target.geno.filter[parRows,]))

## observed hom non-ref, het, and hom-ref counts
oaa<-colSums(target.geno.filter[parRows,] == 2, na.rm=TRUE)
oAa<-colSums(target.geno.filter[parRows,] == 1, na.rm=TRUE)
oAA<-colSums(target.geno.filter[parRows,] == 0, na.rm=TRUE)

## expected hom non-ref, het, and hom-ref counts
eaa<-n*p^2
eAa<-2*n*p*q
eAA<-n*q^2

HWOut<-pchisq((oaa-eaa)^2/eaa + (oAa-eAa)^2/eAa + (oAA-eAA)^2/eAA, df=1)
## conservatively check p < 0.05 (usually it would not be this strict, but you would use a Bonferroni corrected p-value)
## this does not affect many positions
sum(HWOut < 0.05, na.rm = T)
target.geno.filter <- target.geno.filter[,-which(HWOut < 0.05)]   ###remove snp with HWOut < 0.05
str(target.geno.filter)
#####720 individuals,  2566 variants

target.gtdt <- colTDT(target.geno.filter)   ###use this to obtain gTDT for targeted data (no dup with gmkf)
save("target.geno.filter", "target.gtdt", file = "target.gTDT.RData")


###############filter wgs data#############
## which rows of the genotype matrix contain the parents?
parRows<-c(seq(1, nrow(chr8_geno_filt), by=3), seq(2, nrow(chr8_geno_filt), by=3))

## check mafs and fraction of missing calls
mafs<-colSums(chr8_geno_filt[parRows, ], na.rm=TRUE)/(2*colSums(!is.na(chr8_geno_filt[parRows,])))
fracMissing<-colSums(is.na(chr8_geno_filt[parRows,]))/length(parRows)
summary(fracMissing)
sum(fracMissing > 0.05)
gmkf.geno.filter <- chr8_geno_filt[,-which(fracMissing > 0.05)]   ###remove snp with missing calls > 0.05

## Calculate Hardy-Weinberg equlilbrium p-value
## allele frequencies for non-ref (p) and ref (q) alleles
parRows<-c(seq(1, nrow(gmkf.geno.filter), by=3), seq(2, nrow(gmkf.geno.filter), by=3))
mafs<-colSums(gmkf.geno.filter[parRows, ], na.rm=TRUE)/(2*colSums(!is.na(gmkf.geno.filter[parRows,])))
p<-mafs
q<-1-mafs
n<-colSums(!is.na(gmkf.geno.filter[parRows,]))

## observed hom non-ref, het, and hom-ref counts
oaa<-colSums(gmkf.geno.filter[parRows,] == 2, na.rm=TRUE)
oAa<-colSums(gmkf.geno.filter[parRows,] == 1, na.rm=TRUE)
oAA<-colSums(gmkf.geno.filter[parRows,] == 0, na.rm=TRUE)

## expected hom non-ref, het, and hom-ref counts
eaa<-n*p^2
eAa<-2*n*p*q
eAA<-n*q^2

HWOut<-pchisq((oaa-eaa)^2/eaa + (oAa-eAa)^2/eAa + (oAA-eAA)^2/eAA, df=1)
## conservatively check p < 0.05 (usually it would not be this strict, but you would use a Bonferroni corrected p-value)
## this does not affect many positions
sum(HWOut < 0.05, na.rm = T)
gmkf.geno.filter <- gmkf.geno.filter[,-which(HWOut < 0.05)]   ###remove snp with HWOut < 0.05
str(gmkf.geno.filter)
#########981 individuals, 3172 variants
gmkf.gtdt <- colTDT(gmkf.geno.filter)
save("gmkf.geno.filter", "gmkf.gtdt", file = "gmkf.gTDT.RData")




########################combine target geno with gmkf geno######################
inx.target <- which(colnames(target.geno) %in% colnames(chr8_geno_filt))
inx.gmkf <- which(colnames(chr8_geno_filt) %in% colnames(target.geno))
target.geno <- target.geno[,inx.target]
gmkf.geno <- chr8_geno_filt[,inx.gmkf]
combine.geno <- rbind(gmkf.geno, target.geno)   ###1701 individuals, 2769 positions
str(combine.geno)

####################filtering for combined geno file#########################


## check mafs and fraction of missing calls
mafs<-colSums(combine.geno[parRows, ], na.rm=TRUE)/(2*colSums(!is.na(combine.geno[parRows,])))
sum(mafs < 0.01)
combine.geno.filter <- combine.geno[,-which(mafs < 0.01)]
str(combine.geno.filter)

## which rows of the genotype matrix contain the parents?
parRows<-c(seq(1, nrow(combine.geno.filter), by=3), seq(2, nrow(combine.geno.filter), by=3))
fracMissing<-colSums(is.na(combine.geno.filter[parRows,]))/length(parRows)
summary(fracMissing)
sum(fracMissing > 0.05)
combine.geno.filter <- combine.geno.filter[,-which(fracMissing > 0.05)] ###remove snp with missing calls > 0.05
str(combine.geno.filter)

## Calculate Hardy-Weinberg equlilbrium p-value
## allele frequencies for non-ref (p) and ref (q) alleles
parRows<-c(seq(1, nrow(combine.geno.filter), by=3), seq(2, nrow(combine.geno.filter), by=3))   ###re-generate index for parents
mafs<-colSums(combine.geno.filter[parRows, ], na.rm=TRUE)/(2*colSums(!is.na(combine.geno.filter[parRows,])))  ###re-generate maf
p<-mafs
q<-1-mafs
n<-colSums(!is.na(combine.geno.filter[parRows,]))

## observed hom non-ref, het, and hom-ref counts
oaa<-colSums(combine.geno.filter[parRows,] == 2, na.rm=TRUE)
oAa<-colSums(combine.geno.filter[parRows,] == 1, na.rm=TRUE)
oAA<-colSums(combine.geno.filter[parRows,] == 0, na.rm=TRUE)

## expected hom non-ref, het, and hom-ref counts
eaa<-n*p^2
eAa<-2*n*p*q
eAA<-n*q^2

HWOut<-pchisq((oaa-eaa)^2/eaa + (oAa-eAa)^2/eAa + (oAA-eAA)^2/eAA, df=1)
## conservatively check p < 0.05 (usually it would not be this strict, but you would use a Bonferroni corrected p-value)
## this does not affect many positions
sum(HWOut < 0.05, na.rm = T)
combine.geno.filter <- combine.geno.filter[,-which(HWOut < 0.05)]   ###remove snp with HWOut < 0.05
str(combine.geno.filter)
#######1701 individuals, 2646 variants

combine.gtdt <- colTDT(combine.geno.filter)
save("combine.geno.filter", "combine.gtdt", file = "combine.gTDT.RData")


