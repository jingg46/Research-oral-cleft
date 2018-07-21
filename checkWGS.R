## check on HW-equilibrium in variants included in colTDT output from WGS

load("/dcl01/beaty/data/gmkf/euro/vcfs/filtered/chr8_geno_filt_12_06_2017.rda")

## which rows of the genotype matrix contain the parents?
parRows<-c(seq(1, nrow(chr8_geno_filt), by=3), seq(2, nrow(chr8_geno_filt), by=3))

## check mafs and fraction of missing calls
mafs<-colSums(chr8_geno_filt[parRows, ], na.rm=TRUE)/(2*colSums(!is.na(chr8_geno_filt[parRows,])))
fracMissing<-colSums(is.na(chr8_geno_filt[parRows,]))/length(parRows)
summary(fracMissing)
sum(fracMissing > 0.05)

## Calculate Hardy-Weinberg equlilbrium p-value
## allele frequencies for non-ref (p) and ref (q) alleles
p<-mafs
q<-1-mafs
n<-colSums(!is.na(chr8_geno_filt[parRows,]))

## observed hom non-ref, het, and hom-ref counts
oaa<-colSums(chr8_geno_filt[parRows,] == 2, na.rm=TRUE)
oAa<-colSums(chr8_geno_filt[parRows,] == 1, na.rm=TRUE)
oAA<-colSums(chr8_geno_filt[parRows,] == 0, na.rm=TRUE)

## expected hom non-ref, het, and hom-ref counts
eaa<-n*p^2
eAa<-2*n*p*q
eAA<-n*q^2

HWOut<-pchisq((oaa-eaa)^2/eaa + (oAa-eAa)^2/eAa + (oAA-eAA)^2/eAA, df=1)
## conservatively check p < 0.05 (usually it would not be this strict, but you would use a Bonferroni corrected p-value)
## this does not affect many positions
sum(HWOut < 0.05)