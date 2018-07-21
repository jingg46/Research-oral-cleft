setwd("~/Desktop/Research")
##read in ped file
pedi <- read.csv("gmkf_euro_completetrios.csv")
##read in the extracted id from vcf file
id <- read.table("sample.txt", stringsAsFactors = F)
##flag ids that end up with letter "B"
id$b <- substr(id$V1, 19,19)
##flag ids that end up with letter "B"
pedi$b <- substr(pedi$dbGaP.submitted.sample_id, 7,7)
pedi <- pedi[,c(1:3, 14, 4:13)]
##get the short version of id in vcf id for mathcing
id$shor <- substr(id$V1, 6, 11)
##keep only complete trios
a <- id$V1[which(id$shor %in% pedi$Individual.ID)]

write.table(a, file = "sample.with.complete.trios.txt", quote = FALSE, row.names = F, col.names = F)
