setwd("~/Desktop")
pedi <- read.csv("gmkf_euro_completetrios.csv")
id <- read.table("sample.txt", stringsAsFactors = F)
id$V1 <- substr(id$V1, 6, 11)
colnames(id) <- "Individual.ID"
mer <- merge(id, pedi, by = "Individual.ID")

a <- sample[which(!id$Individual.ID %in% pedi$Individual.ID),1]
write.csv(a, "sample with missing pedigree data.csv")