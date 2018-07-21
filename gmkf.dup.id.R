###find dup id in gmkf
setwd("~/Desktop/Research")
combine.ped <- read.csv("gmkf_targeted_merge_peds.csv")
gmkf.id <- combine.ped[which(combine.ped$targeted.status == "both"),2:3]

gmkf.dup.id <- paste0("H_TZ-", gmkf.id[,1], "-", gmkf.id[,2])  ###474

###targeted only have 403 duplicate IDs (some maybe deleted due to mendel error or incomplete trios)
###keep only 403 dup ids in targeted

target.dup.id <- read.table("gmkf.target.dup.id.txt", stringsAsFactors = F)
a <- strsplit(target.dup.id[,1], "-")
a <- sapply(a,function(x) x[2])   ###extract the middle portion

###use the middle portion to find corresponding gmkf id name
gmkf.id <- combine.ped[which(combine.ped$dsiid %in% a),c(2:3, 15)]
gmkf.id$gmkf.id <- paste0("H_TZ-", gmkf.id[,1], "-", gmkf.id[,2])
###match two IDs together
target.dup.id <- cbind(target.dup.id, a)
colnames(target.dup.id)[2] <- "dsiid"
combine.id <- merge(gmkf.id, target.dup.id, by = "dsiid")


###check if all 403 are complete trio
gmkf.complete.trio <- read.table("sample.with.complete.trios.txt", stringsAsFactors = F)
combine.id$gmkf.id[which(!combine.id$gmkf.id %in% gmkf.complete.trio[,1])]  ##1 gmkf.id not found in complete trio data
gmkf.complete.trio$V1[which(grepl(substr(combine.id$gmkf.id[which(!combine.id$gmkf.id %in% gmkf.complete.trio[,1])], 6, 11),gmkf.complete.trio$V1))]
###complete trio data does not have "B" at the end but gmkf.id has "B"
###remove "B"
combine.id$gmkf.id[which(!combine.id$gmkf.id %in% gmkf.complete.trio[,1])] <- "H_TZ-PA1546-PA1546"

###not complete trio, check
target.id <- read.table("targeted.id.txt")

write.table(combine.id$gmkf.id, "gmkf.dup.id.txt", quote = F, row.names = F, col.names = F)
write.table(combine.id$V1, "target.dup.id.txt", quote = F, row.names = F, col.names = F)

combine.id <- combine.id[,c(4,5)]
colnames(combine.id)[2] <- "target.id"
write.csv(combine.id, "combine.id.csv", row.names = F)




