############change format of targeted ped file

target.ped <- read.table("targeted.ped.txt")
target.ped <- target.ped[,c(2,1,3:6)]
colnames(target.ped) <- c("famid", "pid", "fatid", "motid", "sex", "affected")
target.ped$sex <- ifelse(target.ped$sex ==1, "male", "female")
target.ped$affected <- ifelse(target.ped$affected ==1, "affected", "unaffected")

########change format of child id (parents id matched)
target.id <- read.table("targeted.id.txt")
target.id$V1 <- as.character(target.id$V1)
target.ped$pid <- as.character(target.ped$pid)

target.id$sub <- lapply(strsplit(target.id$V1, "-"), function(x) x[2])
inx = which(target.ped$pid %in% target.id$sub)
for (i in 1:length(inx)){
  target.ped$pid[inx[i]] <- target.id$V1[which(target.id$sub == target.ped$pid[inx[i]])]
}

###remove duplicate id from ped file
target.dup.id <- read.table("target.dup.id.txt")
inx2 = which(target.ped$pid %in% target.dup.id$V1)
target.ped <- target.ped[-inx2,]



###check
target.no.dup <- read.table("target.no.dup.id.txt")
which(!target.ped$pid %in% target.no.dup$V1)
target.ped[which(!target.ped$pid %in% target.no.dup$V1),]   ###all have a missing 0, i.e.DS11099_1 not DS11099_01
target.ped$pid[which(!target.ped$pid %in% target.no.dup$V1)] <- sub("(.{8})(.*)", "\\10\\2", target.ped$pid[which(!target.ped$pid %in% target.no.dup$V1)])
inx = which(target.ped$pid %in% target.id$sub)
for (i in 1:length(inx)){
  target.ped$pid[inx[i]] <- target.id$V1[which(target.id$sub == target.ped$pid[inx[i]])]
}

write.csv(target.ped, "target.ped.no.dup.csv")



