setwd("~/Desktop/Research")
sample.8q24 <- read.csv("sample.8q24.txt", header = F, stringsAsFactors = F)
pedi <- read.csv("gmkf_euro_completetrios.csv", stringsAsFactors = F)
pedi <- pedi[,c(1,2,4,5)]
pedi <- pedi[pedi$Father.ID != 0,]

sample.8q24$b <- substr(sample.8q24$V1, 19,19)
temp <- substr(sample.8q24[which(sample.8q24$b == "B"),1], 6,11)

inx <- vector(length=length(temp))
for (j in 1:length(temp)){
  inx[j] <- which(pedi == temp[j])
}

for(i in 2:4){
  pedi[,i] <- paste0("H_TZ-", pedi[,i], "-", pedi[,i])
}

pedi[inx[which(inx <= 664)] -332, 2] <- paste0(pedi[inx[which(inx <= 664)] - 332, 2], "B")
pedi[inx[which(inx > 664 & inx <= 996)] - 664, 3] <- paste0(pedi[inx[which(inx > 664 & inx <= 996)] - 664, 3], "B")
pedi[inx[which(inx > 996)] - 996, 4] <- paste0(pedi[inx[which(inx > 996)] - 996, 4], "B")



write.table(pedi, file = "mendelian.pedi.txt", quote = FALSE, row.names = F, col.names = F, sep = " ")


