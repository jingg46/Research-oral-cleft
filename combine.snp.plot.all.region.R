## need gTDT results for both populations for 8q24
load("/users/jli/combine.gTDT.RData")
tdtEuro <- combine.gtdt

## change rs to position
target.rs.position <- read.csv("target.rs.position.csv")
inx = which(names(tdtEuro$pval) %in% target.rs.position$target.rs.name)
for (i in 1:length(inx)){
  names(tdtEuro$pval)[inx[i]] <- target.rs.position$target.position[which(target.rs.position$target.rs.name %in% names(tdtEuro$pval)[inx[i]])]
}

## SNP.FILE has columns ASSOC (+/-), SNP.NAME, LOC, SS.PVAL
tdtEuroSNPFile<-data.frame(ASSOC=ifelse(tdtEuro$RR>1, "+", "-"), 
                           SNP.NAME = names(tdtEuro$RR),  ###TDT output already has snp name
                           LOC=as.numeric(names(tdtEuro$pval)), 
                           SS.PVAL=tdtEuro$pval, stringsAsFactors=FALSE)
tdtEuroSNPFile[!grepl("^r", tdtEuroSNPFile$SNP.NAME),"SNP.NAME"]<-paste0("x", tdtEuroSNPFile[!grepl("^r", tdtEuroSNPFile$SNP.NAME),"SNP.NAME"])

###remove NAs
table(is.na(tdtEuroSNPFile$SS.PVAL))
tdtEuroSNPFile <- tdtEuroSNPFile[!is.na(tdtEuroSNPFile$SS.PVAL),]

write.table(tdtEuroSNPFile, file="/users/jli/tdtEuroSNPFile.combine.all.region.txt", sep="\t", row.names=FALSE, quote=FALSE)

###count # of variants
tdtEuroSNPFile <- read.table("/users/jli/tdtEuroSNPFile.combine.all.region.txt", header = T)
str(tdtEuroSNPFile)   

makeGenotypeFile<-function(genoDat, famTag="euro"){
  nFams<-nrow(genoDat)/3
  mID<-rep(0, nrow(genoDat))
  mID[seq(3, nrow(genoDat), by=3)]<-rownames(genoDat)[seq(2, nrow(genoDat), by=3)]
  fID<-rep(0, nrow(genoDat))
  fID[seq(3, nrow(genoDat), by=3)]<-rownames(genoDat)[seq(1, nrow(genoDat), by=3)]
  pedInfo<-paste(rep(paste0(famTag, 1:nFams), each=3), rownames(genoDat), fID, mID, rep(1, nrow(genoDat)), rep(0, nrow(genoDat)), sep="\t")
  genoDatRecode<-matrix(NA, nrow=nrow(genoDat), ncol=ncol(genoDat))
  genoDatRecode[genoDat == 0]<-"1\t1"
  genoDatRecode[genoDat == 1]<-"1\t2"
  genoDatRecode[genoDat == 2]<-"2\t2"
  genoDatRecode[is.na(genoDat)]<-"0\t0"
  toWrite<-paste(pedInfo, apply(genoDatRecode, 1, function(x) paste(x, collapse="\t")), sep="\t")
  return(toWrite)
}


toWriteEuro <- makeGenotypeFile(combine.geno.filter)

write(toWriteEuro, file="/users/jli/tdtEuroDATFile.combine.all.region.dat")

library(snp.plotter)
snp.plotter(EVEN.SPACED = FALSE, PVAL.THRESHOLD = 1, USE.GBL.PVAL = TRUE,
            SYMBOLS = NA, SAMPLE.LABELS = NULL, LAB.Y = "log", DISP.HAP = FALSE,
            DISP.SNP = TRUE, DISP.COLOR.BAR = TRUE, DISP.PHYS.DIST = TRUE,
            DISP.LEGEND = FALSE, DISP.MARKER.LINES = TRUE, DISP.LDMAP = FALSE,
            DISP.TYPE = "symbol", DISP.MULT.LAB.X = FALSE, DISP.SNP.NAMES = TRUE,
            DISP.CONNECTING.LINES = FALSE, LD.TYPE = "rsquare",
            LD.COLOR.SCHEME = "heat", USE.COLORS = TRUE, COLOR.LIST = "black",
            PALETTE.FILE = NULL, IMAGE.TITLE = NULL, IMAGE.NAME = "/users/jli/combine.all.region.snp.plot",
            IMAGE.TYPE = "pdf", IMAGE.SIZE = 3.5, CONNECTING.LINES.FACTOR = 1,
            CONNECTING.LINES.ADJ = 0, CONNECTING.LINES.VERT.ADJ = -1,
            CONNECTING.LINES.FLEX = 0, SNP.FILE = "/users/jli/tdtEuroSNPFile.combine.all.region.txt", HAP.FILE = NULL,
            GENOTYPE.FILE = "/users/jli/tdtEuroDATFile.combine.all.region.dat", FONT.FACTOR = NULL, SYMBOL.FACTOR = NULL,
            config.file = NULL)






