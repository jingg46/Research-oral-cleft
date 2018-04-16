## need gTDT results for both populations for 8q24
load("/dcl01/beaty/data/gmkf/euro/vcfs/filtered/tdt_chr8_geno_filt_12_06_2017.rda")
tdtEuro <- tdt_chr8_geno_filt

annoMerge8q24 <- read.table("/dcl01/beaty/data/gmkf/euro/anno/c8.annotation.8124.txt", header = T)
start<-129900000
end<-130000000

## change rs to position
plink.map <- read.table("/dcl01/beaty/data/gmkf/euro/vcfs/filtered/8q24.recode.map", header = F, stringsAsFactors = F)
plink.map <- plink.map[!duplicated(plink.map$V2),]
inx <- which(plink.map$V2 %in% names(tdtEuro$pval))
names(tdtEuro$pval) <- plink.map$V4[inx]

## SNP.FILE has columns ASSOC (+/-), SNP.NAME, LOC, SS.PVAL
tdtEuroSNPFile<-data.frame(ASSOC=ifelse(tdtEuro$RR>1, "+", "-"), 
                           SNP.NAME = names(tdtEuro$RR),  ###TDT output already has snp name
                           LOC=as.numeric(names(tdtEuro$pval)), 
                           SS.PVAL=tdtEuro$pval, stringsAsFactors=FALSE)
tdtEuroSNPFile[!grepl("^r", tdtEuroSNPFile$SNP.NAME),"SNP.NAME"]<-paste0("x", tdtEuroSNPFile[!grepl("^r", tdtEuroSNPFile$SNP.NAME),"SNP.NAME"])

write.table(subset(tdtEuroSNPFile, LOC>=start & LOC <= end), file="/dcl01/beaty/data/gmkf/euro/vcfs/filtered/tdtEuroSNPFile.txt", sep="\t", row.names=FALSE, quote=FALSE)


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

load("/dcl01/beaty/data/gmkf/euro/vcfs/filtered/chr8_geno_filt_12_06_2017.rda")

toWriteEuro <- makeGenotypeFile(chr8_geno_filt)
write(toWriteEuro, file="/dcl01/beaty/data/gmkf/euro/vcfs/filtered/tdtEuroDATFile.dat")

library(snp.plotter)
snp.plotter(EVEN.SPACED = FALSE, PVAL.THRESHOLD = 1, USE.GBL.PVAL = TRUE,
            SYMBOLS = NA, SAMPLE.LABELS = NULL, LAB.Y = "log", DISP.HAP = FALSE,
            DISP.SNP = TRUE, DISP.COLOR.BAR = TRUE, DISP.PHYS.DIST = TRUE,
            DISP.LEGEND = FALSE, DISP.MARKER.LINES = TRUE, DISP.LDMAP = TRUE,
            DISP.TYPE = "symbol", DISP.MULT.LAB.X = FALSE, DISP.SNP.NAMES = TRUE,
            DISP.CONNECTING.LINES = TRUE, LD.TYPE = "rsquare",
            LD.COLOR.SCHEME = "heat", USE.COLORS = TRUE, COLOR.LIST = "black",
            PALETTE.FILE = NULL, IMAGE.TITLE = NULL, IMAGE.NAME = "/dcl01/beaty/data/gmkf/euro/vcfs/filtered/gTDTEuro8q24",
            IMAGE.TYPE = "pdf", IMAGE.SIZE = 3.5, CONNECTING.LINES.FACTOR = 1,
            CONNECTING.LINES.ADJ = 0, CONNECTING.LINES.VERT.ADJ = -1,
            CONNECTING.LINES.FLEX = 0, SNP.FILE = "/dcl01/beaty/data/gmkf/euro/vcfs/filtered/tdtEuroSNPFile.txt", HAP.FILE = NULL,
            GENOTYPE.FILE = "/dcl01/beaty/data/gmkf/euro/vcfs/filtered/tdtEuroDATFile.dat", FONT.FACTOR = NULL, SYMBOL.FACTOR = NULL,
            config.file = NULL)