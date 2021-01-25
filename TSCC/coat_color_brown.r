library(data.table)
############ ALBINO ############
# read in arguments
args <- commandArgs(TRUE) 
ped <- read.table(args[1], strip.white=TRUE, stringsAsFactors=FALSE, header=FALSE)
metadata <- read.table(args[2], strip.white=TRUE, stringsAsFactors=FALSE, header=TRUE, sep=",")
out_path <- args[3]

# merge files
ped$GT_150285633=paste0(ped$V7,ped$V8)
ped$GT_150288295=paste0(ped$V9,ped$V10)
ped$GT_150432118=paste0(ped$V11,ped$V12)
ped$GT_150449245=paste0(ped$V13,ped$V14)
ped$GT_150488934=paste0(ped$V15,ped$V16)
ped$GT_150530733=paste0(ped$V17,ped$V19)
ped$GT_150584522=paste0(ped$V19,ped$V20)
colnames(ped)[1]<-"rfid"
ped$rfid<-toupper(ped$rfid)
ped <- subset(ped, select = c(rfid, GT_150285633,GT_150288295,GT_150432118,GT_150449245,GT_150488934,GT_150530733,GT_150584522)) 
metadata$rfid <- as.character(metadata$rfid)
metadata$rfid <- toupper(metadata$rfid)
out_csv <- merge(ped,metadata,by="rfid",all=F)
out_csv <- subset(out_csv, select = c(rfid, coatcolor, GT_150285633,GT_150288295,GT_150432118,GT_150449245,GT_150488934,GT_150530733,GT_150584522)) 

write.csv(out_csv, paste0(out_path, "/coat_color_brown_GT.csv"), row.names=F, quote=F)
