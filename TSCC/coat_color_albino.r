library("data.table")
############ ALBINO ############
# read in arguments
args <- commandArgs(TRUE) 
ped <- read.table(args[1], strip.white=TRUE, stringsAsFactors=FALSE, header=FALSE)
metadata <- read.table(args[2], strip.white=TRUE, stringsAsFactors=FALSE, header=TRUE, sep=",")
out_path <- args[3]

# merge files
ped$GT=paste0(ped$V7,ped$V8)
ped=ped[,c("V2","GT")]
colnames(ped)[1]<-"rfid"
ped$rfid<-toupper(ped$rfid)
metadata$rfid <- as.character(metadata$rfid)
metadata$rfid <- toupper(metadata$rfid)
all=merge(ped,metadata,by="rfid",all=F)
albino<-all[which(all$coatcolor=="ALBINO"),]

# flag samples
# these should be flagged
hets_which_should_not_be_albino<-albino[which(albino$GT=="TC"),]
non_albino<-all[which(!all$coatcolor=="ALBINO"),]

# these should also be flagged
homo_which_should_be_albino<-non_albino[which(non_albino$GT=="TT"),]
flagged<-do.call("rbind", list(hets_which_should_not_be_albino, homo_which_should_be_albino))
write.table( data.frame(flagged), paste0(out_path, '_albino_outliers.csv'), row.names = FALSE, sep=',' )