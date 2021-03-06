library(plyr)
# read in arguments
args <- commandArgs(TRUE) 
if(file.size(args[2]) > 0){
    metadata_files <- read.table(args[2], strip.white=TRUE, stringsAsFactors=FALSE, header=FALSE)
}else{
    metadata_files <-c()
}
pedigree_files <- read.table(args[3], strip.white=TRUE, stringsAsFactors=FALSE, header=FALSE)
out_path <- args[4]
# read in metadata files
metadata <- vector("list", length(metadata_files)+1)
temp_metadata <- read.table(args[1], stringsAsFactors=FALSE,
                        header=TRUE,  strip.white=TRUE,sep=",")
temp_metadata <- temp_metadata[temp_metadata$strain=="Heterogenous stock",]
metadata[[1]] <- data.frame(temp_metadata)
i<-2
for(file in metadata_files$V1){
    temp_metadata <- read.table(file, stringsAsFactors=FALSE,
                            header=TRUE,  strip.white=TRUE,sep=",")
    temp_metadata <- temp_metadata[temp_metadata$strain=="Heterogenous stock",]
    metadata[[i]] <- data.frame(temp_metadata)
    i <- i+1
}
all_metadata<-do.call('rbind.fill', metadata)
# read in pedigree files
pedigree <- vector("list", length(pedigree_files))
i<-1
for(file in pedigree_files$V1){
    temp_pedigree <- read.table(file, stringsAsFactors=FALSE,
                            header=TRUE,  strip.white=TRUE,sep=",")
    temp_pedigree <- subset(temp_pedigree, select = c(rfid, dames, sires))
    pedigree[[i]] <- data.frame(temp_pedigree)
    i <- i+1
}
all_pedigree<-do.call('rbind.fill', pedigree)
# merge metadata and pedigree
all_info <- merge(all_metadata, all_pedigree, by.x = "rfid", by.y = "rfid", all = TRUE)
metadata_prefix <- paste0("hs_rats_n", length(all_info$rfid), "_", format(Sys.Date(), "%m%d%Y"))
write.table(data.frame(all_info), paste0(out_path, '/', metadata_prefix, '_metadata.csv'), row.names = FALSE, sep=',' )