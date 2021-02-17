library(data.table)
library(ggplot2)
args = commandArgs(trailingOnly=TRUE)

data_path <- args[1]
out_path <- args[2]

info_score_hwe <- read.table(data_path, header=FALSE,  strip.white=TRUE,sep="\t")
colnames(info_score_hwe)<-c("pos","info_score")

png(paste0(out_path, "/stitch_info_score_histo.png"))
hist(info_score_hwe$info_score,  breaks = 100, xlab="Info Score", ylab="Frequency", main="1536 HS Rats Info Score")
dev.off()
