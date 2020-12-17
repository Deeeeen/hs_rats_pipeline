library(data.table)
library(ggplot2)
##################### pairwise error rate #####################
read in arguments
args <- commandArgs(TRUE) 
pairwise_ERR <- read.table(args[1], strip.white=TRUE, stringsAsFactors=FALSE, sep="\t",
                            header=FALSE, colClasses=c("character", "double", "double", "character", "character"))
out_path <- args[2]

# subset the columns and add column names
pairwise_ERR=pairwise_ERR[,2:5]
colnames(pairwise_ERR)<-c("pairwise_error_rate","num_sites_compared","sample_i","sample_j")
pairwise_ERR = subset(pairwise_ERR, select = c("pairwise_error_rate","sample_i","sample_j"))
# construct the full pairwise error rate
pairwise_ERR_rev<-pairwise_ERR
colnames(pairwise_ERR_rev)<-c("pairwise_error_rate","sample_j","sample_i")
pairwise_ERR<-rbind(pairwise_ERR,pairwise_ERR_rev)
pairwise_ERR<-pairwise_ERR[order(pairwise_ERR$sample_i, pairwise_ERR$sample_j),]
rfids<-unique(c(pairwise_ERR$sample_i,pairwise_ERR$sample_j))
same<-data.frame(sample_i=rfids,sample_j=rfids,pairwise_error_rate=rep(NA,length(rfids)),stringsAsFactors = F)
pairwise_ERR<-pairwise_ERR[,colnames(same)]
pairwise_ERR<-rbind(pairwise_ERR,same)
pairwise_ERR$concordance_rate<-1-pairwise_ERR$pairwise_error_rate

# plot histogram 
ggplot(pairwise_ERR, aes(x=concordance_rate)) + 
    geom_histogram(fill="white", color="black", bins=100)+
    xlab("Concordance Rate")+
    ylab("Frequency")
ggsave(paste0(out_path, "/concordance_rate_hstogram.png"))