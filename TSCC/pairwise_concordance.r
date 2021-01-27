library(data.table)
library(ggplot2)
##################### pairwise error rate #####################
# read in arguments
args <- commandArgs(TRUE) 
pairwise_ERR <- read.table(args[1], strip.white=TRUE, stringsAsFactors=FALSE, sep="\t",
                            header=FALSE, colClasses=c("character", "double", "double", "character", "character"))
out_path <- args[3]

# subset the columns and add column names
pairwise_ERR=pairwise_ERR[,2:5]
colnames(pairwise_ERR)<-c("pairwise_error_rate","num_sites_compared","sample_i","sample_j")
pairwise_ERR = subset(pairwise_ERR, select = c("pairwise_error_rate","sample_i","sample_j"))
pairwise_ERR$concordance_rate <- 1-pairwise_ERR$pairwise_error_rate
# plot histogram 
ggplot(pairwise_ERR, aes(x=concordance_rate)) + 
    geom_histogram(fill="white", color="black", bins=100)+
    xlab("Concordance Rate")+
    ylab("Frequency")
ggsave(paste0(out_path, "/concordance_rate_hstogram.png"), width = 5, height = 5, dpi = 300, units = "in")

# concordance rate filter > 0.8
pairwise_ERR_outliers <- pairwise_ERR[pairwise_ERR$concordance_rate > 0.8, ]
row.names(pairwise_ERR_outliers) <- NULL

# keep the sample pairs that are NOT from the same family
metadata <- read.table(args[2], strip.white=TRUE, stringsAsFactors=FALSE, header=TRUE, sep=",")
metadata$family_id <- paste(metadata$dames,metadata$sires, sep="-")
pairwise_ERR_outliers$sample_i_family <- sapply(pairwise_ERR_outliers$sample_i, function(x) metadata[metadata$rfid==x,]$family_id)
pairwise_ERR_outliers$sample_j_family <- sapply(pairwise_ERR_outliers$sample_j, function(x) metadata[metadata$rfid==x,]$family_id)
pairwise_ERR_outliers$same_family <- ifelse(pairwise_ERR_outliers$sample_j_family == pairwise_ERR_outliers$sample_i_family, TRUE, FALSE)
pairwise_ERR_outliers <- pairwise_ERR_outliers[pairwise_ERR_outliers$same_family == FALSE, ]
pairwise_ERR_outliers = subset(pairwise_ERR_outliers, select = c("sample_i", "sample_j", "sample_i_family", "sample_j_family", "concordance_rate"))

#### Output outliers 
write.table(data.frame(pairwise_ERR_outliers), paste0(out_path, "/concordance_rate_outliers.csv"), row.names = FALSE, sep=',' )

######## version 1 ########
# # keep the sample pairs that are NOT from the same family
# metadata <- read.table(args[2], strip.white=TRUE, stringsAsFactors=FALSE, header=TRUE, sep=",")
# metadata$family_id <- paste(metadata$dames,metadata$sires, sep="-")
# pairwise_ERR_outliers$sample_i_family <- sapply(pairwise_ERR_outliers$sample_i, function(x) metadata[metadata$rfid==x,]$family_id)
# pairwise_ERR_outliers$sample_j_family <- sapply(pairwise_ERR_outliers$sample_j, function(x) metadata[metadata$rfid==x,]$family_id)
# pairwise_ERR_outliers$same_family <- ifelse(pairwise_ERR_outliers$sample_j_family == pairwise_ERR_outliers$sample_i_family, TRUE, FALSE)
# pairwise_ERR_outliers <- pairwise_ERR_outliers[pairwise_ERR_outliers$same_family == FALSE, ]
# pairwise_ERR_outliers = subset(pairwise_ERR_outliers, select = c("sample_i", "sample_j", "concordance_rate"))

# # construct the full pairwise error rate
# pairwise_ERR_outliers_rev<-pairwise_ERR_outliers
# colnames(pairwise_ERR_outliers_rev)<-c("sample_j", "sample_i", "concordance_rate")
# unique_rfids<-unique(c(pairwise_ERR_outliers$sample_i,pairwise_ERR_outliers$sample_j))
# same<-data.frame(sample_i=unique_rfids, sample_j=unique_rfids, concordance_rate=rep(NA,length(unique_rfids)), stringsAsFactors = F)
# pairwise_ERR_outliers<-rbind(pairwise_ERR_outliers,pairwise_ERR_outliers_rev,same)
# row.names(pairwise_ERR_outliers) <- NULL
# pairwise_ERR_outliers$concordance_rate <- as.double(pairwise_ERR_outliers$concordance_rate)
# pairwise_ERR_outliers<-pairwise_ERR_outliers[order(pairwise_ERR_outliers$sample_i, pairwise_ERR_outliers$sample_j),]

# ggplot(pairwise_ERR_outliers, aes(sample_i, sample_j, fill=concordance_rate))+
#     geom_tile()+
#     xlab("Samples")+
#     ylab("Samples")
# ggsave(paste0(out_path, "concordance_rate_outlier_heatmap.png"), width = 5, height = 5, dpi = 300, units = "in")


######## version 2 ########
# # keep the sample pairs that are NOT from the same family
# metadata <- read.table(args[2], strip.white=TRUE, stringsAsFactors=FALSE, header=TRUE, sep=",")
# metadata$family_id <- paste(metadata$dames,metadata$sires, sep="-")
# pairwise_ERR_outliers$sample_i_family <- sapply(pairwise_ERR_outliers$sample_i, function(x) metadata[metadata$rfid==x,]$family_id)
# pairwise_ERR_outliers$sample_j_family <- sapply(pairwise_ERR_outliers$sample_j, function(x) metadata[metadata$rfid==x,]$family_id)
# pairwise_ERR_outliers$same_family <- ifelse(pairwise_ERR_outliers$sample_j_family == pairwise_ERR_outliers$sample_i_family, TRUE, FALSE)
# pairwise_ERR_outliers <- pairwise_ERR_outliers[pairwise_ERR_outliers$same_family == FALSE, ]
# unique_rfids<-unique(c(pairwise_ERR_outliers$sample_i,pairwise_ERR_outliers$sample_j))
# pairwise_ERR_outliers <- pairwise_ERR[pairwise_ERR$sample_i %in% unique_rfids & pairwise_ERR$sample_j %in% unique_rfids, ]
# pairwise_ERR_outliers = subset(pairwise_ERR_outliers, select = c("sample_i", "sample_j", "concordance_rate"))

# # construct the full pairwise error rate
# pairwise_ERR_outliers_rev<-pairwise_ERR_outliers
# colnames(pairwise_ERR_outliers_rev)<-c("sample_j", "sample_i", "concordance_rate")
# # unique_rfids<-unique(c(pairwise_ERR_outliers$sample_i,pairwise_ERR_outliers$sample_j))
# same<-data.frame(sample_i=unique_rfids, sample_j=unique_rfids, concordance_rate=rep(NA,length(unique_rfids)), stringsAsFactors = F)
# pairwise_ERR_outliers<-rbind(pairwise_ERR_outliers,pairwise_ERR_outliers_rev,same)
# row.names(pairwise_ERR_outliers) <- NULL
# pairwise_ERR_outliers$concordance_rate <- as.double(pairwise_ERR_outliers$concordance_rate)
# pairwise_ERR_outliers<-pairwise_ERR_outliers[order(pairwise_ERR_outliers$sample_i, pairwise_ERR_outliers$sample_j),]

# # plot heatmap
# ggplot(pairwise_ERR_outliers, aes(sample_i, sample_j, fill=concordance_rate))+
#     geom_tile()+
#     xlab("Samples")+
#     ylab("Samples")
# ggsave(paste0(out_path, "concordance_rate_outlier_heatmap.png"), width = 5, height = 5, dpi = 300, units = "in")



