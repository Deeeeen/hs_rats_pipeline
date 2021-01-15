library(data.table)
library(ggplot2)
library(plyr)
# read in arguments
args <- commandArgs(TRUE)
dir_path <- args[1]
flow_cell_dir_path <- args[2]
all_metadata <- read.table(paste0(flow_cell_dir_path, "/results/", list.files(paste0(flow_cell_dir_path, "/results/"), pattern='*metadata.csv')[1]), stringsAsFactors=FALSE,
                        header=TRUE,  strip.white=TRUE,sep=",")
all_metadata <- all_metadata[all_metadata$strain == "Heterogenous stock",]
# read in missing rate
missing_rate <- read.table(paste0(flow_cell_dir_path, "/results/genotype_result/stitch_result/plink/", list.files(paste0(flow_cell_dir_path, "/results/genotype_result/stitch_result/plink/"), pattern='*_stitch.smiss')[1]), stringsAsFactors=FALSE,
                        header=FALSE,  strip.white=TRUE,sep="\t")
colnames(missing_rate)<-c("rfid", "IID", "MISSING_CT", "OBS_CT", "F_MISS")
# missing_rate$reads <- 0  
missing_rate$reads <- as.integer(NA)

# find the folders of all flow cells
folders <- vector("character", length(unique(all_metadata$runid)))
i<-1
for(run_id in unique(all_metadata$runid)){
    folder <- list.files(dir_path, pattern=paste0('*', run_id, '*'))[1]
    folders[i] <- folder
    i <- i+1
}

# read in all demux metric files
metrics <- vector("list", length(unique(all_metadata$Library_ID)))
i<-1
for(folder in unique(folders)){
    metrics_files<-list.files(paste0(dir_path, '/', folder, '/demux/metrics/'))
    for(file in metrics_files){
        temp_metrics<-read.table(paste0(dir_path, '/', folder, '/demux/metrics/', file), stringsAsFactors=FALSE,
                        header=TRUE,  strip.white=TRUE,sep="\t")
        temp_metrics <- temp_metrics[temp_metrics$barcode_name != "unmatched",]
        temp_metrics <- temp_metrics[c("barcode_name", "templates")]
        metrics[[i]] <- data.frame(temp_metrics)
        i <- i+1 
    }
}

# merge data
all_metrics <- do.call('rbind.fill', metrics)
all_info <- merge(missing_rate, all_metrics, by.x = "rfid", by.y = "barcode_name", all = TRUE)

# plot scatter plot
all_info$reads <- all_info$templates/1000000
ggplot(all_info, aes(x=reads,y=F_MISS)) +
    geom_point(size=1,shape=1) +
    xlab("Reads (million)") +
    ylab("Missing Rate") 
    # +theme(legend.position="none")
ggsave(paste0(flow_cell_dir_path, "/results/genotype_result/stitch_result/plink/missing_vs_reads.png"), width = 5, height = 5, dpi = 300, units = "in")
