library(data.table)
library(ggplot2)
library(scales)
#### read in the directory name
args = commandArgs(trailingOnly=TRUE)
demux_metrics_f <- list.files(path=args[1])
demux_metrics_ls = list()
unmatched_ls = list()
#### load all the demux metrics files
for(i in 1:length(demux_metrics_f)){
    temp <- read.table(paste0(args[1], "/", demux_metrics_f[i]), header=TRUE,  strip.white=TRUE,sep="\t")
    demux_metrics_ls[[i]] <- temp[temp$barcode_name != "unmatched", ]
    unmatched_ls[[i]] <- temp[temp$barcode_name == "unmatched", ]["templates"]/sum(temp$templates)
}
demux_metrics = do.call(rbind, demux_metrics_ls)
unmatched = do.call(rbind, unmatched_ls)
demux_metrics$templates <- demux_metrics$templates/1e6

metadata1 <- read.table(args[2], header=TRUE,  strip.white=TRUE,sep=",")
metadata1 <- metadata1[metadata1$strain=="Heterogenous stock",]
metadata1$library_project <- paste(metadata1$library_name,metadata1$project_name, sep="-")
metadata <- metadata1[order(metadata1$library_project),]
metadata$index <- 1:nrow(metadata)
metadata$templates <- demux_metrics$templates[match(metadata$rfid, demux_metrics$barcode_name)]

#### box plot for # of reads per sample
means <- aggregate(templates ~  library_name, metadata, mean)
ggplot(metadata, aes(x=library_name, y=templates)) +
    geom_boxplot() +
    stat_summary(fun="mean", geom="point", col="red", shape=1, size=2) +
    geom_text(data=means, aes(label=comma(templates), y=templates+0.25), size=3) +
    ylab("# of reads (million)")+
    theme(text=element_text(size=8), axis.text.x=element_text(angle=-45, size=8))
ggsave(paste0(args[3], "/demux_matched_reads_boxplot.png"), width = 5, height = 5, dpi = 300, units = "in")

#### box plot for % of unmatched reads per fastq file
ggplot(unmatched, aes(x="", y=templates)) +
    geom_boxplot() +
    stat_summary(fun="mean", geom="point", col="red", shape=1, size=2) +
    geom_text(label=percent(mean(unmatched$templates), accuracy=0.001), y =mean(unmatched$templates)+0.0002) +
    scale_y_continuous(labels=scales::percent_format(accuracy=0.1))+
    xlab("hs rats") +
    ylab("demux unmatched")
ggsave(paste0(args[3], "/demux_unmatched_reads_boxplot.png"), width = 5, height = 5, dpi = 300, units = "in")

#### scatter plot for # of reads per sample scatter plot
ggplot(metadata,aes(x=index,y=templates,shape=runid,color=library_name)) +
    geom_point(size=2) +
    xlab("Samples") +
    ylab("# of reads (million)") #+ theme(legend.position="none")
ggsave(paste0(args[3], "/demux_matched_reads_scatterplot.png"), width = 7, height = 5, dpi = 300, units = "in")

