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
    temp$group <- as.factor(c(rep(substr(demux_metrics_f[i], 1, nchar(demux_metrics_f[i])-26), length(temp[,1]))))
    demux_metrics_ls[[i]] <- temp[temp$barcode_name != "unmatched", ]
    unmatched_ls[[i]] <- temp[temp$barcode_name == "unmatched", ]["templates"]/sum(temp$templates)
}
demux_metrics = do.call(rbind, demux_metrics_ls)
unmatched = do.call(rbind, unmatched_ls)
#### plot for # of reads per sample
means <- aggregate(templates ~  group, demux_metrics, mean)
ggplot(demux_metrics, aes(x=group, y=templates)) +
    geom_boxplot() +
    stat_summary(fun="mean", geom="point", col="red", shape=1, size=2) +
    geom_text(data=means, aes(label=comma(templates), y=templates+2e5), size=3) +
    ylab("# of reads per sample")+
    theme(text=element_text(size=8), axis.text.x=element_text(angle=-45, size=8))
ggsave(args[2])
#### plot for % of unmatched reads per fastq file
ggplot(unmatched, aes(x="", y=templates)) +
    geom_boxplot() +
    stat_summary(fun="mean", geom="point", col="red", shape=1, size=2) +
    geom_text(label=percent(mean(unmatched$templates), accuracy=0.001), y =mean(unmatched$templates)+0.0002) +
    scale_y_continuous(labels=scales::percent_format(accuracy=0.1))+
    xlab("hs rats") +
    ylab("demux unmatched")
ggsave(args[3])