library(data.table)
library(ggplot2)
library(scales)

args = commandArgs(trailingOnly=TRUE)
#### read in the metadata file
metadata1 <- read.table(args[1], header=TRUE,  strip.white=TRUE,sep=",")
metadata1 <- metadata1[metadata1$strain=="Heterogenous stock",]
metadata1$library_project <- paste(metadata1$library_name,metadata1$project_name, sep="-")
metadata <- metadata1[order(metadata1$library_project),]
metadata$index <- 1:nrow(metadata)

########## markDuplicates metrics##########
mkDup_metrics <- read.table(args[2], header=TRUE,  strip.white=TRUE,sep="\t")
metadata$mapped_read_pairs <- mkDup_metrics$READ_PAIRS_EXAMINED[match(metadata$rfid, mkDup_metrics$rfid)]
metadata$mapped_unpaired_reads <- mkDup_metrics$UNPAIRED_READS_EXAMINED[match(metadata$rfid, mkDup_metrics$rfid)]
metadata$unmapped_reads <- mkDup_metrics$UNMAPPED_READS[match(metadata$rfid, mkDup_metrics$rfid)]
metadata$mapped_examined_reads <- metadata$mapped_read_pairs*2+metadata$mapped_unpaired_reads+metadata$unmapped_reads
metadata$mapped_read_pairs <- metadata$mapped_read_pairs/1e6

metadata$unmapped_proportion <- metadata$unmapped_reads/metadata$mapped_examined_reads
metadata$duplication <- mkDup_metrics$PERCENT_DUPLICATION[match(metadata$rfid, mkDup_metrics$rfid)]

########## uniquelly mapped reads##########
uniq_mapped <- read.table(args[3], header=TRUE,  strip.white=TRUE,sep="\t")
metadata$uniq_mapped <- uniq_mapped$uniq_mapped[match(metadata$rfid, uniq_mapped$rfid)]
metadata$uniq_mapped_proportion <- metadata$uniq_mapped/metadata$mapped_examined_reads

#### plot for # of mapped read pairs per sample scatter plot
ggplot(metadata,aes(x=index,y=mapped_read_pairs,shape=runid,color=library_name)) +
    geom_point(size=2) +
    xlab("Samples") +
    ylab("# of mapped read pairs (million)") #+ theme(legend.position="none")
ggsave(paste0(args[4],"/mapped_read_pairs_scatter.png"), width = 7, height = 5, dpi = 300, units = "in")

#### plot for # of mapped read pairs per sample box plot
means <- aggregate(mapped_read_pairs ~ library_name, metadata, mean)
ggplot(metadata,aes(x=library_name,y=mapped_read_pairs)) +
    geom_boxplot() +
    stat_summary(fun="mean", geom="point", col="red", shape=1, size=2) +
    geom_text(data=means, aes(label=comma(mapped_read_pairs), y=mapped_read_pairs+0.25), size=3) +
    theme(text=element_text(size=8), axis.text.x=element_text(angle=-45, size=8)) +
    xlab("Library names") +
    ylab("# of mapped read pairs (million)") #+ theme(legend.position="none")
ggsave(paste0(args[4],"/mapped_read_pairs_box.png"), width = 5, height = 5, dpi = 300, units = "in")

#### plot for proportion of unmapped reads per sample scatter plot
ggplot(metadata,aes(x=index,y=unmapped_reads,shape=runid,color=library_name)) +
    geom_point(size=2) +
    xlab("Samples") +
    ylab("# of unmapped reads") #+ theme(legend.position="none")
ggsave(paste0(args[4],"/unmapped_scatter.png"), width = 7, height = 5, dpi = 300, units = "in")

#### plot for proportion of unmapped reads per sample box plot
means <- aggregate(unmapped_proportion ~ library_name, metadata, mean)
ggplot(metadata,aes(x=library_name,y=unmapped_proportion)) +
    geom_boxplot() +
    stat_summary(fun="mean", geom="point", col="red", shape=1, size=2) +
    geom_text(data=means, aes(label=percent(unmapped_proportion), y=unmapped_proportion+0.0002), size=3) +
    theme(text=element_text(size=8), axis.text.x=element_text(angle=-45, size=8)) +
    scale_y_continuous(labels=scales::percent_format(accuracy=0.1)) +
    xlab("Library names") +
    ylab("Unmapped / Examined") #+ theme(legend.position="none")
ggsave(paste0(args[4],"/unmapped_box.png"), width = 5, height = 5, dpi = 300, units = "in")

#### plot for duplication scatter plot
ggplot(metadata,aes(x=index,y=duplication,shape=runid,color=library_name)) +
    geom_point(size=2) +
    xlab("Samples") +
    ylab("Duplication") #+ theme(legend.position="none")
ggsave(paste0(args[4],"/duplication_scatter.png"), width = 7, height = 5, dpi = 300, units = "in")

#### plot for duplication box plot
means <- aggregate(duplication ~ library_name, metadata, mean)
ggplot(metadata,aes(x=library_name,y=duplication)) +
    geom_boxplot() +
    stat_summary(fun="mean", geom="point", col="red", shape=1, size=2) +
    geom_text(data=means, aes(label=percent(duplication), y=duplication+0.0015), size=3) +
    theme(text=element_text(size=8), axis.text.x=element_text(angle=-45, size=8)) +
    scale_y_continuous(labels=scales::percent_format(accuracy=0.1)) +
    xlab("library names") +
    ylab("Duplication") #+ theme(legend.position="none")
ggsave(paste0(args[4],"/duplication_box.png"), width = 5, height = 5, dpi = 300, units = "in")

#### plot for uniq_mapped scatter plot
ggplot(metadata,aes(x=index,y=uniq_mapped,shape=runid,color=library_name)) +
    geom_point(size=2) +
    scale_y_continuous(labels=scales::percent_format(accuracy=0.1)) +
    xlab("Samples") +
    ylab("Uniquely mapped") #+ theme(legend.position="none")
ggsave(paste0(args[4],"/uniq_mapped_scatter.png"), width = 7, height = 5, dpi = 300, units = "in")

#### plot for uniq_mapped box plot
means <- aggregate(uniq_mapped_proportion ~ library_name, metadata, mean)
ggplot(metadata,aes(x=library_name,y=uniq_mapped_proportion)) +
    geom_boxplot() +
    stat_summary(fun="mean", geom="point", col="red", shape=1, size=2) +
    geom_text(data=means, aes(label=percent(uniq_mapped_proportion), y=uniq_mapped_proportion+0.0015), size=3) +
    theme(text=element_text(size=8), axis.text.x=element_text(angle=-45, size=8)) +
    scale_y_continuous(labels=scales::percent_format(accuracy=0.1)) +
    xlab("Library names") +
    ylab("Uniquely mapped / Examined") #+ theme(legend.position="none")
ggsave(paste0(args[4],"/uniq_mapped_box.png"), width = 5, height = 5, dpi = 300, units = "in")

#### Output mapped_read_pairs < 1M as outliers 
out_csv <- subset(metadata, select = c(rfid, mapped_read_pairs))
out_csv$QC_reads <- "pass"
out_csv$QC_reads[out_csv$heterozygosity >= 1] <- "reject"
write.table(data.frame(out_csv), paste0(args[4], '/mapped_reads_1M_outliers.csv'), row.names = FALSE, sep=',' )