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

########## RN6 markDuplicates metrics##########
mkDup_metrics_rn6 <- read.table(args[2], header=TRUE,  strip.white=TRUE,sep="\t")
metadata$mapped_read_pairs_rn6 <- mkDup_metrics_rn6$READ_PAIRS_EXAMINED[match(metadata$rfid, mkDup_metrics_rn6$rfid)]
metadata$mapped_unpaired_reads_rn6 <- mkDup_metrics_rn6$UNPAIRED_READS_EXAMINED[match(metadata$rfid, mkDup_metrics_rn6$rfid)]
metadata$unmapped_reads_rn6 <- mkDup_metrics_rn6$UNMAPPED_READS[match(metadata$rfid, mkDup_metrics_rn6$rfid)]
metadata$mapped_examined_reads_rn6 <- metadata$mapped_read_pairs_rn6*2+metadata$mapped_unpaired_reads_rn6+metadata$unmapped_reads_rn6
metadata$mapped_read_pairs_rn6 <- metadata$mapped_read_pairs_rn6/1e6

metadata$unmapped_proportion_rn6 <- metadata$unmapped_reads_rn6/metadata$mapped_examined_reads_rn6
metadata$duplication_rn6 <- mkDup_metrics_rn6$PERCENT_DUPLICATION[match(metadata$rfid, mkDup_metrics_rn6$rfid)]

########## RN6 uniquelly mapped reads##########
uniq_mapped_rn6 <- read.table(args[3], header=TRUE,  strip.white=TRUE,sep="\t")
metadata$uniq_mapped_rn6 <- uniq_mapped_rn6$uniq_mapped[match(metadata$rfid, uniq_mapped_rn6$rfid)]
metadata$uniq_mapped_proportion_rn6 <- metadata$uniq_mapped_rn6/metadata$mapped_examined_reads_rn6

#### plot for # of mapped read pairs per sample scatter plot
ggplot(metadata,aes(x=index,y=mapped_read_pairs_rn6,shape=runid,color=library_name)) +
    geom_point(size=2) +
    xlab("Samples") +
    ylab("# of mapped read pairs rn6 (million)") #+ theme(legend.position="none")
ggsave(paste0(args[4],"/mapped_read_pairs_scatter.png"))

#### plot for # of mapped read pairs per sample box plot
means <- aggregate(mapped_read_pairs_rn6 ~ library_name, metadata, mean)
ggplot(metadata,aes(x=library_name,y=mapped_read_pairs_rn6)) +
    geom_boxplot() +
    stat_summary(fun="mean", geom="point", col="red", shape=1, size=2) +
    geom_text(data=means, aes(label=comma(mapped_read_pairs_rn6), y=mapped_read_pairs_rn6+0.25), size=3) +
    theme(text=element_text(size=8), axis.text.x=element_text(angle=-45, size=8)) +
    xlab("Library names") +
    ylab("# of mapped read pairs rn6 (million)") #+ theme(legend.position="none")
ggsave(paste0(args[4],"/mapped_read_pairs_box.png"))

#### plot for proportion of unmapped reads per sample scatter plot
ggplot(metadata,aes(x=index,y=unmapped_reads_rn6,shape=runid,color=library_name)) +
    geom_point(size=2) +
    xlab("Samples") +
    ylab("# of unmapped reads (rn6)") #+ theme(legend.position="none")
ggsave(paste0(args[4],"/unmapped_scatter.png"))

#### plot for proportion of unmapped reads per sample box plot
means <- aggregate(unmapped_proportion_rn6 ~ library_name, metadata, mean)
ggplot(metadata,aes(x=library_name,y=unmapped_proportion_rn6)) +
    geom_boxplot() +
    stat_summary(fun="mean", geom="point", col="red", shape=1, size=2) +
    geom_text(data=means, aes(label=comma(unmapped_proportion_rn6), y=unmapped_proportion_rn6+0.0002), size=3) +
    theme(text=element_text(size=8), axis.text.x=element_text(angle=-45, size=8)) +
    xlab("Library names") +
    ylab("Unmapped / Examined (rn6)") #+ theme(legend.position="none")
ggsave(paste0(args[4],"/unmapped_box.png"))

#### plot for duplication scatter plot
ggplot(metadata,aes(x=index,y=duplication_rn6,shape=runid,color=library_name)) +
    geom_point(size=2) +
    xlab("Samples") +
    ylab("Duplication (rn6)") #+ theme(legend.position="none")
ggsave(paste0(args[4],"/duplication_scatter.png"))

#### plot for duplication box plot
means <- aggregate(duplication_rn6 ~ library_name, metadata, mean)
ggplot(metadata,aes(x=library_name,y=duplication_rn6)) +
    geom_boxplot() +
    stat_summary(fun="mean", geom="point", col="red", shape=1, size=2) +
    geom_text(data=means, aes(label=comma(duplication_rn6), y=duplication_rn6+0.0007), size=3) +
    theme(text=element_text(size=8), axis.text.x=element_text(angle=-45, size=8)) +
    xlab("library names") +
    ylab("Duplication (rn6)") #+ theme(legend.position="none")
ggsave(paste0(args[4],"/duplication_box.png"))

#### plot for uniq_mapped scatter plot
ggplot(metadata,aes(x=index,y=uniq_mapped_rn6,shape=runid,color=library_name)) +
    geom_point(size=2) +
    xlab("Samples") +
    ylab("Uniquely mapped (rn6)") #+ theme(legend.position="none")
ggsave(paste0(args[4],"/uniq_mapped_scatter.png"))

#### plot for uniq_mapped box plot
means <- aggregate(uniq_mapped_proportion_rn6 ~ library_name, metadata, mean)
ggplot(metadata,aes(x=library_name,y=uniq_mapped_proportion_rn6)) +
    geom_boxplot() +
    stat_summary(fun="mean", geom="point", col="red", shape=1, size=2) +
    geom_text(data=means, aes(label=comma(uniq_mapped_proportion_rn6), y=uniq_mapped_proportion_rn6+0.001), size=3) +
    theme(text=element_text(size=8), axis.text.x=element_text(angle=-45, size=8)) +
    xlab("Library names") +
    ylab("Uniquely mapped / Examined (rn6)") #+ theme(legend.position="none")
ggsave(paste0(args[4],"/uniq_mapped_box.png"))
