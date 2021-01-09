library(data.table)
library(ggplot2)
args = commandArgs(trailingOnly=TRUE)

data_path <- args[1]
metadata_path <- args[2]
out_path <- args[3]

metadata1 <- read.table(metadata_path, header=TRUE,  strip.white=TRUE,sep=",")
metadata1 <- metadata1[metadata1$strain=="Heterogenous stock",]
metadata1$library_project <- paste(metadata1$library_name,metadata1$project_name, sep="-")
metadata <- metadata1[order(metadata1$library_project),]
metadata$index <- 1:nrow(metadata)
### read in chr reads number
metadata$mapped <- 0
for(ch in c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", 
            "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15",
            "chr16", "chr17", "chr18", "chr19", "chr20", "chrX", "chrY")){
    fraction_chr <- read.table(paste0(data_path, '/mapped_', ch), header=FALSE,  strip.white=TRUE,sep=" ")
    metadata[,paste0(ch)] <- fraction_chr$V2[match(metadata$rfid, fraction_chr$V1)]
    metadata$mapped <- metadata$mapped+metadata[,paste0(ch)]
}

max_val_reads = 0
min_val_reads = 0
max_val_fra = 0
min_val_fra = 0

### calculate fraction and max, min
for(ch in c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", 
            "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15",
            "chr16", "chr17", "chr18", "chr19", "chr20", "chrX", "chrY")){
    metadata[,paste0(ch, "_fra")] <- metadata[,paste0(ch)]/metadata$mapped
    if(max(metadata[,paste0(ch, "_fra")]) > max_val_fra){
        max_val_fra = max(metadata[,paste0(ch, "_fra")])
    }
    if(min(metadata[,paste0(ch, "_fra")]) < min_val_fra){
        min_val_fra = min(metadata[,paste0(ch, "_fra")])
    }
    metadata[,paste0(ch)] <- metadata[,paste0(ch)]/1e+6
    if(max(metadata[,paste0(ch)]) > max_val_reads){
        max_val_reads = max(metadata[,paste0(ch)])
    }
    if(min(metadata[,paste0(ch)]) < min_val_reads){
        min_val_reads = min(metadata[,paste0(ch)])
    }
}

### plot number of mapped reads on each chromosome on one plot
library(reshape2)
metadata2 <- melt(metadata,id.vars='library_name', measure.vars=c("chr1",
            "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", 
            "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15",
            "chr16", "chr17", "chr18", "chr19", "chr20", "chrX", "chrY"))
metadata2$chr <- metadata2$variable
ggplot(metadata2) +
    geom_boxplot(aes(x=chr, y=value, color=library_name))+
    ylab("# of mapped reads (million)")
ggsave(paste0(out_path, "/mapped_chrs_lib_reads_box.png"), width = 12, height = 5, dpi = 300, units = "in")

ggplot(metadata2) +
    geom_boxplot(aes(x=library_name, y=value, color=chr))+
    ylab("# of mapped reads (million)")
ggsave(paste0(out_path, "/mapped_lib_chrs_reads_box.png"), width = 12, height = 5, dpi = 300, units = "in")

### plot ratio of mapped reads on each chromosome on one plot
library(reshape2)
metadata2 <- melt(metadata,id.vars='library_name', measure.vars=c("chr1_fra", 
            "chr2_fra", "chr3_fra", "chr4_fra", "chr5_fra", "chr6_fra", "chr7_fra", "chr8_fra", 
            "chr9_fra", "chr10_fra", "chr11_fra", "chr12_fra", "chr13_fra", "chr14_fra", "chr15_fra",
            "chr16_fra", "chr17_fra", "chr18_fra", "chr19_fra", "chr20_fra", "chrX_fra", "chrY_fra"))
metadata2$chr <- metadata2$variable
metadata2$chr <- c("chr1",
            "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", 
            "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15",
            "chr16", "chr17", "chr18", "chr19", "chr20", "chrX", "chrY")
ggplot(metadata2) +
    geom_boxplot(aes(x=chr, y=value, color=library_name))+
    scale_x_discrete(labels=c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", 
            "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15",
            "chr16", "chr17", "chr18", "chr19", "chr20", "chrX", "chrY"))+
    ylab("ratio of mapped reads")
ggsave(paste0(out_path, "/fraction_mapeed_chrs_lib_box.png"), width = 12, height = 5, dpi = 300, units = "in")

ggplot(metadata2) +
    geom_boxplot(aes(x=library_name, y=value, color=chr))+
    ylab("ratio of mapped reads")
ggsave(paste0(out_path, "/fraction_mapeed_lib_chrs_box.png"), width = 12, height = 5, dpi = 300, units = "in")

### plot fraction of mapped reads on chromosome X and Y
ggplot(metadata,aes(x=metadata[,"chrX_fra"],y=metadata[,"chrY_fra"],shape=library_name,color=sex)) +
    geom_point(size=2) +
    xlab("Ratio of mapped reads on chrX") +
    ylab("Ratio of mapped reads on chrY") #+ theme(legend.position="none")
ggsave(paste0(out_path, "/fraction_mapeed_on_chrY_chrX_scatter.png"), width = 6, height = 5, dpi = 300, units = "in")

ggplot(metadata,aes(x=index,y=metadata[,"chrX_fra"],shape=library_name,color=sex)) +
    geom_point(size=2) +
    xlab("Sample") +
    ylab("Ratio of mapped reads on chrX") #+ theme(legend.position="none")
ggsave(paste0(out_path, "/fraction_mapeed_on_chrX_sex_filter_scatter.png"), width = 6, height = 5, dpi = 300, units = "in")

ggplot(metadata,aes(x=index,y=metadata[,"chrY_fra"],shape=library_name,color=sex)) +
    geom_point(size=2) +
    xlab("Sample") +
    ylab("Ratio of mapped reads on chrY") #+ theme(legend.position="none")
ggsave(paste0(out_path, "/fraction_mapeed_on_chrY_sex_filter_scatter.png"), width = 6, height = 5, dpi = 300, units = "in")
