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
    GC_chr <- read.table(paste0(data_path, '/mapped_GC_', ch), header=FALSE,  strip.white=TRUE,sep=" ")
    GC_chr$fraction_GC <- (GC_chr$V2+GC_chr$V3)/GC_chr$V4
    metadata[,paste0(ch)] <- GC_chr$fraction_GC[match(metadata$rfid, GC_chr$V1)]
}

### plot number of mapped reads on each chromosome on one plot
library(reshape2)
metadata2 <- melt(metadata,id.vars='library_name', measure.vars=c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", 
            "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15",
            "chr16", "chr17", "chr18", "chr19", "chr20", "chrX", "chrY"))
metadata2$chr <- metadata2$variable
ggplot(metadata2) +
    geom_boxplot(aes(x=chr, y=value, color=library_name))+
    ylab("GC content")
ggsave(paste0(out_path, "/mapped_GC_chrs_lib_reads_box.png"), width = 12, height = 5, dpi = 300, units = "in")

ggplot(metadata2) +
    geom_boxplot(aes(x=library_name, y=value, color=chr))+
    ylab("GC content")
ggsave(paste0(out_path, "/mapped_GC_lib_chrs_reads_box.png"), width = 12, height = 5, dpi = 300, units = "in")
