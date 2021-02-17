library(data.table)
library(ggplot2)
# read in arguments
args <- commandArgs(TRUE) 
eigenvec <- read.table(args[1], strip.white=TRUE, stringsAsFactors=FALSE, header=FALSE)
eigenval <- read.table(args[2], strip.white=TRUE, stringsAsFactors=FALSE, header=FALSE)
metadata <- read.table(args[3], strip.white=TRUE, stringsAsFactors=FALSE, header=TRUE,sep=",")
out_path <- args[4]

# merge data
eigenvec$rfid<-eigenvec$V1
eigenvec$sex <- metadata$sex[match(eigenvec$V1, metadata$rfid)]
eigenvec$library_name <- metadata$library_name[match(eigenvec$V1, metadata$rfid)]
eigenvec$project_name <- metadata$project_name[match(eigenvec$V1, metadata$rfid)]
metadata$family_id <- paste(metadata$dames,metadata$sires, sep="-")
eigenvec$family_id <- metadata$family_id[match(eigenvec$V1, metadata$rfid)]

# plot PCA plots
ggplot(eigenvec,aes(x=V3,y=V4,color=library_name,shape=sex)) +
    geom_point(size=2) +
    xlab(paste0("PC1 (", sprintf("%.2f%%",eigenval$V1[1]/sum(eigenval$V1)*100), ")")) +
    ylab(paste0("PC2 (", sprintf("%.2f%%",eigenval$V1[2]/sum(eigenval$V1)*100), ")")) 
    # +theme(legend.position="none")
ggsave(paste0(out_path,"/pc1_pc2_library_name.png"), width = 6, height = 5, dpi = 300, units = "in")

ggplot(eigenvec,aes(x=V3,y=V4,color=project_name,shape=sex)) +
    geom_point(size=2) +
    xlab(paste0("PC1 (", sprintf("%.2f%%",eigenval$V1[1]/sum(eigenval$V1)*100), ")")) +
    ylab(paste0("PC2 (", sprintf("%.2f%%",eigenval$V1[2]/sum(eigenval$V1)*100), ")")) 
    # +theme(legend.position="none")
ggsave(paste0(out_path,"/pc1_pc2_project_name.png"), width = 7, height = 5, dpi = 300, units = "in")

ggplot(eigenvec,aes(x=V3,y=V4,color=family_id,shape=sex)) +
    geom_point(size=2) +
    xlab(paste0("PC1 (", sprintf("%.2f%%",eigenval$V1[1]/sum(eigenval$V1)*100), ")")) +
    ylab(paste0("PC2 (", sprintf("%.2f%%",eigenval$V1[2]/sum(eigenval$V1)*100), ")")) +
    theme(legend.position="none")
ggsave(paste0(out_path,"/pc1_pc2_family_id.png"), width = 5, height = 5, dpi = 300, units = "in")

ggplot(eigenvec,aes(x=V3,y=V5,color=family_id,shape=sex)) +
    geom_point(size=2) +
    xlab(paste0("PC1 (", sprintf("%.2f%%",eigenval$V1[1]/sum(eigenval$V1)*100), ")")) +
    ylab(paste0("PC3 (", sprintf("%.2f%%",eigenval$V1[3]/sum(eigenval$V1)*100), ")")) +
    theme(legend.position="none")
ggsave(paste0(out_path,"/pc1_pc3_family_id.png"), width = 5, height = 5, dpi = 300, units = "in")

ggplot(eigenvec,aes(x=V3,y=V6,color=family_id,shape=sex)) +
    geom_point(size=2) +
    xlab(paste0("PC1 (", sprintf("%.2f%%",eigenval$V1[1]/sum(eigenval$V1)*100), ")")) +
    ylab(paste0("PC4 (", sprintf("%.2f%%",eigenval$V1[4]/sum(eigenval$V1)*100), ")")) +
    theme(legend.position="none")
ggsave(paste0(out_path,"/pc1_pc4_family_id.png"), width = 5, height = 5, dpi = 300, units = "in")

