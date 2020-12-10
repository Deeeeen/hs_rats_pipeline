library(data.table)
library(ggplot2)
dir <- '/Users/den/Desktop/PalmerLab/hs_rats/201002_A00953_0161_AHLGVKDSXY'
# read in missing rate
missing_rate <- read.table(paste0(dir, '/hs_rats_stitch_n1536_11042020_sample.smiss'), stringsAsFactors=FALSE,
                        header=FALSE,  strip.white=TRUE,sep="\t")
colnames(missing_rate)<-c("rfid", "IID", "MISSING_CT", "OBS_CT", "F_MISS")
# missing_rate$reads <- 0
missing_rate$reads <- as.integer(NA)

# 201002_A00953_0161_AHLGVKDSXY 
metrics_files<-list.files(paste0(dir, '/metrics'))
for(f in metrics_files){
    temp_metrics<-read.table(paste0(dir, '/metrics/', f), stringsAsFactors=FALSE,
                        header=TRUE,  strip.white=TRUE,sep="\t")
    setDT(missing_rate); setDT(temp_metrics)
    missing_rate[temp_metrics, reads := ifelse(is.na(reads), i.templates, reads), on=c("barcode_name")]
}

# 200915_A00953_0152_AHJGF7DSXY
dir <- '/Users/den/Desktop/PalmerLab/hs_rats/200915_A00953_0152_AHJGF7DSXY'
metrics_files<-list.files(paste0(dir, '/metrics'))
for(f in metrics_files){
    temp_metrics<-read.table(paste0(dir, '/metrics/', f), stringsAsFactors=FALSE,
                        header=TRUE,  strip.white=TRUE,sep="\t")
    setDT(missing_rate); setDT(temp_metrics)
    missing_rate[temp_metrics, reads := ifelse(is.na(reads), i.templates, reads), on=c("barcode_name")]
}

# 200221_A00953_0069_BH5T5LDSXY_and_200326_A00953_0086_BHC2FMDSXY
dir <- '/Users/den/Desktop/PalmerLab/hs_rats/200221_A00953_0069_BH5T5LDSXY_and_200326_A00953_0086_BHC2FMDSXY' 
metrics_files<-list.files(paste0(dir, '/metrics'))
for(f in metrics_files){
    temp_metrics<-read.table(paste0(dir, '/metrics/', f), stringsAsFactors=FALSE,
                        header=TRUE,  strip.white=TRUE,sep="\t")
    setDT(missing_rate); setDT(temp_metrics)
    missing_rate[temp_metrics, reads := ifelse(is.na(reads), i.templates, reads), on=c("barcode_name")]
}

# outliers<-c("933000120138682", "933000120138743", "933000320046375", "933000320047795",
#             "933000320046365", "933000320046491", "933000320046705", "933000320046473",
#             '933000320045989', "933000320045996", "933000320046414", "933000320047537",
#             "933000320046411", "933000320047370", "933000320045999", "933000320045991",
#             "1DCD1762", "933000120138409", "933000120138750", "933000320045791",
#             "933000320046039", "933000320046318", "933000320046510", "933000320046604",
#             "933000320046913", "933000320047011", "933000320047075", "933000320047278")
# missing_rate$QC_outliers <- c(missing_rate$barcode_name %in% outliers)

# missing_rate$reads <- missing_rate$reads/1000000
# p<-ggplot(missing_rate, aes(x=reads,y=F_MISS , color=QC_outliers)) +
#     geom_point(size=1,shape=1) +
#     xlab("Reads (million)") +
#     ylab("Missing Rate") 
#     # +theme(legend.position="none")
# print(p)
# ggsave("/Users/den/Desktop/PalmerLab/hs_rats/201002_A00953_0161_AHLGVKDSXY/results/missing_vs_reads.png")
# write.table( data.frame(temp), '/Users/den/Desktop/PalmerLab/hs_rats/201002_A00953_0161_AHLGVKDSXY/hs_missing_rate_outliers_n1536_20201130.csv', row.names = FALSE, sep=',' )
