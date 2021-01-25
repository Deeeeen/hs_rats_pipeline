library(data.table)

args <- commandArgs(TRUE)
# read in first argument 
dir_path <- args[1]
out_path <- args[2]

########## demux stats ##########
demux_metrics_f <- list.files(path=paste0(dir_path, "/demux/metrics"))
demux_metrics_ls = list()
for(i in 1:length(demux_metrics_f)){
    temp <- read.table(paste0(dir_path, "/demux/metrics/", demux_metrics_f[i]), header=TRUE, strip.white=TRUE, sep="\t", stringsAsFactors = FALSE)
    demux_metrics_ls[[i]] <- temp[temp$barcode_name != "unmatched", ]
}
demux_metrics <- do.call(rbind, demux_metrics_ls)
genotype_log <- subset(demux_metrics, select = c("barcode_name", "templates"))
colnames(genotype_log) <- c("rfid", "demux_reads")

########## mapping stats ##########
mkDup_metrics <- read.table(paste0(dir_path, "/results/mapping_result/mkDup_metrics"), header=TRUE, strip.white=TRUE, sep="\t", stringsAsFactors = FALSE)
# mapped reads
mkDup_metrics$mapped_read_pairs <- mkDup_metrics$READ_PAIRS_EXAMINED
# unmapped ratio
mkDup_metrics$mapped_examined_reads <- mkDup_metrics$READ_PAIRS_EXAMINED*2+mkDup_metrics$UNPAIRED_READS_EXAMINED+mkDup_metrics$UNMAPPED_READS
mkDup_metrics$unmapped_reads_ratio <- mkDup_metrics$UNMAPPED_READS/mkDup_metrics$mapped_examined_reads
# duplication ratio
mkDup_metrics$duplication_ratio <- mkDup_metrics$PERCENT_DUPLICATION
# uniquely mapped ratio
uniq_mapped <- read.table(paste0(dir_path, "/results/mapping_result/uniq_mapped"), header=TRUE, strip.white=TRUE, sep="\t", stringsAsFactors = FALSE)
mkDup_metrics$uniq_mapped <- uniq_mapped$uniq_mapped[match(mkDup_metrics$rfid, uniq_mapped$rfid)]
mkDup_metrics$uniq_mapped_ratio <- mkDup_metrics$uniq_mapped/mkDup_metrics$mapped_examined_reads
# subset mkDup_metrics
mkDup_metrics <- subset(mkDup_metrics, select = c("rfid", "mapped_read_pairs", "unmapped_reads_ratio", "duplication_ratio", "uniq_mapped_ratio"))
# mapped reads QC
mapped_reads_outliers <- read.table(paste0(dir_path, "/results/mapping_result/mapped_reads_1M_outliers.csv"),
                                    stringsAsFactors=FALSE, header=TRUE,  strip.white=TRUE,sep=",")
mapped_reads_outliers <- subset(mapped_reads_outliers, select = c(rfid, QC_reads))
mkDup_metrics <- merge(mkDup_metrics, mapped_reads_outliers, by.x = "rfid", by.y = "rfid", all = TRUE)
# merge genotype log
genotype_log <- merge(genotype_log, mkDup_metrics, by.x = "rfid", by.y = "rfid", all = TRUE)

########## sex QC ##########
sex_outliers_rfid <- c(933000320047998, 933000320187340)
# merge genotype log
genotype_log$QC_sex <- "pass"
genotype_log$QC_sex[genotype_log$rfid %in% sex_outliers_rfid] <- "reject"

########## combine genotype_log ##########
# read in second argument if there is any
library(plyr)
if(length(args)==4){
    previous_genotype_log <- read.table(args[4], header=TRUE, strip.white=TRUE, sep=",", stringsAsFactors = FALSE)
    previous_genotype_log <- subset(previous_genotype_log, select = colnames(genotype_log))
    genotype_log <- merge(previous_genotype_log, genotype_log, all = TRUE)
}

########## add date and file_location ##########
genotype_log$date <- format(Sys.Date(), "%m%d%Y")
genotype_log$file_location <- args[3]

########## missing rate QC ##########
missing_rate_outliers <- read.table(paste0(dir_path, "/results/genotype_result/stitch_result/plink/missing_rate_outliers.csv"),
                                    stringsAsFactors=FALSE, header=TRUE,  strip.white=TRUE,sep=",")
colnames(missing_rate_outliers) <- c("rfid", "missing_rate", "QC_missing")
genotype_log <- merge(genotype_log, missing_rate_outliers, by.x = "rfid", by.y = "rfid")

########## heterozygosity rate QC ##########
het_outliers <- read.table(paste0(dir_path, "/results/genotype_result/stitch_result/plink/het_rate_outliers.csv"),
                            stringsAsFactors=FALSE, header=TRUE,  strip.white=TRUE,sep=",")
colnames(het_outliers) <- c("rfid", "heterozygosity_rate", "QC_heterozygosity")
genotype_log <- merge(genotype_log, het_outliers, by.x = "rfid", by.y = "rfid", all = TRUE)

########## albino coat color QC ##########
albino_outliers <- read.table(paste0(dir_path, "/results/genotype_result/beagle_result/plink/coat_color_albino_outliers.csv"),
                                stringsAsFactors=FALSE, header=TRUE,  strip.white=TRUE,sep=",")
albino_outliers <- subset(albino_outliers, select = c("rfid", "coatcolor", "QC_coat_color_albino"))
genotype_log <- merge(genotype_log, albino_outliers, by.x = "rfid", by.y = "rfid", all = TRUE)


########## brown coat color QC ##########
outliers <- genotype_log[,c("rfid","coatcolor",colnames(genotype_log)[grepl("QC_",colnames(genotype_log))])]
outliers$QC_pass <- apply(outliers[,colnames(outliers)[grepl("QC_",colnames(outliers))]] != "pass",1,any)
outliers$QC_pass <- ifelse(outliers$QC_pass, "reject", "pass")
write.table(data.frame(outliers), "temp", row.names = FALSE, sep=',' )
brown_GTs <- read.table(paste0(dir_path, "/results/genotype_result/beagle_result/plink/coat_color_brown_GT.csv"),
                                stringsAsFactors=FALSE, header=TRUE,  strip.white=TRUE,sep=",")
flagged_by_other<-brown_GTs[which(brown_GTs$rfid %in% outliers$rfid[outliers$QC_pass=="reject"]),]
brown_flagged_by_other<-flagged_by_other[which(flagged_by_other$coatcolor=="BROWN" | flagged_by_other$coatcolor=="BROWNHOOD"),]
​write.csv(brown_flagged_by_other, paste0(dir_path, "/results/genotype_result/beagle_result/plink/hs_rats_n", length(unique(genotype_log$rfid)), '_', format(Sys.Date(), "%m%d%Y"), "_Brown_flagged_by_other.csv"),row.names=FALSE,quote = F)

library(dplyr)
​
GT_counts<-flagged_by_other %>%
  count(coatcolor,GT_150285633,GT_150288295,GT_150432118,GT_150449245,GT_150488934,GT_150530733,GT_150584522,sort = TRUE)​
GT_counts<-GT_counts[which(!GT_counts$coatcolor=="ALBINO"),]
write.csv(GT_counts[which(GT_counts$coatcolor=="BROWN" | GT_counts$coatcolor=="BROWNHOOD"),], paste0(dir_path, "/results/genotype_result/beagle_result/plink/Brown_flagged_by_other.csv"),row.names=FALSE,quote = F)
​write.csv(GT_counts[which(!(GT_counts$coatcolor=="BROWN") | (!GT_counts$coatcolor=="BROWNHOOD")),], paste0(dir_path, "/results/genotype_result/beagle_result/plink/Not_Brown_flagged_by_other.csv"),row.names=FALSE,quote = F)


########## output genotype_log ##########
outliers <- subset(outliers, select = c("rfid", "QC_pass"))
genotype_log <- merge(genotype_log, outliers, by.x = "rfid", by.y = "rfid", all = TRUE)
write.table(data.frame(genotype_log), paste0(out_path, '/hs_rats_n', length(unique(genotype_log$rfid)), '_', format(Sys.Date(), "%m%d%Y"), '_genotype_log.csv'), row.names = FALSE, sep=',',quote = F)