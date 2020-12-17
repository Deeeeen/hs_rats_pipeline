library(data.table)
library(ggplot2)

# read in arguments
args <- commandArgs(TRUE) 
missing_rate <- read.table(args[1], strip.white=TRUE, stringsAsFactors=FALSE, header=FALSE)
het <- read.table(args[2], strip.white=TRUE, stringsAsFactors=FALSE, header=TRUE)
out_path <- args[3]

# read in missing rate
colnames(missing_rate)<-c("FID", "IID", "MISSING_CT", "OBS_CT", "F_MISS")
missing_rate$logF_MISS = log10(missing_rate[, 5])

#read in het
het$meanHet = (het$N.NM. - het$O.HOM.)/het$N.NM
merged <- merge(missing_rate, het, by.x = "FID", by.y = "FID")

# plot
colors <- densCols(merged$logF_MISS, merged$meanHet)
jpeg(paste0(out_path,'/het_vs_missing.png'))
plot(merged$logF_MISS, merged$meanHet, col = colors, xlim = c(-3, 0), ylim = c(0, 
    0.5), pch = 20, xlab = "Proportion of missing genotypes", ylab = "Heterozygosity rate", 
    axes = F)
axis(2, at = c(0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5), tick = T)
axis(1, at = c(-3, -2, -1, 0), labels = c(0.001, 0.01, 0.1, 1))
abline(h = mean(merged$meanHet) - (2 * sd(merged$meanHet)), col = "RED", lty = 2)
abline(h = mean(merged$meanHet) + (2 * sd(merged$meanHet)), col = "RED", lty = 2)
abline(v = log10(0.03), col = "RED", lty = 2)
abline(v = log10(0.07), col = "RED", lty = 2)
dev.off()
