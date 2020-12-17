# REHH - selection scans
# 12.15.2020

# Libraries
library(data.table)
library(foreach)
library(ggplot2)
library(cowplot)
library(rehh)

setwd("C:/Users/Conno/Desktop/Fall2020/Berglandlab/rehh/")

# Chromosomes
chrom <- c('Scaffold_1863_HRSCAF_2081', 'Scaffold_1931_HRSCAF_2197', 'Scaffold_2158_HRSCAF_2565', 
           'Scaffold_2217_HRSCAF_2652', 'Scaffold_2373_HRSCAF_2879', 'Scaffold_6786_HRSCAF_7541', 
           'Scaffold_7757_HRSCAF_8726', 'Scaffold_9197_HRSCAF_10753', 'Scaffold_9198_HRSCAF_10754', 
           'Scaffold_9199_HRSCAF_10755', 'Scaffold_9200_HRSCAF_10757', 'Scaffold_9201_HRSCAF_10758')

# Combine rehh output
y <- foreach(i = c(1:12), .errorhandling="remove") %do% {
    load(paste0("/project/berglandlab/connor/rehh/MapJune2020_ann.hyrbid_strategy.3species.whatshap.shapeit.wildtrust.", 
                chrom[i], ".ihs.Rdata"))
    return(data.table(scan.ihs$ihs))
}

y <- rbindlist(y)
y[,markernum:=c(1:nrow(y))]

# Write csv
# write.csv(y, "MapJune2020_ann.hyrbid_strategy.3species.whatshap.shapeit.ihs.wildtrust.csv", quote = F, row.names = F)

# IHS data
y <- data.table(read.csv("MapJune2020_ann.hyrbid_strategy.3species.whatshap.shapeit.ihs.csv", header = T))

# Peaks
load("gprime_peaks.replicates.250K.05.Rdata")
colnames(peaks)[1] <- "chr"

# Not normally distributed
distribplot(y$IHS, xlab = "iHS")
distribplot(y$IHS, xlab = "iHS", qqplot = TRUE)

peaks$region <- 1

# Find regions
foreach(i = 1:length(peaks$chr)) %do% {
    peaks[i]$region <- y[CHR==peaks[i]$chr][which(y[CHR==peaks[i]$chr]$POSITION %in% 
            c(peaks[i]$start:(peaks[i]$start+peaks[i]$length))==T)][1]$markernum
}

# iHS
ggplot(y, aes(x=markernum, y=IHS, color=CHR)) +
    geom_point() +
    geom_vline(data=peaks, aes(xintercept = region, group=chr), size=1.2, linetype=2) +
    geom_hline(yintercept = c(2,-2), linetype=2, size=1.2) +
    labs(x="", y="iHS", color="") +
    theme_classic() + 
    theme(axis.text.x = element_blank(),
          axis.text.y = element_text(face="bold", size=16),
          axis.title.x = element_text(face="bold", size=18),
          axis.title.y = element_text(face="bold", size=18),
          legend.position = "none") 
