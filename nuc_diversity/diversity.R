# Nucleotide diversity from VCF
# 12.14.2020

# Libraries
library(data.table)
library(tidyverse)
library(foreach)

setwd("C:/Users/Conno/Desktop/Fall2020/Berglandlab/PopGenome")

# Chromosomes
chrom <- c('Scaffold_1863_HRSCAF_2081', 'Scaffold_1931_HRSCAF_2197', 'Scaffold_2158_HRSCAF_2565', 
           'Scaffold_2217_HRSCAF_2652', 'Scaffold_2373_HRSCAF_2879', 'Scaffold_6786_HRSCAF_7541', 
           'Scaffold_7757_HRSCAF_8726', 'Scaffold_9197_HRSCAF_10753', 'Scaffold_9198_HRSCAF_10754', 
           'Scaffold_9199_HRSCAF_10755', 'Scaffold_9200_HRSCAF_10757', 'Scaffold_9201_HRSCAF_10758')

# Empirical diversity data
data <- data.table(read.csv("MapJune2020.recode.finalsnpset.phasedBCF.wildtrust.csv", header = T))

# Diversity long format
data.w <- data.table(data %>% select(colnames(data[,c(20:23)]), chr) %>%
   pivot_longer(colnames(data[,c(20:22)]), names_to = "data", values_to = "stat"))
data.w[, winMBstop:=winMB+250000]
data.w[,region:=c(1:nrow(data.w))]

# Peaks
load("gprime_peaks.replicates.250K.05.Rdata")
colnames(peaks)[1] <- "chr"

# Find regions
foreach(i = 1:length(peaks$chr)) %do% {
   peaks[i]$region <- data.w[chr==peaks[i]$chr][which(data.w[chr==peaks[i]$chr]$winMB %in% 
                   c(peaks[i]$start:(peaks[i]$start+peaks[i]$length))==T)][1]$markernum
}

# Diversity and neutrality stats
ggplot() +
   geom_line(data=data.w[data=="tajima.D.win"], 
             aes(x=region, y=stat, color=chr, group=interaction(data, chr)), size=1) +
   geom_vline(data=peaks, 
              aes(xintercept = start, group=chr), size=1) +
   #facet_grid(rows = vars(data), cols = vars(chr), scales = "free") +
   labs(x="Mb", y="Tajima's D", color="") +
   scale_x_continuous(breaks = seq(0, 12e6, by = 2e6),
                      labels = c(0, 2, 4, 6, 8, 10, 12)) +
   theme_classic() + 
   theme(axis.text.x = element_text(face="bold", size=16),
         axis.text.y = element_text(face="bold", size=16),
         axis.title.x = element_text(face="bold", size=18),
         axis.title.y = element_text(face="bold", size=18),
         legend.position = "none") 
