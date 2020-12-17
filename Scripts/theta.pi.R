# Libraries
library(PopGenome)
library(data.table)
library(foreach)
library(tidyverse)
library(SeqArray)
library(cowplot)

setwd("C:/Users/Conno/Desktop/Fall2020/Berglandlab/")

### Theta pi
data <- data.table(read.delim(file = "maf_filter_0.25_pi_windowed.pi"))

ggplot(data, aes(x=BIN_START, y=PI, color=CHROM)) +
   geom_line() +
   facet_wrap(~CHROM) +
   theme_classic() +
   labs(x="Position (mb)", y="Theta pi") +
   theme(legend.position="none",
         axis.text.x = element_text(face="bold", size=8),
         axis.text.y = element_text(face="bold", size=8),
         axis.title.x = element_text(face="bold", size=12),
         axis.title.y = element_text(face="bold", size=12))

data %>% group_by(CHROM) %>% summarise(pi=mean(PI))
data %>% group_by(CHROM) %>% summarise(pi=median(PI))

# Ne
data %>% summarise(pi=mean(PI, na.rm = TRUE)) / (4*(5.69e-09)) / 20
data %>% summarise(pi=median(PI, na.rm = TRUE)) / (4*(5.69e-09)) / 20

### Tajima's D
data <- data.table(read.delim(file = "maf_filter_0.25_tajima"))

ggplot(data, aes(x=BIN_START, y=TajimaD, color=CHROM)) +
   geom_line() +
   facet_wrap(~CHROM) +
   theme_classic() +
   labs(x="Position (mb)", y="Tajima's D") +
   theme(legend.position="none",
         axis.text.x = element_text(face="bold", size=8),
         axis.text.y = element_text(face="bold", size=8),
         axis.title.x = element_text(face="bold", size=12),
         axis.title.y = element_text(face="bold", size=12))

data %>% group_by(CHROM) %>% summarise(tajima=mean(TajimaD, na.rm = TRUE))
data %>% group_by(CHROM) %>% summarise(tajima=median(TajimaD, na.rm = TRUE))
