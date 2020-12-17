# Using PopGenome to get neutrality and diversity statistics from Daphnia VCF file
# Connor Murray
# 4.28.2020

module load intel/18.0 intelmpi/18.0 R/3.6.0; R

# Libraries
library(PopGenome)
library(data.table)
library(foreach)
library(tidyverse)
library(doParallel)
library(tidyverse)

# Daphnia data 
setwd("/project/berglandlab/connor/PopGenome/")

# Use to verify sample ids
vcf_handle <-.Call("VCF_open", "vcf/MapDec19PulexandObtusaandPulicaria_filtsnps10bpindels_snps_filter_pass_lowGQmiss_ann.vcf.gz")
.Call("VCF_getSampleNames",vcf_handle)

# Get sample ids and clonal associations
samps <- fread("/project/berglandlab/Karen/MappingDec2019/WithPulicaria/Superclones201617182019withObtusaandPulicaria_kingcorr_20200402_wmedrd.txt")

# Good snps and chromosomes
snps <- fread("/project/berglandlab/Karen/MappingDec2019/WithPulicaria/finalsetsnpset01pulex_table_wpulicaria_20200401")
snps.ag <- snps[,list(.N), chr]
snps[,goodChr:=F]
snps[chr%in%snps.ag[N>5000]$chr, goodChr:=T]

# Chromosomes of interest
chrom <- unique(snps[goodChr == TRUE]$chr)

# Clonal lineages of interest - used to define populations
A <- c(samps[SC=="A" & Nonindependent==0]$clone)
B <- c(samps[SC=="B" & Nonindependent==0]$clone)
C <- c(samps[SC=="C" & Nonindependent==0]$clone)
E <- c(samps[SC=="E" & Nonindependent==0]$clone)
f <- c(samps[SC=="F" & Nonindependent==0]$clone)
O <- c(samps[SC=="OO" & Nonindependent==0]$clone)

# Reads in vcf with bootstrap sampling
vcf.fun <- function(sample.n, upper.maf) {
   
   # Goes through each chromosome
   foreach(chr.num = 1:12, .combine="rbind") %do% {
      sample.n = 15; upper.maf = 0.25; i = 2; chr.num = 4
      # sample.n = number of samples
      # upper.maf = upper boundary of minor allele frequency
      # i = bootstrap replicate
      # chr.num = chromosome number
      
      message(paste("Chromosome", chr.num, "replicate #", i, sep=" "))
   
      # Sample clones
      SC.A <- sample(A, size = sample.n)
      SC.B <- sample(B, size = sample.n)
      SC.C <- sample(C, size = sample.n)
      SC.E <- sample(E, size = sample.n)
      SC.F <- sample(f, size = sample.n)
      SC.O <- sample(O, size = sample.n)
   
      # Load tabindex vcf file - load in by chromosome
      vcf <- readVCF(filename = "vcf/MapDec19PulexandObtusaandPulicaria_filtsnps10bpindels_snps_filter_pass_lowGQmiss_ann.vcf.gz", 
         tid = c(chrom[chr.num]), 
         gffpath = "gff/Daphnia.aed.0.6.gff",
         samplenames = c(SC.A, SC.B, SC.C, SC.E, SC.F, SC.O),
         frompos = 1, 
         topos = max(snps[chr== chrom[chr.num]]$pos), 
         numcols = 10000, 
         include.unknown = TRUE,
         approx = FALSE)

      # Filter out sites based on MAF
      genofile <- set.filter(vcf, minor.freqs = TRUE, maf.lower.bound = 0, 
                             maf.upper.bound = upper.maf, missing.freqs = TRUE, 
                             miss.lower.bound = 0, miss.upper.bound = 0.20)
      
      # Populations within data
      genofile <- set.populations(genofile, list(SC.A, SC.B, SC.C, SC.E, SC.F, SC.O), 
                                 diploid = TRUE)

      # Neutrality & diversity statistics
      genofile <- neutrality.stats(genofile, subsites = "included")
      genofile <- diversity.stats(genofile, subsites = "included")
      genofile <- F_ST.stats(genofile, subsites = "included")
      
      # Combine output data
      total <- data.table(Tajima.D = t(genofile@Tajima.D), 
              tW = t(genofile@theta_Watterson), 
              tP = t(genofile@nuc.diversity.within)/(genofile@n.sites),
              eff = t(genofile@nuc.diversity.within)/(genofile@n.sites)/(4*5.69e-09),
              nuc.fst = genofile@nucleotide.F_ST,
              nuc.fst.pair = t(genofile@nuc.F_ST.vs.all),
              chr = chrom[chr.num],
              nSites = genofile@n.sites,
              biall = genofile@n.biallelic.sites, 
              tt.ratio = genofile@trans.transv.ratio,
              Sites.inc = table(genofile@region.data@included)[2],
              n.samp = sample.n,
              up.maf = upper.maf,
              SC = c("A", "B", "C", "E", "F", "OO"),
              rep = i,
              nuc.fst.all.pair = t(genofile@nuc.F_ST.pairwise))
      
      # Finish
      return(total)
   
   }
}

# i = Number of bootstraps; j = Number of chromosomes
# Forloop of vcf function
vcf.out <- foreach(i = 1:2, .combine="rbind") %do% {
                   message(paste("Sample #", i, sep=" "))
                   vcf.fun(sample.n = 15, upper.maf = 0.25)
}

saveRDS(vcf.out, file= "/project/berglandlab/connor/PopGenome/10rep-15samp-0.25maf-A-OO.all.rds")

# Data analyses
data <- readRDS("C:/Users/Conno/Desktop/Spring2020/PopGenome/10rep-15samp-0.25maf-A-OO.all.rds")

avgs <- data.table(data %>% group_by(chr, SC) %>% summarise(avg.tP = mean(tP.V1),
                                                 avg.tW = mean(tW.V1),
                                                 avg.TajimaD = mean(Tajima.D.V1),
                                                 avg.eff = mean(eff.V1),
                                                 avg.nuc.fst = mean(`nuc.fst.FST (Nucleotide)`),
                                                 avg.pair.fst = mean(nuc.fst.pair.V1)))


dt <- avgs[!chr %in% c("Scaffold_2158_HRSCAF_2565", "Scaffold_9201_HRSCAF_10758")]

unique(data[eff.V1 > 20000 & !SC=="OO"]$chr)

ggplot(dt, aes(x=SC, y=avg.eff/26)) +
   geom_boxplot() +
   geom_jitter(alpha=0.4, width = 0.2) +
   labs(x="Clonal lineage", y="Avgerage age (years)") +
   scale_x_discrete(labels=c("A" = "A", "B" = "B", "C" = "C", 
                             "E" = "E", "F" = "F", "OO" = "Outbred")) +
   ggtitle("Clonal age - chromosome averages from bootstraps") +
   theme_classic()

avgs <- data %>% filter(!chr %in% c("Scaffold_7757_HRSCAF_8726", 
                                    "Scaffold_9198_HRSCAF_10754")) %>%
                      group_by(SC) %>% summarise(avg.tP = mean(tP.V1),
                                                 avg.tW = mean(tW.V1),
                                                 avg.TajimaD = mean(Tajima.D.V1),
                                                 avg.eff = mean(eff.V1),
                                                 avg.nuc.fst = mean(`nuc.fst.FST (Nucleotide)`))

nucdiv <-data.table(data %>% group_by(SC, chr, win) %>% summarise(pie = mean(pie)))
                                                             
ggplot(nucdiv, aes(x=win, y=pie, col=SC)) +
   geom_smooth(method = "loess", span=0.05, se = FALSE) +
   geom_hline(yintercept = 0.0002) +
   facet_wrap(~chr) +
   labs(x="Position (Mb)", y="Nucleotide diversity (pie)", col="Clonal lineage") + 
   scale_color_discrete(labels=c("A"="A","OO"="Outbred")) +
   ggtitle("Chromosome nucleotide diversities (10Kb windows)") +
   theme_bw() 

nucdiv.avg <- data.table(nucdiv %>% group_by(SC, chr) %>% summarise(pi = mean(pie)))
nucdiv.tot.avg <- data.table(nucdiv.avg %>% group_by(SC) %>% summarise(pi = mean(pi),
                                                                       eff = mean(pi)/(4*5.69e-09),
                                                                       age = (mean(pi)/(4*5.69e-09))/24))

ggplot(nucdiv[SC %in% c("A", "OO")], aes(x=win, y=pie, col=SC)) +
   geom_smooth(method = "loess", span=0.05, se = FALSE) +
   geom_hline(yintercept = 0.0002) +
   facet_wrap(~chr) +
   labs(x="Position (Mb)", y="Nucleotide diversity (pie)", col="Clonal lineage") + 
   scale_color_discrete(labels=c("A"="A","OO"="Outbred")) +
   ggtitle("Chromosome nucleotide diversities (10Kb windows)") +
   theme_bw() 
