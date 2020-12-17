# Nucleotide diversity from VCF
# 8.18.2020

module load intel/18.0 intelmpi/18.0 R/3.6.3; R

# Libraries
library(PopGenome)
library(data.table)
library(foreach)
library(tidyverse)
library(SeqArray)

setwd("/project/berglandlab/connor/PopGenome/")

# Get sample ids and clonal associations
samps <- fread("/project/berglandlab/Karen/MappingDec2019/WithPulicaria/June2020/Superclones201617182019withObtusaandPulicaria_kingcorr_20200623_wmedrd.txt")
samps <- samps[Nonindependent==0][is.na(LowReadDepth)]

# Good snps and chromosomes
snps <- fread("/project/berglandlab/Karen/MappingDec2019/WithPulicaria/June2020/snpsvarpulexpresentinhalf_table_20200623")

# Chromosomes of interest
chrom <- unique(snps$chr)
chrom.assign <- fread("vcf/5kbchrassignHiCnew.csv")
chrom.assign$Scaffold <- str_replace_all(chrom.assign$Scaffold, c(";" = "_", "=" = "_"))

# Bed files
bed <- rbind(fread("/project/berglandlab/Karen/MappingDec2019/WithPulicaria/NsandDepthandChrEnd.sorted.500merged.bed"), 
             fread("/project/berglandlab/Karen/MappingDec2019/WithPulicaria/RMoutHiCGMgoodscaff.bed"))
names(bed) <- c("chr", "start", "stop")

# Gets denominator for pi
bed <- bed[chr %in% c(chrom)]
bed <- bed[, len:=stop-start]
bed.tot <- bed[, list(total=sum(len)), list(chr)]

# Names of super clones
SC.names <- c("A", "B", "C", "D", "OO")

# Reads in vcf with bootstrap sampling
vcf.fun <- function(sample.n, lower.maf, upper.maf) {
   # sample.n=6; lower.maf=0; upper.maf=0.25; cur.pop=1; chr.num=10; i=1
   
   # Goes through each SC
   foreach(cur.pop = 1:length(SC.names), .combine="rbind") %do% {
      
      # Sample clones
      SC <- sample(unlist(samps[SC==SC.names[cur.pop]]$clone), 
                   size = sample.n, replace = FALSE)
      
   # Goes through each chromosome
   foreach(chr.num = 1:12, .combine="rbind", .errorhandling = "remove") %do% {
      
      message(paste("SC :", SC.names[cur.pop], "Chromosome", chr.num, "replicate #", i, sep=" "))
      
      # Load tabindex vcf file - load in by chromosome
      vcf <- readVCF(filename = "vcf/maf_filter_goodpos.vcf.gz", 
         tid = c(chrom[chr.num]),
         samplenames = c(SC),
         frompos = 1, 
         topos = chrom.assign[Scaffold==chrom[chr.num]]$Length, 
         numcols = 10000,
         approx = FALSE)
      
      # Filter out sites based on MAF
      genofile <- set.filter(vcf, missing.freqs = FALSE, minor.freqs = TRUE, 
                             maf.lower.bound = lower.maf, maf.upper.bound = upper.maf)
      
      # Diversity and neutrality statistics
      genofile <- neutrality.stats(genofile, subsites = "included")
      
      # Denominator for pi
      denom.pi <- (chrom.assign[Scaffold==chrom[chr.num]]$Length - bed.tot[chr==chrom[chr.num]]$total - sum(table(genofile@region.stats@minor.allele.freqs[[1]][1,]>0.25)[2]))
      
      # Some replicates will not have any sites > MAF 0.25
      if(is.na(denom.pi) =="TRUE") {
         denom.pi = genofile@n.sites
      }
      
      # Compile data
      total <- data.table(tW = genofile@theta_Watterson, 
                  tW.stan = genofile@theta_Watterson/denom.pi,
                  tP = genofile@theta_Tajima,
                  tP.stan = genofile@theta_Tajima/denom.pi,
                  tP.eff = (genofile@theta_Tajima/denom.pi)/(4*5.69e-09),
                  tajima.D = genofile@Tajima.D,
                  chr = chrom[chr.num],
                  n.tot.Sites = genofile@n.sites,
                  denom.pi = denom.pi,
                  nMAF.Above = sum(table(genofile@region.stats@minor.allele.freqs[[1]][1,]>0.25)[2]), 
                  nMAF.Below = sum(table(genofile@region.stats@minor.allele.freqs[[1]][1,]<=0.25)[2]),
                  bedSites = bed.tot[chr==chrom[chr.num]]$total,
                  biall.Sites = genofile@n.biallelic.sites, 
                  Sites.inc = table(genofile@region.data@included)[2],
                  n.samp = sample.n,
                  up.maf = upper.maf,
                  low.maf = lower.maf,
                  SC = c(as.character(SC.names[cur.pop])),
                  rep = i)
      
      # Quality control - if replicate does not have any sites > MAF 0.25
      if(is.na(total$nMAF.Below) == "TRUE") {
         total$nMAF.Below = genofile@n.biallelic.sites
         total$Sites.inc = table(genofile@region.data@included)[1]
      }
      
      # Progress message
      print(paste("Iteration - ", i, " : ", round((i/length(filenames))*100, digits = 2), "%", sep=""))
      
      # Finish by chrom
      return(total)
      
      }
   }
}

# Forloop of vcf function
vcf.out <- foreach(i = 1:10, .combine="rbind") %do% {
                   message(paste("Sample #", i, sep=" "))
                   vcf.fun(sample.n = 15, lower.maf=0, upper.maf = 0.25)
}

saveRDS(vcf.out, file="/project/berglandlab/connor/PopGenome/filter-10rep-15samp.rds")
