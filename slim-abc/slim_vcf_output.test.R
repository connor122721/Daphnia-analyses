# Outputs compiled genomic data from SLiM
# 9.6.2020

# Libraries
library(data.table)
library(tidyverse)
library(PopGenome)
library(foreach)

# Executable in command line
arg <- commandArgs(TRUE)
path.name <- "/scratch/csm6hg/slim/test/125340297" # pathway to output
file <- "*125340297_.vcf.gz" # seed for vcf
out.name <- "/scratch/csm6hg/slim/constant.output.csv" # the output file name

# Manual input commands
upper.maf = 0.25
lower.maf = 0

# VCF output folder
setwd(path.name)

# All vcfs
filenames <- list.files(pattern = paste(file, "$", sep=""))

# ReadVCF forloop - calculate diversity statistics.
out <- foreach(i = 1:length(filenames), .combine="rbind", .errorhandling="remove") %do% {

  # Load tabindex vcf file
  vcf <- readVCF(filename = filenames[i],
        tid = c("1"),
        frompos = 1,
        topos = 150000000,
        numcols = 10000,
        approx = FALSE)

  # Filter out sites based on MAF
  genofile <- set.filter(vcf, minor.freqs = TRUE, 
                       maf.lower.bound = lower.maf, maf.upper.bound = upper.maf, 
                       missing.freqs = FALSE)

  # Diversity and neutrality statistics
  genofile <- neutrality.stats(genofile, subsites = "included")

  # Denominator for pi
  denom.pi <- (genofile@n.sites - sum(table(genofile@region.stats@minor.allele.freqs[[1]][1,] > upper.maf)[2]))
  
  # Some replicates will not have any sites > MAF 0.25
  # Included sites will be based off number of biallelic sites
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
            chr = c("1"),
            n.tot.Sites = genofile@n.sites,
            denom.pi = denom.pi,
            nMAF.Above = sum(table(genofile@region.stats@minor.allele.freqs[[1]][1,] > upper.maf)[2]), 
            nMAF.Below = sum(table(genofile@region.stats@minor.allele.freqs[[1]][1,] <= upper.maf)[2]),
            biall.Sites = genofile@n.biallelic.sites, 
            Sites.inc = table(genofile@region.data@included)[2],
	    up.maf = upper.maf,
	    low.maf = lower.maf,
            vcf = filenames[i],
            model = tstrsplit(filenames[i], "_")[[2]],
            EG = tstrsplit(filenames[i], "_")[[3]],
            K = tstrsplit(filenames[i], "_")[[4]],
            seed = tstrsplit(filenames[i], "_")[[5]],
	    rep = tstrsplit(filenames[i], "_")[[6]],
	    nSamp = tstrsplit(filenames[i], "_")[[7]],
	    sim.gen = tstrsplit(filenames[i], "_")[[8]],
	    gen = tstrsplit(filenames[i], "_")[[9]],
	    iteration = i)
  
  # Quality control - if replicate does not have any sites > MAF 0.25
  if(is.na(total$nMAF.Below) == "TRUE") {
    total$nMAF.Below = genofile@n.biallelic.sites
    total$Sites.inc = table(genofile@region.data@included)[1]
  }
  
  # Progress message
  print(paste("Iteration-", i, ": ", round((i/length(filenames))*100, digits=2), "%", " complete", sep=""))

  # Finish
  return(total)
  Sys.sleep(0.5)

}

# Write final output
write.csv(out, file=out.name, append = TRUE, row.names = FALSE, col.names = NA)
