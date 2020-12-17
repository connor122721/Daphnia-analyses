# Nucleotide diversity statistics from VCF
# 10.28.2020

# Libraries
suppressMessages(library(PopGenome))
suppressMessages(library(data.table))
suppressMessages(library(foreach))
suppressMessages(library(tidyverse))
suppressMessages(library(SeqArray))

readData()

# Working directory
setwd("/project/berglandlab/connor/PopGenome/")

# Executable in command line
arg <- commandArgs(TRUE)

# Execute in parallel
SC.names <- as.character(arg[1])
chr.n <- as.numeric(arg[2])
upper.maf <- as.numeric(arg[5])
lower.maf <- as.numeric(arg[6])
out.name <- as.character(arg[7])
tmp.dir <- as.character(arg[8])

SC.names="X";chr.n=1;upper.maf=0.25;lower.maf=0;out.name="test";tmp.dir="test"

# Suppress message function
quiet <- function(x) { 
  sink(tempfile()) 
  on.exit(sink()) 
  invisible(force(x)) 
} 

# Get sample and clone ids
samps <- fread("/project/berglandlab/Karen/MappingDec2019/WithPulicaria/June2020/Superclones201617182019withObtusaandPulicaria_kingcorr_20200623_wmedrd.txt")
samps <- samps[Nonindependent==0][is.na(LowReadDepth)][WildSequenced==1]

# Good snps and chromosomes
snps <- fread("/project/berglandlab/Karen/MappingDec2019/WithPulicaria/June2020/finalsetsnpset01pulex_table_20200623")

# Chromosomes of interest
chrom <- unique(snps$chr)
chrom.assign <- fread("vcf/5kbchrassignHiCnew.csv")
chrom.assign$Scaffold <- str_replace_all(chrom.assign$Scaffold, c(";" = "_", "=" = "_"))

# Bed files
bed <- rbind(fread("/project/berglandlab/Karen/MappingDec2019/WithPulicaria/June2020/NsandDepthandChrEnd.sorted.500merged.bed"), 
             fread("/project/berglandlab/Karen/MappingDec2019/WithPulicaria/June2020/RMoutHiCGMgoodscaff.bed"))
names(bed) <- c("chr", "start", "stop")

# Gets denominator for pi
bed <- bed[chr %in% c(chrom)]
bed <- bed[, len:=stop-start]
bed.tot <- bed[, list(total=sum(len)), list(chr)]

# Goes through each clonal lineage
SC <- unlist(samps[SC==SC.names]$clone)
      
message(paste("SC :", SC.names, "Chromosome", chr.n, "replicate #", rep.i, sep=" "))
      
# Load tabindex vcf file - load in by chromosome
vcf <- quiet(readVCF(filename = "vcf/maf_filter_0.25.recode.vcf.gz", 
         tid = c(chrom[chr.n]),
         samplenames = c(SC),
         frompos = 1, 
         topos = chrom.assign[Scaffold==chrom[chr.n]]$Length, 
         numcols = 10000,
         approx = FALSE,
	 out = tmp.dir))
      
# Filter out sites based on MAF
genofile <- quiet(set.filter(vcf, missing.freqs = FALSE, minor.freqs = TRUE, 
                             maf.lower.bound = lower.maf, maf.upper.bound = upper.maf))
      
# Diversity and neutrality statistics
genofile <- quiet(neutrality.stats(genofile, subsites = "included"))
      
# Denominator for pi
denom.pi <- (chrom.assign[Scaffold==chrom[chr.n]]$Length - bed.tot[chr==chrom[chr.n]]$total - sum(table(genofile@region.stats@minor.allele.freqs[[1]][1,]>0.25)[2]))
      
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
                  chr = chrom[chr.n],
                  n.tot.Sites = genofile@n.sites,
                  denom.pi = denom.pi,
                  nMAF.Above = sum(table(genofile@region.stats@minor.allele.freqs[[1]][1,]>0.25)[2]), 
                  nMAF.Below = sum(table(genofile@region.stats@minor.allele.freqs[[1]][1,]<=0.25)[2]),
                  bedSites = bed.tot[chr==chrom[chr.n]]$total,
                  biall.Sites = genofile@n.biallelic.sites, 
                  Sites.inc = table(genofile@region.data@included)[2],
                  n.samp = length(SC),
                  up.maf = upper.maf,
                  low.maf = lower.maf,
                  SC = c(as.character(SC.names)),
                  rep = "total")
      
# Quality control - if replicate does not have any sites > MAF 0.25
if(is.na(total$nMAF.Below) == "TRUE") {
         total$nMAF.Below = genofile@n.biallelic.sites
         total$Sites.inc = table(genofile@region.data@included)[1]

}

write.csv(total, file=paste("/project/berglandlab/connor/PopGenome/scripts/", out.name, ".csv", sep=""), 
          append = T, quote = F, row.names = F)
