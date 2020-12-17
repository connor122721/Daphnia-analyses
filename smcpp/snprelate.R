# Using SNPRelate for genomic similarity across samples: PCA for fun
# Connor Murray
# 5.1.2020

module load intel/18.0 intelmpi/18.0 R/3.6.0; R

# Libraries
library(data.table)
library(foreach)
library(ggplot2)
library(tidyverse)
library(doParallel)
library(gdsfmt)
library(SNPRelate)
library(PopGenome)

# Daphnia data 
setwd("C:/Users/Conno/Desktop/Spring2020/PopGenome/")

# Use to verify sample ids
vcf_handle <-.Call("VCF_open", "vcf/daphnia.ann.vcf.gz")
.Call("VCF_getSampleNames",vcf_handle)

# Get sample ids and clonal associations
samps <- fread("samples.txt")

# Good snps and chromosomes
snps <- fread("snps.txt")
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

snpgdsSummary("vcf/daphnia.ann.seq.snp.gds")
genofile <- snpgdsOpen("vcf/daphnia.ann.seq.snp.gds")

# LD pruning
# set.seed(1000)
# snpset = snpgdsLDpruning(genofile, ld.threshold = 0.5, autosome.only = FALSE, snp.id = snps[goodChr==TRUE]$variant.ids, sample.id = c(A, B, C, E, f, O))
# snp.id = unlist(snpset)

?snpgdsDiss

#dendogram
dissMatrix <- snpgdsDiss(genofile , sample.id = c(A, B, C, E, f, O), snp.id = NULL, 
                         autosome.only = FALSE, remove.monosnp = TRUE, 
                         maf = NaN, missing.rate = NaN, verbose = TRUE, num.thread = )

snpHCluster <-  snpgdsHCluster(dissMatrix, sample.id = c(A, B, C, E, f, O), 
                               need.mat = TRUE, hang = 0.25)

cutTree <- snpgdsCutTree(snpHCluster, z.threshold=15, outlier.n=5,  
                         n.perm = 5000, samp.group=NULL, col.outlier="red", 
                         col.list=NULL, pch.outlier=4, pch.list=NULL,
                         label.H=FALSE, label.Z=TRUE, verbose=TRUE)

snpgdsDrawTree(cutTree, main = "A - F & O", leaflab="perpendicular")


#pca
pca <- snpgdsPCA(genofile)


tab <- data.frame(sample.id = pca$sample.id,pop = factor(pop_code)[match(pca$sample.id, sample.id), EV1 = pca$eigenvect[,1], EV2 = pca$eigenvect[,2], stringsAsFactors = FALSE)

plot(tab$EV2, tab$EV1, col=as.integer(tab$pop),xlab="eigenvector 2", ylab="eigenvector 1")
legend("topleft", legend=levels(tab$pop), pch="o", col=1:nlevels(tab$pop))