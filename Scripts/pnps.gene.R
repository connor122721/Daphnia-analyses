# Libraries
library(data.table)
library(SeqArray)
library(foreach)
library(tidyverse)
library(SNPRelate)
library(seqinr)

# set working directory
setwd("/project/berglandlab/Karen/MappingDec2019/WithPulicaria/")

# Load meta-data file
samps <- fread("June2020/Superclones201617182019withObtusaandPulicaria_kingcorr_20200623_wmedrd.txt")
samps <- samps[Nonindependent==0][is.na(LowReadDepth)]

# Only OO clones from 2017 DBunk
samps <- samps[!clone %in% c(samps[SC%in%c("OO")][year==2017][population=="DBunk"]$clone)]

# Gene analyses
pro <- read.fasta("/project/berglandlab/Karen/genomefiles/Daphnia.proteins.aed.0.6.fasta", seqtype="AA")
pro <- data.table(gene.name=getName(pro), AA.length=getLength(pro))
pro <- pro[,splice:=tstrsplit(gene.name, "-")[[2]]]

# Load GDS file
genofile <- seqOpen("June2020/MapJune2020_ann.seq.gds")

# Good snps and chromosomes
snpFilter <- fread("June2020/finalsetsnpset01pulex_table_20200623")

# Gene pnps function
pnps.fun <- function(SC.i, sample.n) {
  
  SC.i="A"; sample.n=15;j=1
    
    # Filter
    seqResetFilter(genofile)
    seqSetFilter(genofile, sample.id=sample(samps[population%in%c("D8", "DBunk", "DCat", "Dcat", "D10")][SC==SC.i]$clone, sample.n),
               variant.id=snpFilter$variant.ids)

    # Extract data
    snp.dt <- data.table(variant.id=seqGetData(genofile, "variant.id"),
                         position=seqGetData(genofile, "position"), 
                         chr=seqGetData(genofile, "chromosome"),
                         alleleCount=seqAlleleCount(genofile, ref.allele=1L),
                         af=seqAlleleFreq(genofile, ref.allele=1L))
    
    # Gets annotations
    tmp <- seqGetData(genofile, "annotation/info/ANN") 
    len1 <- tmp$length
    len2 <- tmp$data
    snp.dt1 <- data.table(len=rep(len1, times=len1),
                          ann=len2,
                          variant.id=rep(snp.dt$variant.id, times=len1))

    # Extract SNP class and genes
    snp.dt1[,class:=tstrsplit(snp.dt1$ann,"\\|")[[2]]]
    snp.dt1[,gene:=tstrsplit(snp.dt1$ann,"\\|")[[4]]]
    snp.dt1[,feature.id:=tstrsplit(snp.dt1$ann,"\\|")[[7]]]
    
    # Keep first row - most influential SNP class presented by SeqArray
    snp.dt1 <- data.table(snp.dt1 %>% group_by(variant.id) %>% 
                          filter(row_number() == 1) %>% select(-c(ann)))
  
    # Merge with gene data
    snp.dt1.an <- merge(snp.dt1, snp.dt, by="variant.id")

    # Summarize
    m <- snp.dt1.an[class %in% c("missense_variant", "synonymous_variant")][,
                                     list(pn = sum(class=="missense_variant"), 
                                          ps = sum(class=="synonymous_variant"), 
                                          SC = SC.i, rep = j, sample.n,
                                          variant.id, position, class,
                                          chr, af, MAC = alleleCount), 
                                     list(gene, feature.id)]
    
    # Final merge
    m <- merge(m, pro, by.x="feature.id", by.y="gene.name")
    
    # Return
    return(m)
}

# Clones of interest
clones <- c("A", "B", "C", "D", "E", "F", "OO")

# Run pnps.out function
pnps.out <- foreach(j = 1:10, .combine="rbind") %do% {
  o <- foreach(i = clones, .combine="rbind") %do% {
                 message(paste("Sample", i, "#", j, sep=" "))
                 pnps.fun(SC.i=i, sample.n=15)
  }
}

# Saving pnps output
save(pnps.out, file="/project/berglandlab/connor/pnps/gene.pnps.15samp.10rep.filt.Rdata")
