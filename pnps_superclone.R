# Connor Murray - 3.27.2020
# In rivanna ssh - type: 
# ijob -c1 -p standard -A berglandlab
# module load intel/18.0 intelmpi/18.0 R/3.6.0; R

# Libraries
library(data.table)
library(SeqArray)
library(foreach)
library(tidyverse)
library(doParallel)

# set working directory
setwd("/project/berglandlab/Karen/MappingDec2019/WithPulicaria/")

# Load meta-data file
samps <- fread("/project/berglandlab/Karen/MappingDec2019/WithPulicaria/Superclones201617182019withObtusaandPulicaria_kingcorr_20200402_wmedrd.txt")

# Load GDS file
genofile <- seqOpen("MapDec19PulexandObtusaandPulicaria_filtsnps10bpindels_snps_filter_pass_lowGQmiss_ann.seq.gds")

# SNP filter file
snpFilter <- fread("finalsetsnpset01pulex_table_wpulicaria_20200401")

# Allele frequency pnps function
pnps.fun <- function(SC.i, MAC, sample.n) {
    # SC.i="A"; MAC=1; sample.n=25; j=1

    seqResetFilter(genofile)
    seqSetFilter(genofile, 
                 sample.id=sample(samps[population%in%c("DCat", "D8", "DBunk")][SC==SC.i][Nonindependent==0][year==2017]$clone, sample.n), 
                 variant.id=snpFilter$variant.ids)

    message("Getting allele counts")
    snp.dt <- data.table(variant.id=seqGetData(genofile, "variant.id"),
                         alleleCount=seqAlleleCount(genofile, ref.allele=1L))

    # Annotations
    message("Annotating variants")

    snp.dt <- snp.dt[alleleCount==MAC]
    seqSetFilter(genofile, variant.id=snp.dt$variant.id)

    # Gets annotations
    tmp <- seqGetData(genofile, "annotation/info/ANN") 
    
    # Gets allele frequencies
    tmp.af <- data.table(variant.id= snp.dt$variant.id,
                         af= seqAlleleFreq(genofile, ref.allele=1L))
    
    len1 <- tmp$length
    len2 <- tmp$data

    snp.dt1 <- data.table(len=rep(len1, times=len1),
                          ann=len2,
                          id=rep(snp.dt$variant.id, times=len1))

    # Extracting data between the 2nd and third | symbol
    snp.dt1[,class:=tstrsplit(snp.dt1$ann,"\\|")[[2]]]

    # Collapsing additional annotations to original SNP vector length
    snp.dt1.an <- snp.dt1[,list(n=length(class), col= paste(class, collapse=",")),
                           list(variant.id=id)]
    
    snp.dt1.an[,col:=tstrsplit(snp.dt1.an$col,"\\,")[[1]]]
    
    # Merges allele frequencies
    snp.dt1.an <- merge(snp.dt1.an, tmp.af, by="variant.id")

    # Final merge
    m <- merge(snp.dt, snp.dt1.an, by="variant.id")

    # Summarize
    m.ag <- m[alleleCount==MAC][,list(pn=sum(col=="missense_variant"), 
                                      ps=sum(col=="synonymous_variant"), 
                                      SC=SC.i, rep=j), 
                                 list(MAC=alleleCount, af=af)]
    
    # Return
    return(m.ag)
}

# Run pnps.out function
pnps.out <- foreach(j=1:2, .combine="rbind") %do% {
  o <- foreach(i=unique(samps[superClone.size > 15 & population %in% c("D8", "DBunk", "Dcat") & !SC =="B"|SC=="OO"]$SC), .combine="rbind", 
               .packages = c("SeqArray")) %do% {
                 message(paste("Sample", i, "#", j, sep=" "))
                 pnps.fun(SC.i=i, MAC=1, sample.n=15)
               }
}

# Saving pnps output
save(pnps.out, file="/project/berglandlab/connor/pnps/pnps.15samps.50rep.Rdata")

### Graphing & Analyses ###

# pnps Rdata
load("C:/Users/Conno/Desktop/Spring2020/pnps.daphnia/updated.pnps.out.15samps.50rep.Rdata")

# metadata file
samps <- fread("C:/Users/Conno/Desktop/Spring2020/pnps.daphnia/cloneinfo.txt")

# adds pnps and filters 
pnps.out[,pnps:=pn/ps]
pnps.out <- na.omit(pnps.out)
pnps.out <- pnps.out[!pnps > 100]
pnps.out.ag <- pnps.out[,list(pnps= mean(pnps), 
                              lci= quantile(pnps, .025), uci= quantile(pnps, .975),
                              var = var(pnps)), 
                         list(SC,af)]

# pnps at several mutation classes boxplot
ggplot(data=pnps.out[af < 0.042 & SC %in% c("A", "C", "OO")], aes(x= as.factor(af), y= pnps, fill=SC)) +
  geom_boxplot() +
    scale_x_discrete(labels = c("1/2N", "2/2N", "3/2N", "4/2N")) +
    scale_fill_discrete(name = "Clonal lineage", labels = c("A", "B", "C", "Outbred")) +
    labs(x="Mutation frequencies", y="pN/pS") +
    theme_bw() 

ggplot(data=pnps.out[af < 0.042 & SC %in% c("A","OO")], 
       aes(x= as.factor(af), y= pnps, color=SC, group=SC)) +
  geom_jitter(alpha=0.3, width = 0.2) + 
  geom_smooth(method = "glm", position = "identity", se = FALSE) +
  scale_x_discrete(labels = c("1/2N", "2/2N", "3/2N", "4/2N")) +
  scale_color_discrete(name = "Clonal lineage", labels = c("A", "Outbred")) +
  labs(x="Mutation frequencies", y="pN/pS") +
  theme_bw() 

ggplot(pnps.out[af < 0.042 & SC %in% c("F")], 
aes(x= pnps, color=SC, group=SC)) +
  geom_density()

# pnps at several mutation classes boxplot 1/2N
ggplot(data=pnps.out[af ==min(af) & SC %in% c("A", "B", "C", "OO")], aes(x= SC, y= pnps, fill=SC)) +
  geom_boxplot()  +
  scale_x_discrete(labels = c("A", "B", "C", "Outbred")) +
  labs(x="Clonal lineage", y="pN/pS") +
  theme_bw() +
  ggtitle("1/2N pN/pS") +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5))

# pnps at several mutation classes points
ggplot(data= pnps.out.ag[af < 0.042 & SC %in% c("A", "B", "C", "OO")], aes(x= interaction(SC,af), y= pnps, color=SC)) +
  geom_point() +
  geom_errorbar(aes(ymin= lci, ymax= uci), width= 0.3) +
    scale_x_discrete(labels = c("1/2N", "2/2N", "3/2N", "4/2N")) +
    labs(color="Clonal lineage", x="Mutation frequencies", y="pN/pS") +
    theme_classic() 

# pnps for 1/2N mutations
ggplot(data=pnps.out.ag[af==min(af)]) +
  geom_point(aes(x=SC, y=pnps, color=SC)) +
  geom_errorbar(aes(x=SC, ymin=lci, ymax=uci, color=SC), width=0.3) +
    labs(x="Clonal lineage", y="pN/pS", title= "1/2N Mutation freuquency") +
    theme_classic() +
    theme(legend.position="none")

# pnps lines across all frequencies
ggplot(data= pnps.out.ag, aes(x=interaction(SC,af), y= pnps, color=SC)) +
  geom_point(size=2.4) +
  geom_errorbar(aes(ymin= lci, ymax= uci), width= 0.3) +
  geom_line(aes(group=SC), size=1.2) +
    labs(color="Clonal lineage", x="Mutation frequencies", y="pN/pS") +
    theme_classic() +
    theme(axis.text.x = element_blank())

ggplot(pnps.out.ag[af==min(af)], aes(x=SC, y=var)) +
  geom_boxplot() 
