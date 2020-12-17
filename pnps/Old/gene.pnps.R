# Gene pnps ratios across clonal lineages

# Libraries
library(data.table)
library(tidyverse)
library(moments)
library(cowplot)
library(ggpubr)
library(ggsignif)

# pnps Rdata
load("C:/Users/Conno/Desktop/Fall2020/Berglandlab/gene.pnps.15samp.10rep.filt.Rdata")

# Change name of SC OO
pnps.out[SC=="OO"]$SC <- "Outbred"

# Add pnps and filters
pnps.out[af > 0.5, maf:= 1-af]
pnps.out[is.na(maf), maf:=af]
pnps.out[, nMissing:= 30 - MAC/af]
pnps.out <- pnps.out[nMissing==0 & af < 1]

# Recalculate pnps for folded SFS
pnps.out <- pnps.out[, list(pn=mean(pn), ps=mean(ps)), 
                     list(gene, SC, AA.length)]
pnps.out <- pnps.out[, pnps:=(pn/ps)]

# Removes strange replicates and infinte values
pnps.out <- pnps.out[!pnps > 100] 

# Genes common across all lineages
genes <- data.table(pnps.out %>% count(gene))
pnps.out <- pnps.out[gene %in% genes[n==7]$gene]

# Gene pnps
ggplot(data=pnps.out, aes(x=gene, y=pnps, color=SC)) +
  geom_point() +
  geom_hline(yintercept=1, linetype=2, color="red", size=1.2) +
  facet_wrap(~SC, nrow=1) +
  coord_flip() +
  theme_bw() +
  theme(legend.position = "none",
        axis.text.x = element_text(face="bold", size=8), 
        axis.text.y = element_blank(), 
        axis.ticks.y = element_blank(),
        axis.title.x = element_text(face="bold", size=12),
        axis.title.y = element_text(face="bold", size=12),
        strip.text = element_text(face="bold", size=10)) +
  labs(x="Gene", y="pN/pS")
