library(data.table)
library(ggplot2)
library(foreach)
library(cowplot)

growth <- function(N0, r, t) {
  N0*exp(r*t)
}

max.time <- 100
gens <- 2

o.w <- foreach(g=1:gens, .combine="rbind")%do%{
  r <- c(rnorm(mean=1, sd=.05, n=sample.int(5)[1]),
         1.1,
         rnorm(mean=1, sd=.05, n=sample.int(5)[1]),
         1.105,
         rnorm(mean=1, sd=.05, n=sample.int(5)[1]))
  
  
  o <- foreach(i=1:length(r), .combine="rbind")%do%{
    data.table(t=c(1:max.time),
               clone=i,
               N=growth(N0=100, r=r[i], t=c(1:max.time)))
  }
  o.ag <- o[,list(freq=cumsum(N/sum(N)), clone=clone), list(t)]
  #ggplot(data=o.ag, aes(x=t, y=freq, group=clone)) + geom_line()
  

  
  o.w <- foreach(i=1:max(o.ag$t), .combine="rbind")%do%{
    tmp <- o.ag[t==i]
    
    tmp[,start:=c(0, freq[-max(clone)])]
    tmp[,list(y=seq(from=start, to=freq, by=.001)), list(t, clone)]
  }
  
  o.w[,gen:=g]
  o.w
}

ggplot(data=o.w, aes(x=(t), y=y, color=as.factor(clone))) + 
  geom_tile() + 
  facet_grid(~gen) +
  labs(x="Generation", y="Frequency", color="") +
  theme_classic() + 
  theme(axis.text.x = element_text(face="bold", size=16),
        axis.text.y = element_text(face="bold", size=16),
        axis.title.x = element_text(face="bold", size=18),
        axis.title.y = element_text(face="bold", size=18),
        legend.position = "none") 

ggsave(pl, file="~/clone_cartoon.pdf")