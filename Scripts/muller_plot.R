library(ggmuller)
library(tidyverse)
library(cowplot)
library(data.table)

data <- data.frame(Generation=1:100, 
                   Identity=c("A", "B", "C", "D"),
                   Population=1:100)

# a function for generating exponential growth curves:
pop_seq <- function(gens, lambda, start_gen) c(rep(0, start_gen),
                                               exp(lambda * gens[0:(length(gens) - start_gen)]))
lambda <- 0.1 # baseline fitness
gens <- 0:150 # time points
fitnesses <- c(1, 2, 2.2, 2.5, 3, 3.2, 3.5, 3.5, 3.8) # relative fitnesses of genotypes
pop3 <- data.frame(Generation = rep(gens, 9),
                   Identity = paste0("clone_", LETTERS[rep(1:9, each = length(gens))]),
                   Population = c(1E2 * pop_seq(gens, fitnesses[1]*lambda, 0), 
                                  pop_seq(gens, fitnesses[2]*lambda, 0), 
                                  pop_seq(gens, fitnesses[3]*lambda, 10), 
                                  pop_seq(gens, fitnesses[4]*lambda, 20),
                                  pop_seq(gens, fitnesses[5]*lambda, 30),
                                  pop_seq(gens, fitnesses[6]*lambda, 40),
                                  pop_seq(gens, fitnesses[7]*lambda, 50),
                                  pop_seq(gens, fitnesses[8]*lambda, 50),
                                  pop_seq(gens, fitnesses[9]*lambda, 60)),
                   Fitness = rep(fitnesses, each = length(gens)))

Muller_df3 <- get_Muller_df(edges3, pop3)

ggplot(Muller_df3, aes_string(x = "Generation", y = "Frequency", group = "Group_id", fill = "Identity", colour = "Identity")) + 
  geom_area() +
  theme(legend.position = "right") +
  guides(linetype = FALSE, color = FALSE) + 
  scale_y_continuous(labels = 25 * (0:4), name = "Percentage") +
  theme_classic()

Muller_df1 <- get_Muller_df(example_edges, example_pop_df, cutoff = 0.2)
Muller_df1$Val <- as.numeric(Muller_df1$Identity)
p <- Muller_plot(Muller_df1, colour_by = "Val", conceal_edges = T) +
  theme_classic() +
  labs(title="Year 1") +
  theme(title = element_text(face="bold", size=14),
        axis.text.x = element_text(face="bold", size=10),
        axis.text.y = element_text(face="bold", size=10),
        axis.title.x = element_text(face="bold", size=14),
        axis.title.y = element_text(face="bold", size=14), 
        plot.title = element_text(hjust = 0.5, vjust = -4),
        legend.position = "none") 

Muller_df2 <- get_Muller_df(example_edges, example_pop_df, cutoff = 0.1)
q <- Muller_plot(Muller_df2, ylab = "", conceal_edges = T) +
  theme_classic() +
  labs(title="Year 2") +
  theme(title = element_text(face="bold", size=14),
        axis.text.x = element_text(face="bold", size=10),
        axis.text.y = element_blank(),
        axis.title.x = element_text(face="bold", size=14),
        axis.title.y = element_blank(), 
        axis.line.y.left = element_blank(),
        axis.ticks.y = element_blank(),
        plot.title = element_text(hjust = 0.5, vjust = -4),
        legend.position = "none")

plot_grid(p,q)

old <- data.table(piNS=c(1,1.4,1.2,1.6,1.3,1.8, 0.4, 0.6, 0.8, 0.9, 0.7, 0.7),
                  clone=rep(c("New", "Old"), each=6))

a <- ggplot(old, aes(x=clone, y=piNS)) +
  geom_boxplot() +
  # geom_hline(yintercept = 1, linetype=2, size=1.5) +
  theme_classic() +
  labs(x="Clonal lineage", y="piN/piS", title="Purifying selection stronger in old clones") +
  theme(title = element_text(face="bold", size=12),
        axis.text.x = element_text(face="bold", size=10),
        axis.text.y = element_text(face="bold", size=10),
        axis.title.x = element_blank(),
        axis.title.y = element_text(face="bold", size=14), 
        plot.title = element_text(hjust = 0.5, vjust = -4),
        legend.position = "none")
  
old2 <- data.table(piNS=c(1,1.4,1.2,1.6,1.3,1.8, 1.3,1.4,1.2,1.6,1.3,1.5),
                  clone=rep(c("New", "Old"), each=6))

b <- ggplot(old2, aes(x=clone, y=piNS)) +
  geom_boxplot() +
  # geom_hline(yintercept = 1, linetype=2, size=1.5) +
  theme_classic() +
  ylim(c(0.4, 1.8)) +
  labs(x="Clonal lineage", y="piN/piS", title="Purifying selection strength equal across clones") +
  theme(title = element_text(face="bold", size=12),
        axis.text.x = element_text(face="bold", size=10),
        axis.text.y = element_text(face="bold", size=10),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(), 
        plot.title = element_text(hjust = 0.5, vjust = -4),
        legend.position = "none")

old3 <- data.table(piNS=c(0.4, 0.6, 0.8, 0.9, 0.7, 0.7, 1,1.4,1.2,1.6,1.3,1.8),
                   clone=rep(c("New", "Old"), each=6))

c <- ggplot(old3, aes(x=clone, y=piNS)) +
  geom_boxplot() +
  # geom_hline(yintercept = 1, linetype=2, size=1.5) +
  theme_classic() +
  ylim(c(0.4, 1.8)) +
  labs(x="Clonal lineage", y="piN/piS", title="Purifying selection weaker in old clones") +
  theme(title = element_text(face="bold", size=12),
        axis.text.x = element_text(face="bold", size=10),
        axis.text.y = element_text(face="bold", size=10),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(), 
        plot.title = element_text(hjust = 0.5, vjust = -4),
        legend.position = "none")

plot_grid(a,b,c, nrow = 1)

data <- data.table(frac=c(0.5, 0.3, 0.4, 0.3),
                   sel=c("Low", "Low",  "High", "High"),
                   clone=c("New", "Old", "New", "Old"))

data$sel = factor(data$sel, levels = c("Low", "High"))

ggplot(data, aes(x=as.factor(sel), y=frac, group=clone, fill=clone)) +
  geom_col(position = "dodge") +
  theme_classic() +
  labs(x="Purifying selection strength", y="Fraction of sites",
       fill="Clonal Lineage") +
  theme(title = element_text(face="bold", size=12),
        axis.text.x = element_text(face="bold", size=10),
        axis.text.y = element_blank(),
        axis.title.x = element_text(face="bold", size=14),
        axis.title.y = element_text(face="bold", size=14), 
        axis.ticks.y = element_blank(),
        plot.title = element_text(hjust = 0.5, vjust = -4),
        legend.title = element_text(face="bold", size=14),
        legend.position = c(0.8,0.9))

