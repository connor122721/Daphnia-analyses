# Libraries
library(data.table)
library(tidyverse)
library(ggmap)
library(maps)
library(mapdata)
library(ggrepel)

setwd("C:/Users/Conno/Desktop/Fall2020/Berglandlab/daphnia_phylo")

# Metadata 
dt <- data.table(read.csv("Daphnia_pulex_1703.csv"))

test <- data.table(read.csv("samps.csv", header = T))
samps <- data.table(read.csv("samples.csv", header = T))
samps <- data.table(test %>% left_join(samps, by = "clone"))
unique(samps$population)
write.csv(samps[population %in% c('D8', 'DBunk', 'Dcat', 'DCat', 'DOil', 'Dramp')]$clone, 'wildtrust.samples.csv')

# Extracting latitude and longitude
dt$lat_lon <- gsub("[\\(\\)]", "", regmatches(dt$Isolation_source, 
                   gregexpr("\\(.*?\\)", dt$Isolation_source)))

# Population genomics of Da. pulex
dt[BioProject %in% c("PRJNA351263")]$lat_lon <- c("40.1224,−87.7366")

# PA42 A New Reference Genome Assembly for the Microcrustacean Daphnia pulex
dt[BioProject %in% c("PRJNA307976")]$lat_lon <- c("40.2013,-87.3294")

# Genetic control of male production in Daphnia pulex
dt[BioProject %in% c("PRJNA513203")][Isolate %like% c("KAP")]$lat_lon <- c("40.1224,-87.7366")
dt[BioProject %in% c("PRJNA513203")][Isolate %like% c("POV")]$lat_lon <- c("42.75,-85.35")
dt[BioProject %in% c("PRJNA513203")][Isolate %like% c("TEX")]$lat_lon <- c("42.20,-83.60")
dt[BioProject %in% c("PRJNA513203")][Isolate %like% c("LPB")]$lat_lon <- c("42.6798,-80.4528")
dt[BioProject %in% c("PRJNA513203")][Isolate %like% c("NFL")]$lat_lon <- c("39.9002,-84.9223")

# Spontaneous MA in Da. pulex in selection and selection free environments
dt[BioProject %in% c("PRJNA341529")]$lat_lon <- c("42.12,-82.98")
dt <- dt[!Sample.Name == "Daphnia pulex non-MA isolates"]

# Rapid New Zealand samples
dt[BioProject %in% c("PRJNA573527")][Library.Name %in% "Alexandrina"]$lat_lon <- c("-43.57,170.27")
dt[BioProject %in% c("PRJNA573527")][Library.Name %in% "Tekapo"]$lat_lon <- c("-43.17,171.30")
dt[BioProject %in% c("PRJNA573527")][Library.Name %in% "Hawea"]$lat_lon <- c("-44.30,169.17")
dt[BioProject %in% c("PRJNA573527")][Library.Name %in% "Hayes"]$lat_lon <- c("-44.59,168.48")
dt[BioProject %in% c("PRJNA573527")][Library.Name %in% "Johnson"]$lat_lon <- c("-45,168.43")
dt[BioProject %in% c("PRJNA573527")][Library.Name %in% "Ohau"]$lat_lon <- c("-44.15,169.51")
dt[BioProject %in% c("PRJNA573527")][Library.Name %in% "Pukaki"]$lat_lon <- c("-44.07,170.10")
dt[BioProject %in% c("PRJNA573527")][Library.Name %in% "Sullivans Dam"]$lat_lon <- c("-45.48,170.31")
dt[BioProject %in% c("PRJNA573527")][Library.Name %in% "Moke"]$lat_lon <- c("-44.998117,168.573089")
dt[BioProject %in% c("PRJNA573527")][Library.Name %in% "Tekapo"]$lat_lon <- c("-44.00,170.29")
dt[BioProject %in% c("PRJNA573527")][Library.Name %in% c("Wakatipu", "Von")]$lat_lon <- c("-45.3,168.30")
dt[BioProject %in% c("PRJNA573527")][Library.Name %in% "Wanaka"]$lat_lon <- c("-44.30,169.08")
dt[BioProject %in% c("PRJNA573527")][Library.Name %in% "Kapoai"]$lat_lon <- c("-36.0468,173.8326")

# New columns 
dt$lat <- as.numeric(tstrsplit(dt$lat_lon, ",")[[1]])
dt$lon <- as.numeric(tstrsplit(dt$lat_lon, ",")[[2]])

# All missing data
missing <- dt[is.na(lat)]
dt <- dt[!is.na(lat)][!Isolate=="Mutation accumulation lines"]

# Add counts of samples
dt <- data.table(merge(dt, dt %>% group_by(BioProject) %>% 
                          dplyr::summarize(counts = n()), by="BioProject"))

# All found data
pule <- data.table(read.csv("Dpulex.found.csv"))
pule <- pule[!is.na(lat)][!Isolate=="Mutation accumulation lines"]

# Pulicaria
puli <- data.table(read.csv("pulicaria.csv"))
puli <- data.table(merge(puli, puli %>% group_by(BioProject) %>% 
                         dplyr::summarize(counts = n()), by="BioProject"))

# Adding location
puli[geo_loc_name %like% c("South Center Lake"),lat:=45.3784978][geo_loc_name %like% c("South Center Lake"),lon:=-92.821978]

puli[Library.Name %in% "Alexandrina", lat_lon:=c("-43.57,170.27")]
puli[Library.Name %in% "Tekapo"]$lat_lon <- c("-43.17,171.30")
puli[Library.Name %in% "Hawea"]$lat_lon <- c("-44.30,169.17")
puli[Library.Name %in% "Hayes"]$lat_lon <- c("-44.59,168.48")
puli[Library.Name %in% "Johnson"]$lat_lon <- c("-45,168.43")
puli[Library.Name %in% "Ohau"]$lat_lon <- c("-44.15,169.51")
puli[Library.Name %in% "Pukaki"]$lat_lon <- c("-44.07,170.10")
puli[Library.Name %in% "Sullivans Dam"]$lat_lon <- c("-45.48,170.31")
puli[Library.Name %in% "Moke"]$lat_lon <- c("-44.998117,168.573089")
puli[Library.Name %in% "Tekapo"]$lat_lon <- c("-44.00,170.29")
puli[Library.Name %in% c("Wakatipu", "Von")]$lat_lon <- c("-45.3,168.30")
puli[Library.Name %in% "Wanaka"]$lat_lon <- c("-44.30,169.08")
puli[Library.Name %in% "Kapoai"]$lat_lon <- c("-36.0468,173.8326")

# Obtusa
obtu <- data.table(read.csv("obtusa.csv"))
obtu <- data.table(merge(obtu, obtu %>% group_by(BioProject) %>% 
                         dplyr::summarize(counts = n()), by="BioProject"))

obtu[geo_loc_name %like% c("Ziggy's Head 21"),lat:=33.820472][geo_loc_name %like% c("Ziggy's Head 21"),lon:=80.822361]

# Merge metadata
spp <- data.table(pule %>% select(BioProject, BioSample, Run, Organism, lat, lon, counts) %>% 
                    full_join(obtu %>% select(BioProject, BioSample, Run, Organism, lat, lon, counts))) %>% 
                    full_join(puli %>% select(BioProject, BioSample, Run, Organism, lat, lon, counts))

# World map
canada <- map_data("worldHires", "Canada")
usa <- map_data("worldHires", "usa")
world <- map_data("world")

# Random sample
rand <- data.table(spp[!is.na(lat)] %>% group_by(BioProject, lat, lon) %>% sample_n(1))
write.csv(spp, "Daphnia.SRR.csv")

rand.counts <- data.table(rand[!is.na(lat)] %>% group_by(BioProject) %>% sample_n(1))

# World map
ggplot() +
  geom_polygon(data=world, aes(x = long, y = lat, group=group), fill="white", colour = "black") +
  # geom_polygon(data=usa, aes(x = long, y = lat, group=group), fill="white", colour = "black") +
  # geom_polygon(data=canada, aes(x = long, y = lat, group=group), fill="white", colour = "black") +
  geom_point(data=rand, aes(x=lon, y=lat, color=Organism)) +
  # coord_fixed(xlim = c(-100, -65),  ylim = c(25, 50), ratio = 1.2) +
  labs(x="Longitude", y="Latitude") +
  theme_classic() +
  theme(axis.text.x = element_text(face = "bold", size = 12),
        axis.text.y = element_text(face = "bold", size = 12),
        axis.title.x = element_text(face = "bold", size = 14),
        axis.title.y = element_text(face = "bold", size = 14)) 
