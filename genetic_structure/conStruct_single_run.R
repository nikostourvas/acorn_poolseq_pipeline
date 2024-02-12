# conStruct spatial and non-spatial analysis

library(conStruct)
library(dplyr)

# Set working directory
setwd("")

# Genetic Data input
gen <- read.table(file="",
                   header=TRUE,
                   sep = "\t",
                   check.names=FALSE)

# Remove populations 460 & 500 (outgroups) and Q.robur
gen <- gen[, -c(76:117)]
rownames(gen) <- gen[,1]
gen <- gen[,-1]

# Transpose the data
gen_lfmm <- t(gen)

# Tree data input
treedata <- read.csv("20240102_Tree_data.csv", header=TRUE)
# Use only pools included for Poolseq and in this case only Q. petraea
treedata <- treedata %>% 
  filter(Pool_seq == TRUE) %>% 
  filter(Species_pop_level != "Q.robur")

# Transform the data so that there in only one line per pool
popdata <- treedata[,-c(4,5,6,11:24,26,27,31,32,33,34)]
popdata <- unique(popdata)

popdata <- popdata[,c(1,6:7)]
popdata$Plot_ID <- as.character(popdata$Plot_ID)
CoordRaw <- popdata

Coord=as.matrix(CoordRaw[,3:2])
colnames(Coord)=c("Lon","Lat")

# Load coordinates of each sampling site and geographic distances in kilometers
CoordRaw=read.csv("../../../data/AllPools_woOutgroups_geographic_coords.csv",header=T)
Coord=as.matrix(CoordRaw[,2:3])
colnames(Coord)=c("Lon","Lat")

# run a conStruct analysis

#   you have to specify:
#       the number of layers (K)
#       the allele frequency data (freqs)
#       the sampling coordinates (coords)
#
#   if you're running the nonspatial model, 
#       you do not have to specify 
#       the geographic distance matrix (geoDist)

my.run <- conStruct(spatial = FALSE, 
                    K = 3, 
                    freqs = Freq, 
                    geoDist = NULL, 
                    coords = Coord,
                    prefix = "nspK3",
                    n.chains = 1,
                    n.iter = 1e6,
                    make.figs = TRUE,
                    save.files = TRUE)


