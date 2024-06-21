###############################################
### Running LFMM for ACORN. Credit:         ###
### Benjamin Dauphin and Christian Rellstab ###
### Modified for paralel analysis by        ###
### Lars Littmann                           ###
###############################################

####README####

#This version of the script should be run first.
#This script can be used to establish which K to use in the main analysis.
#The script only looks at a single environmental factor. It is set to 'bio1'. This is arbitrary.
#The script loops for several Ks; from 1 to the maximum K set by the user.
#The script can in principle take on large genomic datasets (Millions),
#but for the sake of speed and RAM usage,it is best to limit it to chunks of 1M SNPs.


#### Load packages ####
#install.packages("BiocManager")
#BiocManager::install("LEA")
library(LEA) # https://bcm-uga.github.io/lea/index.html
#install.packages("viridis")
library(viridis)
#install.packages("scales")
library(scales)
#BiocManager::install("qvalue")
library(qvalue)
rm(list = ls())

#### Accept arguments from the command line ####
args=scan(text=commandArgs(trailingOnly = TRUE),what="")
dir.path <- args[1] #First string to receive is the output directory
dat.gen <- args[2] #Path to genomic data under consideration. Needs to be imputed first.
dat.gen.thin <- args[3] #Path to the thinned genomic dataset. Needs to be imputed first.
dat.env <- args[4] #Path to the environmental data 
max.k <- args [5] #The maximum number of K that the script should analyse.
populations <- args [6] #The code for the subset (e.g. 1GA)
chunk <- args [7] #The chunk (parallelisation) that is currently analysed.
#How to set the species:
#sp <- 1 # if running Pet
#sp <- 2 # if running Pub
#sp <- 3 # if running Rob

#### Set variables ####
labsp <- c("Quercus_petraea","Quercus_pubescens", "Quercus_robur")
codesp <- c("petraea","pubescens", "robur")
cols <- c("#E69F00","#0072B2", "green")

#### Set the working directory ####
print(dir.path)
print(dat.gen)
setwd(dir.path)
getwd()

#### Import the environmental dataset ####
env.data <- read.table(paste(dat.env, sep=""), header=T, sep=",")
#head(env.data)

#decide which variables to test
colnames(env.data)
env <- env.data[c('bio1')]
rownames(env) <- env.data$Plot_ID

#### Import the genetic dataset ####

gen.data <- read.table(paste(dat.gen, sep=""), header=T, sep="\t", row.names = "chrom_pos")

#transpose dataframe
gen <- as.data.frame(t(gen.data))

snp.info <- as.vector(colnames(gen))

#reduce the env set to the gen set
rownames(gen) <- gsub("X","",rownames(gen))
env.reduced <- env[rownames(env) %in% rownames(gen),] 
env <- env.reduced
identical(as.character(rownames(env)),as.character(rownames(gen)))#check NAs

gen.matrix <- as.matrix(gen)
colnames(gen.matrix) <- NULL
rownames(gen.matrix) <- NULL
dim(gen.matrix)
write.table(gen.matrix, (paste(dir.path, "/res/", populations,"/", "genetic_data_", populations, "_Chunk_", chunk, ".lfmm", sep = "")), row.names = F, col.names = F, quote=F)

####import the thinned genetic data ####
gen.data.thin <- read.table(paste(dat.gen.thin, sep=""), header = T, sep = "\t", row.names = "chrom_pos")

#transpose dataframe
gen.thin <- as.data.frame(t(gen.data.thin))
snp.info.thin <- as.vector(colnames(gen.thin))

gen.thin.matrix <- as.matrix(gen.thin)
colnames(gen.thin.matrix) <- NULL
rownames(gen.thin.matrix) <- NULL
write.table(gen.thin.matrix, (paste("./res/", populations, "/", "gen_thinned_matrix", populations, "_Chunk_", chunk, ".lfmm", sep = "")), row.names = F, col.names = F, quote=F)

#prepare
X <- as.matrix(env) #for testing
Y <- gen.matrix #The SNPs we are analysing in this implementation of the script
Z <- gen.thin.matrix #The thinned, genome-wide SNPs that we use to account for structure.
Ks <- max.k #Kmax

dir.create(paste(dir.path, "/res/", populations, sep=""), recursive=F)
dir.create(paste(dir.path, "/res/", populations, "/selectingK", sep=""), recursive=F)
for (i in 1:NCOL(X)) {
  setwd(paste(dir.path, "/res/", populations, "/selectingK", sep=""))
  for (j in 1:Ks) {
    dir.create(paste("K", j, "/", sep=""), recursive=F)
  }
} # delete envX folders before running the script
fdr.thres <- c(0.05,0.01,0.001)
fdr.output <- 0.05

gif <- matrix(0, nrow=1, ncol=as.integer(Ks))
rownames(gif) <- colnames(X)
colnames(gif) <- paste("K", 1:Ks, sep="")

nb.asso.q <- matrix(0,nrow=1, ncol=length(fdr.thres))
rownames(nb.asso.q) <- "value"
colnames(nb.asso.q) <- c("q0.05","q0.01","q0.001")

nb.asso.k <- matrix(0, nrow=1, ncol=as.integer(Ks))
rownames(nb.asso.k) <- colnames(X)
colnames(nb.asso.k) <- paste("K", 1:Ks, sep="")

#### Fit an LFMM based on ridge estimates, i.e, compute B, U, V estimates ####

for (i in 1:Ks) {

  print(paste("Generating Z-scores for K = ", i, sep=""))

  setwd(paste(dir.path, "/res/", populations, "/selectingK", "/K", i, "/", sep=""))
  res <- matrix(nrow=NCOL(Y), ncol=2); rownames(res) <- colnames(Y); colnames(res) <- c("SNPid","zscore")
  mod.lfmm2 <- NULL; stats.lfmm2 <- NULL
    
  #Estimate latent factors and environmental effects using the regularised least-squares problem "ridge estimates"

  mod.lfmm2 <- lfmm2(input=Z, env=X[,1], K=i, lambda=1e-5, effect.sizes=T)
    
  # Statistical tests on genotypic data with imputed missing dat
  stats.lfmm2 <- lfmm2.test(object=mod.lfmm2, input=Y, env=X[,1], full=F, genomic.control=F) 
  res[,"SNPid"] <- snp.info #"SNPid" #res[,"SNPid"] <- snp.info$SNPid
  res[,"zscore"] <- stats.lfmm2$zscores

  # res[,"qvalue"] <- p.adjust(as.vector(stats.lfmm2$pvalues), method="fdr", n=length(stats.lfmm2$pvalues))
  write.table(res, paste("LFMM_Zscores_", populations, "_K", i, "_Chunk_", chunk, ".csv", sep=""), sep=",", row.names=F, col.names=T, quote=F) # save all information per SNP
  }

