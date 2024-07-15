###############################################
### Running LFMM for ACORN. Credit:         ###
### Benjamin Dauphin and Christian Rellstab ###
### Modified for paralel analysis by        ###
### Lars Littmann                           ###
###############################################

####README####

#This script can be used to output raw z-scores from lfmm.
#The script needs to know K beforehand. Another set of scripts is used to determine which value for K is most appropriate.
#The script loops through all environmental factors the user specifies. 
#The script can in principle take on large genomic datasets (Millions),
#but for the sake of speed and RAM usage,it is best to limit it to chunks of 1M SNPs.
#The model for structure is constructed using a thinned dataset, which has to be provided seperately.
#This thinned dataset should be the same across all the chunks of genomic data, if run in parallel.

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
set.k <- as.integer(args[5]) #The K that the model should assume. NEEDS TO BE DETERMINED WITH ANOTHER SCRIPT.
populations <- args [6] #The code for the subset (e.g. 1GA)
chunk <- args [7] #The chunk (parallelisation) that is currently analysed.
selected.env <- gsub(",", " ", args[8]) #A list of the environmental factors that should be analysed. Names seperated by commas WITHOUT SPACES.
                                        #Here the commas are immediately replaced by spaces to make the environmental factors 'legible' as seperate entities.

#### Set the working directory ####
setwd(dir.path)
getwd()

#### Import the environmental dataset ####

env.data <- read.table(paste(dat.env, sep=""), header=T, sep=",")

#extract only the selected environmental variables from the full environmental dataset.
env.variables<-scan(text=selected.env, what= "") #Places the environmental factors in an arrray rather than a continuous string. 
#print(env.variables)

#Make the population names (Plot_ID) the row names of the environmental dataset.
env <- env.data[c(env.variables)]
rownames(env) <- env.data$Plot_ID

#### Import the genetic dataset ####

gen.data <- read.table(paste(dat.gen, sep=""), header=T, sep="\t", row.names = "chrom_pos")

#transpose dataframe
gen <- as.data.frame(t(gen.data))

#Store the names of all the SNPs in a vector.
snp.info <- as.vector(colnames(gen))

#reduce the env set to the gen set. Only the populations that occur in both the environmental and genetic dataset remain. 
rownames(gen) <- gsub("X","",rownames(gen)) #Accounts for a difference in the formating of the population names. 
env.reduced <- env[rownames(env) %in% rownames(gen),] 
env <- env.reduced
identical(as.character(rownames(env)),as.character(rownames(gen)))#check NAs.

#Change from a dataframe to a matrix.
gen.matrix <- as.matrix(gen)
colnames(gen.matrix) <- NULL
rownames(gen.matrix) <- NULL
dim(gen.matrix)

#OPTIONAL: At this stage, output an lfmm file for future use to avoid re-processing the data up to this point. Currently disabled to save disk space.
#write.table(gen.matrix, (paste(dir.path, "/res/", populations,"/", "genetic_data_", populations, "_Chunk_", chunk, ".lfmm", sep = "")), row.names = F, col.names = F, quote=F)

####import the thinned genetic data ####
gen.data.thin <- read.table(paste(dat.gen.thin, sep=""), header = T, sep = "\t", row.names = "chrom_pos")

#transpose dataframe
gen.thin <- as.data.frame(t(gen.data.thin))
snp.info.thin <- as.vector(colnames(gen.thin))

#Store as matrix
gen.thin.matrix <- as.matrix(gen.thin)
colnames(gen.thin.matrix) <- NULL
rownames(gen.thin.matrix) <- NULL

#OPTIONAL: At this stage, output an lfmm file for future use to avoid re-processing the data up to this point. Currently disabled to save disk space.
#write.table(gen.thin.matrix, (paste("./res/", populations, "/", "gen_thinned_matrix", populations, "_Chunk_", chunk, ".lfmm", sep = "")), row.names = F, col.names = F, quote=F)

#Rename datasets for easier shorthand. 
X <- as.matrix(env) #The environmental data
Y <- gen.matrix #The SNPs we are analysing in this implementation of the script
Z <- gen.thin.matrix #The thinned, genome-wide SNPs that we use to account for structure.
Ks <- set.k #The maximum value of K we want the script to loop to (loop goes from 1 to max.k)

print(Ks)

#Prepare the necessary directories. These commands are ignored if the directories already exist.
dir.create(paste(dir.path, "/res/", populations, "/", sep=""), recursive=F)
dir.create(paste(dir.path, "/res/", populations, "/Full_analysis_K_", Ks, "/", sep=""), recursive=F)

for (i in 1:NCOL(X)) {
  dir.create(paste(dir.path, "/res/", populations, "/Full_analysis_K_", Ks, "/environment_", env.variables[i], "/", sep=""), recursive=F)
} 

#### Fit an LFMM based on ridge estimates, i.e, compute B, U, V estimates ####
for (i in 1:NCOL(X)) {
    
  print(paste("Generating Z-scores for environmental factor ", env.variables[i], " and K = ", Ks, sep=""))

  output.dir<-paste(dir.path, "/res/", populations, "/Full_analysis_K_", Ks, "/environment_", env.variables[i], "/", sep="")

  res <- matrix(nrow=NCOL(Y), ncol=2); rownames(res) <- colnames(Y); colnames(res) <- c("SNPid","zscore")
  
  print("First instance of LFMM")
  mod.lfmm2 <- NULL; stats.lfmm2 <- NULL
    
  #Estimate latent factors and environmental effects using the regularised least-squares problem "ridge estimates"

  print("Second instance of LFMM")

  mod.lfmm2 <- lfmm2(input=Z, env=X[,i], K=Ks, lambda=1e-5, effect.sizes=T)
    
  # Statistical tests on genotypic data with imputed missing dat
  print("Third implementation of LFMM")
  stats.lfmm2 <- lfmm2.test(object=mod.lfmm2, input=Y, env=X[,i], full=F, genomic.control=F) 
  res[,"SNPid"] <- snp.info #"SNPid" #res[,"SNPid"] <- snp.info$SNPid
  res[,"zscore"] <- stats.lfmm2$zscores

  # Save a simple table that stores the raw z-score found for every SNP.
  write.table(res, paste(output.dir, "LFMM_Zscores_", populations, "environment_", env.variables[i], "_K", Ks, "_Chunk_", chunk, ".csv", sep=""), sep=",", row.names=F, col.names=T, quote=F)
}
