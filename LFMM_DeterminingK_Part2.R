###############################################
### Running LFMM for ACORN. Credit:         ###
### Benjamin Dauphin and Christian Rellstab ###
### Modified for paralel analysis by        ###
### Lars Littmann                           ###
###############################################

####README####

#This is the second script that is run to determine an appropriate value for K in LFMM.
#This script should only be used to establish which K to use in the main analysis.
#The script can take on all environmental factors, but it is best to use a representative subset to keep runtimes short.
#The script loops for several Ks; from 1 to the maximum K set by the user.
#The script merges data from parallel runs across chunks of the genome (created in part 1) and outputs plots for determining K.

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

args=scan(text=commandArgs(trailingOnly = TRUE),what="") #Turn the arguments string into a string that R can read and split up.
dir.path <- args[1] #First string to receive is the output directory 
populations <- args[2] #Path to the thinned genomic dataset. Needs to be imputed first.
max.k <- as.numeric(args[3]) #Maximum number of Ks.
env.variable <- args[4] #The environmental factor that this instance of the script evaluates. This is what step 2 parallelises across. 

#Create a matrix in which we can store genomic inflation factor values. 
gif.matrix <- matrix(nrow = 1, ncol = max.k)
rownames(gif.matrix) <- c("gif")
colnames(gif.matrix) <- c(1:max.k)

subdir.path<-(paste(dir.path, "/res/", populations, "/selectingK/EnvironmentalFactor_", env.variable, "/", sep="")) #Create a handy string we can use to shorten directory paths. All directory paths are kept explicit.

for(i in 1:max.k) { #Loop for the number of Ks we wish to evaluate.
  merged.zscores<-read.csv(paste(subdir.path, "K", i, "/LFMM_Zscores_", populations, "_EnvironmentalFactor_", env.variable, "_K", i, "_Chunk_1.csv", sep=""), header = T, sep = ",") #Create an initial dataframe from the first chunk that we can append to.
  print(paste("K ", i,", Chunk 1", sep="")) #Provide user with progress update. Not shown during parallel runs.
  for(j in 2:50) { #Loop for as many chunks as the genomic data is split over. 1 is kept outside of the loop to provide a 'basis'.
# This chunk of code is kept in case leading zeros ever become an issue. It was in the past, and this is one of those problems that 'fixed itself' (which I do not trust. That's why the code is here, just in case...)
#    if (j<10) { 
#      zscores.chunk<-read.csv(paste(dir.path, "/res/", populations, "/selectingK/K", i, "/LFMM_Zscores_", populations, "_K", i,"_Chunk_0", j, ".csv", sep="" ), header = T, sep = ",")
#      print(paste("K ", i,", Chunk ", j, sep=""))  
#  } else {
      zscores.chunk<-read.csv(paste(subdir.path, "K", i, "/LFMM_Zscores_", populations, "_EnvironmentalFactor_", env.variable, "_K", i, "_Chunk_", j, ".csv", sep=""), sep=",", header=T) #Read the next chunk into a dataframe.
      print(paste("K ", i,", Chunk ", j, sep="")) #Give the user a progress update.
    #} #This bracket is still part of the 'spare' code meant to tackle leading zeros.
    merged.zscores<-rbind(merged.zscores, zscores.chunk) #Append current chunk dataframe to overall dataframe.
  }
  print(paste("Calculating GIF for K ", i, sep="")) #Give the user a progress update.

  gif <- median((merged.zscores$zscore)^2)*(qchisq(0.5, df = 1, lower.tail = FALSE)) #Calculate the genomic inflation factor
  gif.matrix[1,i] <- gif #Add the genomic inflation factor to the storage matrix.

  print(paste("Calculating P-values for K ", i, sep="")) #Give the user a progress update.
  results.df <- as.data.frame(merged.zscores) #Rename dataframe
  results.df$pvalues <- pchisq(results.df$zscore^2/gif, df = 1, lower.tail = FALSE) #Calculate P-values and add them to the final results dataframe. 
  
  #Plot the distribution of p-values for this particular K and environmental factor.
   png(paste(subdir.path, "PvalueDistribution_", populations, "EnvironmentalFactor_", env.variable, "_K", i, ".png", sep=""), units = "px", width=2500, height=1500)
    par(mfrow=c(1,2), mar=c(5, 5, 4, 1))
    hist(results.df$pvalues, col="red", main="P-value distribution")
    qqplot(rexp(length(results.df$pvalues), rate=log(10)), -log10(results.df$pvalues), xlab="Expected quantile", pch=19, cex=.4)
    abline(coef=c(0,1))
    dev.off()
  
}
write.table(gif.matrix, file = paste(subdir.path, "/EnvironmentalFactor_", env.variable, "_GIFs.txt", sep = ""), sep = ",", quote = F, row.names = F, col.names = T) #Output the GIFs in a seperate file for each environmental variable. 

