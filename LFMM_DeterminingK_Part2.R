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

args=commandArgs(trailingOnly = TRUE)
dir.path <- args[1] #First string to receive is the output directory
populations <- args[2] #Path to the thinned genomic dataset. Needs to be imputed first.
max.k <- as.numeric(args[3]) #Maximum number of Ks.

gif.matrix <- matrix(nrow = 1, ncol = max.k)
rownames(gif.matrix) <- c("gif")
colnames(gif.matrix) <- c(1:max.k)

dir.create(paste(dir.path, "/res/", populations, "/selectingK/Pvalue_distributions/", sep=""), recursive = F)

for(i in 1:max.k) {
  merged.zscores<-read.csv(paste(dir.path, "/res/", populations, "/selectingK/K", i, "/LFMM_Zscores_", populations, "_K", i, "_Chunk_1.csv", sep=""), header = T, sep = ",")
  print(paste("K ", i,", Chunk 1", sep=""))
  for(j in 2:50) {
#    if (j<10) {
#      zscores.chunk<-read.csv(paste(dir.path, "/res/", populations, "/selectingK/K", i, "/LFMM_Zscores_", populations, "_K", i,"_Chunk_0", j, ".csv", sep="" ), header = T, sep = ",")
#      print(paste("K ", i,", Chunk ", j, sep=""))  
#  } else {
      zscores.chunk<-read.csv(paste(dir.path, "/res/", populations, "/selectingK/K", i, "/LFMM_Zscores_", populations, "_K", i,"_Chunk_", j, ".csv", sep=""), sep=",", header=T)
      print(paste("K ", i,", Chunk ", j, sep=""))
    #}
    merged.zscores<-rbind(merged.zscores, zscores.chunk)
  }
  print(paste("Calculating GIF for K ", i, sep=""))

  gif <- median((merged.zscores$zscore)^2)*(qchisq(0.5, df = 1, lower.tail = FALSE))
  gif.matrix[1,i] <- gif

  print(paste("Calculating P-values for K ", i, sep=""))
  results.df <- as.data.frame(merged.zscores)
  results.df$pvalues <- pchisq(results.df$zscore^2/gif, df = 1, lower.tail = FALSE)
  
  #Distribution of p-values 
   png(paste(dir.path, "/res/", populations, "/selectingK/Pvalue_distributions/PvalueDistribution_", populations, "_K", i, ".png", sep=""), units = "px", width=2500, height=1500)
    par(mfrow=c(1,2), mar=c(5, 5, 4, 1))
    hist(results.df$pvalues, col="red", main="P-value distribution")
    qqplot(rexp(length(results.df$pvalues), rate=log(10)), -log10(results.df$pvalues), xlab="Expected quantile", pch=19, cex=.4)
    abline(coef=c(0,1))
    dev.off()
  
}

write.table(gif.matrix, file = paste(dir.path, "/res/", populations, "/selectingK/", populations, "_GIFs.txt", sep = ""), sep = ",", quote = F, row.names = F, col.names = T)

