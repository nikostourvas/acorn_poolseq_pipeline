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
env.variable<-args[4] #The environmental factor this scirpt operates on (parallel to other environmental factors)

gif.matrix <- matrix(nrow = 1, ncol = max.k)
rownames(gif.matrix) <- c("gif")
colnames(gif.matrix) <- c(1:max.k)

fdr.thres <- c(0.05,0.01,0.001)
fdr.output <- 0.05

nb.asso.q <- matrix(nrow=1, ncol=length(fdr.thres))
rownames(nb.asso.q) <- "value"
colnames(nb.asso.q) <- c("q0.05","q0.01","q0.001")

nb.asso.k <- matrix(nrow=1, ncol=as.integer(max.k))
rownames(nb.asso.k) <- colnames(env.variable)
colnames(nb.asso.k) <- paste("K", 1:max.k, sep="")

for(i in 1:max.k) {
  setwd(paste(dir.path, "/res/", populations, "/environment_", env.variable, "/K", i, "/",sep=""))

  merged.zscores<-read.csv(paste("LFMM_Zscores_", populations, "environment_", env.variable, "_K", i, "_Chunk_1.csv", sep=""), header = T, sep = ",")
  print(paste("K ", i,", Chunk 1", sep=""))
  for(j in 2:50) {
#    if (j<10) {
#      zscores.chunk<-read.csv(paste(dir.path, "/res/", populations, "/selectingK/K", i, "/LFMM_Zscores_", populations, "_K", i,"_Chunk_0", j, ".csv", sep="" ), header = T, sep = ",")
#      print(paste("K ", i,", Chunk ", j, sep=""))  
#  } else {
      zscores.chunk<-read.csv(paste("LFMM_Zscores_", populations, "environment_", env.variable, "_K", i,"_Chunk_", j, ".csv", sep=""), sep=",", header=T)
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
  qv.lfmm2 <- qvalue::qvalue(as.vector(results.df$pvalues), fdr.level=fdr.output)
  results.df$qvalues <- qv.lfmm2$qvalues
  #write.table(results.df, paste("LFMM_AllResults_env_", env.variable, "_K", i, ".csv", sep=""), sep=",", row.names=F, col.names=T, quote=F) # save all information per SNP

  #Distribution of p-values 
   png(paste("/PvalueDistribution_", populations, "_environment_", env.variable, "_K", i, ".png", sep=""), units = "px", width=2500, height=1500)
    par(mfrow=c(1,2), mar=c(5, 5, 4, 1))
    hist(results.df$pvalues, col="red", main="P-value distribution")
    qqplot(rexp(length(results.df$pvalues), rate=log(10)), -log10(results.df$pvalues), xlab="Expected quantile", pch=19, cex=1)
    abline(coef=c(0,1))
    dev.off()

   # Summarise results and generate Manhattan plot at different significance thresholds
    for (x in 1:length(fdr.thres)) {
      q <- NULL; w <- NULL; v <- NULL; candidate <- NULL
      q <- fdr.thres[x]
      v <- qv.lfmm2$qvalues
      w <- which(sort(v) <= q)
      
      print(paste(Assessing candidates for K=", i, ", and fdr threshold ", q, sep=""))

      candidate <- order(results.df$qvalues, decreasing=F)[w]
      png(paste("ManhattanPlot_env", env.variable, "_K", i, "_q", q, ".png", sep=""), units="px", width=2500, height=1500)
      par(mfrow=c(1,1), mar=c(5, 5, 4, 1))
      manhattan <- plot(-log10(results.df$pvalues), main=paste("Manhattan plot | env", env.variable, " | K", i, " | q=", q, sep=""), cex.main=1.2, xlab="Locus", ylab="-Log(P-value)", cex=.7, col="grey")
      points(candidate, -log10(results.df$pvalues)[candidate], pch=19, cex=1, col="red")
      dev.off()

      if(q == fdr.output) {
        nb.asso.k[1,i] <- length(candidate)
      }

      nb.asso.q["value",x] <- length(candidate)
      candidate <- results.df[candidate,]
      
      if(NCOL(candidate) == 1) {
        print("In the first possible scenario of the if else statement")
        write.table(t(candidate), paste("CandidatesOrdered_env", env.variable, "_K", i, "_q", q, ".csv", sep=""), sep=",", row.names=F, col.names=T, quote=F)
      } else {
        if(NCOL(candidate) > 1) {
          print("In the second possible scenario of the if else statement")
          write.table(candidate[order(as.numeric(candidate$pvalue), decreasing=F),], paste("CandidatesOrdered_env", env.variable, "_K", i, "_q", q, ".csv", sep=""), sep=",", row.names=F, col.names=T, quote=F)
        } else {
          print("In the third possile scenario of the if else statement")
          write.table(candidate, paste("CandidatesOrdered_env", env.variable, "_K", i, "_q", q, ".csv", sep=""), sep=",", row.names=T, col.names=F, quote=F)
        }
      }
    }
    write.table(nb.asso.q, paste("AssociationsNb_env", env.variable, "_K", i, "_q.csv", sep=""), sep=",", row.names=F, col.names=T, quote=F)
}

write.table(gif.matrix, file = paste(dir.path, "/res/", populations, "/environment_", env.variable, "/", populations, "_GIFs.txt", sep = ""), sep = ",", quote = F, row.names = F, col.names = T)

