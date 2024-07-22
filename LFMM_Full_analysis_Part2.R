###############################################
### Running LFMM for ACORN. Credit:         ###
### Benjamin Dauphin and Christian Rellstab ###
### Modified for paralel analysis by        ###
### Lars Littmann                           ###
###############################################

####README####

#This script can be used to process raw z-scores from lfmm to get final results.
#The script merges z-scores from 50 chunks of the unthinned genomic dataset.
#The script processes only one environmental variable. It should be run in parallel by environmental factor.

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

args=scan(text=commandArgs(trailingOnly = TRUE),what="")
dir.path <- args[1] #First string to receive is the output directory
populations <- args[2] #Path to the thinned genomic dataset. Needs to be imputed first.
set.k <- as.numeric(args[3]) #The value of K that was set for this analysis. 
env.variable<-args[4] #The environmental factor this scirpt operates on (parallel to other environmental factors)
full.gen <-args[5] #Path to the full, unthinned genomic dataset. Needed for plotting towards the end of this script. 
dat.env <-args[6] #Path to the environmental dataset.
fdr<-args[7] #The significance threshold that needs to be reached for a SNP to be included in the output. 

print("checkpoint 0")

gif.matrix <- matrix(nrow = 1, ncol = 1)
rownames(gif.matrix) <- c("gif")
colnames(gif.matrix) <- c("1")

print("checkpoint 1")

fdr.output <- fdr
fdr.thres <- c(0.05, 0.01, 0.001)

nb.asso.q <- matrix(nrow=1, ncol=3)
rownames(nb.asso.q) <- "value"
colnames(nb.asso.q) <- c("q0.05", "q0.01", "q0.001")

print("checkpoint 2")

nb.asso.k <- matrix(nrow=1, ncol=1)
rownames(nb.asso.k) <- env.variable
colnames(nb.asso.k) <- paste("K1")

print("checkpoint 3")

###Merge the Z-scores of all 50 chunks of unthinned genetic data. This is the main output of Part 1. 

setwd(paste(dir.path, "/res/", populations, "/Full_analysis_K_", set.k, "/environment_", env.variable, "/", sep="")) #Set the working directory that contains all the chunk files.

merged.zscores<-read.csv(paste("LFMM_Zscores_", populations, "environment_", env.variable, "_K", set.k, "_Chunk_1.csv", sep=""), header = T, sep = ",") #Start the data frame off with chunk 1.
print(paste("Environmental factor ", env.variable, ", Chunk 1", sep="")) #Give the user some output to track the script's progress. Remains invisible in parallel applications.
for(j in 2:50) {

#This piece of code is defunct for now, but I would still like to keep this solution to the problem of leading zeros in the code, in case it's ever relevant.    
#    if (j<10) {
#      zscores.chunk<-read.csv(paste(dir.path, "/res/", populations, "/Full_analysis_K_", set.k, "/environment_", env.variable, "/LFMM_Zscores_", populations, "_K", set.k,"_Chunk_0", j, ".csv", sep="" ), header = T, sep = ",")
#      print(paste("K ", i,", Chunk ", j, sep=""))  
#  } else {
    
    zscores.chunk<-read.csv(paste("LFMM_Zscores_", populations, "environment_", env.variable, "_K", set.k,"_Chunk_", j, ".csv", sep=""), sep=",", header=T) #Load the z-scores of the chunk we are looping over into a df.
    print(paste("Environmental factor ", env.variable, ", Chunk ", j, sep=""))  #Give the user some output to track the script's progress. Remains invisible in parallel applications.
  
  #} #Still part of the solution for leading zeros. This bracket is commented out and can be ignored. 
  
    merged.zscores<-rbind(merged.zscores, zscores.chunk) #Append the chunk we are looping over at the bottom of the df that we have compiled up to now.
  
  }
  print("Calculating GIF")

  gif <- median((merged.zscores$zscore)^2)/(qchisq(0.5, df = 1, lower.tail = FALSE)) #Calculate the genomic inflation factor based on all the Z-scores.
  gif.matrix[1,1] <- gif #Store the genomic inflation factor for later.

  print("Calculating P-values") #Give the user an indicator of the progress. Not visible during parallel processing. 
  results.df <- as.data.frame(merged.zscores) #Rename the dataframe that we used to merge all the z-scores.

  results.df$pvalues <- pchisq(results.df$zscore^2/gif, df = 1, lower.tail = FALSE) #Calculate all the P-values based on z-scores corrected with the GIF.
  qv.lfmm2 <- qvalue::qvalue(as.vector(results.df$pvalues), fdr.level=fdr.output) #Calculate q-values.
  results.df$qvalues <- qv.lfmm2$qvalues #Add q-values to the dataframe with our results. 

  #write.table(results.df, paste("LFMM_AllResults_env_", env.variable, "_K", set.k, ".csv", sep=""), sep=",", row.names=F, col.names=T, quote=F) # save all information per SNP. Commented out to save disk space.

  #Plot the distribution of p-values. We plot them as PNGs, because the PDFs contained so many points that it took 5 minutes to load each one.
   png(paste(dir.path, "/res/", populations, "/Full_analysis_K_", set.k, "/environment_", env.variable, "/PvalueDistribution_", populations, "_environment_", env.variable, "_K", set.k, ".png", sep=""), units = "px", width=2500, height=1500)
    par(mfrow=c(1,2), mar=c(5, 5, 4, 1))
    hist(results.df$pvalues, col="red", main="P-value distribution")
    qqplot(rexp(length(results.df$pvalues), rate=log(10)), -log10(results.df$pvalues), xlab="Expected quantile", pch=19, cex=1)
    abline(coef=c(0,1))
    dev.off()

   # Summarise results and generate Manhattan plot at different significance thresholds
    for (x in 1:length(fdr.thres)) {
      q <- NULL; w <- NULL; v <- NULL; candidate <- NULL
      q <- fdr.thres[x]
      v <- results.df$qvalues
      w <- which(sort(v) <= q)
      
      print(paste("Assessing candidates for fdr threshold ", q, sep=""))

      candidate <- order(results.df$qvalues, decreasing=F)[w]
      png(paste("ManhattanPlot_env", env.variable, "_K", set.k, "_q", q, ".png", sep=""), units="px", width=2500, height=1500)
      par(mfrow=c(1,1), mar=c(5, 5, 4, 1))
      manhattan <- plot(-log10(results.df$pvalues), main=paste("Manhattan plot | env", env.variable, " | K", set.k, " | q=", q, sep=""), cex.main=1.2, xlab="Locus", ylab="-Log(P-value)", cex=.7, col="grey")
      points(candidate, -log10(results.df$pvalues)[candidate], pch=19, cex=1, col="red")
      dev.off()

      if(q == fdr.output) {
        nb.asso.k[1,1] <- length(candidate)
      }

      nb.asso.q["value",x] <- length(candidate)
      candidate <- results.df[candidate,]
      
      if(NCOL(candidate) == 1) {
        write.table(t(candidate), paste("CandidatesOrdered_env_", env.variable, "_K", set.k, "_q", q, ".csv", sep=""), sep=",", row.names=F, col.names=T, quote=F)
      } else {
        if(NCOL(candidate) > 1) {
          write.table(candidate[order(as.numeric(candidate$pvalue), decreasing=F),], paste("CandidatesOrdered_env", env.variable, "_K", set.k, "_q", q, ".csv", sep=""), sep=",", row.names=F, col.names=T, quote=F)
        } else {
          write.table(candidate, paste("CandidatesOrdered_env_", env.variable, "_K", set.k, "_q", q, ".csv", sep=""), sep=",", row.names=T, col.names=F, quote=F)
        }
      }
    }
    write.table(nb.asso.q, paste("AssociationsNb_env_", env.variable, "_K", set.k, "_q.csv", sep=""), sep=",", row.names=F, col.names=T, quote=F)

write.table(gif.matrix, file = paste("GIF_", env.variable, "_K_", set.k, ".txt", sep = ""), sep = ",", quote = F, row.names = F, col.names = T)

#Remove a bunch of objects to free up memory space.
rm(results.df)
rm(merged.zscores)
rm(qv.lfmm2)

#### Import the full genetic and environmental dataset for plotting purposes ####
gen.data <- read.table(paste(full.gen, sep=""), header=T, sep="\t", row.names = "chrom_pos")

#transpose dataframe
gen <- as.data.frame(t(gen.data))

#Store the names of all the SNPs in a vector.
snp.info <- as.vector(colnames(gen))

env.data <- read.table(paste(dat.env, sep=""), header=T, sep=",")

#Make the population names (Plot_ID) the row names of the environmental dataset.
env <- env.data[c(env.variable)]
rownames(env) <- env.data$Plot_ID

#reduce the env set to the gen set. Only the populations that occur in both the environmental and genetic dataset remain. 
rownames(gen) <- gsub("X","",rownames(gen)) #Accounts for a difference in the formating of the population names. 
env.reduced <- env[rownames(env) %in% rownames(gen),] 
env <- env.reduced
identical(as.character(rownames(env)),as.character(rownames(gen)))#check NAs.

#Rename datasets for easier shorthand. 
X <- as.matrix(env) #The environmental data

#### Plot SNPs significantly associated to environmental conditions ####

for (j in 1:NCOL(X)) {
    candidate.list <- read.table(paste(dir.path, "/res/", populations, "/Full_analysis_K_", set.k, "/environment_", env.variable, "/CandidatesOrdered_env", env.variable, "_K", set.k, "_q", fdr.output, ".csv", sep=""), row.names="SNPid", header=T, sep=",")
    tmp.gen <- gen
    colnames(tmp.gen) <- snp.info # "SNPid" # snp.info$SNPid
    tmp.gen.t <- t(tmp.gen); rownames(tmp.gen.t) <- colnames(tmp.gen); colnames(tmp.gen.t) <- rownames(tmp.gen)
    candidate.gen <- merge(candidate.list, tmp.gen.t, by="row.names"); colnames(candidate.gen)[1] <- "SNPid"
    candidate.gen.ordered <- candidate.gen[order(candidate.gen$pvalue, decreasing=F),]
    if (NROW(candidate.list) > 0) {
      genotypes <- t(candidate.gen.ordered[,5:NCOL(candidate.gen.ordered)])
      pdf(paste(dir.path, "/res/", populations, "/Full_analysis_K_", set.k, "/environment_", env.variable, "/PlotOfSignificantCandidates_", populations, "_env_", env.variable, "_K", set.k, "_q", fdr.output, ".pdf", sep=""), width=10, height=10)
      par(mfrow=c(2,2), mar=c(5, 5, 1, 1))
      for (p in 1:NROW(candidate.list)) {
        plot(X[,j], genotypes[,p], pch=20, cex=1.5, xlab=colnames(X)[j], ylab="Genotype frequency [-]",
             cex.lab=1.2, ylim=c(0, 1), cex.main=1.2, col=alpha("blue",0.2))
        zs <- candidate.list$zscore[p]
        pv <- candidate.list$pvalue[p]
        qv <- candidate.list$qvalue[p]
        legend(min(X[,j]), 1.5, legend=c(paste(rownames(candidate.list)[p]), paste("z-score =", format(zs, digits=3)),
                                           paste("p-value =", format(pv, digits=2)), paste("q-value =", format(qv, digits=4))), bty='n', cex=0.7)
      }
      dev.off()
 } 
}


