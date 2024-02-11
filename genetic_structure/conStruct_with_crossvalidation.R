# conStruct spatial and non-spatial analysis
# script adjusted from https://github.com/capoony/DrosEU_pipeline#d-calculation-of-unbiased-population-genetics-estimators-tajimas-pi-wattersons-theta-and-tajimas-d
library(conStruct)

# Set working directory
setwd("/home/rstudio/working/zurich_workshop/analysis/construct/20240107/")

# Genetic Data input
gen <- read.table(file="../../../data/AllPools_woOutgroups_subsampled0.005_miss0_ADP100_lfmm.txt",
                   header=TRUE,
                   sep = "\t")

gen_lfmm <- t(gen)
poolnames <- rownames(gen_lfmm)
poolnames <- gsub(pattern="^Pool_", replacement = "", x=poolnames)
rownames(gen_lfmm) <- poolnames

Freq <- gen_lfmm

# Load coordinates of each sampling site and geographic distances in kilometers
CoordRaw=read.csv("../../../data/AllPools_woOutgroups_geographic_coords.csv",header=T)
Coord=as.matrix(CoordRaw[,2:3])
colnames(Coord)=c("Lon","Lat")
DistRaw=read.csv("../../../data/AllPools_woOutgroups_geographic_dist_matrix.csv",
                 header = TRUE)
Dist=as.matrix(DistRaw)
Dist <- Dist[ ,-1]
colnames(Dist) <- poolnames
rownames(Dist) <- poolnames

# test values of K ranging from 1 to 10 in 8-fold replication with cross-validation
my.xvals <- x.validation(train.prop = 0.9,
    n.reps = 8,
    K = 1:14,
    freqs = as.matrix(Freq),
    data.partitions = NULL,
    geoDist = Dist,
    coords = Coord,
    prefix = "example2",
    n.iter = 1e3,
    make.figs = TRUE,
    save.files = TRUE,
    parallel = TRUE,
    n.nodes = 12)

# load both the results for the spatial and non-spation models
sp.results <- as.matrix(
    read.table("example2_sp_xval_results.txt",
    header = TRUE,
    stringsAsFactors = FALSE))
nsp.results <- as.matrix(
    read.table("example2_nsp_xval_results.txt",
    header = TRUE,
    stringsAsFactors = FALSE))

# format results from the output list
sp.results <- Reduce("cbind",lapply(my.xvals,function(x){unlist(x$sp)}),init=NULL)
nsp.results <- Reduce("cbind",lapply(my.xvals,function(x){unlist(x$nsp)}),init=NULL)
sp.CIs <- apply(sp.results,1,function(x){mean(x) + c(-1.96,1.96) * sd(x)/length(x)})
nsp.CIs <- apply(nsp.results,1,function(x){mean(x) + c(-1.96,1.96) * sd(x)/length(x)})

# then, plot cross-validation results for K=1:10 with 8 replicates and visualize results with confidence interval bars
pdf("cross-validate-sp.pdf",width=4,height=12)
plot(rowMeans(sp.results),
pch=19,col="blue",
ylab="predictive accuracy",xlab="values of K",
ylim=range(sp.CIs),
main="spatial cross-validation results")
segments(x0 = 1:nrow(sp.results),
y0 = sp.CIs[1,],
x1 = 1:nrow(sp.results),
y1 = sp.CIs[2,],
col = "blue",lwd=2)
dev.off()

# plot all Admixture plots for values of K ranging from 1:10
for (i in seq(1,4,1)){
      my.run <- conStruct(spatial = TRUE,
            K = i,
            freqs = as.matrix(Freq),
            geoDist = Dist,
            coords = Coord,
            prefix = paste("test_",i,seq=""),
            n.chains = 1,
            n.iter = 1e3,
            make.figs = T,
            save.files = T)

      admix.props <- my.run$chain_1$MAP$admix.proportions
      pdf(paste("STRUCTURE_",i,"_1.pdf"),width=8,height=4)
      make.structure.plot(admix.proportions = admix.props,
            sample.names=CoordRaw$X,
            mar = c(6,4,2,2),
            sort.by=NULL)
      dev.off()

      # plot map with pie-charts showing proportion admixture  
      pdf(paste("AdmPIE_",i,"_1.pdf"),width=9,height=8)
      maps::map(xlim = range(Coord[,1]) + c(-5,5), ylim = range(Coord[,2])+c(-2,2), col="gray")
      make.admix.pie.plot(admix.proportions = admix.props,
            coords = Coord,
            add = TRUE)
      box()
      axis(1)
      axis(2)
      dev.off()
}

# load both the results for the spatial and non-spation models
sp.results <- as.matrix(
    read.table("example2_sp_xval_results.txt",
    header = TRUE,
    stringsAsFactors = FALSE))
nsp.results <- as.matrix(
    read.table("example2_nsp_xval_results.txt",
    header = TRUE,
    stringsAsFactors = FALSE))

# format results from the output list
sp.results <- Reduce("cbind",lapply(my.xvals,function(x){unlist(x$sp)}),init=NULL)
nsp.results <- Reduce("cbind",lapply(my.xvals,function(x){unlist(x$nsp)}),init=NULL)
sp.CIs <- apply(sp.results,1,function(x){mean(x) + c(-1.96,1.96) * sd(x)/length(x)})
nsp.CIs <- apply(nsp.results,1,function(x){mean(x) + c(-1.96,1.96) * sd(x)/length(x)})

# then, plot cross-validation results for K=1:10 with 8 replicates and visualize results with confidence interval bars
pdf("non-spatial-cross-validate-sp.pdf",width=4,height=12)
plot(rowMeans(nsp.results),
pch=19,col="blue",
ylab="predictive accuracy",xlab="values of K",
ylim=range(nsp.CIs),
main="spatial cross-validation results")
segments(x0 = 1:nrow(nsp.results),
y0 = nsp.CIs[1,],
x1 = 1:nrow(nsp.results),
y1 = nsp.CIs[2,],
col = "blue",lwd=2)
dev.off()

# plot all Admixture plots for values of K ranging from 1:10
for (i in seq(1,4,1)){
      my.run <- conStruct(spatial = FALSE,
            K = i,
            freqs = as.matrix(Freq),
            geoDist = Dist,
            coords = Coord,
            prefix = paste("test_",i,seq=""),
            n.chains = 1,
            n.iter = 1e3,
            make.figs = T,
            save.files = T)

      admix.props <- my.run$chain_1$MAP$admix.proportions
      pdf(paste("STRUCTURE_nsp_",i,"_1.pdf"),width=8,height=4)
      make.structure.plot(admix.proportions = admix.props,
            sample.names=CoordRaw$X,
            mar = c(6,4,2,2),
            sort.by=NULL)
      dev.off()

      # plot map with pie-charts showing proportion admixture  
      pdf(paste("AdmPIE_nsp_",i,"_1.pdf"),width=9,height=8)
      maps::map(xlim = range(Coord[,1]) + c(-5,5), ylim = range(Coord[,2])+c(-2,2), col="gray")
      make.admix.pie.plot(admix.proportions = admix.props,
            coords = Coord,
            add = TRUE)
      box()
      axis(1)
      axis(2)
      dev.off()
}

