library(dplyr)

# Read the original data
counts <- read.table("/data/genetics_tmp/results/vcf/VCF_AllPools_Outgroups/VCF_AllPools_OutgroupsIndelFilteredSNPs_Biallelic_no-miss_ADP100_thin-15000_counts.txt",
                     sep = "\t", header = TRUE)
original_data <- counts[ ,-c(1:4)]

# List of unique pool IDs (e.g., "Pool_102")
pool_ids <- unique(sub("(.ref.cnt|.alt.cnt)$", "", colnames(original_data)))

# Initialize an empty data frame with the same number of rows as the original data
result <- data.frame(matrix(ncol = length(pool_ids), nrow = nrow(original_data)))

# Name the columns of the result data frame
colnames(result) <- pool_ids

# Loop through each pool ID and merge the .ref and .alt columns
for (pool_id in pool_ids) {
  result[[pool_id]] <- paste(original_data[[paste0(pool_id, ".ref.cnt")]], 
                             original_data[[paste0(pool_id, ".alt.cnt")]], 
                             sep = ",")
}

# Remove prefix from pool IDs
colnames(result) <- sub("P01.....ACORN.BOKU.Pl..", "", colnames(result))

# Write the result to a new file
write.table(result, "/data/genetics_tmp/results/vcf/VCF_AllPools_Outgroups/VCF_AllPools_OutgroupsIndelFilteredSNPs_Biallelic_no-miss_ADP100_thin-15000_treemix.txt", 
            quote = FALSE, row.names = FALSE)

