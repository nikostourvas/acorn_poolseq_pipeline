# Load necessary library
library(fields)
library(tidyverse)

setwd("/home/rstudio/working/zurich_workshop/")

# Read the data from the CSV file
coords <- read.csv("data/20240102_Tree_data.csv")

# Filter according to species and other needed parameters
# coords <- coords %>% 
#   filter(Pool_seq == TRUE) %>% 
#   filter(Species_pop_level == "Q.petraea")

coords <- coords %>% 
  filter(Pool_seq == TRUE)

pools_coords <- coords %>% 
  select(Plot_ID, Plot_Latitude, Plot_Longitude)

coords <- unique(pools_coords)
coords <- coords %>% filter(!Plot_ID %in% c(202,300,304)) # rm pops 202,300,304

# Select the longitude (x) and latitude (y) columns
coordinates <- coords[, c("Plot_Longitude", "Plot_Latitude")]
rownames(coordinates) <- coords$Plot_ID
write.csv(coordinates, file="data/AllPools_woOutgroups_geographic_coords.csv")


# Calculate the distance matrix
dist_matrix_km <- rdist.earth(coordinates, miles=FALSE)

rownames(dist_matrix_km) <- coords$Plot_ID
colnames(dist_matrix_km) <- coords$Plot_ID

# Create a dist object from the lower triangle of the distance matrix
# dist_object <- as.dist(dist_matrix_km)

# Print the dist object
# print(dist_object)

# Export matrix
write.csv(dist_matrix_km, file = "data/AllPools_woOutgroups_geographic_dist_matrix.csv")
