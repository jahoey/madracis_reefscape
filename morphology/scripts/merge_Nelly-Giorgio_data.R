library (here)

giorgio <- read.delim(here("../Giorgio_morphology/Girogio_Corallite_Morphology_Measurements_Aug2024.csv"), header = TRUE, sep = ',') #682 x 25

nelly <- read.delim(here("Nelly_data_May2024_updatedspeciesMPHA_rescaled.csv"), header = TRUE, sep = ',') #691 x 44

gn_data <- merge(nelly, giorgio[,-c(5:7)], by = c('ID', 'Site', 'Depth', 'Date'), all.x = TRUE, all.y = FALSE) #691 x 62

write.csv(gn_data, here("madracis_morphology_GN_merged.csv"), col.names = TRUE, row.names = FALSE, quote = FALSE)

a <- which(rowSums(is.na(gn_data[10:24])) == 15) #septa
b <- which(rowSums(is.na(gn_data[25:39])) == 15) #ic dist
c <- which(rowSums(is.na(gn_data[45:59])) == 15) #c area

no_data <- intersect(intersect(a,b),c) # these are all corals where no data was taken for the 3 different traits

# Filter the merged Nelly-Giorgio dataset and remove rows that have no data at all
gn_data_filtered <- gn_data[-no_data,] # 669 x 62
