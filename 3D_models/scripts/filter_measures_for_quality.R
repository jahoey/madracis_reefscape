library(dplyr)
library(here)

# Read in Madracis master datasheet
mad_master <- read.csv(here("madracis_master_merged_Oct2025.csv"), header = TRUE)

quant_files <- list.files(here("3D_models/output/hab_measures_03132026/"), pattern = "\\.csv$", full.names = TRUE)

# Loop through each smaller file
for (quant_file in quant_files) {
  # Read the smaller file
  quant <- read.csv(quant_file, header = TRUE)
  
  # Merge with specific columns from the large file
  merged <- merge(quant, mad_master[, c(5, 9, 21, 24)], by.x = 'id', by.y = 'ID', all.x = TRUE)
  
  # Subset to annotations actually on colony
  merged.yes <- merged[merged$Annotation_point_on_actual_coral == 'yes',]
  
  # Extract the basename of the quant file
  output_filename <- paste0(tools::file_path_sans_ext(basename(quant_file)), "_quality.csv")
  
  # Construct the full path to the output file
  output_file_path <- file.path(here("3D_models/output/hab_measures_03132026/", output_filename))
    
  # Write high quality annotations to a new CSV
  write.csv(merged.yes, row.names = FALSE, file = output_file_path)
}

