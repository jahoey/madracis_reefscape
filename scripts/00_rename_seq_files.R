library(here)
library(stringr)

#### Read in metadata ####
metadata_for_names <- read.csv(here("raw_seq_data/metadata_for_renaming.csv"), header = TRUE) # metadata for renaming raw seq files

#### Read in samples names contained in library ####
libs_original_names <- read.csv(here("raw_seq_data/from_Davis/lib1-lib7_original_names.csv"), header = FALSE)

#### Determine which extracts are in library ####
libs_original_names_split <- str_split_fixed(libs_original_names$V1, '_', 5)
libs_ext <- unique(libs_original_names_split[,1])

#### Subsample metadata to sample extracts that are in library ####
libs_metadata <- metadata_for_names[metadata_for_names$extraction_id %in% libs_ext,]
libs_metadata_ordered <- libs_metadata[order(libs_metadata$extraction_id),]

#### Create new labels for the appropriate samples ####
libs_new_labels <- paste(libs_metadata_ordered$species, libs_metadata_ordered$site, libs_metadata_ordered$depth, libs_metadata_ordered$sample_id_standardized, libs_metadata_ordered$extraction_id, sep = "_")
libs_new_labels_R1 <- paste(libs_new_labels, "_R1.fastq.gz", sep = '')
libs_new_labels_R2 <- paste(libs_new_labels, "_R2.fastq.gz", sep = '')

# Join R1 and R2 labels together
libs_new_labels_all <- c(rbind(libs_new_labels_R1, libs_new_labels_R2))

# Now combine with original labels
libs_new_labels_all2 <- cbind(libs_original_names[1:1344,], libs_new_labels_all) #there are two files that couldn't be split out by barcode

#### Write old and new labels as a text file ####
write.table(libs_new_labels_all2, here("raw_seq_data/from_Davis", "libs_all_label_convert.txt"), row.names = F, col.names = F, quote = F)


