library(here)
library(stringr)

#### Read in metadata ####
metadata_for_names <- read.csv(here("raw_seq_data/metadata_for_renaming.csv"), header = TRUE) # metadata for renaming raw seq files

#### Read in samples names contained in library ####
lib1_original_names <- read.csv(here("raw_seq_data/from_Davis/lib1_original_names.csv"), header = FALSE)

#### Determine which extracts are in library ####
lib1_original_names_split <- str_split_fixed(lib1_original_names$V1, '_', 5)
lib1_ext <- unique(lib1_original_names_split[,1])

#### Subsample metadata to sample extracts that are in library ####
lib1_metadata <- metadata_for_names[metadata_for_names$extraction_id %in% lib1_ext,]
lib1_metadata_ordered <- lib1_metadata[order(lib1_metadata$extraction_id),]

#### Create new labels for the appropriate samples ####
lib1_new_labels <- paste(lib1_metadata_ordered$species, lib1_metadata_ordered$site, lib1_metadata_ordered$depth, lib1_metadata_ordered$sample_id_standardized, lib1_metadata_ordered$extraction_id, sep = "_")
lib1_new_labels_R1 <- paste(lib1_new_labels, "_R1.fastq.gz", sep = '')
lib1_new_labels_R2 <- paste(lib1_new_labels, "_R2.fastq.gz", sep = '')

# Join R1 and R2 labels together
lib1_new_labels_all <- c(rbind(lib1_new_labels_R1, lib1_new_labels_R2))

# Now combine with original labels
lib1_new_labels_all2 <- cbind(lib1_original_names, lib1_new_labels_all)

#### Write old and new labels as a text file ####
write.table(lib1_new_labels_all2, here("raw_seq_data/from_Davis", "lib1_label_convert.txt"), row.names = F, col.names = F, quote = F)


