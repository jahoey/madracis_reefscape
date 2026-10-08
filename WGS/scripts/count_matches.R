library(here)
library(tidyr) 

### BLAST ###
### Calculate number of reads matching to each genome from a BLAST analysis using a local database of Pocillopora, B. minutum and M. mirabilus mito genomes ###
# 384000 lines/192000 reads in lib12blast.fasta  (2000 reads for 96 individuals each; 1000 from R1 and 1000 from R2)
matches <- read.table(here("blast_results", "matches_localdb.txt"), header = FALSE, sep = '\t')

dim(matches) # 71998 x 2
table(matches$V2)

species_matches <- separate(matches, 2, c('genus', 'species', 'rest'), sep = '\\s', extra = 'merge') 
species_read_counts <- table(species_matches$species)
species_read_prop <- species_read_counts/192000

total_reads <- 952004170 # 1904008340 lines/2 lines per read

estimated_species_reads <- total_reads * species_read_prop # 356991647
non_mapped_reads <- total_reads - sum(estimated_species_reads) # 595012523

# Explore BLAST results against NCBI nt database for 2000 queries for 96 library 1 Madracis (192000) #
total_queries <- 192000
nt_matches <- read.table(here("blast_results", "library1_madracis_ncbi_hits_summary.txt"), header = FALSE, sep = '\t', col.names = c('count', 'phylum'))

sum(nt_matches$count)
nt_matches[which(nt_matches$phylum == 'Cnidaria'),1]/sum(nt_matches$count)
nt_matches[which(nt_matches$phylum == 'Arthropoda'),1]/sum(nt_matches$count)
nt_matches[which(nt_matches$phylum == 'Chordata'),1]/sum(nt_matches$count)

nt_matches[which(nt_matches$phylum == 'Cnidaria'),1]/total_queries


### Fastq screen ###
### Calculate number of reads matching to each genome from a fastq screen analysis for selected coral colonies
# Number of lines counted using wc -l
fastq_number_of_lines <- read.table(here("fastq_screen_results", "fastq.line.counts.txt"), col.names = c("line_count", "coral_id"))
fastq_number_of_reads <- cbind.data.frame(fastq_number_of_lines$line_count/4, fastq_number_of_lines$coral_id)
  
# Read in fastq screen results
fastq_names <- list.files(here("fastq_screen_results"), pattern = "*screen.txt", full.names = TRUE)
fastq_files <- lapply(fastq_names, function(x)read.table(x, header=F, fill=T))
genomes <- fastq_files[[1]]$V1[-c(1,10)]
fastq_files_array <- array(as.numeric(unlist(fastq_files)), dim=c(10, 12, 7))
fastq_files_array_cleaned <- fastq_files_array[-c(1,10),-1,] # remove rows 1 and 10 and column 1
rownames(fastq_files_array_cleaned) <- genomes
colnames(fastq_files_array_cleaned) <- c('reads_processed', 'unmapped', 'percent_unmapped', 'one_hit_one_genome', 'percent_one_hit_one_genome', 'multi_hit_one_genome', 'percent_multi_hit_one_genome', 'one_hit_multi_genome', 'percent_one_hit_multi_genome', 'multi_hits_multi_genome', 'percent_multi_hits_multi_genome')

fastq_files_array_cleaned_percent <- fastq_files_array_cleaned[,-c(1:4,6,8,10),] # keep percent columns only. Percents based off of screening a subset of the total number of reads.

hit_prop_sum <- data.frame(matrix(NA, nrow = 7, ncol = 8))
no_hit_prop <- data.frame(matrix(NA, nrow = 7, ncol = 8))
for (i in 1:7) {
  hit_prop_sum[i,] <- (rowSums(fastq_files_array_cleaned_percent[,,i]))/100 # Calculate total percent, then convert to proportion
  no_hit_prop[i,] <- (fastq_files_array_cleaned[,3,i])/100
}
colnames(hit_prop_sum) <- genomes
rownames(hit_prop_sum) <- fastq_number_of_lines$coral_id
colnames(no_hit_prop) <- genomes
rownames(no_hit_prop) <- fastq_number_of_lines$coral_id

estimated_fastq_reads <- hit_prop_sum * fastq_number_of_reads$`fastq_number_of_lines$line_count/4` # mapped reads for an individual
estimated_fastq_unmapped_reads <- no_hit_prop* fastq_number_of_reads$`fastq_number_of_lines$line_count/4` # unmapped reads for an individual

# Sanity check
total_reads <- estimated_fastq_reads + estimated_fastq_unmapped_reads # Check that unmapped reads and mapped reads to a given genome is equal to the total number of reads for an individual
check <- vector()
for (j in 1:ncol(estimated_fastq_unmapped_reads)) {
  check[j] <- all.equal(total_reads[,j], fastq_number_of_reads$`fastq_number_of_lines$line_count/4`)
}

# Now write file
all_reads <- cbind(estimated_fastq_reads, estimated_fastq_unmapped_reads)
all_reads <- all_reads[,c(matrix(1:ncol(all_reads), nrow = 2, byrow = T))] #alternating mapped then unmapped for each genome
all_reads_plus_total <- cbind.data.frame(all_reads, fastq_number_of_reads$`fastq_number_of_lines$line_count/4`)
write.table(all_reads_plus_total, here("fastq_screen_results", "fastq_summary_numbers.txt"), row.names = T, col.names = T)
