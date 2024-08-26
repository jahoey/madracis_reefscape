#!/bin/bash

# This script takes the unlinked SNP IDs resulting from prune_graph as input. It filters the angsd global snp list so that only the unlinked SNPs remain.

# Define the directories containing the input files, where to store the output, and the angsd global snp list:
input_dir="$HOME/analyses/06_Madracis_672/angsd_MMIRmapped_606noclones_results/ld/output_chromosomes_unlinked_flor"
output_dir="$HOME/analyses/06_Madracis_672/angsd_MMIRmapped_606noclones_results/ld/output_chromosomes_unlinked"
global_snps_file="$HOME/analyses/06_Madracis_672/angsd_MMIRmapped_606noclones_results/global_snps/madracis606_global_snps.txt"

# Loop over each .ld file in the input directory
for file in $input_dir/*.ld; do
  # Extract the base name of the file (e.g., OZ060712_snps_unlinked)
  base_name=$(basename "$file" .ld)
  
  # Extract the chromosome ID (e.g., OZ060712)
  chrom_id=$(echo "$base_name" | cut -d'_' -f1)
  
  # Construct the output file name with just the chromosome ID
  output_file="$output_dir/madracis606_${chrom_id}_unlinked_snps.txt"
  
  # Run the awk command
  awk 'NR==FNR {a[$0]; next} ($1"\t"$2) in a' <(sed 's/:/\t/' "$file") "$global_snps_file" > "$output_file"
  
  echo "Processed $file and saved output to $output_file"
done

