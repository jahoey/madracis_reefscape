This directory contains steps for handling raw WGS data and analyses for 672 *Madracis* samples

`/angsd_results`: output from angsd analyses

`/scripts`: scripts for bioinformatics and genomic analyses

`/metadata`: associated metadata needed to rename files

# Overview

## Renaming raw sequencing files with more meaningful names

1. After downloading the demultiplexed reads, I used this [file](https://github.com/jahoey/madracis_reefscape/blob/main/WGS/metadata/libs_all_label_convert.txt) to rename them using example code like this. Lanes were kept separate. 

```
while read line ; do
mv $line
done < libs_all_label_convert.txt
```

## Data processing
2. Download the *Madracis auretenra* aka (*Madracis mirabilis*) chromosome-level genome from [here](https://www.ebi.ac.uk/ena/browser/view/GCA_964059495).
3. Trim the reads for the demultiplexed reads and run `fastqc`.
   * [00_trim.sh](https://github.com/jahoey/madracis_reefscape/blob/main/WGS/scripts/00_trim.sh)
4. Map each of the 672 libraries sequenced across three sequencing lanes to the *Madracis mirabilis* genome using `bwa`.
   * [01a_bwa_mem_pe_lane1.sh](https://github.com/jahoey/madracis_reefscape/blob/main/WGS/scripts/01a_bwa_mem_pe_lane1.sh)
   * [01b_bwa_mem_pe_lane2.sh](https://github.com/jahoey/madracis_reefscape/blob/main/WGS/scripts/01b_bwa_mem_pe_lane2.sh)
   * [01c_bwa_mem_pe_lane3.sh](https://github.com/jahoey/madracis_reefscape/blob/main/WGS/scripts/01c_bwa_mem_pe_lane3.sh)
5. Sort, quality filter, and index the aligned reads with `samtools`.
   * [02a_sam2bam_lane1.sh](https://github.com/jahoey/madracis_reefscape/blob/main/WGS/scripts/02a_sam2bam_lane1.sh)
   * [02b_sam2bam_lane2.sh](https://github.com/jahoey/madracis_reefscape/blob/main/WGS/scripts/02b_sam2bam_lane2.sh)
   * [02c_sam2bam_lane3.sh](https://github.com/jahoey/madracis_reefscape/blob/main/WGS/scripts/02c_sam2bam_lane3.sh)
6. Run `qualimap` to understand lane level mapping quality (optional).
   * [03a_quality_reports.sh](https://github.com/jahoey/madracis_reefscape/blob/main/WGS/scripts/03a_quality_reports.sh)
   * [03b_cat_coverage.sh](https://github.com/jahoey/madracis_reefscape/blob/main/WGS/scripts/03b_cat_coverage.sh)
7. Merge the sample bam files across the three sequencing lanes.
   * [04_merge_bams.sh](https://github.com/jahoey/madracis_reefscape/blob/main/WGS/scripts/04_merge_bams.sh)
8. Deduplicate (`MARKDUPLICATES`), clip overlapping read pairs (`clipOverlap`), and reindex the bams before rerunning `qualimap` to determine library level mapping quality and coverage.
   * [05a_dedup.sh](https://github.com/jahoey/madracis_reefscape/blob/main/WGS/scripts/05a_dedup.sh)
   * [05b_quality_reports.sh](https://github.com/jahoey/madracis_reefscape/blob/main/WGS/scripts/05b_quality_reports.sh)
   * [05c_cat_coverage.sh](https://github.com/jahoey/madracis_reefscape/blob/main/WGS/scripts/05c_cat_coverage.sh)
   * [05d_coverage.R](https://github.com/jahoey/madracis_reefscape/blob/main/WGS/scripts/05d_coverage.R)

## Genotype likelihoods and SNP calling for 672 coral libraries
9. Use angsd to do some basic filtering and then calculate genotype likelihoods, determine allele frequencies, and call SNPs across all 672 libraries. This requires a [bamlist](https://github.com/jahoey/madracis_reefscape/blob/main/WGS/angsd_results/bam_lists/angsd_bam_list_madracis672.txt) that will be used by `angsd`.
```
BASEDIR=~/analyses/06_Madracis_672/angsd_MMIRmapped_results

ls /home/jhoey/references/asg/June2024/reference_split/* | parallel -j 15 "angsd -b $BASEDIR/angsd_bam_list.txt -ref ~/references/asg/June2024/GCA_964059495-chromosomes.fasta -out $BASEDIR/temp_snps/{/} -rf {} -P 1 -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 -minMapQ 20 -minQ 20 -minInd 336 -setMinDepth 20 -setMaxDepth 67200 -doCounts 1 -GL 1 -doGlf 2 -doMajorMinor 1 -doMaf 2 -SNP_pval 1e-6 -skipTriallelic 1"
```

10. Prepare global SNP list across 672 libraries, index, and then generate covariance and IBS matrix.
```
BASEDIR=~/analyses/06_Madracis_672/angsd_MMIRmapped_results
REFERENCE=~/references/asg/June2024/GCA_964059495-chromosomes.fasta

# Keep first four columns of each SNP subfile
(angsd) jhoey@flor:~/analyses/06_Madracis_672/angsd_MMIRmapped_results$ ls $BASEDIR/temp_snps/*.mafs.gz | parallel -j15 "gunzip -c {} | cut -f 1,2,3,4 | tail -n +2 > {.}.txt"

# Move SNP subfiles to snps_global directory
(angsd) jhoey@flor:~/analyses/06_Madracis_672/angsd_MMIRmapped_results/temp_snps$ mv *.txt ../global_snps/

# Join all the SNP subfiles together into a single file
(angsd) jhoey@flor:~/analyses/06_Madracis_672/angsd_MMIRmapped_results/global_snps$ cat *.txt > madracis672_global_snps.txt

# Index the SNP list
(angsd) jhoey@flor:~/analyses/06_Madracis_672/angsd_MMIRmapped_results/global_snps$ angsd sites index ./madracis672_global_snps.txt

# Generate `angsd` covariance and IBS matrix
(angsd) jhoey@flor:~/analyses/06_Madracis_672/angsd_MMIRmapped_results/global_snps$ angsd -b $BASEDIR/angsd_bam_list.txt -anc $REFERENCE -out $BASEDIR/global_snps/all_snps -GL 1 -doGlf 2 -doMajorMinor 3 -doMAF 1 -doPost 1 -doIBS 1 -doCounts 1 -doCov 1 -makeMatrix 1 -sites $BASEDIR/global_snps/madracis672_global_snps.txt
```

## Clone identification and removal
11. Plot clustering dendrogram & distribution of 1-IBS for each species separately to determine 95% quartile cutoff for known technical and sample replicates. Use this to determine putative clone clusters.
    * [06_clone_id.R](https://github.com/jahoey/madracis_reefscape/blob/main/WGS/scripts/06_clone_id.R)
12. Merge the technical, sample, and putative clone group bams. Remove the intermediate bams. This reduces the dataset to 606 *Madracis* colonies.
    * [06a_merge_sample_bams.sh](https://github.com/jahoey/madracis_reefscape/blob/main/WGS/scripts/06a_merge_sample_bams.sh)
    * [06b_remove_clone_bams_for_angsd_bamlist.sh](https://github.com/jahoey/madracis_reefscape/blob/main/WGS/scripts/06b_remove_clone_bams_for_angsd_bamlist.sh)
13. Create a [bamlist](https://github.com/jahoey/madracis_reefscape/blob/main/WGS/angsd_results/bam_lists/angsd_no_clone_bam_list.txt) that does not contain putative clones, reindex bams, and rerun `angsd` again for 606 *Madracis* to calculate genotype likelihoods and call SNPs. Again prepare a global SNP list, index, and generate a new covariance matrix. Replot the resulting clustering dendrogram to ensure that all technical, sample, and putative clone groups were sucessfully dealt with.

```
# Rerun angsd
BASEDIR=~/analyses/06_Madracis_672/angsd_MMIRmapped_606noclones_results

ls /home/jhoey/references/asg/June2024/reference_split/* | parallel -j 15 "angsd -b $BASEDIR/angsd_no_clone_bam_list.txt -ref ~/references/asg/June2024/GCA_964059495-chromosomes.fasta -out $BASEDIR/temp_snps/{/} -rf {} -P 1 -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 -minMapQ 20 -minQ 20 -minInd 303 -setMinDepth 20 -setMaxDepth 60600 -doCounts 1 -GL 1 -doGlf 2 -doMajorMinor 1 -doMaf 2 -SNP_pval 1e-6 -skipTriallelic 1"

# Prepare new global SNP list
BASEDIR=~/analyses/06_Madracis_672/angsd_MMIRmapped_606noclones_results
REFERENCE=~/references/asg/June2024/GCA_964059495-chromosomes.fasta

# Keep first four columns of each SNP subfile
(angsd) jhoey@flor:~/analyses/06_Madracis_672/map2mmir_asg_June2024/bams_merged_dedup_clones_removed$ ls $BASEDIR/temp_snps/*.mafs.gz | parallel -j15 "gunzip -c {} | cut -f 1,2,3,4 | tail -n +2 > {.}.txt"

# Move SNP subfiles to snps_global directory
(angsd) jhoey@flor:~/analyses/06_Madracis_672/map2mmir_asg_June2024/bams_merged_dedup_clones_removed$ mv *.txt ../global_snps/

# Join all the SNP subfiles together into a single file
(angsd) jhoey@flor:~/analyses/06_Madracis_672/angsd_MMIRmapped_606noclones_results/global_snps$ cat *.txt > madracis606_global_snps.txt

# Index the new snps list
(angsd) jhoey@flor:~/analyses/06_Madracis_672/angsd_MMIRmapped_606noclones_results/global_snps$ angsd sites index ./madracis606_global_snps.txt

# Generate new `angsd` covariance matrix
(angsd) jhoey@flor:~/analyses/06_Madracis_672/angsd_MMIRmapped_606noclones_results/global_snps$ angsd -b $BASEDIR/angsd_no_clone_bam_list.txt -anc $REFERENCE -out $BASEDIR/global_snps/all_snps -GL 1 -doGlf 2 -doMajorMinor 3 -doMAF 1 -doPost 1 -doIBS 1 -doCounts 1 -doCov 1 -makeMatrix 1 -sites $BASEDIR/global_snps/madracis606_global_snps.txt
```

## Initial clustering analysis
14. Plot PCA for 606 *Madracis* using the global SNP list.
    * [08_pca_coral.R](https://github.com/jahoey/madracis_reefscape/blob/main/WGS/scripts/08_pca_coral.R)

## Estimate linkeage disequilibrium (LD) and LD prune the global SNP list
15. Estimate LD using a subset of sites and `ngsLD`.
```
# Prepare input files for ngsLD
BASEDIR=~/analyses/06_Madracis_672/angsd_MMIRmapped_606noclones_results

# Make input files using one SNP in every 10 SNPs. This will reduce the computation issues with estimating LD: goes from 27582048 SNPs to 2758204 for LD estimation
# Remove first line and first 3 columns for every SNP
(base) jhoey@flor:~/analyses/06_Madracis_672/angsd_MMIRmapped_606noclones_results$ zcat $BASEDIR/global_snps/all_snps.beagle.gz | awk 'NR % 10 == 0' | cut -f 4- | gzip > $BASEDIR/ld/all_snps.beagle.ld.gz

# Remove first line and keep first 2 columns for every SNP
(base) jhoey@flor:~/analyses/06_Madracis_672/angsd_MMIRmapped_606noclones_results$ zcat $BASEDIR/global_snps/all_snps.mafs.gz | cut -f 1,2 | awk 'NR % 10 == 0' |sed 's/:/_/g' | gzip > $BASEDIR/ld/all_snps.mafs.ld.pos.gz

# Estimate LD
~/programs/ngsLD/ngsLD --geno $BASEDIR/ld/all_snps.beagle.ld.gz --pos $BASEDIR/ld/all_snps.mafs.ld.pos.gz --probs --n_ind 606 --n_sites 2758204 --max_kb_dist 100 --n_threads 60 --rnd_sample 0.05 --out $BASEDIR/ld/all_snps.ld

# Plot LD decay curve
BASEDIR=~/analyses/06_Madracis_672/angsd_MMIRmapped_606noclones_results

(base) jhoey@tdobz:~$ ls $BASEDIR/ld/all_snps.ld | Rscript --vanilla --slave ~/programs/ngsLD/scripts/fit_LDdecay.R --out $BASEDIR/ld/ld_plot.pdf --fit_level 100
```

16. Prepare input files for calcuating pairwise LD. Next, divide these files up by chromosome. Based on estimated LD decay, finally calculate LD for SNPs on each chromosome using `ngsLD`.
```
# Prepare the input files for calculating pairwise LD
BASEDIR=~/analyses/06_Madracis_672/angsd_MMIRmapped_606noclones_results

# Make input files for all SNPs: 27582048 SNPs
# Remove first line and first 3 columns for every SNP
(base) jhoey@flor:~/analyses/06_Madracis_672/angsd_MMIRmapped_606noclones_results$ zcat $BASEDIR/global_snps/all_snps.beagle.gz | awk 'NR>1' | cut -f 4- | gzip > $BASEDIR/ld/all_snps.beagle.ld.gz

# Remove first line and keep first 2 columns for every SNP
(base) jhoey@flor:~/analyses/06_Madracis_672/angsd_MMIRmapped_606noclones_results$ zcat $BASEDIR/global_snps/all_snps.mafs.gz | cut -f 1,2 | awk 'NR>1' | sed 's/:/_/g' | gzip > $BASEDIR/ld/all_snps.mafs.ld.pos.gz
```
  * [07a_divide_mafs_files.sh](https://github.com/jahoey/madracis_reefscape/blob/main/WGS/scripts/07a_divide_mafs_files.sh)
  * [07b_divide_beagle_files.sh](https://github.com/jahoey/madracis_reefscape/blob/main/WGS/scripts/07b_divide_beagle_files.sh)
  * [07c_run_ngsLD_chromosomes.sh](https://github.com/jahoey/madracis_reefscape/blob/main/WGS/scripts/07c_run_ngsLD_chromosomes.sh)

17. Use `prune_graph` to prune SNPs for each chromosome. This generates a list of independent SNPs for each chromosome. For an example of how to run `prune_graph`, see below
```
BASEDIR=~/analyses/06_Madracis_672/angsd_MMIRmapped_606noclones_results

cat $BASEDIR/ld/output_chromosomes/OZ060699_snps.ld | ./target/release/prune_graph --weight-field "column_7" --weight-filter "column_3 <= 12000 && column_7 >= 0.5" -n 4 --out $BASEDIR/ld/output_chromosomes_unlinked_flor/OZ060699_snps_unlinked.ld
```

19. Merge the chromosome-level SNPs lists and use this to filter the `angsd` output.
20. Generate `angsd` covariance matrix using LD-pruned SNPs.
21. Generate `PCangsd` covariance matrix for LD-pruned SNPs.

## Clustering analysis
21. Using the LD-pruned SNPs, plot PCA for 606 *Madracis*
22. Run `ngsAdmix` for 10 replicates of K = 6-12

