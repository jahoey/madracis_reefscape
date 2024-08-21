This directory contains steps for handling raw WGS data and analyses for 672 *Madracis* samples

**/angsd_results**: output from angsd analyses

**/scripts**: scripts for bioinformatics and genomic analyses

**/metadata**: associated metadata needed to rename files

# Overview

## Renaming raw sequencing files with more meaningful names

1. After downloading the demultiplexed reads, I used this [file](https://github.com/jahoey/madracis_reefscape/blob/main/WGS/metadata/libs_all_label_convert.txt) to rename them using example code like this. Lanes were kept separate. 

> while read line ; do
> mv $line
> done < libs_all_label_convert.txt

## Bioinformatics

2. Download the *Madracis auretenra* aka (*Madracis mirabilis*) chromosome-level genome from [here](https://www.ebi.ac.uk/ena/browser/view/GCA_964059495)
3. Trim the reads for the demultiplexed reads
   * 
4. Align each of the 672 libraries sequenced across three sequencing lanes to the *Madracis mirabilis* genome
   * [01a_bwa_mem_pe_lane1.sh](https://github.com/jahoey/madracis_reefscape/blob/main/WGS/scripts/01a_bwa_mem_pe_lane1.sh)
   * [01b_bwa_mem_pe_lane2.sh](https://github.com/jahoey/madracis_reefscape/blob/main/WGS/scripts/01b_bwa_mem_pe_lane2.sh)
   * [01c_bwa_mem_pe_lane3.sh](https://github.com/jahoey/madracis_reefscape/blob/main/WGS/scripts/01c_bwa_mem_pe_lane3.sh)
