#!/bin/bash

for f in $(ls *.bam/genome_results.txt)
do

awk '/bam file/' $f | sed -e 's/^[ \s]*//' | cut -d ' ' -f 4 | cut -c 3- >> bam_names.txt
awk '/mean coverageData/' $f | sed -e 's/^[ \s]*//' -e 's/X//' | cut -d ' ' -f 4 >> coverage_mean.txt
awk '/std coverageData/' $f | sed -e 's/^[ \s]*//' -e 's/X//' | cut -d ' ' -f 4 >> coverage_sd.txt
awk '/number of reads/' $f | sed -e 's/^[ \s]*//' -e 's/,//g' | cut -d ' ' -f 5 >> total_reads.txt
awk '/number of mapped reads/' $f | sed -e 's/^[ \s]*//' -e 's/,//g' | cut -d ' ' -f 6 >> mapped_reads.txt
awk '/mean mapping quality/' $f | sed -e 's/^[ \s]*//' | cut -d ' ' -f 5 >> mapping_quality.txt

#grep -f string_file.txt $f >> msen_coverage.txt
#awk '/mean coverageData|std coverageData|number of reads|number of mapped reads|mean mapping quality/' ORS="\t" $f  >> msen_coverage.txt # all one line, separated by tabs
#echo -e "\n" >> msen_coverage.txt
done

# Essentially cbind each of the files to each other
paste bam_names.txt coverage_mean.txt coverage_sd.txt total_reads.txt mapped_reads.txt mapping_quality.txt >> mapping2mmir_asg_June2024_lane_stats.txt

# Adds header names
cat mapping_names.txt mapping2mmir_asg_June2024_lane_stats.txt >> mapping2mmir_asg_June2024_lane_stats_headers.txt

# Remove intermediate files
rm bam_names.txt
rm coverage_mean.txt
rm coverage_sd.txt
rm total_reads.txt
rm mapped_reads.txt
rm mapping_quality.txt
rm mapping2mmir_asg_June2024_lane_stats.txt
