#!/bin/bash

#Requires bwa conda environment

ls ./lane1*.sam | parallel -j20 "samtools sort -n {} > {.}.bam" #sorts the reads in the sam file in the order they apppear along the chromosome and outputs a bam file; -n sorts by name as required by samtools fixmate

ls ./lane1*.bam | parallel -j20 "samtools view -q 20 -b {} | samtools sort -n - > {.}.minq20.bam" #filter bam files to remove reads mapping to multiple places

ls ./lane1*.minq20.bam | parallel -j20 "samtools fixmate -m {} {.}.fixmate.bam" #samtools fixmate fills in coordinates and flags from name sorted alignment

ls ./lane1*.minq20.fixmate.bam | parallel -j20 "samtools sort {} > {.}.sorted.bam" # then sorts again by genomic coordinates (no -n)

#ls ./*.fixmate.sorted.bam | parallel -j20 "samtools markdup {} {.}.markdup.bam" # removes pcr duplicates

ls ./lane1*.fixmate.sorted.bam | parallel -j20 "samtools index {}" #indexing results in .bai file

ls ./lane1*.fixmate.sorted.bam | parallel -j20 "samtools flagstat {} > {.}.mappingstats.txt"

#cleanup intermediate files, but leave initial bam and final bam/bam.bai files
# hold up on this for now to make sure I have what I need
# rm *.sam
# rm *.fixmate.bam
# rm *.fixmate.sorted.bam

