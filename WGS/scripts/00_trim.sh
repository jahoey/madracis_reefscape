#!/bin/bash

# Requires trim_galore and parallel
ls *.fastq.gz | sed -e 's/_R1.fastq.gz//' -e 's/_R2.fastq.gz//' | sort -u | parallel -j 24 'trim_galore -q 20 --length 30 --stringency 5 --phred33 --paired --retain_unpaired --fastqc {}_R1.fastq.gz {}_R2.fastq.gz'
