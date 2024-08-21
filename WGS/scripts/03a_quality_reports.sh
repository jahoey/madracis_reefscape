#!/bin/bash

# Requires bwa conda environment

mkdir ./quality_results

ls ./merged*.overlapclipped.bam | parallel -j12 "qualimap bamqc -nt 1 --java-mem-size=3000M -bam {} -outdir quality_results/{}/"
