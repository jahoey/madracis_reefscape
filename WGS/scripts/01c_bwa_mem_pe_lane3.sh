#!/bin/bash

ls lane3*val*.fq.gz | sed -e 's/_R1_val_1.fq.gz//' -e 's/_R2_val_2.fq.gz//' | sort -u | parallel -j 30 'bwa mem -t 2 /home/jhoey/references/asg/June2024/GCA_964059495-chromosomes.bwaidx {}_R1_val_1.fq.gz {}_R2_val_2.fq.gz > /home/jhoey/analyses/06_Madracis_672/map2mmir_asg_June2024/{}.sam'
