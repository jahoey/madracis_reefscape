#!/bin/bash

#Requires dedup conda environment: Picard Tools MarkDuplicates, BamUtil clipOverlap, parallel

BASEDIR=~/analyses/06_Madracis_672/map2mmir_asg_June2024
OUTPUT_DIR=$BASEDIR/bams_merged_dedup

ls $BASEDIR/bams_merged/merged*.bam > $BASEDIR/bams_merged/angsd_bam_list.txt
mkdir -p "$OUTPUT_DIR"

# First task: MarkDuplicates
parallel -j 30 --bar '
    bam_file={}  
    file_without_path=$(basename "$bam_file")
    sample_name=${file_without_path%.minq20.fixmate.sorted.bam}
    output_file='"$OUTPUT_DIR"'/${sample_name}.dedup.bam

    java -jar ~/programs/miniconda3/pkgs/picard-2.26.10-hdfd78af_0/share/picard-2.26.10-0/picard.jar MarkDuplicates \
        -I "$bam_file" \
        -O "$output_file" \
        -M '"$OUTPUT_DIR"'/${sample_name}.dupstat.txt \
        --VALIDATION_STRINGENCY SILENT \
        --REMOVE_DUPLICATES true' ::: "$BASEDIR"/bams_merged/merged*.bam

# Second task: bam clipOverlap
parallel -j 30 --bar '
    bam_file={}  
    file_without_path=$(basename "$bam_file")
    sample_name=${file_without_path%.dedup.bam}
    output_file='"$OUTPUT_DIR"'/${sample_name}.dedup.overlapclipped.bam

    ~/programs/miniconda3/pkgs/bamutil-1.0.15-h43eeafb_4/bin/bam clipOverlap \
        --in "$bam_file" \
        --out "$output_file" \
        --stats' ::: "$OUTPUT_DIR"/*dedup.bam

# Third task: index the merged, deduplicated, clipped bams
ls '"$OUTPUT_DIR"'/*overlapclipped.bam | parallel -j15 "samtools index {}" #indexing results in .bai file

#########
# This worked by running line by line
#BASEDIR=~/analyses/02_madracis_wgs/map90_2mdec_spades/bams

#mkdir -p $BASEDIR\_dedup

#for bam_file in $BASEDIR/merged*.bam
#do
#    file_without_path=$( basename "$bam_file" )
#    sample_name=${file_without_path%.fixmate.sorted.markdup.bam}
#    output_file=${BASEDIR}_dedup/${sample_name}.dedup.bam

#    java -jar ~/programs/miniconda3/pkgs/picard-2.26.10-hdfd78af_0/share/picard-2.26.10-0/picard.jar MarkDuplicates \
#	-I $bam_file \
#	-O $output_file \
#	-M ${BASEDIR}_dedup/${sample_name}.dupstat.txt
#	--VALIDATION_STRINGENCY SILENT \
#	--REMOVE_DUPLICATES true \

#done

#sleep 10

#for bam_file in ~/analyses/02_madracis_wgs/map90_2mdec_spades/bams_dedup/*bam
#do
#   file_without_path=$( basename "$bam_file" )
#   sample_name=${file_without_path%.dedup.bam}
#   output_file=${BASEDIR}_dedup/${sample_name}.dedup.overlapclipped.bam

#   ~/programs/miniconda3/pkgs/bamutil-1.0.15-h43eeafb_4/bin/bam clipOverlap \
#	--in $bam_file \
#	--out $output_file \
#	--stats

#done
