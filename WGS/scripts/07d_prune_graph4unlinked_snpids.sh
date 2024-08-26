#!/bin/bash

# prune_graph takes the output of ngsLD for each chromosome as input. This input file is the pairwise LD calculations given user specifications. prune_graph determines which linked SNPs to remove, resulting in an output file containing a list of unlinked SNPs that then needs to be applied to filter the global SNP list from angsd.
# This step is very time intensive, so it is advised to run prune_graph for each chromosome as a separate job, which still takes quite a bit of time.

BASEDIR=~/analyses/06_Madracis_672/angsd_MMIRmapped_606noclones_results

cat $BASEDIR/ld/output_chromosomes/OZ060699_snps.ld | ./target/release/prune_graph --weight-field "column_7" --weight-filter "column_3 <= 12000 && column_7 >= 0.5" -n 4 --out $BASEDIR/ld/output_chromosomes_unlinked_flor/OZ060699_snps_unlinked.ld

cat $BASEDIR/ld/output_chromosomes/OZ060700_snps.ld | ./target/release/prune_graph --weight-field "column_7" --weight-filter "column_3 <= 12000 && column_7 >= 0.5" -n 4 --out $BASEDIR/ld/output_chromosomes_unlinked_flor/OZ060700_snps_unlinked.ld

cat $BASEDIR/ld/output_chromosomes/OZ060701_snps.ld | ./target/release/prune_graph --weight-field "column_7" --weight-filter "column_3 <= 12000 && column_7 >= 0.5" -n 4 --out $BASEDIR/ld/output_chromosomes_unlinked_flor/OZ060701_snps_unlinked.ld

cat $BASEDIR/ld/output_chromosomes/OZ060702_snps.ld | ./target/release/prune_graph --weight-field "column_7" --weight-filter "column_3 <= 12000 && column_7 >= 0.5" -n 4 --out $BASEDIR/ld/output_chromosomes_unlinked_flor/OZ060702_snps_unlinked.ld

cat $BASEDIR/ld/output_chromosomes/OZ060703_snps.ld | ./target/release/prune_graph --weight-field "column_7" --weight-filter "column_3 <= 12000 && column_7 >= 0.5" -n 4 --out $BASEDIR/ld/output_chromosomes_unlinked_flor/OZ060703_snps_unlinked.ld

cat $BASEDIR/ld/output_chromosomes/OZ060704_snps.ld | ./target/release/prune_graph --weight-field "column_7" --weight-filter "column_3 <= 12000 && column_7 >= 0.5" -n 4 --out $BASEDIR/ld/output_chromosomes_unlinked_flor/OZ060704_snps_unlinked.ld

cat $BASEDIR/ld/output_chromosomes/OZ060705_snps.ld | ./target/release/prune_graph --weight-field "column_7" --weight-filter "column_3 <= 12000 && column_7 >= 0.5" -n 4 --out $BASEDIR/ld/output_chromosomes_unlinked_flor/OZ060705_snps_unlinked.ld

cat $BASEDIR/ld/output_chromosomes/OZ060706_snps.ld | ./target/release/prune_graph --weight-field "column_7" --weight-filter "column_3 <= 12000 && column_7 >= 0.5" -n 4 --out $BASEDIR/ld/output_chromosomes_unlinked_flor/OZ060706_snps_unlinked.ld

cat $BASEDIR/ld/output_chromosomes/OZ060707_snps.ld | ./target/release/prune_graph --weight-field "column_7" --weight-filter "column_3 <= 12000 && column_7 >= 0.5" -n 4 --out $BASEDIR/ld/output_chromosomes_unlinked_flor/OZ060707_snps_unlinked.ld

cat $BASEDIR/ld/output_chromosomes/OZ060708_snps.ld | ./target/release/prune_graph --weight-field "column_7" --weight-filter "column_3 <= 12000 && column_7 >= 0.5" -n 4 --out $BASEDIR/ld/output_chromosomes_unlinked_flor/OZ060708_snps_unlinked.ld

cat $BASEDIR/ld/output_chromosomes/OZ060709_snps.ld | ./target/release/prune_graph --weight-field "column_7" --weight-filter "column_3 <= 12000 && column_7 >= 0.5" -n 4 --out $BASEDIR/ld/output_chromosomes_unlinked_flor/OZ060709_snps_unlinked.ld

cat $BASEDIR/ld/output_chromosomes/OZ060710_snps.ld | ./target/release/prune_graph --weight-field "column_7" --weight-filter "column_3 <= 12000 && column_7 >= 0.5" -n 6 --out $BASEDIR/ld/output_chromosomes_unlinked_flor/OZ060710_snps_unlinked.ld

cat $BASEDIR/ld/output_chromosomes/OZ060711_snps.ld | ./target/release/prune_graph --weight-field "column_7" --weight-filter "column_3 <= 12000 && column_7 >= 0.5" -n 6 --out $BASEDIR/ld/output_chromosomes_unlinked_flor/OZ060711_snps_unlinked.ld

cat $BASEDIR/ld/output_chromosomes/OZ060712_snps.ld | ./target/release/prune_graph --weight-field "column_7" --weight-filter "column_3 <= 12000 && column_7 >= 0.5" -n 6 --out $BASEDIR/ld/output_chromosomes_unlinked_flor/OZ060712_snps_unlinked.ld

cat $BASEDIR/ld/output_chromosomes/OZ060713_snps.ld | ./target/release/prune_graph --weight-field "column_7" --weight-filter "column_3 <= 12000 && column_7 >= 0.5" --out $BASEDIR/ld/output_chromosomes_unlinked_flor/OZ060713_snps_unlinked.ld
