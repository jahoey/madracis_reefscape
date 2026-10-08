#!/bin/bash
#SBATCH --job-name=merge_saf
#SBATCH --array=0-7
#SBATCH --cpus-per-task=4
#SBATCH --mem=20G
#SBATCH --time=04:00:00
#SBATCH --output=logs/merge_saf_%A_%a.out
#SBATCH --error=logs/merge_saf_%A_%a.err

set -euo pipefail

source /home/jhoey/programs/miniconda3/etc/profile.d/conda.sh
conda activate angsd

BASEDIR=~/analyses/06_Madracis_672/angsd_MMIRmapped_606noclones_results/fst

OUTDIR=$BASEDIR/merged_saf

TAXA=(MDEC MFOR MMIR MPHA1 MPHA2 MPHA3 MSEN1 MSEN2)

GROUP=${TAXA[$SLURM_ARRAY_TASK_ID]}

echo "Merging SAF files for $GROUP"

realSFS cat \
    "$BASEDIR"/saf/"$GROUP"/chrm_*.saf.idx \
    -outnames "$OUTDIR/$GROUP.saf" \
    -P 4

echo "Finished $GROUP"
