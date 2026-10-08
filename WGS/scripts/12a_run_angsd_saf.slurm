#!/bin/bash
#SBATCH --job-name=madracis_saf
#SBATCH --array=0-7%2
#SBATCH --cpus-per-task=15
#SBATCH --mem=100G
#SBATCH --time=48:00:00
#SBATCH --output=logs/saf_%A_%a.out
#SBATCH --error=logs/saf_%A_%a.err

set -euo pipefail

# --------------------------------------------------
# Activate ANGSD conda environment
# --------------------------------------------------

source /home/jhoey/programs/miniconda3/etc/profile.d/conda.sh
conda activate angsd

echo "Conda environment: $CONDA_DEFAULT_ENV"
echo "ANGSD: $(which angsd)"
echo "GNU parallel: $(which parallel)"

# --------------------------------------------------
# Paths
# --------------------------------------------------

BASEDIR=~/analyses/06_Madracis_672/angsd_MMIRmapped_606noclones_results/fst

REF=~/references/asg/June2024/GCA_964059495-chromosomes.fasta

REFDIR=~/references/asg/June2024/reference_split

BAMDIR=$BASEDIR/bamlists
OUTDIR=$BASEDIR/saf

mkdir -p "$OUTDIR"
mkdir -p "$BASEDIR/logs"

# --------------------------------------------------
# Taxon corresponding to SLURM array task
# --------------------------------------------------

TAXA=(MDEC MFOR MMIR MPHA1 MPHA2 MPHA3 MSEN1 MSEN2)

GROUP=${TAXA[$SLURM_ARRAY_TASK_ID]}

BAMLIST=$BAMDIR/${GROUP}.bamlist

GROUPOUT=$OUTDIR/$GROUP

mkdir -p "$GROUPOUT"

# --------------------------------------------------
# Print information
# --------------------------------------------------

echo "=============================================="
echo "Madracis ANGSD SAF"
echo "=============================================="
echo "Taxon:       $GROUP"
echo "BAM list:    $BAMLIST"
echo "Reference:   $REF"
echo "Reference regions: $REFDIR"
echo "Output:      $GROUPOUT"
echo "=============================================="

echo "Number of BAMs:"
wc -l "$BAMLIST"

# --------------------------------------------------
# Run ANGSD on each reference chunk
# --------------------------------------------------

ls "$REFDIR"/chrm_* | parallel -j 15 \
    "angsd \
    -b $BAMLIST \
    -ref $REF \
    -anc $REF \
    -rf {} \
    -out $GROUPOUT/{/} \
    -P 1 \
    -uniqueOnly 1 \
    -remove_bads 1 \
    -only_proper_pairs 1 \
    -trim 0 \
    -C 50 \
    -baq 1 \
    -minMapQ 20 \
    -minQ 20 \
    -doCounts 1 \
    -GL 1 \
    -doSaf 1"

echo "=============================================="
echo "Finished: $GROUP"
echo "=============================================="
