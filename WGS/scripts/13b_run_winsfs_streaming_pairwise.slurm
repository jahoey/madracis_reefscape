#!/bin/bash
#SBATCH --job-name=winsfs_fst
#SBATCH --array=0-27%2
#SBATCH --cpus-per-task=24
#SBATCH --mem=32G
#SBATCH --time=35-00:00:00
#SBATCH --output=logs/fst_%A_%a.out
#SBATCH --error=logs/fst_%A_%a.err

set -euo pipefail

# --------------------------------------------------
# Activate ANGSD environment
# --------------------------------------------------

source /home/jhoey/programs/miniconda3/etc/profile.d/conda.sh
conda activate angsd

# Add winsfs installed via Cargo to PATH
export PATH="/home/jhoey/.cargo/bin:$PATH"

# --------------------------------------------------
# Paths
# --------------------------------------------------

BASEDIR=~/analyses/06_Madracis_672/angsd_MMIRmapped_606noclones_results/fst

SAFDIR=$BASEDIR/merged_saf
OUTDIR=$BASEDIR/pairwise

# --------------------------------------------------
# Define all 28 pairwise comparisons
# --------------------------------------------------

PAIRS=(
    "MDEC MFOR"
    "MDEC MMIR"
    "MDEC MPHA1"
    "MDEC MPHA2"
    "MDEC MPHA3"
    "MDEC MSEN1"
    "MDEC MSEN2"

    "MFOR MMIR"
    "MFOR MPHA1"
    "MFOR MPHA2"
    "MFOR MPHA3"
    "MFOR MSEN1"
    "MFOR MSEN2"

    "MMIR MPHA1"
    "MMIR MPHA2"
    "MMIR MPHA3"
    "MMIR MSEN1"
    "MMIR MSEN2"

    "MPHA1 MPHA2"
    "MPHA1 MPHA3"
    "MPHA1 MSEN1"
    "MPHA1 MSEN2"

    "MPHA2 MPHA3"
    "MPHA2 MSEN1"
    "MPHA2 MSEN2"

    "MPHA3 MSEN1"
    "MPHA3 MSEN2"

    "MSEN1 MSEN2"
)

# --------------------------------------------------
# Select pair for this SLURM task
# --------------------------------------------------

PAIR=${PAIRS[$SLURM_ARRAY_TASK_ID]}

read POP1 POP2 <<< "$PAIR"

PREFIX=${POP1}_${POP2}

OUT=$OUTDIR/$PREFIX

mkdir -p "$OUT"

SAF1=$SAFDIR/${POP1}.saf.idx
SAF2=$SAFDIR/${POP2}.saf.idx

# --------------------------------------------------
# Print job information
# --------------------------------------------------

echo "=============================================="
echo "Pairwise FST analysis"
echo "=============================================="
echo "Job ID:       $SLURM_JOB_ID"
echo "Array task:   $SLURM_ARRAY_TASK_ID"
echo "Population 1: $POP1"
echo "Population 2: $POP2"
echo "SAF 1:        $SAF1"
echo "SAF 2:        $SAF2"
echo "Output:       $OUT"
echo "CPUs:         $SLURM_CPUS_PER_TASK"
echo "Memory:       32G"
echo "=============================================="

echo ""
echo "Software:"
echo "ANGSD:   $(which angsd)"
echo "realSFS: $(which realSFS)"
echo "winsfs:  $(which winsfs)"
echo ""

# --------------------------------------------------
# Check required programs
# --------------------------------------------------

command -v winsfs >/dev/null 2>&1 || {
    echo "ERROR: winsfs not found in PATH"
    exit 1
}

command -v realSFS >/dev/null 2>&1 || {
    echo "ERROR: realSFS not found in PATH"
    exit 1
}

# --------------------------------------------------
# Check input SAF files
# --------------------------------------------------

if [ ! -f "$SAF1" ]; then
    echo "ERROR: Missing SAF index: $SAF1"
    exit 1
fi

if [ ! -f "$SAF2" ]; then
    echo "ERROR: Missing SAF index: $SAF2"
    exit 1
fi

# --------------------------------------------------
# Define output files
# --------------------------------------------------

SHUF_OUT=$OUT/${PREFIX}.saf.shuf
WINSFS_OUT=$OUT/${PREFIX}.2dsfs.winsfs
SFS_OUT=$OUT/${PREFIX}.2dsfs
FST_PREFIX=$OUT/${PREFIX}

# --------------------------------------------------
# Step 1: Shuffle intersected SAF data
# --------------------------------------------------

echo "=============================================="
echo "Step 1: Shuffling SAF data"
echo "Started: $(date)"
echo "=============================================="

winsfs shuffle \
    --output "$SHUF_OUT" \
    "$SAF1" \
    "$SAF2"

echo "SAF shuffling finished: $(date)"

# Check shuffled file

if [ ! -s "$SHUF_OUT" ]; then
    echo "ERROR: Shuffled SAF file is missing or empty:"
    echo "$SHUF_OUT"
    exit 1
fi

echo "Shuffled SAF created:"
ls -lh "$SHUF_OUT"

# --------------------------------------------------
# Step 2: Estimate 2D SFS in streaming mode
# --------------------------------------------------

echo "=============================================="
echo "Step 2: Estimating 2D SFS with winsfs"
echo "Streaming mode"
echo "Started: $(date)"
echo "=============================================="

# Streaming mode is automatically activated when
# winsfs receives a .saf.shuf file.
#
# Streaming mode can only use one thread.

winsfs \
    -v \
    "$SHUF_OUT" \
    > "$WINSFS_OUT"

echo "winsfs finished: $(date)"

# --------------------------------------------------
# Extract numeric SFS for realSFS
# --------------------------------------------------

echo "Extracting numeric 2D SFS..."

grep -v '^#' \
    "$WINSFS_OUT" \
    > "$SFS_OUT"

# Check that the SFS file is not empty

if [ ! -s "$SFS_OUT" ]; then
    echo "ERROR: Extracted 2D SFS is empty"
    exit 1
fi

# Check that exactly one numeric SFS line was produced

N_LINES=$(wc -l < "$SFS_OUT")

if [ "$N_LINES" -ne 1 ]; then
    echo "ERROR: Expected exactly 1 SFS line, found $N_LINES"
    echo "Contents of SFS output:"
    cat "$SFS_OUT"
    exit 1
fi

echo "2D SFS successfully created:"
echo "$SFS_OUT"

# --------------------------------------------------
# Step 3: Calculate FST index
# --------------------------------------------------

echo "=============================================="
echo "Step 3: Calculating FST index"
echo "Started: $(date)"
echo "=============================================="

realSFS fst index \
    "$SAF1" \
    "$SAF2" \
    -sfs "$SFS_OUT" \
    -fstout "$FST_PREFIX" \
    -whichFst 1 \
    -P "$SLURM_CPUS_PER_TASK"

echo "FST index finished: $(date)"

# --------------------------------------------------
# Step 4: Calculate global FST
# --------------------------------------------------

echo "=============================================="
echo "Step 4: Calculating global FST"
echo "Started: $(date)"
echo "=============================================="

realSFS fst stats \
    "${FST_PREFIX}.fst.idx" \
    > "${FST_PREFIX}.fst.stats"

echo "Global FST finished: $(date)"

# --------------------------------------------------
# Print final results
# --------------------------------------------------

echo "=============================================="
echo "Finished $PREFIX"
echo "Completed: $(date)"
echo "=============================================="

echo ""
echo "Global FST results:"
cat "${FST_PREFIX}.fst.stats"

# --------------------------------------------------
# Optional cleanup
# --------------------------------------------------

# Uncomment the following lines once you are confident
# that the workflow is working and you no longer need
# to retain the shuffled intermediate file.

# echo "Removing shuffled intermediate..."
# rm -f "$SHUF_OUT"
