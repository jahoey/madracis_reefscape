#!/bin/bash

# I recommend breaking up the Ks across different computer clusters to make the runs more efficient

# Define your variables
BASEDIR=~/analyses/06_Madracis_672/angsd_MMIRmapped_606noclones_results
LIKES_FILE=$BASEDIR/global_snps/unlinked_snps.beagle.gz
MIN_MAF=0.05
THREADS=20
MAX_JOBS=10  # Define the maximum number of jobs to run simultaneously

# Function to execute NGSadmix
run_ngsadmix() {
    K=$1
    i=$2
    OUTPUT_FILE="$BASEDIR/admixture/madracis606_ngsAdmix_K${K}_${i}_out"
    NGSadmix -likes $LIKES_FILE -K $K -o $OUTPUT_FILE -P $THREADS -minMaf $MIN_MAF -minInd 425
    echo "NGSadmix finished for K=$K, run $i. Output saved as $OUTPUT_FILE"
}

# Loop through K values from 6 to 12
for K in {6..16}; do
    # Loop through 10 runs for each value of K
    for i in {1..10}; do
        # Run NGSadmix in the background
        run_ngsadmix $K $i &

        # Limit the number of simultaneous jobs
        if (( $(jobs | wc -l) >= MAX_JOBS )); then
            wait -n  # Wait for any background job to finish
        fi
    done
done

# Wait for all background jobs to finish
wait

