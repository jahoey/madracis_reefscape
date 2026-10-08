#!/bin/bash
#SBATCH --job-name=cat_fst
#SBATCH --cpus-per-task=1
#SBATCH --mem=5G
#SBATCH --time=1:00:00
#SBATCH --output=logs/cat_fst_%A_%a.out
#SBATCH --error=logs/cat_fst_%A_%a.err

BASEDIR=~/analyses/06_Madracis_672/angsd_MMIRmapped_606noclones_results/fst/pairwise_folded

TAXA=(
    MDEC
    MFOR
    MMIR
    MPHA1
    MPHA2
    MPHA3
    MSEN1
    MSEN2
)

# Output file
OUT=$BASEDIR/foldedsfs_pairwise_fst.txt

# Header
printf "Taxon"
for TAXON in "${TAXA[@]}"; do
    printf "\t%s" "$TAXON"
done
printf "\n" > "$OUT"

# Matrix
for TAXON1 in "${TAXA[@]}"; do

    printf "%s" "$TAXON1" >> "$OUT"

    for TAXON2 in "${TAXA[@]}"; do

        # Diagonal
        if [[ "$TAXON1" == "$TAXON2" ]]; then
            FST="NA"

        else
            # Determine alphabetical/order-independent pair
            if [[ "$TAXON1" < "$TAXON2" ]]; then
                POP1=$TAXON1
                POP2=$TAXON2
            else
                POP1=$TAXON2
                POP2=$TAXON1
            fi

            PREFIX="${POP1}_${POP2}"
            FILE="$BASEDIR/$PREFIX/${PREFIX}.folded.fst.stats"

            if [[ ! -f "$FILE" ]]; then
                echo "ERROR: Missing file: $FILE" >&2
                exit 1
            fi

            # Extract the second number = weighted Fst
            FST=$(awk '{print $2}' "$FILE")
        fi

        printf "\t%s" "$FST" >> "$OUT"

    done

    printf "\n" >> "$OUT"
done

echo "Matrix written to:"
echo "$OUT"

cat "$OUT"