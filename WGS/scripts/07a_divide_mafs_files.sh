#!/bin/bash

input_file="all_snps.mafs.ld.pos.gz"
suffix="mafs.ld.pos.gz"
output_dir="output_chromosomes"

mkdir -p "$output_dir"

# Create an associative array to store the first occurrence line numbers
declare -A first_occurrence

# Process the input file
zcat "$input_file" | awk -v suffix="$suffix" -v output_dir="$output_dir" '
{
    # Extract the first column value
    split($1, arr, "|")
    key = arr[2]

    # Append the line to the corresponding file based on the first column value
    file = output_dir "/" key "." suffix
    print $0 >> file

    # Check and store the line number of the first occurrence of the key
    if (!(key in first_occurrence)) {
        first_occurrence[key] = NR
    }
}
END {
    for (key in first_occurrence) {
        print key ": " first_occurrence[key]
    }
}
' > first_occurrences.txt

# Print the first occurrences line numbers
cat first_occurrences.txt

