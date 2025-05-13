#!/bin/bash

# Input directory containing BAM files
bam_directory="path/to/bamfiles"

# Output directory for processed BAM files
output_directory="./"

# Check if the input directory exists
if [ ! -d "$bam_directory" ]; then
    echo "Input directory $bam_directory does not exist"
    exit 1
fi

# Create the output directory if it doesn't exist
mkdir -p "$output_directory"

# Iterate over each .bam file in the input directory
for bam in "$bam_directory"*.bam; do
    if [ -f "$bam" ]; then
        base_name=$(basename "${bam%.*}")
        echo "Processing ${bam}"

        java -Xmx10G -jar picard.jar MarkDuplicates \
                    I="${bam}" \
                    O="${output_directory}${base_name}_markdup.bam" \
                    M="${output_directory}${base_name}_markdup_metrics.txt" \
                    OPTICAL_DUPLICATE_PIXEL_DISTANCE=100 \
                    CREATE_INDEX=true \
                    VALIDATION_STRINGENCY=LENIENT \
                    REMOVE_DUPLICATES=true \
                    ASSUME_SORT_ORDER=coordinate

      else
        echo "No .bam files found in $bam_directory"
        exit 1
    fi

done

echo "Completed marking duplicates"


