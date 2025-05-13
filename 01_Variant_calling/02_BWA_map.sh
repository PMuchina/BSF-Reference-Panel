#!/bin/bash

#Index reference genome
bwa index reference_genome.fasta

# Path to the reference genome indexed with BWA
reference_genome="path/to/ref"

# Loop through each forward-read file
# Change the file extension if its not R1.fastq.gz || R2.fastq.gz 

for forward_file in ../Data/*.R1.fastq.gz; do
     # Get the base filename without the extension
     base_filename=$(basename -- "$forward_file" | sed 's/.R1.fq.gz//')

     # Generate the paths for the forward and reverse-read files
     reverse_file="../1_data/${base_filename}.R2.fq.gz"

     echo "Processing: ${base_filename}"

     # Perform the alignment using bwa mem, convert the sam to bam, and sort the bam
     # Change the read group according to read description and the sequencing platform used.
     bwa mem -M -t 8 -R "@RG\tID:${base_filename}\tPL:BGISEQ-500\tSM:${base_filename}" \
     "$reference_genome" "$forward_file" "$reverse_file" | \
     samtools view -bS | \
     samtools sort --output-fmt BAM -@ 8 -o "${base_filename}.sorted.bam" -

     # Index the sorted BAM file
     samtools index ${base_filename}.sorted.bam


    echo "Done Processing: ${base_filename}"

done