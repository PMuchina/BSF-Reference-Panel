#!/bin/bash

#Running freebayes in parrallel
reference_genome="path/to/ref"

# Create a list with all your bam files
ls -l "path/to/markdup.bam" | awk '{print $NF}' | xargs -I{} readlink -f {} > bamlist

# Run freebayes
# Here the genome is chunked into regions of 100000 (1mb). 
# Ran on 20 processor (20 jobs will be ran in parallel).
# This can be adjusted based on the size of the genome and the computational resources

freebayes-parallel \
   <(fasta_generate_regions.py ${ref}.fai 100000) 20 \
   --fasta-reference ${reference_genome}  \
   --bam-list bamlist  > output.vcf