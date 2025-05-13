#Step 1: Quality control

cd 2_fastqc

# Create soft links for all the files in 1_data
for fastq in ../Data/*.fastq.gz; do
  ln -s $fastq
done

# Run fastqc on the files
fastqc -t 4 *.fastq

