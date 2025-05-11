# Peter Muchina October 2024
# This script builds the binary reference panel for QUILT for a single chromosome
#More information on QUILT: https://github.com/rwdavies/QUILT
# It builds reference bins based on the chunks
# Note: I use QUILT2_prepare_reference.R since it takes the VCF file directly for QUILT1, its different
# to generate reasonable chunks use: dat <- QUILT::quilt_chunk_map("Chr1", "genetic_map.txt") i.e Chromosome and the genetic map
# save the output: write.table(dat, file = "chunks_output.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE) 

set -e  # Exit immediately if a command exits with a non-zero status

# Define the chromosome to analyze (update per chromosome if needed)
# You can pass the chromosome as an argument via sbatch instead of hardcoding
Chr="Chr2"

# Define paths, VCF file, genetic map, chunks, output directory
VCF_file="/${Chr}/QUILTv1.0.5/${Chr}_reference_panel.bcf"
Genetic_map="/${Chr}/QUILTv1.0.5/genetic_map.txt" # (optional)
chunks_file="/${Chr}/QUILTv1.0.5/chunks.${Chr}.txt"
Out_dir="/Quilt_${Chr}/"                         

# Ensure output directory exists
mkdir -p "$Out_dir"

# Get the chunk from the array task ID
chunk_info=$(tail -n +2 "$chunks_file" | sed -n "$((SLURM_ARRAY_TASK_ID + 1))p")
chunk=$(echo "$chunk_info" | cut -f1)
chr=$(echo "$chunk_info" | cut -f2)
region=$(echo "$chunk_info" | cut -f3)

# Extract region start and end
regionStart=$(echo "$region" | cut -d":" -f2 | cut -d"-" -f1)
regionEnd=$(echo "$region" | cut -d"-" -f2)

# Check if output exists, skip if found
if [ ! -f "$Out_dir/RData/QUILT_prepared_reference.${Chr}.$regionStart.$regionEnd.RData" ]; then
  echo "Processing chunk $chunk for $Chr (region $regionStart-$regionEnd)"

  # Run QUILT2_prepare_reference.R
  QUILT2_prepare_reference.R \
  --outputdir="$Out_dir" \
  --chr="$Chr" \
  --nGen=100 \
  --regionStart="$regionStart" \
  --regionEnd="$regionEnd" \
  --buffer=500000 \
  --reference_vcf_file="$VCF_file" \
  --genetic_map_file="$Genetic_map" || { echo "Failed on chunk $chunk"; exit 1; }

else
  echo "Skipping chunk $chunk for $Chr as output already exists."
fi

echo "Job completed at $(date '+%d_%m_%y_%H_%M_%S')"
