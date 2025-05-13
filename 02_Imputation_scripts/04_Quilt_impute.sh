# Peter Muchina October 2024
#More information on QUILT: https://github.com/rwdavies/QUILT
#This code uses Quilt to impute across low coverages (0.5x, 1x, 3x)
#Note: I use Quilt1.05 for imputation 
# Pass the chromosome as an argument via sbatch instead of hardcoding
# Check if a chromosome argument has been provided
if [ -z "$1" ]; then
  echo "Error: No chromosome argument provided."
  echo "Usage: sbatch script_name.sh <chromosome>"
  exit 1
fi

# Set the chromosome from the input argument
Chr=$1

set -e   # Exit immediately if a command exits with a non-zero status

echo "Job started at $(date '+%d_%m_%y_%H_%M_%S')"

# Define coverages and paths
coverages=("1x" )
num_coverages=${#coverages[@]}
chunks_file="/${Chr}/QUILTv1.0.5/chunks.${Chr}.txt"

# Calculate chunk and coverage index based on the array task ID
chunk_index=$((SLURM_ARRAY_TASK_ID / num_coverages))
coverage_index=$((SLURM_ARRAY_TASK_ID % num_coverages))

# Read chunk and region based on chunk_index
chunk=$(awk "NR==$((chunk_index + 2))" "$chunks_file" | cut -f1)  # +2 to skip header
region=$(awk "NR==$((chunk_index + 2))" "$chunks_file" | cut -f3)  # Use 3rd column for region
Coverage=${coverages[$coverage_index]}

# Extract start and end positions from region
if [[ "$region" =~ ^[^:]+:([0-9]+)-([0-9]+)$ ]]; then
    regionStart="${BASH_REMATCH[1]}"
    regionEnd="${BASH_REMATCH[2]}"
else
    echo "Error: Invalid region format in $chunks_file on line $((chunk_index + 2))"
    exit 1
fi

# Debug output to verify values
echo "Chunk: $chunk, Region: $region, Start: $regionStart, End: $regionEnd, Coverage: $Coverage"

# Define paths (build reference panel, bamlist of samples to be imputed, and output directory)
Ref="/Quilt_${Chr}/RData/"
Bamlist="bamlist"
Out_dir="/${Chr}/Quilt/"

# Ensure output directory exists
mkdir -p "$Out_dir"

# Check if the output file already exists to avoid reprocessing
if [ ! -f "$Out_dir/quilt.${Chr}.$regionStart.$regionEnd.vcf.gz" ]; then
    echo "Processing $Chr, Coverage $Coverage, Chunk $chunk, Region ${regionStart}-${regionEnd}"

    # Run Quilt2 Imputation
    QUILT.R \
      --prepared_reference_filename="${Ref}QUILT_prepared_reference.${Chr}.${regionStart}.${regionEnd}.RData" \
      --bamlist="${Bamlist}" \
      --method=diploid \
      --chr="$Chr" \
      --regionStart="${regionStart}" \
      --regionEnd="${regionEnd}" \
      --buffer=500000 \
      --output_filename="${Out_dir}/quilt.${Chr}.${regionStart}.${regionEnd}.vcf.gz" \
      --nCores=1
else
    echo "Skipping $Chr, Coverage $Coverage, Chunk $chunk: Output already exists."
fi

echo "Job completed for ${Chr} at $(date '+%d_%m_%y_%H_%M_%S')"