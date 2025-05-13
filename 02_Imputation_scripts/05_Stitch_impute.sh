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
coverages=("0.5x" "1x" "3x")
num_coverages=${#coverages[@]}
chunks_file="/QUILTv1.0.5/chunks.${Chr}.txt"

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

# Define position file and output directory dynamically
posfile="/${Chr}/positions.txt"
Bamlist="/${Coverage}/bamlist"
Out_dir="/${Chr}/${Coverage}/STITCH/"

# Create the output directory if it doesn't exist
mkdir -p ${Out_dir}

echo "Job started for ${Chr} with $Coverage at $(date '+%d_%m_%y_%H_%M_%S')"

# Run the STITCH command for the current chromosome
STITCH.R --chr=${Chr} \
         --bamlist=${Bamlist} \
         --posfile=${posfile} \
         --regionStart="${regionStart}" \
         --regionEnd="${regionEnd}" \
         --buffer=500000 \
         --outputdir=${Out_dir} \
         --output_filename="${Out_dir}stitch.${Chr}.${regionStart}.${regionEnd}.vcf.gz" \
         --K=10 \
         --nGen=100 \
         --nCores=1

echo "Job completed for ${Chr} with $Coverage at $(date '+%d_%m_%y_%H_%M_%S')"