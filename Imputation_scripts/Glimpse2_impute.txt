# Peter Muchina October 2024
# For more information on Glimpse usage, here is the link: https://odelaneau.github.io/GLIMPSE/
# This is the imputation step for Glimpse2
# Pass the chromosome as an argument via sbatch instead of hardcoding
# Check if a chromosome argument has been provided
if [ -z "$1" ]; then
  echo "Error: No chromosome argument provided."
  echo "Usage: sbatch script_name.sh <chromosome>"
  exit 1
fi

# Set the chromosome from the input argument
Chr=$1

set -e

echo "Job started at $(date '+%d_%m_%y_%H_%M_%S')"

# Define coverage levels and chunks file
coverages=("1x")
Chunks="${Chr}/GLIMPSE2/chunks.${Chr}.txt"
num_coverages=${#coverages[@]}
num_chunks=$(wc -l < "$Chunks")

# Calculate chunk and coverage index based on SLURM_ARRAY_TASK_ID
chunk_index=$((SLURM_ARRAY_TASK_ID / num_coverages))
coverage_index=$((SLURM_ARRAY_TASK_ID % num_coverages))

# Get coverage level for this task
Coverage=${coverages[$coverage_index]}

# Read the line corresponding to chunk_index + 2 (to skip the header) from Chunks file
# In this case since chunks from Glimpse2 do not have a header, I use + 1
LINE=$(awk "NR==$((chunk_index + 1))" "$Chunks")
if [ -z "$LINE" ]; then
    echo "Error: No data in line $((chunk_index + 2)) of $Chunks"
    exit 1
fi

# Extract information from line
# Refer to Glimpse2 tutorial https://odelaneau.github.io/GLIMPSE/
ID=$(echo $LINE | cut -d" " -f1)
CHR=$(echo $LINE | cut -d" " -f2)
IRG=$(echo $LINE | cut -d" " -f3)
REGS=$(echo ${IRG} | cut -d":" -f2 | cut -d"-" -f1)
REGE=$(echo ${IRG} | cut -d"-" -f2)
       
# Set paths for this coverage level and chunk
REF="/Glimpse2_${Chr}/Reference_panel_Glimpse2_${CHR}_${REGS}_${REGE}.bin"
Bamlist="/bamlist"
Output_dir="/${Chr}/Glimpse2"

# Ensure output directory exists
mkdir -p "${Output_dir}"

# Check if output already exists to avoid reprocessing
Output_file="${Output_dir}/${Coverage}_imputed_${CHR}_${REGS}_${REGE}.bcf"
if [ ! -f "$Output_file" ]; then
    echo "Processing $CHR, Coverage $Coverage, Chunk $chunk_index, Region ${REGS}-${REGE}"

    # Run GLIMPSE2 phase command
    GLIMPSE2_phase \
        --bam-list ${Bamlist} \
        --reference ${REF} \
        --output ${Output_file} \
        --threads 1
else
    echo "Skipping $CHR, Coverage $Coverage, Chunk $chunk_index: Output already exists."
fi

echo "Job completed for ${CHR}, Coverage ${Coverage}, Chunk ${chunk_index} at $(date '+%d_%m_%y_%H_%M_%S')"
