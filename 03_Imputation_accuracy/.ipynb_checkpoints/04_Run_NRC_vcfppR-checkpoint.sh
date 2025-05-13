# Exit immediately if a command exits with a non-zero status
set -e

# Define the tools and coverages
tools=("Glimpse2" "Quilt" "Stitch_Ref" "Stitch")
coverages=("0.5x" "1x" "3x")

# Calculate the tool and coverage for this job
tool_idx=$((SLURM_ARRAY_TASK_ID / 3))
coverage_idx=$((SLURM_ARRAY_TASK_ID % 3))

tool=${tools[$tool_idx]}
coverage=${coverages[$coverage_idx]}

echo "Job started at $(date '+%d_%m_%y_%H_%M_%S')"

# Run the R script with the current tool and coverage
Rscript NRC_vcfppR.R  $tool $coverage

echo "Job completed at $(date '+%d_%m_%y_%H_%M_%S')"
