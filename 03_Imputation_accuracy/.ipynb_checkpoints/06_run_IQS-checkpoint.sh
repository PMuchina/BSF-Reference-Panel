# Exit immediately if a command exits with a non-zero status
set -e

echo "Job started at $(date '+%d_%m_%y_%H_%M_%S')"

python3.7 IQS.py

echo "Job completed at $(date '+%d_%m_%y_%H_%M_%S')"
