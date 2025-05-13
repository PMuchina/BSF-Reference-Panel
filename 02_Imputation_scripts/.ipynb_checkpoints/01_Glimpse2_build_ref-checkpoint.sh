# Peter Muchina October 2024
# For more information on Glimpse usage, here is the link: https://odelaneau.github.io/GLIMPSE/
#This script builds the binary reference panel
#Create a binary reference panel. The reference panel is converted into GLIMPSE2’s binary file format
#The output is a binary reference panel file for each imputed chunk.
# Pass the chromosome as an argument via sbatch instead of hardcoding
# Check if a chromosome argument has been provided
if [ -z "$1" ]; then
  echo "Error: No chromosome argument provided."
  echo "Usage: sbatch script_name.sh <chromosome>"
  exit 1
fi

# Set the chromosome from the input argument
Chr=$1

echo "Job started at $(date '+%d_%m_%y_%H_%M_%S')"

REF="/${Chr}/GLIMPSE2/${Chr}_reference_panel.bcf"
MAP="/${Chr}/GLIMPSE2/genetic_map.txt.gz"
Chunks="/${Chr}/GLIMPSE2/chunks.${Chr}.txt"    
Ref_panel="/Glimpse2_${Chr}/"

# Ensure output directory exists
mkdir -p "$Ref_panel"

while IFS="" read -r LINE || [ -n "$LINE" ];
do
  printf -v ID "%02d" $(echo $LINE | cut -d" " -f1)
  IRG=$(echo $LINE | cut -d" " -f3)
  ORG=$(echo $LINE | cut -d" " -f4)

 GLIMPSE2_split_reference --reference ${REF} --map ${MAP} --input-region ${IRG} --output-region ${ORG} --output ${Ref_panel}Reference_panel_Glimpse2 --threads 30
done < $Chunks

echo "Job completed at $(date '+%d_%m_%y_%H_%M_%S')"


