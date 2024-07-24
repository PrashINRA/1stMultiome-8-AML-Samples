#!/bin/bash
#SBATCH --job-name=cellranger_arc_mkfastq
#SBATCH --output=cellranger_arc_mkfastq_%j.out
#SBATCH --error=cellranger_arc_mkfastq_%j.err
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=20
#SBATCH --mem=200G  # Increased memory allocation 
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=prashnitrr@gmail.com

# Check the exit status of the last command
if [ $? -ne 0 ]; then
    echo "cellranger-arc mkfastq failed with an error."
else
    echo "cellranger-arc mkfastq completed successfully."
fi

module load cellranger/arc-2.0.0

export PATH=/trinity/home/psingh/usr/local/bin:$PATH

BCL_PATH="/net/beegfs/scratch/psingh/240514_VH00563_167_AACHTM5HV/"
SAMPLESHEET_PATH="/trinity/home/psingh/SampleSheets/atac.csv"
OUTPUT_DIR="/net/beegfs/scratch/psingh/240514_VH00563_167_AACHTM5HV/Data/atac-fastqs"

cellranger-arc mkfastq --run="$BCL_PATH" --csv="$SAMPLESHEET_PATH" --output-dir="$OUTPUT_DIR" --localcores=20 --localmem=200

