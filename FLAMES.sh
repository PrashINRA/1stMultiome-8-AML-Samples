##This is the script to Run FLAMES pipeline (https://changqingw.github.io/IgniteRNAseq/articles/FLAMESWorkflow.html),
#within the singularity container for Rstudio.

#!/bin/bash
#SBATCH --job-name=flames_sc_long
#SBATCH --output=flames_pipelie__%j.out
#SBATCH --error=flames_pipeline_%j.err
#SBATCH --ntasks=1
#SBATCH --nodelist=node004
#SBATCH --partition=defq
#SBATCH --cpus-per-task=20
#SBATCH --mem=350G
#SBATCH --time=168:00:00
#SBATCH --mail-user=p.singh2@amsterdamumc.nl
#SBATCH --mail-type=END,FAIL

echo "Job started on $(date)"
module load singularity
module load minimap2
module load CUDA/11.8.0

BIND_PATH="/net/beegfs:/net/beegfs,/trinity/shared:/trinity/shared"
CONTAINER_IMAGE="docker://bioconductor/bioconductor_docker:latest"
R_SCRIPT_PATH="/trinity/home/psingh/Sandbox/1stMultiome/sc_long.R"

# Check if minimap2 is available
MINIMAP2_PATH="/trinity/shared/apps/easybuild/software/minimap2/2.20-GCCcore-10.3.0/bin/minimap2"
if [[ ! -x "$MINIMAP2_PATH" ]]; then
    echo "Error: minimap2 not found or not executable at $MINIMAP2_PATH"
    exit 1
else
    echo "Minimap2 found: $MINIMAP2_PATH"
fi

# Check if minimap2 is accessible in the Singularity container
singularity exec --bind "$BIND_PATH" "$CONTAINER_IMAGE" "$MINIMAP2_PATH" --version
if [[ $? -ne 0 ]]; then
    echo "Error: minimap2 not accessible in the Singularity container"
    exit 1
fi

# Check if the R script exists
if [[ ! -f "$R_SCRIPT_PATH" ]]; then
    echo "Error: R script not found at $R_SCRIPT_PATH"
    exit 1
fi

# Run the R script inside the container
echo "Starting FLAMES workflow..."
singularity exec --bind "$BIND_PATH" "$CONTAINER_IMAGE" Rscript "$R_SCRIPT_PATH"

# Capture the exit status of the R script
STATUS=$?
if [[ $STATUS -ne 0 ]]; then
    echo "Error: FLAMES workflow failed with exit status $STATUS"
    exit $STATUS
else
    echo "FLAMES workflow completed successfully"
fi

echo "Job completed on $(date)"

