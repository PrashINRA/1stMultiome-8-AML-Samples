#!/bin/bash

#This is the shell script to generate demultiplexed fastq files using bcl files with cellranger-arc 2.0

#SBATCH --job-name=cellranger_arc_mkfastq
#SBATCH --output=cellranger_arc_mkfastq_%j.out
#SBATCH --error=cellranger_arc_mkfastq_%j.err
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=20
#SBATCH --mem=200G  
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=p.singh2@amsterdamumc.nl


#Load modules
module load cellranger/arc-2.0.0

export PATH=/trinity/home/psingh/usr/local/bin:$PATH #Require to invoke illumina's bcl2fastq, module load bcl2fastq could also work.

BCL_PATH="/net/beegfs/scratch/psingh/240514_VH00563_167_AACHTM5HV/" #Path to .bcl files
SAMPLESHEET_PATH="/trinity/home/psingh/SampleSheets/atac.csv"  #Prepare the samplesheet 1st
OUTPUT_DIR="/net/beegfs/scratch/psingh/240514_VH00563_167_AACHTM5HV/Data/atac-fastqs"

cellranger-arc mkfastq --run="$BCL_PATH" --csv="$SAMPLESHEET_PATH" 
--output-dir="$OUTPUT_DIR" 
--localcores=40 
--localmem=200
