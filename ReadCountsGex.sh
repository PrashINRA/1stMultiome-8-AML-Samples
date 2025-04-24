#!/bin/bash

#Scripts to count reads per samples using fastq files from GEX data
#SBATCH --job-name=CountReads
#SBATCH --output=CountReads.out
#SBATCH --error=CountReads.err
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=20
#SBATCH --mem=200G   
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=p.singh2@amsterdamumc.nl


# Define the path to your FASTQ folder
fastq_path="/net/beegfs/scratch/psingh/Multiome10x_June2024/GEX/data/FR34156623_10x-GEX-Library/B22KLCKLT3/"

# Define your samples
samples=("6108-DN" "6108-FU2" "6108-FU3" "6279-DN" "6279-FU2" "6905-DN" "6905-FU2" "6905-FU3")

# Initialize an associative array to store read counts
declare -A read_counts

# Iterate over each sample
for sample in "${samples[@]}"; do
  # Initialize count for the current sample
  read_counts[$sample]=0
  
  # Find all FASTQ files for the current sample and count reads
  for file in $(ls ${fastq_path}${sample}*.fastq.gz); do
    # Count the number of lines in the FASTQ file and divide by 4
    reads=$(zcat $file | wc -l)
    reads=$((reads / 4))
    
    # Add to the total count for the current sample
    read_counts[$sample]=$((read_counts[$sample] + reads))
  done
done

# Print the read counts for each sample
for sample in "${!read_counts[@]}"; do
  echo "Sample: $sample, Reads: ${read_counts[$sample]}"
done
