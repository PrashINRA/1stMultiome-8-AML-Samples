#!/bin/bash

#This is the shell script to generate a count matrix for RNA and ATAC using cellranger-arc 2.0

#SBATCH --job-name=cellranger-arc-count
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --time=72:00:00
#SBATCH --cpus-per-task=40
#SBATCH --mem=300G
#SBATCH --partition=defq
# Email notifications
#SBATCH --mail-type=END,FAIL         
#SBATCH --mail-user=p.singh2@amsterdamumc.nl


# Load the cellranger-arc module
module load cellranger/arc-2.0.0

# Define the base paths
gex_path="/net/beegfs/scratch/psingh/Multiome10x_June2024/GEX/data/FR34156623_10x-GEX-Library/B22KLCKLT3/"
atac_path="/net/beegfs/scratch/psingh/Multiome10x_June2024/ATAC/Data/atac-fastqs/"
output_dir_base="/trinity/home/psingh/OUTS"
csv_dir="/trinity/home/psingh/SampleSheets"

# Define samples
samples=("6108-DN" "6108-FU2" "6108-FU3" "6279-DN" "6279-FU2" "6905-DN" "6905-FU2" "6905-FU3")
atac_samples=("6108_DN" "6108_FU2" "6108_FU3" "6279_DN" "6279_FU2" "6905_DN" "6905_FU2" "6905_FU3")

# Create libraries.csv for each sample and run cellranger-arc count
for i in "${!samples[@]}"; do
  sample_gex=${samples[$i]}
  sample_atac=${atac_samples[$i]}

  csv_file="$csv_dir/libraries_${sample_gex}.csv"
  output_dir="$output_dir_base/${sample_gex}_output"

  # Write the CSV content
  echo "fastqs,sample,library_type" > $csv_file
  echo "$gex_path,$sample_gex,Gene Expression" >> $csv_file
  echo "$atac_path,$sample_atac,Chromatin Accessibility" >> $csv_file

  # Ensure output directory exists
  mkdir -p $output_dir

  # Change to output directory
  cd $output_dir

  # Run cellranger-arc count for each sample
  cellranger-arc count --id="${sample_gex}_count_output" \
                       --libraries=$csv_file \
                       --reference="/trinity/home/psingh/Refs/refdata-cellranger-arc-GRCh38-2020-A-2.0.0" \
                       --no-bam \
                       --localmem=300 \
                       --localcores=40
done
