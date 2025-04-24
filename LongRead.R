##This is the sc_long.R scripts required to run FLAMES.sh
# Load FLAMES library
library('FLAMES')
library('SingleCellExperiment')


# Environment variables to enforce UTF-8 encoding
Sys.setenv(PYTHONIOENCODING = "utf-8")
Sys.setenv(LC_ALL = "C.UTF-8")
Sys.setenv(LANG = "C.UTF-8")

basilisk::basiliskRun(env = FLAMES:::flames_env, fun = function(){})

basilisk::basiliskRun(env = FLAMES::: bins_env, fun = function(){})

# # Run the FLAMES pipeline
# sce_6108FU3 <- sc_long_pipeline(
#   annotation = '/trinity/home/psingh/ensembl/Homo_sapiens.GRCh38.113.gtf',
#   fastq = '/net/beegfs/scratch/psingh/LongRead_AMC/merged_fastqs_longRead/P6108/6108_FU3_merged.fastq.gz',
#   outdir = '/net/beegfs/scratch/psingh/LongRead_AMC/merged_fastqs_longRead/P6108_outputFlames/FU3',
#   genome_fa = '/trinity/home/psingh/ensembl/Homo_sapiens.GRCh38.dna.primary_assembly.fa',
#   minimap2 = '/trinity/shared/apps/easybuild/software/minimap2/2.20-GCCcore-10.3.0/bin/minimap2',
#   config_file = '/trinity/home/psingh/ensembl/sclong.json',
#   barcodes_file ='/trinity/home/psingh/OUTS/6108-FU3_output/6108-FU3_count_output/outs/filtered_feature_bc_matrix/6108_FU3_barcodes.csv'
# )

##To genrate compatible barcode file go to the output folder of countmatrix by 10X in terminal then
#  zcat barcodes.tsv.gz | cut -f1 -d'-' > 6108_FU3_barcodes.csv

#Then use realpath 6108_FU3_barcodes.csv as input for barcodes_file above.


# Save the resulting R environment(sce object)
save.image('/trinity/home/psingh/Sandbox/1stMultiome/sce_6108FU3.RData')

gc()

##You have to repeat this for each Sample

