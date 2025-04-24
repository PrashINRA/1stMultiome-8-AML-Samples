# This is the R-code for generating a merged seurat object from cell-ranger count matrices of multiom data(RNA+ATAC)
# Load pkgs
pkgs <- c( 'Signac','Seurat','dplyr','rstatix','ggpubr','EnsDb.Hsapiens.v86','BSgenome.Hsapiens.UCSC.hg38' )
sapply(pkgs, library, character.only = T)
theme_set(theme_bw())


#Create Seurat object for each sample at once

base_path <- "/trinity/home/psingh/OUTs"
samples <- c("6108-DN", "6108-FU2", "6108-FU3", "6279-DN", "6279-FU2", "6905-DN", "6905-FU2", "6905-FU3")

# Initialize lists to store individual Seurat objects for RNA and ATAC
rna_seurat_objects <- list()
atac_seurat_objects <- list()

annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevels(annotation) <- paste0('chr', seqlevels(annotation))


# Loop through each sample and create Seurat objects
for (sample in samples) {
  # Define paths to filtered data
  data_path <- file.path(base_path, paste0(sample, "_output"), paste0(sample, "_count_output/outs/filtered_feature_bc_matrix.h5"))
  fragpath <- file.path(base_path, paste0(sample, "_output"), paste0(sample, "_count_output/outs/atac_fragments.tsv.gz"))
  
  # Read gene expression and ATAC data
  counts <- Read10X_h5(data_path)
  
  # Create a Seurat object containing the RNA data
  rna_seurat <- CreateSeuratObject(
    counts = counts$`Gene Expression`,
    assay = "RNA",
    project = sample
  )
  
  # Create ATAC assay
  atac_assay <- CreateChromatinAssay(
    counts = counts$Peaks,
    sep = c(":", "-"),
    fragments = fragpath,
    annotation = annotation
  )
  
  # Create a Seurat object for ATAC data
  atac_seurat <- CreateSeuratObject(
    counts = atac_assay,
    assay = "ATAC",
    project = sample
  )
  
  # Store the Seurat objects in their respective lists
  rna_seurat_objects[[sample]] <- rna_seurat
  atac_seurat_objects[[sample]] <- atac_seurat
}

# Merge RNA Seurat objects
combined_rna <- merge(rna_seurat_objects[[1]], y = rna_seurat_objects[-1], add.cell.ids = samples)

# Merge ATAC Seurat objects

library(GenomicRanges)

# Create a unified set of peaks to quantify in each dataset
peak_list <- lapply(atac_seurat_objects, function(x) {
  granges(x@assays$ATAC)
})

grl <- GRangesList(peak_list)

# Reduce to a single set of non-overlapping peaks
combined_peaks <- reduce(unlist(grl))

# Filter out bad peaks based on length
peakwidths <- width(combined_peaks)
combined_peaks <- combined_peaks[peakwidths < 10000 & peakwidths > 20]

# Quantify peaks in each dataset

quantified_matrices <- lapply(names(atac_seurat_objects), function(sample_name) {
  seurat_obj <- atac_seurat_objects[[sample_name]]
  fragments_path <- Fragments(seurat_obj@assays$ATAC)[[1]]@path
  fragments <- CreateFragmentObject(path = fragments_path)
  FeatureMatrix(
    fragments = fragments,
    features = combined_peaks,
    cells = colnames(seurat_obj)
  )
})


# Create Seurat objects with the quantified peaks

for (i in seq_along(quantified_matrices)) {
  sample_name <- names(atac_seurat_objects)[i]
  fragments_path <- Fragments(atac_seurat_objects[[sample_name]]@assays$ATAC)[[1]]@path
  fragments <- CreateFragmentObject(path = fragments_path)
  
  assay <- CreateChromatinAssay(
    counts = quantified_matrices[[i]],
    fragments = fragments
  )
  
  atac_seurat_objects[[sample_name]] <- CreateSeuratObject(
    counts = assay,
    assay = "ATAC"
  )
}

# Add dataset information for merging
for (i in seq_along(atac_seurat_objects)) {
  atac_seurat_objects[[i]]$dataset <- names(atac_seurat_objects)[i]
}

# Merge all atac datasets
combined_atac <- merge(atac_seurat_objects[[1]], y = atac_seurat_objects[-1], add.cell.ids = samples)

# Get gene annotations for hg38
annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevels(annotation) <- paste0('chr', seqlevels(annotation))


# Set the default assay to the combined peaks assay
DefaultAssay(combined_atac) <- "ATAC"

# Combine RNA and ATAC assays into one Seurat object
combined_rna <- JoinLayers(combined_rna, assay = 'RNA')
mdata <- combined_rna$orig.ident
Seu <- CreateSeuratObject(counts = GetAssayData(combined_rna, layer = "counts", assay = "RNA"), min.cell=50, meta.data = mdata)

# Add ATAC assay with fragment files to the combined Seurat object
Seu[["ATAC"]] <- CreateChromatinAssay(counts = GetAssayData(combined_atac, layer = 'counts', assay = 'ATAC'), 
                                      fragments = Fragments(combined_atac))

# Ensure the annotations are set
Annotation(Seu) <- annotation



##Process RNA data for QC

counts <- GetAssayData(Seu, assay = 'RNA', layer = 'counts')
counts[1:10,1:3]
genes_per_cell <- Matrix::colSums(counts>0) # count a gene only if it has non-zero reads mapped.
counts_per_cell <- Matrix::colSums(counts)
plot(sort(genes_per_cell), xlab='cell', log='y', main='genes per cell (ordered)')

MIN_GENES_PER_CELL <- 500
MAX_GENES_PER_CELL <- 2500  

# now replot with the thresholds being shown:
plot(sort(genes_per_cell), xlab='cell', log='y', main='genes per cell (ordered)')
abline(h=MIN_GENES_PER_CELL, col='magenta')  # lower threshold
abline(h=MAX_GENES_PER_CELL, col='gold') # upper threshold

Seu[["percent.mt"]] <- PercentageFeatureSet(Seu, pattern = "^MT-")

seu_rnaqcd <-  subset(Seu, subset = nFeature_RNA >500 & nFeature_RNA <2500 & percent.mt <2 )
DefaultAssay(seu_rnaqcd) <- 'RNA'

counts <- GetAssayData(seu_rnaqcd, assay = 'RNA', layer = 'counts')
mdata <- seu_rnaqcd$orig.ident
rm(seu_rnaqcd)
seu_rnaqcd <- CreateSeuratObject(counts = counts, assay = 'RNA', min.cells = 50)

seu_rnaqcd <-  AddMetaData(seu_rnaqcd, metadata = mdata, col.name = 'Samples')


ALs <- grep('^AL[0-9]',rownames(seu_rnaqcd),value = TRUE)
ACs <- grep('^AC[0-9]',rownames(seu_rnaqcd),value = TRUE)
ADs <- grep('^AD[0-9]',rownames(seu_rnaqcd),value = TRUE)
MTs <- grep('^MT',rownames(seu_rnaqcd),value = TRUE)

APs <- grep('^AP[0-9]',rownames(seu_rnaqcd),value = TRUE)
EFs <-  grep('^EEF',rownames(seu_rnaqcd),value = TRUE)
LINCs <-  grep('^LINC',rownames(seu_rnaqcd),value = TRUE)


gtr <- c('MALAT1','HBB','IGKC', ALs, ACs, ADs,MTs,APs,EFs,LINCs)

idx <- which(rownames(seu_rnaqcd)%in%gtr)
seu_rnaqcd <- seu_rnaqcd[-idx,]
dim(seu_rnaqcd)

#compare 3 different normalization only on RNA data
TCs <- seu_rnaqcd
library(patchwork)

# Perform SCTransform normalization
TCs_sct <- SCTransform(TCs, verbose = TRUE)
normalized_depth_sct <- colSums(GetAssayData(TCs_sct, slot = "data", assay = "SCT"))
TCs_sct <- AddMetaData(TCs_sct, metadata = normalized_depth_sct, col.name = "normalized_depth_sct")

# Perform LogNormalize normalization
TCs_log <- NormalizeData(TCs, normalization.method = "LogNormalize", scale.factor = median(TCs$nCount_RNA))
normalized_depth_log <- colSums(GetAssayData(TCs_log, slot = "data", assay = "RNA"))
TCs_log <- AddMetaData(TCs_log, metadata = normalized_depth_log, col.name = "normalized_depth_log")

# Perform CLR normalization
TCs_clr <- NormalizeData(TCs, normalization.method = "CLR", scale.factor = median(TCs$nCount_RNA))
normalized_depth_clr <- colSums(GetAssayData(TCs_clr, slot = "data", assay = "RNA"))
TCs_clr <- AddMetaData(TCs_clr, metadata = normalized_depth_clr, col.name = "normalized_depth_clr")

# Prepare data for plotting
plot_data_sct <- data.frame(Samples = TCs_sct@meta.data$orig.ident, NormalizedDepth = TCs_sct@meta.data$normalized_depth_sct)
plot_data_log <- data.frame(Samples = TCs_log@meta.data$orig.ident, NormalizedDepth = TCs_log@meta.data$normalized_depth_log)
plot_data_clr <- data.frame(Samples = TCs_clr@meta.data$orig.ident, NormalizedDepth = TCs_clr@meta.data$normalized_depth_clr)

# Generate individual boxplots
plot_sct <- ggplot(plot_data_sct, aes(x = Samples, y = NormalizedDepth)) +
  geom_boxplot() +
  stat_summary(fun = median, geom = "text", aes(label = format(..y.., digits = 2, nsmall = 2)), vjust = -0.5, color = 'cyan3', size = 3.5) +
  theme_minimal() +
  labs(title = "Normalized Sequencing Depth: SCTransform", x = "", y = "") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

plot_log <- ggplot(plot_data_log, aes(x = Samples, y = NormalizedDepth)) +
  geom_boxplot() +
  stat_summary(fun = median, geom = "text", aes(label = format(..y.., digits = 2, nsmall = 2)), vjust = -0.5, color = 'cyan3', size = 3.5) +
  theme_minimal() +
  labs(title = "Normalized Sequencing Depth: LogNormalize", x = "", y = "Normalized Sequencing Depth (UMI Count per Cell)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

plot_clr <- ggplot(plot_data_clr, aes(x = Samples, y = NormalizedDepth)) +
  geom_boxplot() +
  stat_summary(fun = median, geom = "text", aes(label = format(..y.., digits = 2, nsmall = 2)), vjust = -0.5, color = 'cyan3', size = 3.5) +
  theme_minimal() +
  labs(title = "Normalized Sequencing Depth: CLR", x = "", y = "") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Combine plots using patchwork
combined_plot <- plot_sct + plot_log + plot_clr + plot_layout(ncol = 1)

# Display the combined plot
print(combined_plot)

#Go with SCT

##Process ATAC data for QC-

DefaultAssay(Seu) <- 'ATAC'
# Now run NucleosomeSignal and TSSEnrichment
Seu <- NucleosomeSignal(Seu)
Seu <- TSSEnrichment(Seu)

a1 <- DensityScatter(Seu, x = 'nCount_ATAC', y = 'TSS.enrichment', log_x = TRUE, quantiles = TRUE)
a2 <- DensityScatter(Seu, x = 'nucleosome_signal', y = 'TSS.enrichment', log_x = TRUE, quantiles = TRUE)

a1 | a2

dittoSeq::dittoPlot(Seu, c('nCount_ATAC', 'nFeature_ATAC'), group.by = "Samples", plots = c("jitter", "vlnplot"))+xlab('')

dittoSeq::dittoPlot(Seu, c('TSS.enrichment', 'nucleosome_signal'), group.by = "Samples", plots = c("jitter", "vlnplot"))+xlab('')

# ....Filtering poor quality cells --------------------

counts <- GetAssayData(Seu, assay = 'ATAC', layer = 'counts')
peaks_per_cell <- Matrix::colSums(counts>0)
peaks_per_cell <- Matrix::colSums(counts)
plot(sort(peaks_per_cell), xlab='cell', log='y', main='sorted(peaks_per_cell)')

min_peaks_perCell <- 1000
max_peaks_perCell <- 15000  

# now replot with the thresholds being shown:
plot(sort(peaks_per_cell), xlab='cell', log='y', main='peaks per cell (ordered)')
abline(h=min_peaks_perCell, col='purple')  # lower threshold
abline(h=max_peaks_perCell, col='gold') # upper threshold


# ....Filtering poor quality cells ------------------
Seu$Samples <- Seu$orig.ident

samp_cols <- c("grey", "#f56505", "#dec400", "#006630", "#0223c7","#5b02c7","#00b0e6", "#c40080")
VlnPlot(
  Seu,
  features = c( "nCount_RNA","nCount_ATAC", "TSS.enrichment", "nucleosome_signal"),
  ncol = 4,group.by = 'Samples',cols = samp_cols,
  pt.size = 0
)

Seu <- subset(Seu,
              subset = nCount_ATAC >1000 &
                nCount_ATAC < 15000 &
                nucleosome_signal < 4 &
                TSS.enrichment > 3)
dim(pbmc)


rm(list=setdiff(ls(), c('Seu', 'samp_cols')))
gc()
