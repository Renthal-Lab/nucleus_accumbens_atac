library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(EnsDb.Mmusculus.v79)
library(ggplot2)
library(patchwork)

counts <- Read10X_h5(".../outs/filtered_peak_bc_matrix.h5")
metadata <- read.csv(
  file = "../outs/singlecell.csv",
  header = TRUE,
  row.names = 1
)

nac_assay <- CreateChromatinAssay(
  counts = counts,
  sep = c(":", "-"),
  genome = "mm10",
  fragments = '../outs/fragments.tsv.gz',
  min.cells = 1
)

## Computing hash

nac.obj <- CreateSeuratObject(
  counts = nac_assay,
  assay = 'peaks',
  project = 'ATAC',
  meta.data = metadata
)



# extract gene annotations from EnsDb
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)

# change to UCSC style since the data was mapped to hg19
seqlevelsStyle(annotations) <- 'UCSC'

# add the gene information to the object
Annotation(nac.obj) <- annotations

nac.obj <- NucleosomeSignal(object = nac.obj)
nac.obj$nucleosome_group <- ifelse(nac.obj$nucleosome_signal > 4, 'NS > 4', 'NS < 4')

nac.obj <- TSSEnrichment(nac.obj, fast = FALSE)
nac.obj$high.tss <- ifelse(nac.obj$TSS.enrichment > 2, 'High', 'Low')
TSSPlot(nac.obj, group.by = 'high.tss') + NoLegend()

nac.obj$pct_reads_in_peaks <- nac.obj$peak_region_fragments / nac.obj$passed_filters * 100
nac.obj$blacklist_ratio <- nac.obj$blacklist_region_fragments / nac.obj$peak_region_fragments

VlnPlot(
  object = nac.obj,
  features = c('pct_reads_in_peaks', 'peak_region_fragments',
               'TSS.enrichment', 'blacklist_ratio', 'nucleosome_signal'),
  pt.size = 0.1,
  ncol = 5
)

nac.obj <- subset(
  x = nac.obj,
  subset = peak_region_fragments > 3000 &
    peak_region_fragments < 100000 &
    pct_reads_in_peaks > 40 &
    blacklist_ratio < 0.025 &
    nucleosome_signal < 4 &
    TSS.enrichment > 2
)

saveRDS(nac.obj, "NAc.QCFiltered.SeuratObj.rds")

nac.obj <- RunTFIDF(nac.obj)
nac.obj <- FindTopFeatures(nac.obj, min.cutoff = 'q0')
nac.obj <- RunSVD(object = nac.obj)

# We exclude the first dimension as this is typically correlated with sequencing depth

nac.obj <- RunUMAP(
  object = nac.obj,
  reduction = 'lsi',
  dims = 2:30
)
nac.obj <- FindNeighbors(
  object = nac.obj,
  reduction = 'lsi',
  dims = 2:30
)
nac.obj <- FindClusters(
  object = nac.obj,
  algorithm = 3,
  resolution = 0.4,
  verbose = FALSE
)

DimPlot(object = nac.obj, label = TRUE) + NoLegend()

# compute gene activities
gene.activities <- GeneActivity(nac.obj)

# add the gene activity matrix to the Seurat object as a new assay
nac.obj[['RNA']] <- CreateAssayObject(counts = gene.activities)
nac.obj <- NormalizeData(						# normalize gene activities
  object = nac.obj,
  assay = 'ACTIVITY',
  normalization.method = 'LogNormalize',
  scale.factor = median(nac.obj$nCount_RNA)
)



nac.rna <- readRDS("../NAc.Cleaned.Clustered.Labeled.scRNA.rds")

transfer.anchors <- FindTransferAnchors(			# Identify anchors
  reference = nac.rna,
  query = nac.obj,
  reduction = 'cca',
  dims = 1:40
)

predicted.labels <- TransferData(
  anchorset = transfer.anchors,
  refdata = nac.rna$Identities,
  weight.reduction = nac.obj[['lsi']],
  dims = 2:30
)

nac.obj <- AddMetaData(object = nac.obj, metadata = predicted.labels)

saveRDS(nac.obj, "NAc.snATAC.Anchored.FullObj.rds")

atac.sub <- subset(nac.obj, subset = prediction.score.max <= 0.6)
table(atac.sub$predicted.id)
atac.sub <- subset(atac.sub, subset = predicted.id %in% c("Endothelial","IN_1","IN_2","IN_3","Mural","NB_cells")) #removing anchored clusters with less than 35

saveRDS(atac.sub, "NAc.Filtered.snATAC.NoINs.rds")

# note that we restrict the imputation to variable genes from scRNA-seq, but could impute the
# full transcriptome if we wanted to
genes.use <- VariableFeatures(nac.rna)
refdata <- GetAssayData(nac.rna, assay = "RNA", slot = "data")[genes.use, ]

# refdata (input) contains a scRNA-seq expression matrix for the scRNA-seq cells.  imputation
# (output) will contain an imputed scRNA-seq matrix for each of the ATAC cells
imputation <- TransferData(anchorset = transfer.anchors, refdata = refdata, weight.reduction = nac.obj[["lsi"]],
    dims = 2:30)
nac.obj[["RNA"]] <- imputation

coembed <- merge(x = nac.rna, y = nac.obj)

# Finally, we run PCA and UMAP on this combined object, to visualize the co-embedding of both
# datasets
coembed <- ScaleData(coembed, features = genes.use, do.scale = FALSE)
coembed <- RunPCA(coembed, features = genes.use, verbose = FALSE)
coembed <- RunUMAP(coembed, dims = 1:30)

saveRDS(coembed, "NAc.Coembedded.Cleaned.FinalObj.rds")