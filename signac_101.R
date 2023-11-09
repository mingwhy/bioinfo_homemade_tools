#data access from `PBMC scATAC-seq Vignette`, https://satijalab.org/seurat/archive/v3.1/atacseq_integration_vignette

# my lab iMac is only compatible with '0.2.5' verion
# but 0.2.5 version didn't contain several functions
# I solved this problem via `Tutorial:Using R in Conda`: https://www.biostars.org/p/450316/

# if you want to create a new environment, Use 
$ conda create -n R_signac #to create new environment.
$ conda install -c r r-essentials
$ conda install -c r rstudio.
#install signac, https://stuartlab.org/signac/articles/install
$ conda install -c bioconda r-signac

$ conda install -c conda-forge r-biocmanager   
$ conda install -c bioconda r-seurat
# if (!requireNamespace("EnsDb.Hsapiens.v75", quietly = TRUE)) BiocManager::install("EnsDb.Hsapiens.v75")
$ conda install -c bioconda bioconductor-ensdb.hsapiens.v75 
$ rstudio

setwd("~/Downloads")
library(Signac) 
packageVersion('Signac') #‘1.10.0’
library(ggplot2)
library(patchwork)
library(Seurat)
library(EnsDb.Hsapiens.v75)


#download data: https://satijalab.org/seurat/archive/v3.0/atacseq_integration_vignette.html
# save link as
#scATAC-seq, scATAC-seq metadata, scRNA-seq
#http://cf.10xgenomics.com/samples/cell-atac/1.0.1/atac_v1_pbmc_10k/atac_v1_pbmc_10k_filtered_peak_bc_matrix.h5
#http://cf.10xgenomics.com/samples/cell-atac/1.0.1/atac_v1_pbmc_10k/atac_v1_pbmc_10k_singlecell.csv
#http://cf.10xgenomics.com/samples/cell-exp/3.0.0/pbmc_10k_v3/pbmc_10k_v3_filtered_feature_bc_matrix.h5

#also check: https://github.com/stuart-lab/signac/issues/2
#http://cf.10xgenomics.com/samples/cell-atac/1.1.0/atac_v1_pbmc_10k/atac_v1_pbmc_10k_fragments.tsv.gz
#http://cf.10xgenomics.com/samples/cell-atac/1.1.0/atac_v1_pbmc_10k/atac_v1_pbmc_10k_fragments.tsv.gz.tbi

#need to install hdf5r to read HDF5 file
#$ brew install hdf5
$ conda install -c conda-forge r-hdf5r

counts <- Read10X_h5(filename = "atac_v1_pbmc_10k_filtered_peak_bc_matrix.h5")
metadata <- read.csv(
  file = "atac_v1_pbmc_10k_singlecell.csv",
  header = TRUE,
  row.names = 1
)

chrom_assay <- CreateChromatinAssay(
  counts = counts,
  sep = c(":", "-"),
  fragments = 'atac_v1_pbmc_10k_fragments.tsv.gz',
  min.cells = 10,
  min.features = 200
)

pbmc <- CreateSeuratObject(
  counts = chrom_assay,
  assay = "peaks",
  meta.data = metadata
)

pbmc
## An object of class Seurat 
## 87561 features across 8728 samples within 1 assay 
## Active assay: peaks (87561 features, 0 variable features)
##  2 layers present: counts, data


#########
# extract gene annotations from EnsDb
#need to install biovizBase
$ conda install -c bioconda bioconductor-biovizbase

annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v75)

# change to UCSC style since the data was mapped to hg19
seqlevels(annotations) <- paste0('chr', seqlevels(annotations))
genome(annotations) <- "hg19"

# add the gene information to the object
Annotation(pbmc) <- annotations

#########
#Computing QC Metrics
# compute nucleosome signal score per cell
pbmc <- NucleosomeSignal(object = pbmc)

# compute TSS enrichment score per cell
pbmc <- TSSEnrichment(object = pbmc, fast = FALSE)

# add blacklist ratio and fraction of reads in peaks
pbmc$pct_reads_in_peaks <- pbmc$peak_region_fragments / pbmc$passed_filters * 100
pbmc$blacklist_ratio <- pbmc$blacklist_region_fragments / pbmc$peak_region_fragments

DensityScatter(pbmc, x = 'nCount_peaks', y = 'TSS.enrichment', log_x = TRUE, quantiles = TRUE)

pbmc$high.tss <- ifelse(pbmc$TSS.enrichment > 3, 'High', 'Low')
#TSSPlot(pbmc, group.by = 'high.tss') + NoLegend()

pbmc$nucleosome_group <- ifelse(pbmc$nucleosome_signal > 4, 'NS > 4', 'NS < 4')
FragmentHistogram(object = pbmc, group.by = 'nucleosome_group')

#Finally we remove cells that are outliers for these QC metrics.
pbmc <- subset(
  x = pbmc,
  subset = nCount_peaks > 3000 &
    nCount_peaks < 30000 &
    pct_reads_in_peaks > 15 &
    blacklist_ratio < 0.05 &
    nucleosome_signal < 4 &
    TSS.enrichment > 3
)
pbmc
#saveRDS(pbmc,'pbmc_small.rds')

pbmc=readRDS('pbmc_small.rds')
################
#Normalization and linear dimensional reduction
pbmc <- RunTFIDF(pbmc)
pbmc <- FindTopFeatures(pbmc, min.cutoff = 'q0')

#need to install irlba (check if it's installed or not before run conda)
#$ conda install -c conda-forge r-irlba
#when running $pbmc <- RunSVD(pbmc) #https://stuartlab.org/signac/reference/runsvd
#error with Matrix package: https://github.com/bwlewis/irlba/issues/70
#github issue: Changes in Matrix 1.6-2 breaking irlba: function 'as_cholmod_sparse' not provided by package 'Matrix'
# I reinstall r-matrix to 1.6_1 verion via conda
#https://stackoverflow.com/questions/38411942/anaconda-conda-install-a-specific-package-version
$ conda install 'r-matrix=1.6_1'

pbmc <- RunSVD(pbmc)
DepthCor(pbmc)

#Non-linear dimension reduction and clustering
pbmc <- RunUMAP(object = pbmc, reduction = 'lsi', dims = 2:30)
pbmc <- FindNeighbors(object = pbmc, reduction = 'lsi', dims = 2:30)
pbmc <- FindClusters(object = pbmc, verbose = FALSE, algorithm = 3)
DimPlot(object = pbmc, label = TRUE) + NoLegend()


################
#Create a gene activity matrix
gene.activities <- GeneActivity(pbmc)
# add the gene activity matrix to the Seurat object as a new assay and normalize it
pbmc[['RNA']] <- CreateAssayObject(counts = gene.activities)
pbmc <- NormalizeData(
  object = pbmc,
  assay = 'RNA',
  normalization.method = 'LogNormalize',
  scale.factor = median(pbmc$nCount_RNA)
)
DefaultAssay(pbmc) <- 'RNA'

FeaturePlot(
  object = pbmc,
  features = c('MS4A1', 'CD3D', 'LEF1', 'NKG7', 'TREM1', 'LYZ'),
  pt.size = 0.1,
  max.cutoff = 'q95',
  ncol = 3
)

################
# can download the pre-processed Seurat object here.
# https://signac-objects.s3.amazonaws.com/pbmc_10k_v3.rds
# Integrating with scRNA-seq data
# Load the pre-processed scRNA-seq data for PBMCs
pbmc_rna <- readRDS("pbmc_10k_v3.rds")
pbmc_rna <- UpdateSeuratObject(pbmc_rna)

transfer.anchors <- FindTransferAnchors(
  reference = pbmc_rna,
  query = pbmc,
  reduction = 'cca'
)

predicted.labels <- TransferData(
  anchorset = transfer.anchors,
  refdata = pbmc_rna$celltype,
  weight.reduction = pbmc[['lsi']],
  dims = 2:30
)

pbmc <- AddMetaData(object = pbmc, metadata = predicted.labels)

plot1 <- DimPlot(
  object = pbmc_rna,
  group.by = 'celltype',
  label = TRUE,
  repel = TRUE) + NoLegend() + ggtitle('scRNA-seq')

plot2 <- DimPlot(
  object = pbmc,
  group.by = 'predicted.id',
  label = TRUE,
  repel = TRUE) + NoLegend() + ggtitle('scATAC-seq')

plot1 + plot2

# replace each label with its most likely prediction
for(i in levels(pbmc)) {
  cells_to_reid <- WhichCells(pbmc, idents = i)
  newid <- names(which.max(table(pbmc$predicted.id[cells_to_reid])))
  Idents(pbmc, cells = cells_to_reid) <- newid
}

#################
#Find differentially accessible peaks between cell types
# change back to working with peaks instead of gene activities
DefaultAssay(pbmc) <- 'peaks'

da_peaks <- FindMarkers(
  object = pbmc,
  ident.1 = "CD4 Naive",
  ident.2 = "CD14+ Monocytes",
  test.use = 'LR',
  latent.vars = 'nCount_peaks'
)
head(da_peaks)

plot1 <- VlnPlot(
  object = pbmc,
  features = rownames(da_peaks)[1],
  pt.size = 0.1,
  idents = c("CD4 Naive","CD14+ Monocytes")
)
plot2 <- FeaturePlot(
  object = pbmc,
  features = rownames(da_peaks)[1],
  pt.size = 0.1
)

plot1 | plot2

fc <- FoldChange(pbmc, ident.1 = "CD4 Naive", ident.2 = "CD14+ Monocytes")
# order by fold change
fc <- fc[order(fc$avg_log2FC, decreasing = TRUE), ]
head(fc)

open_cd4naive <- rownames(da_peaks[da_peaks$avg_log2FC > 3, ])
open_cd14mono <- rownames(da_peaks[da_peaks$avg_log2FC < -3, ])

closest_genes_cd4naive <- ClosestFeature(pbmc, regions = open_cd4naive)
closest_genes_cd14mono <- ClosestFeature(pbmc, regions = open_cd14mono)
head(closest_genes_cd4naive)
head(closest_genes_cd14mono)
######################
#Plotting genomic regions
# set plotting order
levels(pbmc) <- c("CD4 Naive","CD4 Memory","CD8 Naive","CD8 effector","Double negative T cell","NK dim", "NK bright", "pre-B cell",'B cell progenitor',"pDC","CD14+ Monocytes",'CD16+ Monocytes')

CoveragePlot(
  object = pbmc,
  region = rownames(da_peaks)[1],
  extend.upstream = 40000,
  extend.downstream = 20000
)
