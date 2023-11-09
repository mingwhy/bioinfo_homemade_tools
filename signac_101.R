#data access from `PBMC scATAC-seq Vignette`, https://satijalab.org/seurat/archive/v3.1/atacseq_integration_vignette
#install signac, https://stuartlab.org/signac/articles/install

# my lab iMac is only compatible with '0.2.5' verion
# replace "0.2.5" with the version that you want to install
#devtools::install_version(package = 'Signac', version = package_version('0.2.5'))

# if want to use the latest version 1.12.0, you can try install via conda
# $ conda install -c bioconda r-signac
# then open Rstudio and library(Signac)

# get started: Analyzing PBMC scATAC-seq
#https://stuartlab.org/signac/articles/pbmc_vignette
#if (!requireNamespace("EnsDb.Hsapiens.v75", quietly = TRUE)) BiocManager::install("EnsDb.Hsapiens.v75")
library(Signac) 
packageVersion('Signac') #0.2.5"
library(Seurat)
library(EnsDb.Hsapiens.v75)
library(ggplot2)
library(patchwork)

#download data: https://satijalab.org/seurat/archive/v3.0/atacseq_integration_vignette.html
# save link as
#scATAC-seq, scATAC-seq metadata, scRNA-seq
#http://cf.10xgenomics.com/samples/cell-atac/1.0.1/atac_v1_pbmc_10k/atac_v1_pbmc_10k_filtered_peak_bc_matrix.h5
#http://cf.10xgenomics.com/samples/cell-atac/1.0.1/atac_v1_pbmc_10k/atac_v1_pbmc_10k_singlecell.csv
#http://cf.10xgenomics.com/samples/cell-exp/3.0.0/pbmc_10k_v3/pbmc_10k_v3_filtered_feature_bc_matrix.h5

#or check: https://github.com/stuart-lab/signac/issues/2

#need to install hdf5r to read HDF5 file
#$ brew install hdf5
#$ conda install -c conda-forge r-hdf5r
#install.packages("hdf5r")
library(hdf5r)
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
