#https://github.com/neurorestore/Augur
devtools::install_github("Bioconductor/MatrixGenerics")
devtools::install_github("const-ae/sparseMatrixStats")
devtools::install_github("neurorestore/Augur")
BiocManager::install('metap') #for `FindConservedMarkers`

#https://www.nature.com/articles/s41596-021-00561-x
library(tidyverse)
library(Seurat)
library(Augur)

# input dataset from Zenodo: https://zenodo.org/record/4473025#.YjThF5rMLm8

input_dir = "augur_protocol_zenodo/rnaseq/raw"
mat = readRDS(file.path(input_dir, "Kang2018_mat.rds"))
meta = readRDS(file.path(input_dir, "Kang2018_meta.rds"))
# create the Seurat object
sc = CreateSeuratObject(mat, min.cells = 3, min.features = 0,
                        meta.data = meta)
# confirm dimensions of the object
dim(sc)
# [1] 15706 24673
# print the experimental conditions
unique(sc$label)
# [1] "ctrl" "stim"


sc = sc %>%
  # Split the object into a list for input the Seurat integration
  SplitObject(split.by = 'label') %>%
  # Normalize the data using regularized negative binomial models
  map(~ SCTransform(.)) %>%
  # Use Seurat to find anchors across the conditions (baseline/stim)
  PrepSCTIntegration(
    anchor.features = SelectIntegrationFeatures(.)) %>%
  FindIntegrationAnchors(
    anchor.features = SelectIntegrationFeatures(.),
    normalization.method = 'SCT') %>%
  # Integrate data
  IntegrateData(normalization.method = 'SCT')

sc = sc %>%
  # Run principal component analysis
  RunPCA(npcs = 30, verbose = F) %>%
  # Embed in two dimensions
  RunUMAP(dims = 1:20, do.fast = T)

sc = sc %>%
  FindNeighbors(dims = 1:10) %>%
  FindClusters(resolution = 0.4)

# Find marker genes for cluster 7
markers = FindConservedMarkers(sc, ident.1 = 7, grouping.var = "label")
# Print the top 10 markers for cluster 7
rownames(markers)[1:10]
# [1] "FCGR3A" "FAM26F" "VMO1" "GBP5" "TNFSF10" "C3AR1" "ATP1B3" "MS4A7" "CFD" "SERPINA1"
# Repeat these lines for each cluster in turn
# Alternatively, inspect some known marker genes
FeaturePlot(sc, features = c("FCGR3A", "LYZ", "MS4A1", "NKG7"))
## See Fig. 3c

# Annotate the clusters
sc = RenameIdents(sc,
                  `0` = "CD4 T cells",
                  `1` = "CD14+ Monocytes",
                  `2` = "CD4 T cells",
                  `3` = "B cells",
                  `4` = "NK cells",
                  `5` = "CD8 T cells",
                  `6` = "CD4 T cells",
                  `7` = "FCGR3A+ Monocytes",
                  `8` = "Dendritic cells",
                  `9` = "B cells",
                  `10` = "Megakaryocytes",
                  `11` = "CD4 T cells",
                  `12` = "Dendritic cells")
# Add cell type annotations into the metadata of the Seurat object
sc$cell_type = Idents(sc)


DimPlot(sc, label = TRUE)
## See Fig. 3a

DefaultAssay(sc) = "RNA"

# Cell-type prioritization
#augur = calculate_auc(sc)

#Alternatively, if the input data are not in a Seurat, monocle or SingleCellExperiment object, the user can pass in the expression matrix and the accompanying metadata data frame separately, as described above in ‘Overview of the procedure’. For instance, these can be extracted from the Seurat object and provided as input into Augur directly, using the following commands:
expr = GetAssayData(sc)
meta = sc@meta.data
dim(expr); #15706 24673
dim(meta); #24673    14
colnames(meta)
table(meta$label) 
#ctrl  stim 
#12315 12358 
#https://rdrr.io/github/neurorestore/Augur/man/calculate_auc.html
Sys.time()
system.time(
  augur <- calculate_auc(input = expr, meta = meta,
                         n_threads = 4,
                         cell_type_col='cell_type',label_col='label')
)
Sys.time()
#user  system elapsed 
#881.340  10.679 259.279 
#elapsed 259.279 seconds

