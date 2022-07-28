#https://satijalab.org/seurat/articles/integration_mapping.html
#Mapping and annotating query datasets

library(Seurat)
#devtools::install_github('satijalab/seurat-data')
library(SeuratData)

#InstallData("panc8")

#Dataset preprocessing
#To construct a reference, we will identify ‘anchors’ between the individual datasets. First, we split the combined object into a list, with each dataset as an element (this is only necessary because the data was bundled together for easy distribution).
data("panc8")
pancreas.list <- SplitObject(panc8, split.by = "tech")
pancreas.list <- pancreas.list[c("celseq", "celseq2", "fluidigmc1", "smartseq2")]

#Prior to finding anchors, we perform standard preprocessing (log-normalization), and identify variable features individually for each. Note that Seurat implements an improved method for variable feature selection based on a variance stabilizing transformation ("vst")
for (i in 1:length(pancreas.list)) {
  pancreas.list[[i]] <- NormalizeData(pancreas.list[[i]], verbose = FALSE)
  pancreas.list[[i]] <- FindVariableFeatures(pancreas.list[[i]], selection.method = "vst", nfeatures = 2000,
                                             verbose = FALSE)
}
dim(pancreas.list[[1]]@assays$RNA@counts) #34363  1004
dim(pancreas.list[[2]]@assays$RNA@counts) #34363  2285
dim(pancreas.list[[3]]@assays$RNA@counts) #34363   638

#Integration of 3 pancreatic islet cell datasets
reference.list <- pancreas.list[c("celseq", "celseq2", "smartseq2")]
reference.list
pancreas.anchors <- FindIntegrationAnchors(object.list = reference.list, dims = 1:30)
pancreas.integrated <- IntegrateData(anchorset = pancreas.anchors, dims = 1:30)

#The returned object will contain a new Assay, which holds an integrated (or ‘batch-corrected’) expression matrix for all cells, enabling them to be jointly analyzed.

#downstream after integration
library(ggplot2)
library(cowplot)
library(patchwork)
# switch to integrated assay. The variable features of this assay are automatically set during
# IntegrateData
DefaultAssay(pancreas.integrated) <- "integrated"

# Run the standard workflow for visualization and clustering
pancreas.integrated <- ScaleData(pancreas.integrated, verbose = FALSE)
pancreas.integrated <- RunPCA(pancreas.integrated, npcs = 30, verbose = FALSE)
pancreas.integrated <- RunUMAP(pancreas.integrated, reduction = "pca", dims = 1:30, verbose = FALSE)
p1 <- DimPlot(pancreas.integrated, reduction = "umap", group.by = "tech")
p2 <- DimPlot(pancreas.integrated, reduction = "umap", group.by = "celltype", label = TRUE, repel = TRUE) +
  NoLegend()
p1 + p2

#Cell type classification using an integrated reference
pancreas.query <- pancreas.list[["fluidigmc1"]]
pancreas.anchors <- FindTransferAnchors(reference = pancreas.integrated, query = pancreas.query,
                                        dims = 1:30, reference.reduction = "pca")
predictions <- TransferData(anchorset = pancreas.anchors, refdata = pancreas.integrated$celltype,
                            dims = 1:30)
pancreas.query <- AddMetaData(pancreas.query, metadata = predictions)
pancreas.query$predicted.id

#Because we have the original label annotations from our full integrated analysis, we can evaluate how well our predicted cell type annotations match the full reference. In this example, we find that there is a high agreement in cell type classification, with over 96% of cells being labeled correctly.
pancreas.query$prediction.match <- pancreas.query$predicted.id == pancreas.query$celltype
table(pancreas.query$prediction.match)
## 
## FALSE  TRUE 
##    21   617
table(pancreas.query$predicted.id)

#Unimodal UMAP Projection
pancreas.integrated <- RunUMAP(pancreas.integrated, dims = 1:30, reduction = "pca", return.model = TRUE)
pancreas.query <- MapQuery(anchorset = pancreas.anchors, reference = pancreas.integrated, query = pancreas.query,
                           refdata = list(celltype = "celltype"), reference.reduction = "pca", reduction.model = "umap")

p1 <- DimPlot(pancreas.integrated, reduction = "umap", group.by = "celltype", label = TRUE, label.size = 3,
              repel = TRUE) + NoLegend() + ggtitle("Reference annotations")
p2 <- DimPlot(pancreas.query, reduction = "ref.umap", group.by = "predicted.celltype", label = TRUE,
              label.size = 3, repel = TRUE) + NoLegend() + ggtitle("Query transferred labels")
p1 + p2

