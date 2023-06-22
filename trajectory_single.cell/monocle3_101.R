
#https://cole-trapnell-lab.github.io/monocle3/docs/installation/

library(BiocGenerics)
library(DelayedArray)
library(DelayedMatrixStats)
library(limma)
library(S4Vectors)
library(SingleCellExperiment)
library(SummarizedExperiment)
library(batchelor)
library(Matrix.utils)
#devtools::install_github('cole-trapnell-lab/leidenbase')
library(leidenbase)
#devtools::install_github('cole-trapnell-lab/monocle3')
library(monocle3)

# The tutorial shown below and on subsequent pages uses two additional packages:
library(ggplot2)
library(dplyr)

#https://cole-trapnell-lab.github.io/monocle3/docs/trajectories/
# Load the data
expression_matrix <- readRDS(url("http://staff.washington.edu/hpliner/data/packer_embryo_expression.rds"))
cell_metadata <- readRDS(url("http://staff.washington.edu/hpliner/data/packer_embryo_colData.rds"))
gene_annotation <- readRDS(url("http://staff.washington.edu/hpliner/data/packer_embryo_rowData.rds"))

class(expression_matrix) #sparse matrix
class(cell_metadata) #data.frame
class(gene_annotation) #data.frame

expression_matrix[1:3,1:3] #umi count matrix

cds <- new_cell_data_set(expression_matrix,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)

# pre-process the data
cds <- preprocess_cds(cds, num_dim = 50)
cds <- align_cds(cds, alignment_group = "batch", residual_model_formula_str = "~ bg.300.loading + bg.400.loading + bg.500.1.loading + bg.500.2.loading + bg.r17.loading + bg.b01.loading + bg.b02.loading")
# Reduce dimensionality and visualize the results
cds <- reduce_dimension(cds)
plot_cells(cds, label_groups_by_cluster=FALSE,  color_cells_by = "cell.type")
# cluster cells
cds <- cluster_cells(cds)
plot_cells(cds, color_cells_by = "partition")

# Learn the trajectory graph
cds <- learn_graph(cds)
plot_cells(cds,
           color_cells_by = "cell.type",
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE)
# Order the cells in pseudotime
plot_cells(cds,
           color_cells_by = "embryo.time.bin",
           label_cell_groups=FALSE,
           label_leaves=TRUE,
           label_branch_points=TRUE,
           graph_label_size=1.5)
# launch a graphical user interface for selecting one or more root nodes.
cds <- order_cells(cds)
plot_cells(cds,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5)
#Note that some of the cells are gray. This means they have infinite pseudotime, because they were not reachable from the root nodes that were picked. 


## specify the root programmatically
# a helper function to identify the root principal points:
get_earliest_principal_node <- function(cds, time_bin="130-170"){
  cell_ids <- which(colData(cds)[, "embryo.time.bin"] == time_bin)
  
  closest_vertex <-
    cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  root_pr_nodes <-
    igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names
                                                              (which.max(table(closest_vertex[cell_ids,]))))]
  
  root_pr_nodes
}
cds <- order_cells(cds, root_pr_nodes=get_earliest_principal_node(cds))
#Passing the programatically selected root node to order_cells() via the root_pr_nodeargument yields:
plot_cells(cds,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5)


######################################################
expression_matrix <- readRDS(url("http://staff.washington.edu/hpliner/data/cao_l2_expression.rds"))
cell_metadata <- readRDS(url("http://staff.washington.edu/hpliner/data/cao_l2_colData.rds"))
gene_annotation <- readRDS(url("http://staff.washington.edu/hpliner/data/cao_l2_rowData.rds"))
dim(expression_matrix)# 20271 42035

# Make the CDS object
cds <- new_cell_data_set(expression_matrix,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)
# preprocess, pick PC1-100 by PCA
cds <- preprocess_cds(cds, num_dim = 100)
plot_pc_variance_explained(cds)

# reduce dimension
cds <- reduce_dimension(cds)
plot_cells(cds)
plot_cells(cds, color_cells_by="cao_cell_type")
plot_cells(cds, genes=c("cpna-2", "egl-21", "ram-2", "inos-1"))

# group cells into cluster
cds = cluster_cells(cds, resolution=1e-5)
plot_cells(cds)
plot_cells(cds, color_cells_by="partition", group_cells_by="partition")

## https://cole-trapnell-lab.github.io/monocle3/docs/trajectories/
#Constructing single-cell trajectories
cds <- preprocess_cds(cds, num_dim = 50)
cds <- align_cds(cds, alignment_group = "batch", residual_model_formula_str = "~ bg.300.loading + bg.400.loading + bg.500.1.loading + bg.500.2.loading + bg.r17.loading + bg.b01.loading + bg.b02.loading")

cds <- reduce_dimension(cds)
plot_cells(cds, label_groups_by_cluster=FALSE,  color_cells_by = "cell.type")

