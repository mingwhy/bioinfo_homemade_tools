
setwd("~/Documents/sc_transcriptome.index/pseudotime_analysis/")

#BiocManager::install(c('BiocGenerics', 'DelayedArray', 'DelayedMatrixStats',
#                       'limma', 'S4Vectors', 'SingleCellExperiment',
#                       'SummarizedExperiment', 'batchelor', 'Matrix.utils'))

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

library(ggplot2)
library(dplyr)

#install.packages('Seurat')
library(Seurat)

########################################################
## read in dataset
#BiocManager::install("zellkonverter")
library(zellkonverter)
library(SummarizedExperiment)

## read in data
sce<-readH5AD('~/Documents/Data_mouse_aging_atlas//TMS.gene.data_final/tabula-muris-senis-facs-official-raw-obj.h5ad');       
class(assay(sce))
sce


dim(assay(sce))
df.expr=assay(sce)
df.expr[1:10,1:10] #umi count
dim(df.expr) #22966 110824

## sample meta information
cell.meta=colData(sce)
head(cell.meta)
dim(cell.meta) #110824     13
colnames(cell.meta)


table(cell.meta$age) #6 age group: 1, 3, 18, 21, 24,30month. fac: 3,18,21,24 month
table(cell.meta$tissue) #23 tissue
length(table(cell.meta$cell_ontology_class)) #120 cell types
unique(paste(cell.meta$tissue,cell.meta$cell_ontology_class)) #207 unique #consistent with `elife-62293-supp1-v2.xlsx`
x=as.data.frame(cell.meta) %>% group_by(tissue) %>% summarise(n.cell.type=length(unique(cell_ontology_class)))
x$n.cell.type

tissue_cell.type=paste(cell.meta$tissue,cell.meta$cell_ontology_class,sep=':')
cell.meta$tissue_cell.type=tissue_cell.type

table(cell.meta[cell.meta$age=='21m',]$tissue_cell.type) #remove these cells
i=cell.meta$age!='21m'
cell.meta=cell.meta[i,]
df.expr=df.expr[,i]
dim(cell.meta) # 110096     14
dim(df.expr) # 22966 110096


#######################################################################
## https://cole-trapnell-lab.github.io/monocle3/docs/starting/
tc=sort(unique(cell.meta$tissue_cell.type))[7]
tc

expression_matrix=df.expr[,cell.meta$tissue_cell.type==tc]
class(expression_matrix)
dim(expression_matrix)

#dat=subset(dat,annotation=="spermatocyte 3")
gene.names <- rownames(expression_matrix)
cell.names <- colnames(expression_matrix)

cell_metadata<-cell.meta[cell.meta$tissue_cell.type==tc,]
rownames(cell_metadata)=cell.names
gene_annotation<-data.frame('gene_short_name'=gene.names)
rownames(gene_annotation)=gene.names

cds <- new_cell_data_set(expression_matrix,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)

head(cds@colData)

# pre-process the data
cds <- preprocess_cds(cds, num_dim = 100)

# Reduce dimensionality and visualize the results
cds <- reduce_dimension(cds,reduction_method ='UMAP')

plot_cells(cds, label_groups_by_cluster=FALSE)

# cluster cells, each cell is assigned to a cluster and a partition
cds <- cluster_cells(cds,cluster_method='leiden')
cds@clusters
x=(cds@clusters$UMAP)
names(x) #"cluster_result" "partitions"     "clusters"      

plot_cells(cds, color_cells_by = "partition")

# Learn the trajectory graph
cds <- learn_graph(cds)
plot_cells(cds,
           #color_cells_by = "partition",
           color_cells_by = "age",
           show_trajectory_graph = T,
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE)

# check one gene
grep('Rpl13a',gene.names)
aging_genes <- c("Rpl13a")

plot_cells(cds,
           genes=aging_genes,
           label_cell_groups=FALSE,
           show_trajectory_graph=FALSE)


# launch a graphical user interface for selecting one or more root nodes.

cds <- order_cells(cds)
x=pseudotime(cds)
length(x) #number of cells
sum(is.infinite(x))

plot_cells(cds,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5)
#Note that some of the cells are gray. This means they have infinite pseudotime, 
#because they were not reachable from the root nodes that were picked. 

#In general, any cell on a parition that lacks a root node will be assigned an infinite pseudotime. 
#In general, you should choose at least one root per partition.
x=pseudotime(cds)
length(x)


class(cds@colData) #"DFrame"
#cds@colData$test=1
cds@colData$test=pseudotime(cds)

plot_cells(cds,
           color_cells_by = "test",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5)

