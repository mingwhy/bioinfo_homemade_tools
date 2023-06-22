
library(scMerge)

########################################################
## read in gene chro info
df.gene.table=data.table::fread('../../single.cell_datasets/embryo/gene.meta_embryo.txt',header=T,sep='\t')
head(df.gene.table)
df.gene.table=df.gene.table[df.gene.table$LOCATION_ARM %in% c('2R','3R','3L','4','2L','X'),]
table(df.gene.table$LOCATION_ARM)

df.gene.table=as.data.frame(df.gene.table)
rownames(df.gene.table)=df.gene.table$current_symbol

########################################################
## read in dataset
tissue='embryo'
dat=readRDS('../../single.cell_datasets/embryo/whole_embryo_filtered_valid.rds')

mat <- dat@assays$RNA@counts
gene.names <- rownames(mat)
cell.names <- colnames(mat)
sex='unknown'
cell.type='unknown';

##
dim(mat) #9710 gene x 1297 cells
length(gene.names)
length(cell.names)
rownames(mat)=gene.names
colnames(mat)=cell.names
mat[1:3,1:3]

umi.mat=as.matrix(mat)

## create SingleCellExperiment object
if(!file.exists('embryo_sce.norm.rds')){
  library(scran)
  sce <- SingleCellExperiment(list(counts=umi.mat),
                              colData=DataFrame(cell.type=rep(cell.type,ncol(mat)),sex=sex),
                              rowData=DataFrame(gene=rownames(umi.mat)) )
  sce
  
  clusters <- quickCluster(sce)
  sce <- computeSumFactors(sce, clusters=clusters)
  summary(sizeFactors(sce))
  
  #sce.count <- normalizeCounts(sce, log=FALSE)
  #sce.count <- normalizeCounts(sce,transform="log",pseudo.count=1)
  assayNames(sce)
  logcounts(sce) <- log2(t(t(counts(sce))/sizeFactors(sce)) + 1)
  assayNames(sce)
  saveRDS(sce,file='embryo_sce.obj.norm.rds')
}
sce.count=readRDS('embryo_sce.obj.norm.rds')

## subsetting genes to illustrate usage.
dim(sce.count)
exprs_mat = SummarizedExperiment::assay(sce.count, 'logcounts')
set.seed(1)
result1 = scSEGIndex(exprs_mat = exprs_mat)
x=result1[order(result1$segIdx,decreasing = T),]
head(x)
saveRDS(result1,'embryo_scSEGIndex.rds')

