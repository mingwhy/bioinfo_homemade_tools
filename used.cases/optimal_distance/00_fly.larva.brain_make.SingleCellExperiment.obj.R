#https://bioconductor.org/packages/devel/bioc/vignettes/SingleCellExperiment/inst/doc/intro.html
#http://bioconductor.org/packages/release/bioc/vignettes/scuttle/inst/doc/qc.html
library(Matrix);library(Seurat)
library(ggplot2);library(gridExtra);
library(patchwork);library(dplyr)
library(SingleCellExperiment)
library(SummarizedExperiment)
library(scuttle);
library(tidyverse);
library(scater);library(scran)
#options(future.globals.maxSize = 3.2 * 1024^3)

###############
## raw expr.mat
(mtx.files=Sys.glob('./GSE134722/dataset/*final*mtx.gz'))
(features.files=Sys.glob('./GSE134722/dataset/*final*genes.tsv.gz'))
(cells.files=Sys.glob('./GSE134722/dataset/*final*barcodes.tsv.gz'))

sample.names=c('Normal','Starvation')

two.samples<-lapply(1:2,function(i){
  expr.mtx=readMM(mtx.files[i])
  dim(expr.mtx) #14393 gene x 6k~25k cell
  
  barcode=read.table(cells.files[i],as.is = T)
  dim(barcode) 
  
  geneId=read.table(features.files[i],as.is = T)
  head(geneId)
  colnames(geneId)=c("flybase",'symbol')
  
  rownames(expr.mtx)<-geneId$symbol
  colnames(expr.mtx)<-barcode[,1]
  #sum(Matrix::rowSums(expr.mtx>0)==0)
  #sum(Matrix::colSums(expr.mtx>0)<100) #expr at least 100 genes
  keep.gene=Matrix::rowSums(expr.mtx>0)!=0
  #geneId=geneId[keep.gene,]
  #expr.mtx=expr.mtx[keep.gene,]
  expr.mtx=expr.mtx[,Matrix::colSums(expr.mtx>0)>=100]
  
  sce <- SingleCellExperiment(list(counts=expr.mtx),
                              rowData=DataFrame(geneId))
  sce
  
  is.mito <- grep("^mt:", rownames(sce))
  per.cell <- perCellQCMetrics(sce, subsets=list(Mito=is.mito))
  head(per.cell)
  summary(per.cell$subsets_Mito_percent)
  
  colData(sce) <- cbind(colData(sce), per.cell)
  sce
})
names(two.samples)<-sample.names

sapply(two.samples,dim) 

saveRDS(two.samples,file='two.samples.sceObjs.rds')
#######################################################################
## quality control
two.samples=readRDS('two.samples.sceObjs.rds')
sample.names=names(two.samples)

two.samples2=lapply(1:2,function(i){
  one.sample=two.samples[[i]]
  # filter cell
  one.sample=one.sample[,one.sample$detected>=100]
  one.sample$condition=sample.names[i]
  # add gene stats
  one.sample <- scater::addPerFeatureQC(one.sample, exprs_values = "counts") ## Remove genes with zero total counts across all cells
  rowData(one.sample)$detected_cells <- rowData(one.sample)$detected * ncol(one.sample) / 100
  # cell-lib size normalization
  clusters <- quickCluster(one.sample)
  one.sample <- computeSumFactors(one.sample, clusters=clusters)
  summary(sizeFactors(one.sample))
  one.sample <- logNormCounts(one.sample,log = FALSE)
  assayNames(one.sample) #"counts"     "normcounts"
  one.sample
})
names(two.samples2)=sample.names
saveRDS(two.samples2,file='two.samples.filterCell.sceObjs.rds')

sapply(two.samples2,dim)
