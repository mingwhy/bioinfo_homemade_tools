
library(zellkonverter)
library(SingleCellExperiment)
library(scater);library(scran)
library(ggplot2);library(gridExtra);
library(tidyverse);
library(zellkonverter)
library(SummarizedExperiment)
library(ggpointdensity)
library(viridis)
library(patchwork); #for plot_annotation
library(lme4)
library(variancePartition)
library(ggpubr)

## read in data
sce=readH5AD('../1120_TMS_male_analysis/select.tc.h5ad') # 22966 31001 
unique(sce$age) #3m 18m 24m
tc.names=sort(unique(sce$tissue_cell.type))
tc.names #38
sce.shared=lapply(tc.names,function(tc){
  sce[,sce$tissue_cell.type==tc]
})
names(sce.shared)<-tc.names
###############################################################################################
## reading in gene id.mapping and dn.ds infomation
id.mapping=data.table::fread('~/Documents/Data_mouse_aging_atlas/fac_20449genes_id.mapping.txt')  #readin_h5ad.R
mouse_geneOrigin=data.table::fread('../gene.properties/gene.age/Mus_musculus.csv')
head(mouse_geneOrigin)
dim(mouse_geneOrigin) #22832 genes
dim(id.mapping) #20449
sum(mouse_geneOrigin$ensembl_id %in% id.mapping$ensembl_gene_id) #18252


gene.meta=merge(mouse_geneOrigin,id.mapping,by.x='ensembl_id',by.y='ensembl_gene_id')
dim(gene.meta) #18252

unique(gene.meta$gene_age)
table(gene.meta$gene_age)
gene.meta[gene.meta$gene_age==">4290",]$gene_age='5000'

###############################################################
## for each cell type, commonly expressed gene age distribution
sapply(sce.shared,dim)
(tc.names=names(sce.shared))

mouse_tcs_expr.genes<-lapply(tc.names,function(tc){
  sce_naive=sce.shared[[tc]]
  assayNames(sce_naive)<-'counts'
  mat=assay(sce_naive,'counts')
  mat=mat[,sce_naive$age=='3m']
  
  #ncell=ncol(mat)
  #i=Matrix::rowSums(mat>0)>=ncell*0.2
  #rownames(mat)[i]
  i=Matrix::rowSums(mat)
  i=i[i!=0]
  df.i=data.frame(gene=names(i),expr=log10(i));
  
  df=merge(df.i,gene.meta,by.x='gene',by.y='mgi_symbol')
  df=df[df$gene_age!='5000',]
  tai=sum(df$expr * as.numeric(df$gene_age))/sum(df$expr)
})
names(mouse_tcs_expr.genes)<-tc.names
df=data.frame(cell.type=tc.names,TAI=unlist(mouse_tcs_expr.genes))

