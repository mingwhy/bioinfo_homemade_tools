
library(EasyMMD)
library(SingleCellExperiment)
library(scater);library(scran)
library(ggplot2);library(gridExtra)
library(tidyverse);
library(SummarizedExperiment)

sce.minCellCounts=readRDS('~/Documents/aging_cell.turnover/1026_TMS_male_analysis/sce_minCellCounts.rds')
tc.names=names(sce.minCellCounts) #39 tc
tc.names=tc.names[tc.names!='Marrow:NK cell']

n.sample.cell=50;
i.rep=10;
DM.out<-lapply(tc.names,function(tc){
  sce_naive=sce.minCellCounts[[tc]]
  assayNames(sce_naive)
  expr.m=assay(sce_naive,'counts')
  
  #n.expr.gene=Matrix::colSums(expr.m>0)
  n.expr.prop.cell=Matrix::rowSums(expr.m>0)/ncol(expr.m)
  expr.m=expr.m[n.expr.prop.cell>0.2,] #remove gene which expr in less than 20% cells
  
  young.cells=expr.m[,sce_naive$age=='3m']
  old.cells=expr.m[,sce_naive$age=='24m']
  n.young.cell=ncol(young.cells)
  n.old.cell=ncol(old.cells)
  dist.mmds<-lapply(1:i.rep,function(i){
    sample1=young.cells[,sample(1:n.young.cell,n.sample.cell,replace = F)]
    sample2=old.cells[,sample(1:n.old.cell,n.sample.cell,replace = F)]
    x=as.matrix(sample1);y=as.matrix(sample2)
    dist.mmd=MMD(y, x,var = diag(rep(0.5, ncol(y))), bias = TRUE) #https://rdrr.io/github/AnthonyEbert/EasyMMD/man/MMD.html
    return(dist.mmd)
  })
  out=data.frame(tc=tc,mmd=unlist(dist.mmds))
  return(out)
})
saveRDS(DM.out,'MMD_sample50rep10.rds')
