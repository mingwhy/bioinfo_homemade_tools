
library(zellkonverter)
library(SingleCellExperiment)
library(tidyverse)
library(ggExtra) #
library(scRNAseq)
library(ggplot2);theme_set(theme_classic())
library(scran)
library(ggpubr)
library(Seurat)
library(grDevices);library(RColorBrewer)
library(SingleCellExperiment)
plotcol <- brewer.pal(8,'Dark2')
library(glmpca)
library(scry)

out=readRDS('sce.filtered_glm_PCA.rds')
######################################################################
## install python package: https://gitlab.com/olgaibanez/decibel
# git clone https://gitlab.com/olgaibanez/decibel.git

# https://gitlab.com/olgaibanez/scallop
# check https://gitlab.com/olgaibanez/scallop/requirements.txt
# make sure all required python packages are installed
# then: pip install scallop

######################################################################
library(reticulate)
Sys.setenv(RETICULATE_PYTHON = "/Users/mingyang/anaconda3/envs/RCFGL/bin/python")
py_config()
pd=import('pandas')
np=import('numpy')
sc=import('scanpy') #py_install('scanpy')
sce=import('scanpy.external')
tqdm=import('tqdm')
sys=import('sys')

# https://gitlab.com/olgaibanez/decibel
sys$path=c(sys$path,'./decibel/module/') #sys.path.append('./decibel/module/')
sys$path 
dcb=import('decibel') 
dcb
dcb$enge_euclidean_dist #make sure this function works

#https://gitlab.com/olgaibanez/scallop
if(F){
  # check out https://gitlab.com/olgaibanez/scallop/requirements.txt
  py_install('dask',pip=TRUE)
  py_install('leidenalg',pip=TRUE)
  py_install('louvain',pip=TRUE)
  py_install('phate',pip=TRUE)
  py_install('cython',pip=TRUE)
  py_install('pytest',pip=TRUE)
  py_install('ray',pip=TRUE)
  py_install('psutil',pip=TRUE)
  # then install scallop
  py_install('scallop',pip=TRUE)
}
sl=import('scallop') #py_install('scallop') 

# begin scallop
#test1.obs['condition']=1  #must have a 'condition' column, if already per cell type per age, set it to 1
#dcb.scallop_pipeline(test1, res_vals=None) #need neighbor information, so sc.pp.neighbors is a must
#out=test1.obs
tmp=out[[1]]
reducedDimNames(tmp)
reducedDims(tmp) <- list(X_PCA=reducedDims(tmp)[[1]]) #https://bioconductor.org/packages/devel/bioc/vignettes/SingleCellExperiment/inst/doc/intro.html#3_Adding_low-dimensional_representations
reducedDimNames(tmp)

test1 = SCE2AnnData(tmp)  #http://www.bioconductor.org/packages/devel/bioc/vignettes/zellkonverter/inst/doc/zellkonverter.html

#Then process data (normalization, log-transformation, QC filtering on cells and genes).
sc$pp$normalize_total(test1, target_sum=1e4)
sc$pp$log1p(test1)
sc$pp$filter_genes(test1, min_cells=5) #test
sc$pp$filter_cells(test1, min_genes=100) #test
sc$pp$pca(test1) #test.obsm
sc$pp$neighbors(test1,n_neighbors=15) #test.obsp['distances']

dcb$scallop_pipeline(adata) #https://decibel.readthedocs.io/en/latest/modules.html
