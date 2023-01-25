## read in data
library(anndata) #https://anndata.dynverse.org/articles/getting_started.html
library(reticulate)
#Sys.setenv(RETICULATE_PYTHON = "/Users/mingyang/anaconda3/envs/RCFGL/bin/python")
#py_config()
sc=import('scanpy')  #py_install('scanpy')
ad=sc$read_h5ad('/Users/mingyang/Documents/Data_worm_aing/ad_worm_aging_umi.h5ad') #https://scanpy.readthedocs.io/en/stable/generated/scanpy.read_h5ad.html
ad$obs
ad$var
ad$X$shape #47423   20305
class(ad$X)
mat=ad$X$toarray()
class(mat) #"matrix" "array" 
dim(mat) #  47423 20305