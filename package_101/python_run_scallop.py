
######################################################################
## install python package: https://gitlab.com/olgaibanez/decibel
# git clone https://gitlab.com/olgaibanez/decibel.git

# https://gitlab.com/olgaibanez/scallop
# check https://gitlab.com/olgaibanez/scallop/requirements.txt
# make sure all required python packages are installed
# then: pip install scallop

######################################################################

import pandas as pd
import numpy as np
import scanpy as sc
import scanpy.external as sce
import tqdm as tqdm
import scallop as sl
import triku as tk #pip install triku
import sys
sys.path.append('./decibel/module/')
import decibel as dcb
dcb
dcb.enge_euclidean_dist


inp_sce=sc.read('/Users/mingyang/Documents/Data_mouse_aging_atlas/TMS.gene.data_final/tabula-muris-senis-facs-official-raw-obj.h5ad');       
inp_sce # cell x gene =110824 × 22966
inp_sce.obs['tissue_ct']=inp_sce.obs['tissue'].str.cat(inp_sce.obs['cell_ontology_class'],sep="_")
set(inp_sce.obs['tissue_ct'])

pick=(inp_sce.obs['tissue_ct']=='Brain_Non-Myeloid_neuron') & (inp_sce.obs['age']=='3m') #'Brain_Non-Myeloid_oligodendrocyte'
pick=(inp_sce.obs['tissue_ct']=='Brain_Non-Myeloid_neuron') & (inp_sce.obs['age']=='24m')
#pick=(inp_sce.obs['cell_ontology_class']=='skeletal muscle satellite cell') & (inp_sce.obs['age']=='3m')
#pick=(inp_sce.obs['cell_ontology_class']=='skeletal muscle satellite cell') & (inp_sce.obs['age']=='24m')
test = inp_sce[pick]

#Then process data (normalization, log-transformation, QC filtering on cells and genes).
sc.pp.normalize_total(test, target_sum=1e4)
sc.pp.log1p(test)
sc.pp.filter_genes(test, min_cells=3) #test
sc.pp.filter_cells(test, min_genes=100) #test
sc.pp.pca(test) #test.obsm
sc.pp.neighbors(test) #test.obsp['distances']

#Run PCA, feature selection (triku)
tk.tl.triku(test) #test.var
sc.pp.pca(test) #test.obsm['X_pca']

test.obs['condition']=1  #must have a 'condition' column, if already per cell type per age, set it to 1
dcb.scallop_pipeline(test, res_vals=None) #need neighbor information, so sc.pp.neighbors is a must
test.obs

young=test.obs['scallop_noise']
old=test.obs['scallop_noise']

print(young.mean(),old.mean())
print(young.median(),old.median())


