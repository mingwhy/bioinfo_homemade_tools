
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

#inp_sce=sc.read('/Users/mingyang/Documents/Data_mouse_aging_atlas/TMS.gene.data_final/tabula-muris-senis-facs-official-raw-obj.h5ad');       
inp_sce=sc.read('/Users/mingyang/Documents/aging_cell.turnover/DV_0619_normcounts_lifespan/select.tc_subsampled.h5ad');
inp_sce # cell x gene =110824 × 22966, 22958 × 22966
inp_sce.obs['tissue_ct']=inp_sce.obs['tissue'].str.cat(inp_sce.obs['cell_ontology_class'],sep=":")
set(inp_sce.obs['tissue_ct'])


tc='Brain_Non-Myeloid:neuron'
tc='Brain_Non-Myeloid:oligodendrocyte'
tc='Marrow:NK cell'
tc="Heart:atrial myocyte"   
pick1=(inp_sce.obs['tissue_ct']==tc) & (inp_sce.obs['age']=='3m') #'Brain_Non-Myeloid_oligodendrocyte'
pick2=(inp_sce.obs['tissue_ct']==tc) & (inp_sce.obs['age']=='24m')
#pick=(inp_sce.obs['cell_ontology_class']=='skeletal muscle satellite cell') & (inp_sce.obs['age']=='3m')
#pick=(inp_sce.obs['cell_ontology_class']=='skeletal muscle satellite cell') & (inp_sce.obs['age']=='24m')
test1 = inp_sce[pick1]
test2 = inp_sce[pick2]

#Then process data (normalization, log-transformation, QC filtering on cells and genes).
sc.pp.normalize_total(test1, target_sum=1e4)
sc.pp.log1p(test1)
sc.pp.filter_genes(test1, min_cells=5) #test
sc.pp.filter_cells(test1, min_genes=100) #test
sc.pp.pca(test1) #test.obsm
sc.pp.neighbors(test1) #test.obsp['distances']
tk.tl.triku(test1) #Run PCA, feature selection (triku), test.var
sc.pp.pca(test1) #test.obsm['X_pca']
test1.obs['condition']=1  #must have a 'condition' column, if already per cell type per age, set it to 1
dcb.scallop_pipeline(test1, res_vals=None) #need neighbor information, so sc.pp.neighbors is a must
young=test1.obs['scallop_noise']

sc.pp.normalize_total(test2, target_sum=1e4)
sc.pp.log1p(test2)
sc.pp.filter_genes(test2, min_cells=5) #test
sc.pp.filter_cells(test2, min_genes=100) #test
sc.pp.pca(test2) #test.obsm
sc.pp.neighbors(test2) #test.obsp['distances']
tk.tl.triku(test2) #Run PCA, feature selection (triku), test.var
sc.pp.pca(test2) #test.obsm['X_pca']
test2.obs['condition']=1  #must have a 'condition' column, if already per cell type per age, set it to 1
dcb.scallop_pipeline(test2, res_vals=None) #need neighbor information, so sc.pp.neighbors is a must
old=test2.obs['scallop_noise']

print(young.mean(),old.mean())
print(young.median(),old.median())


