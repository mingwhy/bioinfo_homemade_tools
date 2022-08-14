# https://www.biorxiv.org/content/10.1101/2022.05.18.492432v1.full

######################################################################
## install python package
# https://gitlab.com/olgaibanez/decibel
git clone https://gitlab.com/olgaibanez/decibel.git

# https://gitlab.com/olgaibanez/scallop
# check https://gitlab.com/olgaibanez/scallop/requirements.txt
# make sure all required python packages are installed
# then: pip install scallop


######################################################################
## download datasets
# generate h5ad of pbmc3k data following: https://scanpy-tutorials.readthedocs.io/en/latest/pbmc3k.html
# or download pbmc.h5ad directly: https://github.com/david758/DESC_Impute/blob/master/desc/datasets/pbmc.h5ad

# download from 20k Human PBMCs, 3' HT v3.1, Chromium X used in `2022-Lack of evidence for increased transcriptional noise in aged tissues` manuscript: https://www.10xgenomics.com/resources/datasets/20-k-human-pbm-cs-3-ht-v-3-1-chromium-x-3-1-high-6-1-0
`Feature / cell matrix HDF5 (raw)`


######################################################################
## test decibel
# following: https://gitlab.com/olgaibanez/decibel
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


##############
import h5py
filename = "20k_PBMC_3p_HT_nextgem_Chromium_X_filtered_feature_bc_matrix.h5"
#filename = '20k_PBMC_3p_HT_nextgem_Chromium_X_raw_feature_bc_matrix.h5'
h5=sc.read_10x_h5(filename)
h5
#AnnData object with n_obs × n_vars = 23837 × 36601
#    var: 'gene_ids', 'feature_types', 'genome'
adata=h5

##############
adata = sc.read('pbmc.h5ad')
adata 
#AnnData object with n_obs × n_vars = 2700 × 32738
#    obs: 'barcode'
#    var: 'gene_ids', 'gene_symbols'

sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.filter_genes(adata, min_cells=3)
sc.pp.filter_cells(adata, min_genes=100)

# pip install harmonypy
sc.pp.pca(adata)
sc.pp.neighbors(adata)
tk.tl.triku(adata)
sc.pp.pca(adata)
sce.pp.harmony_integrate(adata, 'batch')
sc.pp.neighbors(adata, use_rep='X_pca_harmony')
sc.tl.umap(adata)

# https://decibel.readthedocs.io/en/latest/modules.html
#Compute transcriptional noise as in Enge et al, (2017):
dcb.enge_transcriptional_noise(adata, 'batch')
dcb.distance_to_celltype_mean(adata,'batch') 

# batch could be cell.type or cell.cluster or other meta.information
#module.decibel.gcl(adata, num_divisions)
#module.decibel.gcl_per_cell_type_and_batch(adata, num_divisions, batch)
gcl_out=dcb.gcl(adata,50) #Return type pd.DataFrame
gcl_out #length=50 array

inp_sce=sc.read('/Users/mingyang/Documents/Data_mouse_aging_atlas/TMS.gene.data_final/tabula-muris-senis-facs-official-raw-obj.h5ad');       
inp_sce # cell x gene =110824 × 22966
inp_sce.obs['tissue_ct']=inp_sce.obs['tissue'].str.cat(inp_sce.obs['cell_ontology_class'],sep=":")
set(inp_sce.obs['tissue_ct'])

tc='Brain_Non-Myeloid:neuron'
tc='Brain_Non-Myeloid:oligodendrocyte'
tc='Marrow:NK cell'
tc="Heart:atrial myocyte"   
pick=(inp_sce.obs['tissue_ct']==tc) & (inp_sce.obs ['age']=='3m')
test1 = inp_sce[pick]
pick=(inp_sce.obs['tissue_ct']==tc) & (inp_sce.obs ['age']=='24m')
test2 = inp_sce[pick]

sc.pp.normalize_total(test1, target_sum=1e4)
sc.pp.log1p(test1)
sc.pp.filter_genes(test1, min_cells=5)
sc.pp.filter_cells(test1, min_genes=100)

sc.pp.normalize_total(test2, target_sum=1e4)
sc.pp.log1p(test2)
sc.pp.filter_genes(test2, min_cells=5)
sc.pp.filter_cells(test2, min_genes=100)

young_gcl=dcb.gcl(test1,100)
old_gcl=dcb.gcl(test2,100)
print(young_gcl.mean(),old_gcl.mean())

######################################################################
#https://decibel.readthedocs.io/en/latest/modules.html#module.decibel.scallop_pipeline: module.decibel.scallop_pipeline(adata, res_vals=None)

inp_sce=sc.read('/Users/mingyang/Documents/Data_mouse_aging_atlas/TMS.gene.data_final/tabula-muris-senis-facs-official-raw-obj.h5ad');       
inp_sce # cell x gene =110824 × 22966
inp_sce.obs['cell_ontology_class']

#pick=(inp_sce.obs['cell_ontology_class']=='skeletal muscle satellite cell')
pick=(inp_sce.obs['cell_ontology_class']=='skeletal muscle satellite cell') & (inp_sce.obs['age']=='3m')
pick=(inp_sce.obs['cell_ontology_class']=='skeletal muscle satellite cell') & (inp_sce.obs['age']=='24m')
test = inp_sce[pick]

sc.pp.normalize_total(test, target_sum=1e4)
sc.pp.log1p(test)
sc.pp.filter_genes(test, min_cells=3)
sc.pp.filter_cells(test, min_genes=100)
sc.pp.pca(test)
sc.pp.neighbors(test)

#tk.tl.triku(test)
#sc.pp.pca(test)
#test.obs['leiden']
#test.obs['condition']=test.obs['leiden'] #must have a 'condition' column

#test.obs['condition']=test.obs['age'] #must have a 'condition' column
test.obs['condition']=1  #must have a 'condition' column, if already per cell type per age, set it to 1
dcb.scallop_pipeline(test, res_vals=None) #need neighbor information, so sc.pp.neighbors is a must
test.obs

young=test.obs['scallop_noise']
old=test.obs['scallop_noise']

print(young.mean(),old.mean())
print(young.median(),old.median())

######################################################################
## test scallop
import scanpy as sc
import scallop as sl

adata = sc.read('pbmc.h5ad')
adata #2700 × 32738

#Initialize scallop object:
scal = sl.Scallop(adata)

#Run scallop using on 95% of the cells in each iteration (30 iterations) and giving the resolution parameter a value of 1.2.
sl.tl.getScore(scal, res=1.2, n_trials=30, frac_cells=0.95)
adata.obs


