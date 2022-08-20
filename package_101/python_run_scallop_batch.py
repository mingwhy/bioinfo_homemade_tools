
######################################################################
## install python package: https://gitlab.com/olgaibanez/decibel
# git clone https://gitlab.com/olgaibanez/decibel.git

# https://gitlab.com/olgaibanez/scallop
# check https://gitlab.com/olgaibanez/scallop/requirements.txt
# make sure all required python packages are installed
# then: pip install scallop

######################################################################

## conda activate RCFGL or conda activate myenv beforehand

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
dcb.enge_euclidean_dist #make sure this function works

#inp_sce=sc.read('/Users/mingyang/Documents/Data_mouse_aging_atlas/TMS.gene.data_final/tabula-muris-senis-facs-official-raw-obj.h5ad');       
inp_sce=sc.read('/Users/mingyang/Documents/aging_cell.turnover/0619_DV_normcounts/select.tc_subsampled.h5ad');
inp_sce # cell x gene =110824 × 22966, 22958 × 22966
inp_sce.obs['tissue_ct']=inp_sce.obs['tissue'].str.cat(inp_sce.obs['cell_ontology_class'],sep=":")

all_tc=list(set(inp_sce.obs['tissue_ct']))
scallop_out=pd.DataFrame()

for tc in all_tc:
	print(tc)
	inp_sce_age=inp_sce[inp_sce.obs['tissue_ct']==tc]
	ages=list(set(inp_sce_age.obs['age']))

	for age in ages:
		pick= (inp_sce_age.obs['age']==age) 
		test1 = inp_sce_age[pick]
		#Then process data (normalization, log-transformation, QC filtering on cells and genes).
		sc.pp.normalize_total(test1, target_sum=1e4)
		sc.pp.log1p(test1)
		sc.pp.filter_genes(test1, min_cells=5) #test
		sc.pp.filter_cells(test1, min_genes=100) #test
		sc.pp.pca(test1) #test.obsm
		sc.pp.neighbors(test1,n_neighbors=15) #test.obsp['distances']
		
		tk.tl.triku(test1) #need to run PCA, before feature selection (triku), test.var
		sc.pp.pca(test1) #test.obsm['X_pca']		
		sc.pp.neighbors(test1,n_neighbors=15) #test.obsp['distances'] #have must at least 2 clusters
		# change `#if adata_sample.shape[0] > 50:` in decibel.py into `adata_sample.shape[0] >= 40:` as my tc contain>=50cells

		test1.obs['condition']=1  #must have a 'condition' column, if already per cell type per age, set it to 1
		dcb.scallop_pipeline(test1, res_vals=None) #need neighbor information, so sc.pp.neighbors is a must
		out=test1.obs
		#noise=test1.obs['scallop_noise']
		#print(noise.mean(),noise.mean())
		#print(noise.median(),noise.median())
		#scallop_out.columns.values
		#scallop_out.groupby(['age'])['scallop_noise'].mean()
		scallop_out=scallop_out.append(out)

scallop_out.to_pickle("scallop_out.pkl")  

## conda deactivate 
