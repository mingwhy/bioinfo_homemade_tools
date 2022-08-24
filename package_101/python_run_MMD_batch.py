
# https://github.com/calico/scmmd
######################################################################
import os
#os.environ['KMP_DUPLICATE_LIB_OK']='True'
import scmmd
import anndata

import pandas as pd
import numpy as np
import scanpy as sc
import scanpy.external as sce
import tqdm as tqdm
import triku as tk #pip install triku
import sys
import warnings
warnings.filterwarnings("ignore")

#inp_sce=sc.read('/Users/mingyang/Documents/aging_cell.turnover/0619_DV_normcounts/select.tc_subsampled.h5ad');
inp_sce=sc.read('select.tc_subsampled_minCellCounts.h5ad');
inp_sce # cell x gene =22958 × 22909

all_tc=list(set(inp_sce.obs['tissue_cell.type']))
len(all_tc) #48

scallop_out=pd.DataFrame()

for tc in all_tc:
	print(tc)
	sce=inp_sce[inp_sce.obs['tissue_cell.type']==tc]
	
	#Then process data (normalization, log-transformation, QC filtering on cells and genes).
	#sc.pp.normalize_total(test1, target_sum=1e4)
	#sc.pp.log1p(test1)
	#sc.pp.filter_genes(test1, min_cells=5) #test
	#sc.pp.filter_cells(test1, min_genes=100) #test
	sc.pp.pca(sce) #test.obsm
	
	#sc.pp.neighbors(test1,n_neighbors=15) #test.obsp['distances'] #this step
	#tk.tl.triku(test1) #need to run PCA, before feature selection (triku), test.var
	#sc.pp.pca(test1) #test.obsm['X_pca']		
	#sc.pp.neighbors(test1,n_neighbors=15) #test.obsp['distances'] #have must at least 2 clusters
	# change `#if adata_sample.shape[0] > 50:` in decibel.py into `adata_sample.shape[0] >= 40:` as my tc contain>=50cells

	ages=sce.obs['age'].unique()
	for i in ages[0:1]:
		for j in ages[1:2]:
			print(i,j,tc)
			test=sce[sce.obs['age'].isin([i,j])]			
			distances, p_values = scmmd.compute_mmd_contrast(
			    adata=test, # [Cells, Genes] object
			    representation='X_pca', # representation to use, "X" or key in `adata.obsm`.
			    groupby='tissue_cell.type', # a categorical grouping variable in `adata.obs`
			    contrast='age', # a binary contrast in `adata.obs`
			    n_iters=100, # number of random sampling iterations
			    sample_size=20, # sample size for random samples
			    n_permutations=1000, # permutations for p-val calculations
			)
			#type(distances) #numpy.ndarray
			#distances.shape #(1, 100, 3)
			#p_values.shape #(1, 100, 3)
			arr = np.concatenate([distances[0,:,:], p_values[0,:,:]], axis=1)
			col=np.column_stack((np.repeat(age1,100),np.repeat(age2,100),np.repeat(tc,100))) #n_iters
			col.shape
			distances = np.append(arr, col, axis=1)
			distances.shape #add age1, age2, tc info
			df=pd.DataFrame(distances)
			scallop_out=scallop_out.append(df)


scallop_out.to_pickle("MMD_out.pkl")  


