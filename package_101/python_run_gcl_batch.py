
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
dcb.enge_euclidean_dist #make sure this function works

#inp_sce=sc.read('/Users/mingyang/Documents/Data_mouse_aging_atlas/TMS.gene.data_final/tabula-muris-senis-facs-official-raw-obj.h5ad');       
inp_sce=sc.read('/Users/mingyang/Documents/aging_cell.turnover/DV_0619_normcounts_lifespan/select.tc_subsampled.h5ad');
inp_sce # cell x gene =110824 × 22966, 22958 × 22966
inp_sce.obs['tissue_ct']=inp_sce.obs['tissue'].str.cat(inp_sce.obs['cell_ontology_class'],sep=":")

all_tc=list(set(inp_sce.obs['tissue_ct']))

gcl_out=pd.DataFrame()
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
		
		out=dcb.gcl(test1,100)
		data = {'tc':np.repeat(tc,100), 'age': np.repeat(age,100), 'gcl':out}
		df=pd.DataFrame(data)		
		gcl_out=gcl_out.append(df)	

gcl_out.to_pickle("gcl_out.pkl")  


