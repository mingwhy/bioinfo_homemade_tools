#2024-Combining LIANA and Tensor-cell2cell to decipher cell-cell communication across multiple samples
#python: https://ccc-protocols.readthedocs.io/en/latest/notebooks/ccc_python/QuickStart.html
#R: https://ccc-protocols.readthedocs.io/en/latest/notebooks/ccc_R/QuickStart.html
use_gpu=False
import cell2cell as c2c
import liana as li
import pandas as pd
import decoupler as dc
import scanpy as sc
#import matplotlib.pylot as plt
import matplotlib as plt
%matplotlib inline
import plotnine as p9
import seaborn as sns

data_folder='./data/'
output_folder='./data/outputs/'
c2c.io.directories.create_directory(data_folder)
c2c.io.directories.create_directory(output_folder)

# data download follow:#https://ccc-protocols.readthedocs.io/en/latest/notebooks/ccc_R/QuickStart.html
adata=c2c.datasets.balf_covid(data_folder+'BALF-COVID19-Liao_et_al-NatMed-2020.h5ad')

# data preprocessing
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)

adata.var['mt'] = adata.var_names.str.startswith('MT-')  # annotate the group of mitochondrial genes as 'mt'
sc.pp.calculate_qc_metrics(adata,
                           qc_vars=['mt'],
                           percent_top=None,
                           log1p=False,
                           inplace=True)
adata = adata[adata.obs.pct_counts_mt < 15, :]

# save the raw counts to a layer
adata.layers["counts"] = adata.X.copy()

# Normalize the data
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

adata.obs.head()

#Deciphering Cell-Cell Communication
li.method.show_methods()

liana_resources = li.resource.show_resources()
print(*liana_resources, sep = ', ')

lr_pairs = li.resource.select_resource('consensus')

# adata.obs.columns
# adata.obs.celltype
li.mt.rank_aggregate.by_sample(adata,
                               sample_key='sample_new',
                               groupby='celltype',
                               resource_name = 'consensus',
                               expr_prop=0.1, # must be expressed in expr_prop fraction of cells
                               min_cells = 5,
                               n_perms = 100,
                               use_raw = False, # run on log- and library-normalized counts
                               verbose = True,
                               inplace = True
                              )


