#https://earmingol.github.io/cell2cell/tutorials/ASD/01-Tensor-Factorization-ASD/
#conda info --env
#conda activate python37

#ipython
use_gpu = False
import cell2cell as c2c
import scanpy as sc

import numpy as np
import pandas as pd

from tqdm import tqdm

#import matplotlib.pyplot as plt
import matplotlib as plt
import seaborn as sns
%matplotlib inline

import warnings
warnings.filterwarnings('ignore')
c2c.__version__

# load data
import os
data_folder = './data/'
directory = os.fsencode(data_folder)

output_folder = './results/'
if not os.path.isdir(output_folder):
    os.mkdir(output_folder)

# Load data into an AnnData object
rnaseq = sc.read_text(data_folder + '/exprMatrix.tsv.gz')
rnaseq = rnaseq.transpose()

#Load metadata
meta = pd.read_csv(data_folder + '/meta.tsv', sep='\t', index_col=0)

# Add metadata to the AnnData object
rnaseq.obs = rnaseq.obs.join(meta)
rnaseq
rnaseq.obs.columns
rnaseq.obs.age.value_counts()

meta.head()
meta.columns
meta.age[:,].describe()

#In this analysis, we only included samples from the prefrontal cortex (PFC)
brain_area = 'PFC'
rnaseq = rnaseq[rnaseq.obs.region == brain_area]
rnaseq

#Data for Ligand-Receptor pairs
#lr_pairs = pd.read_csv('https://raw.githubusercontent.com/LewisLabUCSD/Ligand-Receptor-Pairs/master/Human/Human-2020-Jin-LR-pairs.csv')
lr_pairs = pd.read_csv('data/Human-2020-Jin-LR-pairs.csv')
lr_pairs = lr_pairs.astype(str)
lr_pairs.head(2)
int_columns = ('ligand_symbol', 'receptor_symbol')

#Data Preprocessing
context_dict = dict()

for diag, df in rnaseq.obs.groupby('age'):    
#for diag, df in rnaseq.obs.groupby('diagnosis'):
    for donor in df['sample'].unique():
        context_dict[donor] = diag
context_dict
context_names = list(context_dict.keys())

#Generate list of RNA-seq data aggregated into cell types
rnaseq_matrices = []

# Iteraty by sample/context
for context in tqdm(context_names):
    # Obtain metadata for context
    meta_context = rnaseq.obs.loc[rnaseq.obs['sample'] == context].copy()
    # Single cells in the context
    cells = list(meta_context.index)

    # Rename index name to identify the barcodes when aggregating expression
    meta_context.index.name = 'barcode'

    # Subset RNAseq data by the single cells in the sample/context
    tmp_data = rnaseq[cells]

    # Keep genes in each sample with at least 4 single cells expressing it
    genes = sc.pp.filter_genes(tmp_data, min_cells=4, inplace=False)[0]
    tmp_data = tmp_data.to_df().loc[:, genes]

    # Aggregate gene expression of single cells into cell types
    exp_df = c2c.preprocessing.aggregate_single_cells(rnaseq_data=tmp_data,
                                                      metadata=meta_context,
                                                      barcode_col='barcode',
                                                      celltype_col='cluster',
                                                      method='nn_cell_fraction',
                                                     )

    rnaseq_matrices.append(exp_df)


# Change gene names to ensembl (here they are annotated as ENSEMBL|SYMBOL)
matrices = []
for rna in rnaseq_matrices:
    tmp = rna.copy()
    tmp.index = [idx.split('|')[0] for idx in rna.index]
    matrices.append(tmp)

# LR pairs
lr_pairs = c2c.preprocessing.ppi.remove_ppi_bidirectionality(ppi_data=lr_pairs, 
                                                             interaction_columns=int_columns
                                                             )
lr_pairs.shape

# Generate a dictionary with function info for each LR pairs.
ppi_functions = dict()

for idx, row in lr_pairs.iterrows():
    ppi_label = row[int_columns[0]] + '^' + row[int_columns[1]]
    ppi_functions[ppi_label] = row['annotation']

ensembl_symbol = dict()

for idx, row in lr_pairs.iterrows():
    ensembl_symbol[row['interaction_ensembl']] = row['interaction_symbol']

# Tensor-cell2cell Analysis
tensor = c2c.tensor.InteractionTensor(rnaseq_matrices=matrices,
                                      ppi_data=lr_pairs,
                                      context_names=context_names,
                                      how='inner',
                                      complex_sep='&',
                                      interaction_columns=('ligand_ensembl', 'receptor_ensembl'),
                                      communication_score='expression_mean',
                                     )
len(context_names) #23
tensor.tensor.shape
# (23, 749, 16, 16)

# If using a GPU, convert tensor & mask into a GPU-manipulable object.
if use_gpu:
    tensor.tensor = tl.tensor(tensor.tensor, device='cuda:0')
    if tensor.mask is not None:
        tensor.mask = tl.tensor(tensor.mask, device='cuda:0')

# Put LR pair names from ensembl to symbol
tensor.order_names[1] = [ensembl_symbol[lr] for lr in tensor.order_names[1]]

meta_tf = c2c.tensor.generate_tensor_metadata(interaction_tensor=tensor,
                                              metadata_dicts=[context_dict, ppi_functions, None, None],
                                              fill_with_order_elements=True
                                             )
meta_tf[0]

#Run Analysis
elbow, error = tensor.elbow_rank_selection(upper_rank=25,
                                           runs=1, # This can be increased for more robust results
                                           init='random',
                                           automatic_elbow=False,
                                           random_state=888,
                                          )

# If automatic_elbow=True, remove these two lines. To save the figure in that case,
# add the parameter filename=output_folder + 'Elbow.svg' in the previous function.
# The number of factors will be saved in tensor.rank
import matplotlib as matplotlib
_ = matplotlib.pyplot.plot(error) 
_ = matplotlib.pyplot.plot(*error[5], 'ro') # Here we selected a number of 6 factors.
matplotlib.pyplot.savefig(output_folder + 'Elbow.svg', dpi=300, bbox_inches='tight')


#Perform tensor factorization
tensor.compute_tensor_factorization(rank=8,
                                    init='random',
                                    #init='svd', 
                                    random_state=888
                                   )
# init='svd' helps to get an global-optimal solution.
# Replace by 'random' if a memory error is obtained.

# Results
# Color palettes for each of the tensor dimensions.
cmaps = ['viridis', 'Dark2_r', 'tab20', 'tab20']
factors, axes = c2c.plotting.tensor_factors_plot(interaction_tensor=tensor,
                                                 metadata = meta_tf, # This is the metadata for each dimension
                                                 sample_col='Element',
                                                 group_col='Category',
                                                 meta_cmaps=cmaps,
                                                 fontsize=14,
                                                 filename=output_folder + 'Tensor-Factorization.svg'
                                                )

#Top-10 LR pairs from their factor-specific loadings
for i in range(tensor.rank):
        print(tensor.get_top_factor_elements(order_name='Ligand-Receptor Pairs', 
                                             factor_name='Factor {}'.format(i+1), 
                                             top_number=10))
        print('')

#Export Loadings for all dimensions of the 4D-communication tensor.
tensor.export_factor_loadings(output_folder + 'Loadings.xlsx')

