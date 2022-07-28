#https://xuranw.github.io/locCSN/docs/vignettes.html
# Import packages
import locCSN
import os
import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


# Set path to data
#os.chdir('yourpathtodata/Chutype/')

# read in Chutype dataset
data = sc.read_text('logChumarker.txt')
data.shape
data = data.transpose() # 1018 cells * 51 genes

cell_type = pd.read_csv('chutypectname.txt', sep = ' ')
data.obs = cell_type   # The observation are labeled with cell types.

# Plot the Heatmap of gene expression
sc.pl.heatmap(data, data.var.index, groupby= "cell_type", dendrogram = False, swap_axes = True, 
              show_gene_labels= True, cmap='Wistia', figsize=(8,6))

#Calculate Pearson’s Correlation
data_dec = data.X[data.obs.cell_type == "DEC", ]
X_dec = data_dec.transpose()
data_npc = data.X[data.obs.cell_type == 'NPC', ]
X_npc = data_npc.transpose()

corr_dec = np.corrcoef(X_dec)
corr_npc = np.corrcoef(X_npc)


np.fill_diagonal(corr_dec, 0)
np.fill_diagonal(corr_npc, 0)


plt.subplot(1, 2, 1)
plt.imshow(abs(corr_dec), vmin=0, vmax=0.7, cmap='RdPu')
plt.title('DEC', fontweight ="bold")
plt.subplot(1, 2, 2)
plt.imshow(abs(corr_npc), vmin=0, vmax=0.7, cmap='RdPu')
plt.title("NPC", fontweight = "bold")
plt.suptitle("Absolute Pearson`s Correlation", fontsize = 14, fontweight = "bold")
plt.show()

#Calculate CSN test statistics
import time
start = time.time()
csn_dec = locCSN.csn(X_dec, dev = True)
end = time.time()
print(end - start) 
start = time.time()
csn_npc = locCSN.csn(X_npc, dev = True)
end = time.time()
print(end - start) 
#1275 pairs need calculation
#14.49136209487915
#903 pairs need calculation
#6.92561411857605

#Now we show what function csn produces. For a specific cell, we compute each pair of genes and store test statistics in an upper diagnol matrix.
type(csn_dec) 
# list
len(csn_dec) # 138 cells
# Let's see the test statistics for the first cell in DEC
plt.imshow(csn_dec[0].toarray(), vmin = -6, vmax = 6, cmap = 'coolwarm')
plt.title('DEC one cell', fontweight = "bold")
plt.colorbar()
#plt.savefig('dec_one_cell.png')
plt.show()

#For CSN test statistics matrices within a cell group, we first threshold test statistics and averaged the adjacency matrices with the cell group. The averaged CSN is the estimate of gene co-expression of the cell group. In this example, we thresholded at α=0.05.
from scipy.stats import norm

# Cutoff at norm(0.95)
csn_mat = [(item > norm.ppf(0.95)).astype(int) for item in csn_dec]
avgcsn_dec = sum(csn_mat).toarray()/len(csn_mat) + np.transpose(sum(csn_mat).toarray()/len(csn_mat))
csn_mat = [(item > norm.ppf(0.95)).astype(int) for item in csn_npc]
avgcsn_npc = sum(csn_mat).toarray()/len(csn_mat) + np.transpose(sum(csn_mat).toarray()/len(csn_mat))

plt.subplot(1, 2, 1)
plt.imshow(avgcsn_dec, cmap = "Greens", vmin = 0, vmax = 0.7)
plt.title('DEC', fontweight ="bold")
plt.subplot(1, 2, 2)
plt.imshow(avgcsn_npc, cmap = "Greens", vmin = 0, vmax = 0.7)
plt.title('NPC', fontweight = 'bold')
plt.suptitle("Averaged CSN, cut at alpha = 0.05", fontsize=14, fontweight = "bold")
plt.show()

#Comparing networks between two groups of cells (ASD dataset)
# import scanpy as sc
# load data
os.chdir('./Velme/')
data = sc.read_text('Velme_log_mc_cpm_L.txt')
data = data.transpose() 
data.shape # 1778 metacells * 942 genes
meta_L = pd.read_csv('Velme_meta_mc_L.txt', sep = ' ') 
meta_L.columns
# Index(['sampleID', 'broad.cluster', 'cluster', 'diagnosis'], dtype='object')
data.obs = meta_L
sc.pl.heatmap(data, data.var.index, groupby= ["cluster", "diagnosis"], dendrogram = False, swap_axes = True, 
              cmap='Wistia', figsize=(8,4))

data.obs['cluster'].value_counts()
data.obs['diagnosis'].value_counts()
grouped_data = data.obs.groupby(['cluster', 'diagnosis'])
grouped_data.describe()['sampleID']['count']


#CSN Construction of L4 cell-group
ct_name = "L4"
data_L4 = data[data.obs.cluster == ct_name, :]
data_L4.shape # 449 metacell * 942 genes
mcknn = pd.read_csv('mcknn100_' + ct_name + '.txt', sep = ' ')
mcknn = mcknn.to_numpy()
X_L4 = data_L4.X.transpose()

#The runtime of 942 genes is longer than 1 hour. Therefore we provide a toy example only use a subset of 20 genes and the runtime is approximately 20-40s.

start = time.time()
csn_L4_sub = locCSN.csn_loc(X_L4[0:20, :], mcknn)
end = time.time()
print(end-start)
# 13.172302722930908


csn_L4_sub_flat = locCSN.csntoflat(csn_L4_sub)  # 20 genes
csn_L4_sub_flat.shape  #190 gene pairs * 449 cells
# np.savetxt('csn_'+ct_name+'_sub_flat.txt', csn_L4_sub_flat, delimiter = '\t')

#For analysis and visualization, we threshold the CSN test statistics at α=0.01 and average CSN within ASD and Control group respectively.

csn_mat = [(item > norm.ppf(0.99)).astype(int) for item in csn_L4_sub]
meta_L4 = meta_L[meta_L['cluster'] == ct_name]

c_index = (meta_L4['diagnosis'].to_numpy() == 'Control').tolist()
csn_mat_L4_control =  [x for x, y in zip(csn_mat_L4, c_index) if y]
a_index = (meta_L4['diagnosis'].to_numpy() == 'ASD').tolist()
csn_mat_L4_asd =  [x for x, y in zip(csn_mat_L4, a_index) if y]

avgcsn_L4_sub_control = sum(csn_mat_L4_control).toarray()/len(csn_mat_L4_control)
avgcsn_L4_sub_control = + np.transpose(avgcsn_L4_sub_control)
avgcsn_L4_sub_asd = sum(csn_mat_L4_asd).toarray()/len(csn_mat_L4_asd)
avgcsn_L4_sub_asd = + np.transpose(avgcsn_L4_sub_asd)

my_dpi = 50
fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(10, 5), dpi=my_dpi)
print(fig)
print(axes)
fig.suptitle("Avg CSN for L4, 20 genes, cut at alpha = 0.01", fontsize = 14, fontweight = "bold")
axes[0].set_title('Control', fontweight = "bold")
axes[0].imshow(avgcsn_L4_sub_control, cmap = "Greens", vmin = 0, vmax = 0.5)
axes[1].set_title('ASD', fontweight = "bold")
axes[1].imshow(avgcsn_L4_sub_asd, cmap = "Greens", vmin = 0, vmax = 0.5)
#fig.savefig('Velme_Avg_csn_L4_sub.png')


