#https://zhuanlan.zhihu.com/p/341941329

# source code downloaded from:
# https://github.com/cssmillie/ulcerative_colitis
# there is a 'ulcerative_colitis-master' folder which contains a number of r code scripts

## load required libraries (analysis.r code includes all the libraries for running all analyses,
## but as here only shows the Cellular Composition difference test, only libraries listed below need to be loaded.
library(Seurat)
library(RColorBrewer) #for brewer.pal
library(Matrix) #for Matrix
library(DirichletReg)
library(data.table)
library(tidyverse)
library(cowplot)

## this function is extracted from analysis.r 
dirichlet_regression = function(counts, covariates, formula){  
  # Dirichlet multinomial regression to detect changes in cell frequencies
  # formula is not quoted, example: counts ~ condition
  # counts is a [samples x cell types] matrix
  # covariates holds additional data to use in the regression
  #
  # Example:
  # counts = do.call(cbind, tapply(seur@data.info$orig.ident, seur@ident, table))
  # covariates = data.frame(condition=gsub('[12].*', '', rownames(counts)))
  # res = dirichlet_regression(counts, covariates, counts ~ condition)
  
  #ep.pvals = dirichlet_regression(counts=ep.freq, covariates=ep.cov, formula=counts ~ condition)$pvals
  
  # Calculate regression
  counts = as.data.frame(counts)
  counts$counts = DR_data(counts)
  data = cbind(counts, covariates)
  fit = DirichReg(counts ~ condition, data)
  
  # Get p-values
  u = summary(fit)
  #compared with healthy condition, 15 vars. noninflame and inflame, 30pvalues
  pvals = u$coef.mat[grep('Intercept', rownames(u$coef.mat), invert=T), 4] 
  v = names(pvals)
  pvals = matrix(pvals, ncol=length(u$varnames))
  rownames(pvals) = gsub('condition', '', v[1:nrow(pvals)])
  colnames(pvals) = u$varnames
  fit$pvals = pvals
  
  fit
}



## not all **.r in analysis.r need to be 'source'.
## only these three below are required for cellular composition difference test.
source('mtx.r') #for sparse_cbind
source('plot.r') #for matrix_barplot
source('colors.r') #for set.colors

## Load metadata for discovery and validation cohorts
## data downloaded from: 
## https://portals.broadinstitute.org/single_cell/study/SCP259
meta = read.table('all.meta2.txt', sep='\t', header=T, row.names=1, stringsAsFactors=F)

# Read a list of cell subsets, including the group that each one belongs to
# Groups include: Epithelial, Endothelial, Fibroblasts, Glia, Myeloid, B, and T cells
cell_subsets = read.table('cell_subsets.txt', sep='\t', header=F, stringsAsFactors=F)

# load seurat objects
# data downloaded from:
# https://www.dropbox.com/sh/dn4gwdww8pmfebf/AACXYu8rda5LoLwuCZ8aZXfma?dl=0

epi.seur = readRDS('UC_Data/train.Epi.seur.rds')
fib.seur = readRDS('UC_Data/train.Fib.seur.rds')
imm.seur = readRDS('UC_Data/train.Imm.seur.rds')

# set counts matrices
epi.counts = epi.seur@assays[['RNA']]@counts
fib.counts = fib.seur@assays[['RNA']]@counts
imm.counts = imm.seur@assays[['RNA']]@counts


# Use dirichlet-multinomial regression to find significant changes in cell frequencies during disease
# ---------------------------------------------------------------------------------------------------

# Count each cell subset in every sample
epi.freq = as.matrix(as.data.frame.matrix(table(epi.seur@meta.data$Sample, epi.seur@meta.data$Cluster)))
fib.freq = as.matrix(as.data.frame.matrix(table(fib.seur@meta.data$Sample, fib.seur@meta.data$Cluster)))
imm.freq = as.matrix(as.data.frame.matrix(table(imm.seur@meta.data$Sample, imm.seur@meta.data$Cluster)))

# Combine counts into a single matrix
all.freq = sparse_cbind(list(epi.freq, fib.freq, imm.freq))

# For the validation cohort, we need to combine the replicate samples (because they are not independent)
# To construct a list of replicates, we remove "1", "2", "a", and "b" from the sample IDs
reps = gsub('[12]*[ab]*$', '', rownames(all.freq))
temp = as.matrix(data.frame(aggregate(as.matrix(all.freq), list(reps), sum), row.names=1))
colnames(temp) = colnames(all.freq)
all.freq = temp[,colSums(temp) > 0]

# Split matrix into "epithelial" and "lamina propria" cell subsets and samples
ep.ident = levels(epi.seur@active.ident)
lp.ident = c(levels(fib.seur@active.ident), levels(imm.seur@active.ident))
ep.freq = all.freq[grep('Epi', rownames(all.freq)), ep.ident]
lp.freq = all.freq[grep('LP', rownames(all.freq)), lp.ident]

# For the dirichlet-multinomial regression, we need to know the disease state for each sample
# We can get this from the metadata table as follows:
sample2health = data.frame(unique(data.frame(sample=gsub('[12]*[ab]*$', '', meta[,'Sample']), health=meta[,'Health'])), row.names=1)
ep.cov = data.frame(condition=factor(sample2health[rownames(ep.freq),1], levels=c('Healthy', 'Non-inflamed', 'Inflamed')), row.names=rownames(ep.freq))
lp.cov = data.frame(condition=factor(sample2health[rownames(lp.freq),1], levels=c('Healthy', 'Non-inflamed', 'Inflamed')), row.names=rownames(lp.freq))

# Calculate significant changes using dirichlet multinomial regression
# This returns a matrix of p-values for each cell type / disease state
# 3 condition
dim(ep.freq) #34 sample, 15 cell type
dim(ep.cov) #34 1
head(ep.cov) #condition of each sample
table(ep.cov) #20, 7, 7
ep.pvals = dirichlet_regression(counts=ep.freq, covariates=ep.cov, 
                                formula=counts ~ condition)$pvals
colnames(ep.pvals) = colnames(ep.freq)
lp.pvals = dirichlet_regression(counts=lp.freq, covariates=lp.cov, formula=counts ~ condition)$pvals
colnames(lp.pvals) = colnames(lp.freq)

# Plot epithelial cell proportions
ep.pct = 100*ep.freq/rowSums(ep.freq)
p1 = matrix_barplot(ep.pct, group_by=ep.cov$condition, pvals=ep.pvals, colors=set.colors)
save_plot(p1, file='train.Fig2A.epi_freqs.pdf', nrow=1, ncol=2.5)

# Plot lamina propria cell proportions
lp.pct = 100*lp.freq/rowSums(lp.freq)
p2 = matrix_barplot(lp.pct, group_by=lp.cov$condition, pvals=lp.pvals, colors=set.colors)
save_plot(p2, file='train.Fig2A.lp_freqs.pdf', nrow=1, ncol=2.5)