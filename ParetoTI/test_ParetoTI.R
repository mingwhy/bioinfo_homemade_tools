
library(Seurat)
library(ggplot2);library(gridExtra);
options(stringsAsFactors = F)

out.dir1='brain_scRNA-seq_n15c0.005/';
if(!dir.exists(out.dir1)) dir.create(out.dir1)

file="../single.cell_datasets/fly.brain.atlas/wholebrain_filtered_valid.rds";

dat=readRDS(file);
#dat=NormalizeData(dat);
#df.umi=dat@assays$RNA@counts #12094 features across 56192 samples within 1 assay 
df.expr=dat@assays$RNA@data #logNormal

unique(dat@meta.data$sex)
table(dat@meta.data$annotation)

## separate female and male, select cell clusters that contain>=500cells in both sexes
min.cell=200;max.cell=Inf;
i=apply(table(dat$sex,dat$annotation),2,function(i)sum(i>=min.cell & i<max.cell))
#pick.cell.clusters=names(which(i==2)) #79
pick.cell.clusters=names(which(i>=1)) #79
pick.cell.clusters #min200: sn,40. sc,38. min100: sn,54; sc,60.

dat=subset(dat,annotation %in% pick.cell.clusters)
df.expr=dat@assays$RNA@counts;
dim(df.expr); #sn:12602 50588

sparsity=sum(df.expr==0)/nrow(df.expr)/ncol(df.expr)
cat(out.dir1,sparsity,'\n') #min100:sn=0.9537,sc=0.9055
#min200: sc=0.9051, sn=0.9544

#i.cluster='Ensheathing_glia'
#i.cluster='Astrocyte-like'
i.cluster=c('Pm1/Pm2','Ensheathing_glia');

mat.meta=dat@meta.data[dat$annotation %in% i.cluster,]
dim(mat.meta) #1508 cells
table(mat.meta$annotation,mat.meta$Age)

mat=df.expr[,dat$annotation %in% i.cluster]
gene.filter=Matrix::rowSums(mat>0) >= max(15,ncol(mat)*0.005)
mat=mat[gene.filter,]
dim(mat); #6533 x 1508


########################################################################################################
# check that py_pcha library is successfully installed and discoverable
reticulate::py_discover_config("py_pcha")
# To make sure R uses the correct conda enviroment you can run this when you start R:
reticulate::use_condaenv("reticulate_PCHA", conda = "auto",
                         required = TRUE) # set TRUE to force R to use reticulate_PCHA
library(ParetoTI)
library(SingleCellExperiment)
library(ggplot2)
library(cowplot)
library(Matrix)

# convert to single cell experiment
mat[1:3,1:3]
sce = SingleCellExperiment(assays = list(counts = mat),
                                   colData = mat.meta)
# remove genes with too many zeros (> 95% cells)
sce = sce[rowMeans(counts(sce) > 0) > 0.05,]
# remove cells with too many zeros (> 85%)
sce = sce[,colMeans(counts(sce) > 0) > 0.15]
sce #3562 1166 


# Normalise gene expression by cell sum factors and log-transform
sce = scran::computeSumFactors(sce)#https://github.com/LTLA/scuttle/issues/12
sce = scater::logNormCounts(sce,transform='none')
#assays(2): counts normcounts
sce = scuttle::logNormCounts(sce) #default: log2. https://rdrr.io/github/LTLA/scuttle/man/logNormCounts.html
sce #assays(3): counts logcounts normcounts

#Plot below shows first 3 PCs colored by batch.
# Find principal components
sce = scater::runPCA(sce, ncomponents = 7,scale = T, exprs_values = "logcounts")
# Plot PCA colored by batch
colnames(sce@colData)
scater::plotReducedDim(sce, ncomponents = 3, dimred = "PCA", colour_by = "Age")

# extract PCs (centered at 0 with runPCA())
PCs4arch = t(reducedDim(sce, "PCA"))
sce #3562 1166 
dim(PCs4arch) #7PC X 1166 cells


#Fit k=2:8 polytopes to Hepatocytes to find which k best describes the data
# find archetypes
arc_ks = k_fit_pch(PCs4arch, ks = 2:8, check_installed = T,
                   bootstrap = T, bootstrap_N = 200, maxiter = 1000,
                   bootstrap_type = "m", seed = 2543, 
                   volume_ratio = "t_ratio", # set to "none" if too slow
                   delta=0, conv_crit = 1e-04, order_type = "align",
                   sample_prop = 0.75)
arc_ks$summary

# Show variance explained by a polytope with each k (cumulative)
plot_arc_var(arc_ks, type = "varexpl", point_size = 2, line_size = 1.5) + theme_bw()

# Show variance explained by k-vertex model on top of k-1 model (each k separately)
plot_arc_var(arc_ks, type = "res_varexpl", point_size = 2, line_size = 1.5) + theme_bw()

# Show variance in position of vertices obtained using bootstraping 
# - use this to find largest k that has low variance
plot_arc_var(arc_ks, type = "total_var", point_size = 2, line_size = 1.5) +
  theme_bw() +
  ylab("Mean variance in position of vertices")


# Show t-ratio
plot_arc_var(arc_ks, type = "t_ratio", point_size = 2, line_size = 1.5) + theme_bw()


# find archetypes on all data (allows using archetype weights to describe cells)
arc_1 = fit_pch(PCs4arch, volume_ratio = "t_ratio", maxiter = 500,
                noc = 4, delta = 0,
                conv_crit = 1e-04)
dim(PCs4arch) #7 1166
dim(arc_1$S) #4 1166, 4 coordinates of 1166 cells

# 2D plots
p1 <- plot_arc(arc_data = arc_1, data = PCs4arch,
               data_lab = apply(arc_1$S, 2, max), #adds color to vertex by max archetype score
               data_alpha = 0.75,
               data_size = 2,
               which_dimensions = 1:2) + 
  theme_classic(base_size=18)

#BiocManager::install('wesanderson')
library(wesanderson)
COLS <- wes_palette("Darjeeling2", 9, type = c("continuous"))
#COLS
CO <- rank(unique(sce$Age))
p2 <- plot_arc(arc_data = arc_1, data = PCs4arch,
               data_lab = sce$Age,
               data_alpha = 1, 
               data_size = 2,
               #colors = c(COLS[CO], "red"),
               which_dimensions = 1:2)+
               #text_size = 1) + 
  theme_classic(base_size = 28)

p3 <- plot_arc(arc_data = arc_1, data = PCs4arch,
               data_lab = sce$Age,
               data_alpha = 1, 
               data_size = 2,
               #colors = c(COLS[CO], "red"),
               which_dimensions = 2:3)+
  #text_size = 1) + 
  theme_classic(base_size = 28)

grid.arrange(p1,p2,ncol=2)
grid.arrange(p2,p3,ncol=2)


#Find genes and gene sets enriched near vertices
# Map GO annotations and measure activities: https://github.com/vitkl/ParetoTI/blob/master/R/map_go_annot.R
sce; #3562 1166 
expr_mat<-assay(sce,'logcounts')
dim(expr_mat) #3562 1166
activ = measure_activity(expr_mat, # row names are assumed to be gene identifiers
                         which = "BP", return_as_matrix = F,
                         taxonomy_id = 7227, keytype = "ALIAS",
                         #lower = 20, upper = 1000,
                         aucell_options = list(aucMaxRank = nrow(sce) * 0.1,
                                               binary = F, nCores = 3,
                                               plotStats = FALSE))
dim(activ) #1166 6051
activ[1:3,1:3]


# Merge distances, gene expression and gene set activity into one matrix
data_attr = merge_arch_dist(arc_data = arc_1, data = PCs4arch, 
                            feature_data = as.matrix(logcounts(sce)),
                            colData = activ,
                            dist_metric = c("euclidean", "arch_weights")[1],
                            colData_id = "cells", rank = F) 

# Use Wilcox test to find genes maximally expressed in 10% closest to each vertex
enriched_genes = find_decreasing_wilcox(data_attr$data, data_attr$arc_col,
                                        features = data_attr$features_col,
                                        bin_prop = 0.1, method = "BioQC")
dim(enriched_genes) #14248     8


enriched_sets = find_decreasing_wilcox(data_attr$data, data_attr$arc_col,
                                       features = data_attr$colData_col,
                                       bin_prop = 0.1, method = "BioQC")
length(enriched_sets) #8

# Take a look at top genes and functions for each archetype
labs = get_top_decreasing(summary_genes = enriched_genes, summary_sets = enriched_sets,
                          cutoff_genes = 0.01, cutoff_sets = 0.05, 
                          cutoff_metric = "wilcoxon_p_val", 
                          p.adjust.method = "fdr",
                          order_by = "mean_diff", order_decreasing = T,
                          min_max_diff_cutoff_g = 0.4, min_max_diff_cutoff_f = 0.03)

p_pca = plot_arc(arc_data = arc_1, data = PCs4arch,
                 which_dimensions = 1:3, line_size = 1.5,
                 data_lab = activ$ribosomal_large_subunit_biogenesis,
                 text_size = 60, data_size = 6)
plotly::layout(p_pca, title = "ribosomal_large_subunit_biogenesis activity")


grid.arrange(p2,p3,ncol=2)
labs$lab

if(F){
#4. Randomise variables to measure goodness of observed fit
# use permutations within each dimension - this is only possible for less than 8 vertices because computing convex hull gets exponentially slower with more dimensions
start = Sys.time()
pch_rand = randomise_fit_pch(PCs4arch, arc_data = arc_1,
                             n_rand = 100,
                             replace = FALSE, bootstrap_N = NA,
                             volume_ratio = "t_ratio",
                             maxiter = 500, delta = 0, conv_crit = 1e-4,
                             type = "m", clust_options = list(cores = 3))
# use type m to run on a single machine or cloud
# type = "m", clust_options = list(cores = 3))
# use clustermq (type cmq) to run as jobs on a computing cluster (higher parallelisation)
# type = "cmq", clust_options = list(njobs = 10)) 

# This analysis took:
Sys.time() - start #2min
# plot background distribution of t-ratio and show p-value
#plot(pch_rand, type = c("t_ratio"), nudge_y = 5)

pch_rand
summary(pch_rand$rand_dist$obs_vs_rand)
}


###########################################################################################
#https://github.com/U54Bioinformatics/03A_scRNA_Archetype_Mutitasklearning_Analysis
#5. Use archetype scores to train a multi-task model (group lasso penalty) with hallmark pahtway enrichment scores or gene expression as predictors. The coefficients for the pathways can be used to identify core phenotypes associated with each archetype.
if(F){
library(glmnet)

X <- as.matrix(ssgsea.scores) # ssgsea.scores is matrix of Hallmark ssGSEA pathway enrichment scores calculated using GSVA. Cell IDs are in rownames and pathways in colnames

X <- X[intersect(rownames(ssgsea.scores), colnames(arc_data_t$S)), ]
Y <- t(arc_data_t$S[, intersect(rownames(ssgsea.scores), colnames(arc_data_t$S))])

# Perform cross-validation analysis and fit a group-lasso penalized multitask model

cv.fit1 = cv.glmnet(X, Y, family = "mgaussian", alpha=1)
l1se <- match(cv.fit1$lambda.1se, cv.fit1$lambda) #selects a lambda penalty parameter  

write.csv(data.frame("A1"=cv.fit1$glmnet.fit$beta$y1[, l1se], 
                     "A2"=cv.fit1$glmnet.fit$beta$y2[, l1se],
                     "A3"=cv.fit1$glmnet.fit$beta$y3[, l1se]),
          file="./Integrated_10PCs_mgauss_coefs.csv") # save the coefficients from the cross-validation analysis at the selected lambda

}


