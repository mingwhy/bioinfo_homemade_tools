########################################################################################################
## install: set up both python and R env
#https://github.com/vitkl/ParetoTI
if(F){
  # You can also install these python modules directly in terminal,
  # To help ParetoTI find the right python create conda environment named "reticulate_PCHA" (uncomment):
  # (other enviroment containing py_pcha will be used if reticulate_PCHA doesn't exist)    
  conda create -n reticulate_PCHA python=3.7.3 pip    
  # Light install:    
  source activate reticulate_PCHA && pip install --upgrade py_pcha numpy scipy datetime geosketch umap-learn    
  # To use more features:    
  source activate reticulate_PCHA && pip install --upgrade py_pcha numpy scipy datetime tensorflow tensorflow-probability pandas keras h5py geosketch pydot sklearn umap-learn    
  
  # Sometimes on some platforms R sees only the "base" conda enviroment 
  # (like when RStudio Server is setup incorrectly)    
  # In that case use:
  source activate base && pip install --upgrade py_pcha numpy scipy datetime geosketch umap-learn
}

## install in Rstudio
#BiocManager::install("vitkl/ParetoTI", dependencies = c("Depends", "Imports", "LinkingTo"))

library('ParetoTI')

## Finally, check that py_pcha library is successfully installed and discoverable
reticulate::py_discover_config("py_pcha")

# To make sure R uses the correct conda enviroment you can run this when you start R:
reticulate::use_condaenv("reticulate_PCHA", conda = "auto",
                         required = TRUE) # set TRUE to force R to use reticulate_PCHA

library(ggplot2)
library(cowplot)
########################################################################################################
## https://vitkl.github.io/ParetoTI/articles/Comparison_to_kmeans.html#find-archetypes-with-pcha-and-cluster-centers-with-k-means
# 1, Comparison of achetypal analysis and k-means clustering

# set random seed
set.seed(4355)
# generate archetype positions
archetypes = generate_arc(arc_coord = list(c(5, 0, 4), c(-10, 15, 0), c(-30, -20, -5)),
                          mean = 0, sd = 1)
archetypes$XC #3X3

# generate 1000 data points that are a convex combination (weighted sum) of archetypes
data = generate_data(archetypes$XC, N_examples = 1e3, jiiter = 0.04, size = 0.99)
colnames(data) = paste0("cell", seq_len(ncol(data)))
dim(data) #3,1000
data[1:3,1:10]

# plot
plot_arc(arc_data = archetypes, data = data,which_dimensions = 1:2) + ylim(-19, 17) +
  ggtitle("Ground truth archetypes")
plot_arc(arc_data = archetypes, data = data,which_dimensions = 2:3) + ylim(-19, 17) +
  ggtitle("Ground truth archetypes")

# 2, Find archetypes with PCHA and cluster centers with k-means and Louvain methods
# find archetypes
arc = fit_pch(data, noc = 3)
# find k-means clusters
clusters = fit_pch(data, noc = 3, method = "kmeans")
# find Louvain clusters
lou_clusters = fit_pch(data, noc = 3, method = "louvain") 

plot_grid(plot_arc(arc_data = arc, data = data,
                   which_dimensions = 1:2) +
            ylim(-18, 17) + ggtitle("Detected archetypes (PCHA)"),
          plot_arc(arc_data = clusters, data = data,
                   which_dimensions = 1:2,
                   data_lab = as.character(apply(clusters$S, 2, which.max))) +
            ylim(-18, 17) + ggtitle("K-means clusters"),
          plot_arc(arc_data = lou_clusters, data = data,
                   which_dimensions = 1:2,
                   data_lab = as.character(apply(lou_clusters$S, 2, which.max))) +
            ylim(-18, 17) + ggtitle("Louvain clusters"),
          align = "vh", nrow = 1)

# 3, Use bootstrapping to find variability in positions of archetypes and cluster centers
# bootstrap archetypes
arc_rob = fit_pch_bootstrap(data, n = 200, sample_prop = 0.65, seed = 2543,
                            noc = 3)
# bootstrap kmeans
clusters_rob = fit_pch_bootstrap(data, n = 200, sample_prop = 0.65, seed = 2543,
                                 noc = 3, method = "kmeans")
# bootstrap Louvain 
# specific noc is achieved by optimising clustering resolution to get noc clusters
lou_clusters_rob = fit_pch_bootstrap(data, n = 200, sample_prop = 0.9, seed = 2543,
                                     noc = 3, method = "louvain",
                                     method_options = list(resolution = 0.1))

# show both side-by-side
plot_grid(plot_arc(arc_data = arc_rob, data = data,
                   which_dimensions = 1:2) +
            ylim(-18, 17) + ggtitle("Detected archetypes (PCHA)"),
          plot_arc(arc_data = clusters_rob, data = data,
                   which_dimensions = 1:2,
                   data_lab = as.character(apply(clusters$S, 2, which.max))) +
            ylim(-18, 17) + ggtitle("K-means clusters"),
          plot_arc(arc_data = lou_clusters_rob, data = data,
                   which_dimensions = 1:2,
                   data_lab = as.character(apply(lou_clusters$S, 2, which.max))) +
            ylim(-18, 17) + ggtitle("Louvain clusters"),
          align = "vh")

# 4,Too many archetypes and cluster centers => increased variability in positions?
# trying different number of archetypes
arc_ks = k_fit_pch(data, ks = 2:5,
                   bootstrap = T, bootstrap_N = 200, maxiter = 500,
                   bootstrap_type = "m", clust_options = list(cores = 3),
                   seed = 2543, replace = FALSE,
                   volume_ratio = "none", # set to "none" if too slow
                   order_type = "align", sample_prop = 0.65, reference = T)
# trying different number of clusters
cluster_ks = k_fit_pch(data, ks = 2:5,
                       bootstrap = T, bootstrap_N = 200, maxiter = 500,
                       bootstrap_type = "m", clust_options = list(cores = 3),
                       seed = 2543, replace = FALSE,
                       volume_ratio = "none", # set to "none" if too slow
                       order_type = "align", sample_prop = 0.65, reference = T, method = "kmeans")

lou_cluster_ks = k_fit_pch(data, ks = 2:5,
                           bootstrap = T, bootstrap_N = 200, maxiter = 500,
                           bootstrap_type = "m", clust_options = list(cores = 3),
                           seed = 2543, replace = FALSE,
                           volume_ratio = "none", # set to "none" if too slow
                           sample_prop = 0.95, method = "louvain",
                           method_options = list(resolution = 0.1,
                                                 noc_optim_iter = 500)) # try resolutions for more iterations

# Show variance explained by k-vertex model on top of k-1 model (each k separately)
plot_grid(plot_arc_var(arc_ks, type = "res_varexpl",
                       point_size = 2, line_size = 1.5) +
            theme_bw() + ggtitle("Archetypes (PCHA)"),
          plot_arc_var(cluster_ks, type = "res_varexpl",
                       point_size = 2, line_size = 1.5) +
            theme_bw() + ggtitle("K-means clusters"),
          plot_arc_var(lou_cluster_ks, type = "res_varexpl",
                       point_size = 2, line_size = 1.5) +
            theme_bw() + ggtitle("Louvain clusters"),
          align = "vh")

# Show variance in position of vertices obtained using bootstraping
# - use this to find largest k that has low variance
plot_grid(plot_arc_var(arc_ks, type = "total_var",
                       point_size = 2, line_size = 1.5) +
            theme_bw() + ylab("Mean variance in position of vertices") +
            ggtitle("Archetypes (PCHA)"),
          plot_arc_var(cluster_ks, type = "total_var",
                       point_size = 2, line_size = 1.5) +
            theme_bw() + ylab("Mean variance in position of vertices") +
            ggtitle("K-means clusters"),
          plot_arc_var(lou_cluster_ks, type = "total_var",
                       point_size = 2, line_size = 1.5) +
            theme_bw() + ggtitle("Louvain clusters"),
          align = "vh")

#############################################################################
#Example of using archetypal analysis to find representative cells & describe 
#heterogeniety in hepatocyte population between those archetypes
#https://vitkl.github.io/ParetoTI/articles/Hepatocyte_example.html

library(ParetoTI)
library(SingleCellExperiment)
library(ggplot2)
library(cowplot)
library(Matrix)
#1. Load data from GEO and filter as described in the paper, normalise and PCs for finding polytopes
# uncomment to load data -------------------------------------------------------
#BiocManager::install("GEOquery")
if(F){
gse = GEOquery::getGEO("GSE84498", GSEMatrix = TRUE)#dir.create first before run the line below
filePaths = GEOquery::getGEOSuppFiles("GSE84498", fetch_files = T, baseDir = "./processed_data/")
}

filePaths = c("./processed_data/GSE84498/GSE84498_experimental_design.txt.gz",
              "./processed_data/GSE84498/GSE84498_umitab.txt.gz")
design = fread(filePaths[1], stringsAsFactors = F)
data = fread(filePaths[2], stringsAsFactors = F, header = T)

data = as.matrix(data, rownames = "gene")
dim(data) #27389gene x 1736cell
data[1:3,1:3]

# convert to single cell experiment
hepatocytes = SingleCellExperiment(assays = list(counts = data),
                                   colData = design)

# look at mitochondrial-encoded MT genes
mito.genes = grep(pattern = "^mt-",
                  x = rownames(data), 
                  value = TRUE)
hepatocytes$perc.mito = colSums(counts(hepatocytes[mito.genes, ])) / colSums(counts(hepatocytes))
#qplot(hepatocytes$perc.mito, geom = "histogram")

# look at nuclear-encoded MT genes (find those genes using GO annotations)
go_annot = map_go_annot(taxonomy_id = 10090, keys = rownames(hepatocytes),
                        columns = c("GOALL"), keytype = "ALIAS",
                        ontology_type = c("CC"))

mitochondria_located_genes = unique(go_annot$annot_dt[GOALL == "GO:0005739", ALIAS])

hepatocytes$all_mito_genes = colSums(counts(hepatocytes[mitochondria_located_genes, ])) / colSums(counts(hepatocytes))
#qplot(hepatocytes$perc.mito, hepatocytes$all_mito_genes, geom = "bin2d")

## Filtering
# remove batches of different cells (probably non-hepatocytes)
hepatocytes = hepatocytes[, !hepatocytes$batch %in% c("AB630", "AB631")]

# remove cells with more less than 1000 or more than 30000 UMI
hepatocytes = hepatocytes[, colSums(counts(hepatocytes)) > 1000 &
                            colSums(counts(hepatocytes)) < 30000]
# remove cells that express less than 1% of albumine
alb_perc = counts(hepatocytes)["Alb",] / colSums(counts(hepatocytes))
hepatocytes = hepatocytes[, alb_perc > 0.01]
# remove genes with too many zeros (> 95% cells)
hepatocytes = hepatocytes[rowMeans(counts(hepatocytes) > 0) > 0.05,]
# remove cells with too many zeros (> 85%)
hepatocytes = hepatocytes[,colMeans(counts(hepatocytes) > 0) > 0.15]

hepatocytes #7269 1224

# Normalise gene expression by cell sum factors and log-transform
hepatocytes = scran::computeSumFactors(hepatocytes)#https://github.com/LTLA/scuttle/issues/12
#hepatocytes = scater::normalize(hepatocytes)
#hepatocytes = scater::normalize(hepatocytes, return_log = FALSE) # just normalise
hepatocytes = scater::logNormCounts(hepatocytes,transform='none')
#assays(2): counts normcounts
hepatocytes = scuttle::logNormCounts(hepatocytes) #default: log2. https://rdrr.io/github/LTLA/scuttle/man/logNormCounts.html
hepatocytes #assays(3): counts logcounts normcounts

#Plot below shows first 3 PCs colored by batch.
# Find principal components
hepatocytes = scater::runPCA(hepatocytes, ncomponents = 7,scale = T, exprs_values = "logcounts")
#hepatocytes = scater::runPCA(hepatocytes, ncomponents = 7,scale = T, exprs_values = "normcounts")
# Plot PCA colored by batch
scater::plotReducedDim(hepatocytes, ncomponents = 3, dimred = "PCA",
                       colour_by = "batch")
hepatocytes
#reducedDimNames(1): PCA

# extract PCs (centered at 0 with runPCA())
PCs4arch = t(reducedDim(hepatocytes, "PCA"))
hepatocytes # 7269 1224 
dim(PCs4arch) #7PC X 1224 cells

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

#Examine the polytope with best k & look at known markers of subpopulations
# fit a polytope with bootstraping of cells to see stability of positions
arc = fit_pch_bootstrap(PCs4arch, n = 200, sample_prop = 0.75, seed = 235,
                        noc = 4, delta = 0, conv_crit = 1e-04, type = "m")

p_pca = plot_arc(arc_data = arc, data = PCs4arch, 
                 which_dimensions = 1:3, line_size = 1.5,
                 data_lab = as.numeric(logcounts(hepatocytes["Alb",])),
                 text_size = 60, data_size = 6) 
plotly::layout(p_pca, title = "Hepatocytes colored by Alb (Albumine)")

p_pca = plot_arc(arc_data = arc, data = PCs4arch, 
                 which_dimensions = 1:3, line_size = 1.5,
                 data_lab = as.numeric(logcounts(hepatocytes["Cyp2e1",])),
                 text_size = 60, data_size = 6) 
plotly::layout(p_pca, title = "Hepatocytes colored by Cyp2e1")

p_pca = plot_arc(arc_data = arc, data = PCs4arch, 
                 which_dimensions = 1:3, line_size = 1.5,
                 data_lab = as.numeric(logcounts(hepatocytes["Gpx1",])),
                 text_size = 60, data_size = 6) 
plotly::layout(p_pca, title = "Hepatocytes colored by Gpx1")

# You can also check which cells have high entropy of logistic regression predictions when classifying all cells in a tissue into cell types. These could have been misclassified by the method and wrongly assigned to Hepatocytes, or these could be doublets.

# find archetypes on all data (allows using archetype weights to describe cells)
arc_1 = fit_pch(PCs4arch, volume_ratio = "t_ratio", maxiter = 500,
                noc = 4, delta = 0,
                conv_crit = 1e-04)

# check that positions are similar to bootstrapping average from above
p_pca = plot_arc(arc_data = arc_1, data = PCs4arch, 
                 which_dimensions = 1:3, line_size = 1.5, 
                 data_lab = as.numeric(logcounts(hepatocytes["Alb",])),
                 text_size = 60, data_size = 6) 
plotly::layout(p_pca, title = "Hepatocytes colored by Alb")

#Find genes and gene sets enriched near vertices
# Map GO annotations and measure activities
activ = measure_activity(hepatocytes, # row names are assumed to be gene identifiers
                         which = "BP", return_as_matrix = F,
                         taxonomy_id = 10090, keytype = "ALIAS",
                         lower = 20, upper = 1000,
                         aucell_options = list(aucMaxRank = nrow(hepatocytes) * 0.1,
                                               binary = F, nCores = 3,
                                               plotStats = FALSE))
dim(activ) #1224, 2729
activ[1:3,1:3]

# Merge distances, gene expression and gene set activity into one matrix
data_attr = merge_arch_dist(arc_data = arc_1, data = PCs4arch, 
                            feature_data = as.matrix(logcounts(hepatocytes)),
                            colData = activ,
                            dist_metric = c("euclidean", "arch_weights")[1],
                            colData_id = "cells", rank = F) 

# Use Wilcox test to find genes maximally expressed in 10% closest to each vertex
enriched_genes = find_decreasing_wilcox(data_attr$data, data_attr$arc_col,
                                        features = data_attr$features_col,
                                        bin_prop = 0.1, method = "BioQC")
dim(enriched_genes) #29076 x 8

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

p_pca = plot_arc(arc_data = arc, data = PCs4arch,
                 which_dimensions = 1:3, line_size = 1.5,
                 data_lab = activ$ribosomal_large_subunit_biogenesis,
                 text_size = 60, data_size = 6)
plotly::layout(p_pca, title = "ribosomal_large_subunit_biogenesis activity")

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
plot(pch_rand, type = c("t_ratio"), nudge_y = 5)

pch_rand
summary(pch_rand$rand_dist$obs_vs_rand)





