#https://bioconductor.org/packages/devel/bioc/vignettes/CellMixS/inst/doc/CellMixS.html
library(CellMixS)

# Load required packages
suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(cowplot)
  library(limma)
  library(magrittr)
  library(dplyr)
  library(purrr)
  library(ggplot2)
  library(scater)
})

# 1) Load sim_list example data
sim_list <- readRDS(system.file(file.path("extdata", "sim50.rds"), 
                                package = "CellMixS"))
names(sim_list)
#"batch0"  "batch20" "batch50"

sce50 <- sim_list[["batch50"]]
class(sce50)
#[1] "SingleCellExperiment"
#attr(,"package")
#[1] "SingleCellExperiment"

table(sce50[["batch"]])

# 2） visualize batch effect
visGroup(sce50, group = "batch")
# Visualize batch distribution in other elements of sim_list 
batch_names <- c("batch0", "batch20")

vis_batch <- lapply(batch_names, function(name){
  sce <- sim_list[[name]]
  visGroup(sce, "batch") + ggtitle(paste0("sim_", name))
})

plot_grid(plotlist = vis_batch, ncol = 2)

# 3) quantify batch effects
# Call cell-specific mixing score (cms) for sce50
# Note that cell_min is set to 4 due to the low number of cells and small k.
# Usually default parameters are recommended!! 
sce50 <- cms(sce50, k = 30, group = "batch", res_name = "unaligned", 
             n_dim = 3, cell_min = 4)
head(colData(sce50))

hist(sce50$cms_smooth.unaligned)
hist(sce50$cms.unaligned)

# Call cell-specific mixing score for all datasets
sim_list <- lapply(batch_names, function(name){
  sce <- sim_list[[name]]
  sce <- cms(sce, k = 30, group = "batch", res_name = "unaligned", 
             n_dim = 3, cell_min = 4)
}) %>% set_names(batch_names)

# Append cms50
sim_list[["batch50"]] <- sce50

# 4) visualize the cell mixing score
visHist(sce50)
# p-value histogram sim30
# Combine cms results in one matrix
batch_names <- names(sim_list)
cms_mat <- batch_names %>% 
  map(function(name) sim_list[[name]]$cms.unaligned) %>% 
  bind_cols() %>% set_colnames(batch_names)
#> New names:
#> * NA -> ...1
#> * NA -> ...2
#> * NA -> ...3

visHist(cms_mat, n_col = 3)

# cms only cms20
sce20 <- sim_list[["batch20"]]
metric_plot <- visMetric(sce20, metric_var = "cms_smooth.unaligned")

# group only cms20
group_plot <- visGroup(sce20, group = "batch")

plot_grid(metric_plot, group_plot, ncol = 2)

# 5) Evaluate data integration
# MNN - embeddings are stored in the reducedDims slot of sce
reducedDimNames(sce20)
#> [1] "TSNE" "PCA"  "MNN"
sce20 <- cms(sce20, k = 30, group = "batch", 
             dim_red = "MNN", res_name = "MNN", n_dim = 3, cell_min = 4)

# Run limma
sce20 <- scater::logNormCounts(sce20)
limma_corrected <- removeBatchEffect(logcounts(sce20), batch = sce20$batch)
# Add corrected counts to sce
assay(sce20, "lim_corrected") <- limma_corrected 
# Run cms
sce20 <- cms(sce20, k = 30, group = "batch", 
             assay_name = "lim_corrected", res_name = "limma", n_dim = 3, 
             cell_min = 4)

names(colData(sce20))

# Compare data integration methods
visHist(sce20, metric = "cms.",  n_col = 3)
