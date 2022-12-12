
# elife brain data download from: https://cells.ucsc.edu/?gene=ENSG00000132854&org=Fruit+fly+(D.+melanogaster)&bp=brain&ds=dros-brain
# some help info: https://cellbrowser.readthedocs.io/en/master/load.html

#https://bioconductor.org/packages/devel/bioc/vignettes/SingleCellExperiment/inst/doc/intro.html
#http://bioconductor.org/packages/release/bioc/vignettes/scuttle/inst/doc/qc.html
library(Matrix);library(Seurat)
library(ggplot2);library(gridExtra);
library(patchwork);library(dplyr)
library(SingleCellExperiment)
library(SummarizedExperiment)
library(scuttle);
library(tidyverse);
library(scater);library(scran)
#options(future.globals.maxSize = 3.2 * 1024^3)

###########################################################################
## processing elife larva brain with cell type annotation data
if(!file.exists('larva.brain_cell.type_sce.rds')){
  #merged.mat <- data.table::fread("ucsc_cellBrowser/merged_normal_starved/exprMatrix.tsv.gz")
  #merged.meta <- data.frame(data.table::fread("ucsc_cellBrowser/merged_normal_starved/meta.tsv"), row.names=1)
  #umap.df=data.table::fread("ucsc_cellBrowser/merged_normal_starved/Seurat_umap.coords.tsv.gz")
  
  merged.mat <- data.table::fread("https://cells.ucsc.edu/dros-brain/merge/exprMatrix.tsv.gz")
  merged.meta <- data.frame(data.table::fread("https://cells.ucsc.edu/dros-brain/merge/meta.tsv"), row.names=1)
  umap.df=data.table::fread("https://cells.ucsc.edu/dros-brain/merge/Seurat_umap.coords.tsv.gz")
  dim(merged.mat);dim(merged.meta)
  head(merged.meta)
  sum(rownames(merged.meta)==colnames(merged.mat)[-1])
  
  genes = merged.mat[,1][[1]]
  genes = gsub(".+[|]", "", genes)
  expr.mat = as.matrix(merged.mat[,-1]) #Values in matrix are: Seurat normalized (counts per 10,000)
  rownames(expr.mat)=genes;
  
  # add normal cell labels
  x=stringr::str_extract(rownames(merged.meta),'_\\d')
  table(x)
  # _1   _2 
  #4349 4347 
  # _1: normal
  # _2: starved
  merged.meta$condition='starved';
  merged.meta[x=='_1',]$condition='normal'
  merged.meta$cell_barcode=rownames(merged.meta)
  
  
  dim(umap.df)
  sum(umap.df$V1==rownames(merged.meta))
  umap.df$condition=merged.meta$condition;
  
  elife.plot<-ggplot(umap.df,aes(x=V2,y=V3,col=condition))+geom_point(pch=16,size=0.1)+theme_classic()+
    scale_color_manual(values=c('darkgreen','pink'))
  
  sum(colnames(expr.mat)==rownames(merged.meta))
  sce <- SingleCellExperiment(list(counts=expr.mat),
                              colData=merged.meta)
  sce
  table(sce$condition)
  saveRDS(sce,'larva.brain_cell.type_sce.rds')
}

###########################################################################
## create ref with elife brain
library(Seurat)
library(symphony)
library(ggplot2);library(gridExtra)
library(zellkonverter)
library(SummarizedExperiment)
library(SingleCellExperiment)

# Other packages for this tutorial
suppressPackageStartupMessages({
  # Analysis
  library(harmony)
  library(irlba)
  library(data.table)
  library(dplyr)
  # Plotting
  library(ggplot2)
  library(ggthemes)
  library(ggrastr)
  library(RColorBrewer)
})
plotBasic = function(umap_labels,                # metadata, with UMAP labels in UMAP1 and UMAP2 slots
                     title = 'Query',         # Plot title
                     color.by = 'annotation',  # metadata column name for coloring
                     facet.by = NULL,         # (optional) metadata column name for faceting
                     color.mapping = NULL,    # custom color mapping
                     legend.position = 'right') {  # Show cell type legend
  umap_labels=as.data.frame(umap_labels)
  p = umap_labels %>%
    #dplyr::sample_frac(1L) %>% # permute rows randomly
    ggplot(aes(x = UMAP1, y = UMAP2)) + 
    geom_point_rast(aes(col = get(color.by)), size = 1, stroke = 0.4, shape = 16)
  if (!is.null(color.mapping)) { p = p + scale_color_manual(values = color.mapping) }
  
  # Default formatting
  p = p + theme_bw() +
    labs(title = title, color = color.by) + 
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(legend.position=legend.position) +
    theme(legend.text = element_text(size=8), legend.title=element_text(size=12)) + 
    guides(colour = guide_legend(override.aes = list(size = 4))) + guides(alpha = 'none')
  
  if(!is.null(facet.by)) {
    p = p + facet_wrap(~get(facet.by)) +
      theme(strip.text.x = element_text(size = 12)) }    
  return(p)
}


#################################################################################
# https://rdrr.io/cran/symphony/f/inst/doc/quickstart_tutorial.Rmd
# build reference 
if(!file.exists('larva.brain_cell.type_ref.rds')){

  sce=readRDS('larva.brain_cell.type_sce.rds')
  ref_metadata=colData(sce)
  
  # build ref from harmony object
  assayNames(sce)
  #ref_exp_full=assay(sce,'raw_counts')
  ref_exp_full=assay(sce,'counts')
  ref_exp_full <- as(ref_exp_full, "sparseMatrix")
  class(ref_exp_full)
  
  #Select variable genes and subset reference expression by variable genes (the command below will select the top 1,000 genes per batch, then pool them)
  var_genes = vargenes_vst(ref_exp_full,topn = 2000)
  ref_exp = ref_exp_full[var_genes, ]
  dim(ref_exp)
  
  #Calculate and save the mean and standard deviations for each gene
  vargenes_means_sds = tibble(symbol = var_genes, mean = Matrix::rowMeans(ref_exp))
  vargenes_means_sds$stddev = rowSDs(ref_exp, vargenes_means_sds$mean)
  head(vargenes_means_sds)
  
  #Scale data using calculated gene means and standard deviations
  ref_exp_scaled = scaleDataWithStats(ref_exp, vargenes_means_sds$mean, vargenes_means_sds$stddev, 1)
  
  #Run PCA (using SVD), save gene loadings (s$u)
  set.seed(0)
  s = irlba(ref_exp_scaled, nv = 100)
  Z_pca_ref = diag(s$d) %*% t(s$v) # [pcs by cells]
  loadings = s$u
  
  #Run Harmony integration. It is important to set return_object = TRUE.
  head(ref_metadata)
  ref_metadata$fake=1
  
  set.seed(0)
  ref_harmObj = harmony::HarmonyMatrix(
    data_mat = t(Z_pca_ref),  ## PCA embedding matrix of cells
    meta_data = ref_metadata, ## dataframe with cell labels
    theta = c(2),             ## cluster diversity enforcement
    vars_use = c('fake'),    ## variable to integrate out
    nclust = 30,             ## number of clusters in Harmony model
    max.iter.harmony = 20,    ## max number of iterations
    return_object = TRUE,     ## return the full Harmony model object
    do_pca = FALSE            ## don't recompute PCs
  )
  
  #To run the next function buildReferenceFromHarmonyObj(), you need to input the saved gene loadings (loadings) and vargenes_means_sds.
  # Compress a Harmony object into a Symphony reference
  reference = buildReferenceFromHarmonyObj(
    ref_harmObj,            # output object from HarmonyMatrix()
    ref_metadata,           # reference cell metadata
    vargenes_means_sds,     # gene names, means, and std devs for scaling
    loadings,               # genes x PCs matrix
    verbose = TRUE,         # verbose output
    do_umap = TRUE,         # set to TRUE to run UMAP
    save_uwot_path = './larva.brain_cell.type_ref_uwot_model') # file path to save uwot model
  
  #Save Symphony reference for later mapping (modify with your desired output path)
  saveRDS(reference, 'larva.brain_cell.type_ref.rds')
}


#################################################################################
## plot reference
reference=readRDS('larva.brain_cell.type_ref.rds')
sce=readRDS('larva.brain_cell.type_sce.rds')
ref_metadata=colData(sce)
str(reference)
#The harmonized embedding is located in the Z_corr slot of the reference object.
dim(reference$Z_corr) # PC by cell matrix

#Visualize reference UMAP
umap_labels = cbind(ref_metadata, reference$umap$embedding)
umap_labels=as.data.frame(umap_labels)

ggplot(umap_labels,aes(x=UMAP1,y=UMAP2,col=condition))+
  geom_jitter(size=0.5)+theme_classic()+
  #scale_color_manual(values=my_colors)+
  guides(colour = guide_legend(override.aes = list(size=5)))


#################################################################################
## map query
reference=readRDS('larva.brain_cell.type_ref.rds')
four.samples=readRDS('../RNA3-042-nova_data/four.samples.filterCell.sceObjs.rds')
names(four.samples)

if(!file.exists('map.to.LarvaBrainCellType_out.rds')){
  query.out=list();
  for(i in 1:4){
    #i=1
    one.sample=four.samples[[i]]
    one.sample.name=names(four.samples)[[i]]
    assayNames(one.sample)
    #query_exp=assay(one.sample,'counts')
    query_exp=assay(one.sample,"normcounts")
    query_metadata=colData(one.sample)
    rownames(query_exp) #fly gene symbols
    
    query = mapQuery(query_exp,             # query gene expression (genes x cells)
                     query_metadata,        # query metadata (cells x attributes)
                     reference,             # Symphony reference object
                     do_normalize = FALSE,  # perform log(CP10k) normalization on query
                     do_umap = TRUE)        # project query cells into reference UMAP
    
    str(query)
    query = knnPredict(query,       # query object
                       reference,   # reference object
                       #reference$meta_data$condition, # reference cell labels for training
                       reference$meta_data$Cluster,
                       k = 10,       # number of reference neighbors to use for prediction
                       confidence = TRUE)
    query.out[[i]]<-query
  }
  names(query.out)<-names(four.samples)
  saveRDS(query.out,'map.to.LarvaBrainCellType_out.rds')
}

#####################################################################
query.out=readRDS('map.to.LarvaBrainCellType_out.rds')

query.result<-lapply(names(query.out),function(i){
  x=query.out[[i]]
  df=x$meta_data
  #summary(df$cell_type_pred_knn_prob)
  #df=df[df$cell_type_pred_knn_prob>0.8,]
  df$condition=i
  df
})
df=as.data.frame(Reduce(`rbind`,query.result))
df$condition=factor(df$condition,levels=names(query.out))
df1=df %>% dplyr::count(condition,cell_type_pred_knn)
df2=df1 %>%  tidyr::spread(key=cell_type_pred_knn,value=n)
df2
dim(df) #31772    10
df.elife=df

#####################################################################
## compare with previously mapped results
query.out=readRDS('../RNA3-042-nova_data/map.to.L1_out.rds')
query.result<-lapply(names(query.out),function(i){
  x=query.out[[i]]
  df=x$meta_data
  #df=df[df$cell_type_pred_knn_prob>0.8,]
  df$condition=i
  df
})
df=as.data.frame(Reduce(`rbind`,query.result))
df$condition=factor(df$condition,levels=names(query.out))
df1=df %>% dplyr::count(condition,cell_type_pred_knn)
df2=df1 %>%  tidyr::spread(key=cell_type_pred_knn,value=n)
df2
dim(df) #31772    10
df.L1=df

sum(rownames(df.L1)==rownames(df.elife)) #31772
mat=table(df.L1$cell_type_pred_knn,df.elife$cell_type_pred_knn)
library(ComplexHeatmap);
max(mat)
#https://bookdown.org/ndphillips/YaRrr/more-colors.html#colorramp2
col_fun = circlize::colorRamp2(c(0, 4910, 2), c("blue", "green", "red"),transparency = .3)

pdf("compare_ref_elife_L1.pdf",useDingbats = T,width = 12)
elife.plot+geom_point(size=0.3)+guides(color = guide_legend(override.aes = list(size = 8)),
                  shape = guide_legend(override.aes = list(size = 8)))
Heatmap(mat, name = "Tabulate mapped cell annotation between two refs", col = col_fun,
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(sprintf("%.0f", mat[i, j]), x, y, gp = gpar(fontsize = 10))
        })
dev.off()


