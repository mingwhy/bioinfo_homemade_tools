# analyze all data together, do this on server
library(Seurat)
library(ggplot2);library(gridExtra);
options(stringsAsFactors = F)

## read in processed wholebrain data
file="../data/wholebrain_filtered.rds";
dat0=readRDS(file); 
dat0; #12616 features across 100527 samples

# for each cell type, the count or proportion of male, female and mix cell numbers
# remove 'mix' cell types
table(dat0@meta.data$sex)
#female   male    mix 
#49105  47409   4013 
dat=subset(dat0,sex!='mix')
dat # 12616 features across 96514 samples within 1 assay 

table(dat$batch) #0~12, 13 batches in total
dat.test=dat;
  
  dat.test <- NormalizeData(dat.test, normalization.method = "LogNormalize", scale.factor = 10000)
  # for PCA, PCA for UMAP
  dat.test <- FindVariableFeatures(dat.test, selection.method = "vst", nfeatures = 2000)
  all.genes <- rownames(dat.test)
  dat.test <- ScaleData(dat.test, features = all.genes)
  dat.test <- RunPCA(dat.test, features = VariableFeatures(object = dat.test))
  dat.test <- RunUMAP(dat.test, dims = 1:10)
  
  #DimPlot(dat.test, reduction = "umap") #just plot all dots

if(F){
  pdf('test.pdf',useDingbats=T)
  # color use sex labels
  dat.test@meta.data$orig.ident=dat.test@meta.data$sex
  print(DimPlot(object = dat.test, reduction = 'umap', 
                group.by = 'orig.ident', 
                cols = c('red','blue'))+ggtitle('color by sex')
  )
  # color use batch labels  
  dat.test@meta.data$orig.ident=dat.test@meta.data$batch
  print(DimPlot(object = dat.test, reduction = 'umap', 
                group.by = 'orig.ident')+ggtitle('color by batch')
  )
  dev.off()

# looks like there is no obvious batch effect
# use cmx store to quantify batch effect per cell
# Seurat -> SingleCellExperiment object
#https://satijalab.org/seurat/archive/v3.1/conversion_vignette.html
# BiocManager::install("almutlue/CellMixS")
# BiocManager::install("SingleCellExperiment")
library(SingleCellExperiment)
library(CellMixS)
dat.sce <- as.SingleCellExperiment(dat.test)

pdf('test_cms.pdf',useDingbats=T)
print(visGroup(dat.sce, group = "batch") )

dat.sce <- cms(dat.sce, k = 30, group = "batch", res_name = "unaligned");
             #n_dim = 3, cell_min = 4)
head(colData(dat.sce))
# p-value histogram of dat.sce
print(visHist(dat.sce))
dev.off()

}

## use limma to correct batch effect which returns a corrected data matrix
# Run limma, use dat.sce, a SingleCellExperiment object
library(SingleCellExperiment)
library(CellMixS)
dat.sce <- as.SingleCellExperiment(dat.test)
table(dat.sce$batch)
library(limma)
dat.sce <- scater::logNormCounts(dat.sce)
limma_corrected <- removeBatchEffect(logcounts(dat.sce), batch = dat.sce$batch)
# Add corrected counts to sce
assay(dat.sce, "lim_corrected") <- limma_corrected 

# Run cms
dat.sce <- cms(dat.sce, k = 30, group = "batch", 
             assay_name = "lim_corrected", res_name = "limma");

names(colData(dat.sce))

pdf('test_cms_limma.pdf',useDingbats=T)
# As pvalue histograms
par(mfrow=c(2,1))
hist(dat.sce$cms_smooth.limma)
hist(dat.sce$cms.limma)
#print(visHist(dat.sce, metric = "cms.",  n_col = 1))
dev.off();

# based on the plot, batch correction didn't help that much

if(F){
  ## Harnomy outputs genex new.embedding coords, which is not what I need
  #https://github.com/immunogenomics/harmony
  library(harmony)
  dat.test <- RunHarmony(dat.test, "batch")
  #object,Pipeline object. Must have PCA computed.
  #group.by.vars, Which variable(s) to remove (character vector).
  dat.test <- RunUMAP(dat.test, reduction = "harmony")
  #https://portals.broadinstitute.org/harmony/SeuratV3.html
  harmony_embeddings <- Embeddings(dat.test, 'harmony')
  harmony_embeddings[1:5, 1:5]

  options(repr.plot.height = 5, repr.plot.width = 12)
  p1 <- DimPlot(object = dat.test, reduction = "harmony", pt.size = .1, group.by = "batch", do.return = TRUE)
  p2 <- VlnPlot(object = dat.test, features = "harmony_1", group.by = "batch", do.return = TRUE, pt.size = .1)
  plot_grid(p1,p2)
}
