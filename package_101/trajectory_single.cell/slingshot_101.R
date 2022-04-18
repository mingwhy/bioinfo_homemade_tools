
library(SummarizedExperiment)
library(SingleCellExperiment)

# https://bioconductor.org/packages/devel/bioc/vignettes/slingshot/inst/doc/vignette.html
# generate synthetic count data representing a single lineage
means <- rbind(
  # non-DE genes
  matrix(rep(rep(c(0.1,0.5,1,2,3), each = 300),100),
         ncol = 300, byrow = TRUE),
  # early deactivation
  matrix(rep(exp(atan( ((300:1)-200)/50 )),50), ncol = 300, byrow = TRUE),
  # late deactivation
  matrix(rep(exp(atan( ((300:1)-100)/50 )),50), ncol = 300, byrow = TRUE),
  # early activation
  matrix(rep(exp(atan( ((1:300)-100)/50 )),50), ncol = 300, byrow = TRUE),
  # late activation
  matrix(rep(exp(atan( ((1:300)-200)/50 )),50), ncol = 300, byrow = TRUE),
  # transient
  matrix(rep(exp(atan( c((1:100)/33, rep(3,100), (100:1)/33) )),50), 
         ncol = 300, byrow = TRUE)
)
counts <- apply(means,2,function(cell_means){
  total <- rnbinom(1, mu = 7500, size = 4)
  rmultinom(1, total, cell_means)
})
rownames(counts) <- paste0('G',1:750) #750 genes
colnames(counts) <- paste0('c',1:300) #300 cells
sce <- SingleCellExperiment(assays = List(counts = counts))

library(slingshot, quietly = FALSE)
data("slingshotExample")
rd <- slingshotExample$rd
cl <- slingshotExample$cl
dim(rd) # data representing cells in a reduced dimensional space
length(cl)  # vector of cluster labels

# filter genes down to potential cell-type markers
# at least M (15) reads in at least N (15) cells
geneFilter <- apply(assays(sce)$counts,1,function(x){
  sum(x >= 3) >= 10
})
sce <- sce[geneFilter, ]

#full quantile normalization, a well-established method which forces each cell to have the same distribution of expression values.
FQnorm <- function(counts){
  rk <- apply(counts,2,rank,ties.method='min')
  counts.sort <- apply(counts,2,sort)
  refdist <- apply(counts.sort,1,median)
  norm <- apply(rk,2,function(r){ refdist[r] })
  rownames(norm) <- rownames(counts)
  return(norm)
}
assays(sce)$norm <- FQnorm(assays(sce)$counts)

#Dimensionality Reduction
# pca
pca <- prcomp(t(log1p(assays(sce)$norm)), scale. = FALSE)
rd1 <- pca$x[,1:2]
plot(rd1, col = rgb(0,0,0,.5), pch=16, asp = 1)

# umap
library(uwot)
rd2 <- uwot::umap(t(log1p(assays(sce)$norm)))
colnames(rd2) <- c('UMAP1', 'UMAP2')
plot(rd2, col = rgb(0,0,0,.5), pch=16, asp = 1)

#add both dimensionality reductions to the SingleCellExperiment object
reducedDims(sce) <- SimpleList(PCA = rd1, UMAP = rd2)

#Clustering Cells
# mclust
library(mclust, quietly = TRUE)
cl1 <- Mclust(rd1)$classification
colData(sce)$GMM <- cl1
length(cl1) #300 cells

library(RColorBrewer)
plot(rd1, col = brewer.pal(9,"Set1")[cl1], pch=16, asp = 1)

# kmeans
cl2 <- kmeans(rd1, centers = 4)$cluster #300cells
colData(sce)$kmeans <- cl2
plot(rd1, col = brewer.pal(9,"Set1")[cl2], pch=16, asp = 1)

#Using Slingshot
sce <- slingshot(sce, clusterLabels = 'GMM', reducedDim = 'PCA')
summary(sce$slingPseudotime_1)

library(grDevices)
colors <- colorRampPalette(brewer.pal(11,'Spectral')[-6])(100)
plotcol <- colors[cut(sce$slingPseudotime_1, breaks=100)]

plot(reducedDims(sce)$PCA, col = plotcol, pch=16, asp = 1)
lines(SlingshotDataSet(sce), lwd=2, col='black')

plot(reducedDims(sce)$PCA, col = brewer.pal(9,'Set1')[sce$GMM], pch=16, asp = 1)
lines(SlingshotDataSet(sce), lwd=2, type = 'lineages', col = 'black')

# 4,Downstream Analysis
# 4.1, Identifying temporally dynamic genes
library(tradeSeq)

# fit negative binomial GAM
sce <- fitGAM(sce)

# test for dynamic expression
ATres <- associationTest(sce)

topgenes <- rownames(ATres[order(ATres$pvalue), ])[1:250]
pst.ord <- order(sce$slingPseudotime_1, na.last = NA)
heatdata <- assays(sce)$counts[topgenes, pst.ord]
heatclus <- sce$GMM[pst.ord]

heatmap(log1p(heatdata), Colv = NA,
        ColSideColors = brewer.pal(9,"Set1")[heatclus])

#5, Detailed Slingshot Functionality
lin1 <- getLineages(rd, cl, start.clus = '1')
lin1
plot(rd, col = brewer.pal(9,"Set1")[cl], asp = 1, pch = 16)
lines(SlingshotDataSet(lin1), lwd = 3, col = 'black')

lin2 <- getLineages(rd, cl, start.clus= '1', end.clus = '3')
plot(rd, col = brewer.pal(9,"Set1")[cl], asp = 1, pch = 16)
lines(SlingshotDataSet(lin2), lwd = 3, col = 'black', show.constraints = TRUE)
#5.2 Constructing smooth curves and ordering cells
crv1 <- getCurves(lin1)
crv1
plot(rd, col = brewer.pal(9,"Set1")[cl], asp = 1, pch = 16)
lines(SlingshotDataSet(crv1), lwd = 3, col = 'black')

# multiple trajectories
rd2 <- rbind(rd, cbind(rd[,2]-12, rd[,1]-6))
cl2 <- c(cl, cl + 10)
pto2 <- slingshot(rd2, cl2, omega = TRUE, start.clus = c(1,11))

plot(rd2, pch=16, asp = 1,
     col = c(brewer.pal(9,"Set1"), brewer.pal(8,"Set2"))[cl2])
lines(SlingshotDataSet(pto2), type = 'l', lwd=2, col='black')

plot(rd2, pch=16, asp = 1,
     col = c(brewer.pal(9,"Set1"), brewer.pal(8,"Set2"))[cl2])
lines(SlingshotDataSet(pto2), lwd=2, col='black')

