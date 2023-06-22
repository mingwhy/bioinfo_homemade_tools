
#paper: 2019-Feature selection and dimension reduction for single-cell RNA-Seq based on a multinomial model
#https://bioconductor.org/packages/release/bioc/html/scry.html
#https://bioconductor.org/packages/release/bioc/vignettes/scry/inst/doc/scry.html

library(ggplot2); theme_set(theme_bw())
#BiocManager::install('DuoClustering2018')
library(DuoClustering2018)
#BiocManager::install("scry")
library(scry)
library(SingleCellExperiment)

sce<-sce_full_Zhengmix4eq()
sce #15568 x 3994
colnames(colData(sce)) #cell meta
colnames(rowData(sce)) #gene meta

#Feature Selection with Deviance: `devianceFeatureSelection` function
sce<-devianceFeatureSelection(sce, assay="counts", sorted=TRUE)
summary(rowData(sce)$binomial_deviance)
plot(rowData(sce)$binomial_deviance, type="l", xlab="ranked genes",
     ylab="binomial deviance", main="Feature Selection with Deviance")
abline(v=2000, lty=2, col="red")

sce2<-sce[1:1000, ]

#Dimension Reduction with GLM-PCA
#GLM-PCA can reduce the dimensionality of UMI counts to facilitate visualization and/or clustering without needing any normalization.
set.seed(101)
sce2<-GLMPCA(sce2, 2, assay="counts")
fit<-metadata(sce2)$glmpca
pd<-cbind(as.data.frame(colData(sce2)), fit$factors)
ggplot(pd, aes(x=dim1, y=dim2, colour=phenoid)) + geom_point(size=.8) +
  ggtitle("GLM-PCA applied to high deviance genes")

#Dimension Reduction with Null Residuals
#GLM-PCA can be slow for large datasets. A fast approximation is to fit a null model of constant expression for each gene across cells, then fit standard PCA to either the Pearson or deviance residuals from the null model.
assayNames(sce)
sce<-nullResiduals(sce, assay="counts", type="deviance")
sce<-nullResiduals(sce, assay="counts", type="pearson")
assayNames(sce)
# two new assays are added: "binomial_deviance_residuals" "binomial_pearson_residuals" 
dim(assay(sce,'binomial_deviance_residuals')) #15568  3994
dim(assay(sce,'binomial_pearson_residuals')) #15568  3994

sce2<-sce[1:1000, ] #use only the high deviance genes

pca<-function(Y, L=2, center=TRUE, scale=TRUE){
  #assumes features=rows, observations=cols
  res<-prcomp(as.matrix(t(Y)), center=center, scale.=scale, rank.=L)
  factors<-as.data.frame(res$x)
  colnames(factors)<-paste0("dim", 1:L)
  factors
}

pca_d<-pca(assay(sce2, "binomial_deviance_residuals"))
pca_d$resid_type<-"deviance_residuals"
head(pca_d)

pca_p<-pca(assay(sce2, "binomial_pearson_residuals"))
pca_p$resid_type<-"pearson_residuals"
head(pca_p)

cm<-as.data.frame(colData(sce2))
pd<-rbind(cbind(cm, pca_d), cbind(cm, pca_p))
dim(sce2)
dim(pd)

ggplot(pd, aes(x=dim1, y=dim2, colour=phenoid)) + geom_point() +
  facet_wrap(~resid_type, scales="free", nrow=2) +
  ggtitle("PCA applied to null residuals of high deviance genes")





