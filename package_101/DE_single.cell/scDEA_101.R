#https://github.com/Zhangxf-ccnu/scDEA

# BPSC
devtools::install_github("nghiavtr/BPSC")

# DEsingle
BiocManager::install('numDeriv',force=T);
BiocManager::install('miscTools',force=T)
BiocManager::install("DEsingle",force=T)

# DESeq2
library(DESeq2)
#BiocManager::install("DESeq2")

# edgeR (http://www.bioconductor.org/packages/release/bioc/html/edgeR.html)
library(edgeR)
#BiocManager::install("edgeR")

#MAST (http://www.bioconductor.org/packages/release/bioc/html/MAST.html)
library(MAST)
#BiocManager::install("MAST")

#monocle (http://www.bioconductor.org/packages/release/bioc/html/monocle.html)
BiocManager::install('combinat',force=T,update = F)
BiocManager::install('FNN',force=T,update = F)
BiocManager::install('fastICA',force=T,update = F)
BiocManager::install('qlcMatrix',force=T,update = F)
BiocManager::install('sparsesvd',force=T,update = F)
BiocManager::install('docopt',force=T,update = F)
BiocManager::install("monocle",force=T,update = F)
library(monocle)

#scDD (http://www.bioconductor.org/packages/release/bioc/html/scDD.html)
BiocManager::install("scDD")
library(scDD)

#limma (http://www.bioconductor.org/packages/release/bioc/html/limma.html)
#BiocManager::install("limma")
library(limma)

#Seurat (http://www.bioconductor.org/packages/release/bioc/html/Seurat.html)
#BiocManager::install("Seurat")
library(Seurat)

#zingeR (https://github.com/statOmics/zingeR/)
devtools::install_github("statOmics/zingeR")
library(zingeR)

#SingleCellExperiment (http://www.bioconductor.org/packages/release/bioc/html/SingleCellExperiment.html)
library(SingleCellExperiment)
#BiocManager::install("SingleCellExperiment")

#scater (http://www.bioconductor.org/packages/release/bioc/html/scater.html)
library(scater)
#BiocManager::install("scater")

#aggregation (https://cran.r-project.org/web/packages/aggregation/)
install.packages("aggregation")
library(aggregation)

## install scDEA from Github
devtools::install_github("Zhangxf-ccnu/scDEA")

# start
library(scDEA)
#Simply run the scDEA on the Grun datasets.
data("Grun.counts.hvg")
data("Grun.group.information")

#load("./scDEA-main/data/Grun.counts.hvg.RData")
#load("./scDEA-main/data/Grun.group.information.RData")
ls()

dim(Grun.counts.hvg) #2000 x 716
length(Grun.group.information) #716, two groups
Grun.counts.hvg[10:20,10:20] #a UMI count matrix

Pvals <- scDEA_individual_methods(raw.count = Grun.counts.hvg, cell.label = Grun.group.information)
dim(Pvals) #2000 genes by 12 pvalues

combination.Pvals <- lancaster.combination(Pvals, weight = TRUE, trimmed = 0.2)
length(combination.Pvals) #2000 combined.pvalues

adjusted.Pvals <- scDEA.p.adjust(combination.Pvals, adjusted.method = "bonferroni")
length(adjusted.Pvals)  #2000 adj pvalues

par(mfrow=c(2,1))
hist(combination.Pvals)
hist(adjusted.Pvals)
