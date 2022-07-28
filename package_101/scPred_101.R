
# install (https://github.com/powellgenomicslab/scPred)
if(F){
  install.packages("harmony")
  devtools::install_github("powellgenomicslab/scPred")
}

## used case (https://powellgenomicslab.github.io/scPred/articles/introduction.html)

library("scPred")

library("Seurat")
library("magrittr")

reference <- scPred::pbmc_1
query <- scPred::pbmc_2
dim(reference); #32838  3500, a seurat obj
dim(query); #32838  3000, gene by cell

# preprocess
system.time( 
  reference <- reference %>% 
    NormalizeData() %>% 
    FindVariableFeatures() %>% 
    ScaleData() %>% 
    RunPCA() %>% 
    RunUMAP(dims = 1:30)
)
# 15sec

DimPlot(reference,group.by='cell_type',label=TRUE,repel=TRUE)

# start training: Training classifiers with scPred
##getFeatureSpace will create a `scPred` object stored in the @misc slot. 
reference <- getFeatureSpace(reference, "cell_type")
names(reference@misc)
reference@misc$scPred

## Secondly, we train the classifiers for each cell using the trainModel function. 
## By default, scPred will use a support vector machine with a radial kernel.
system.time( reference<-trainModel(reference) ) #80secs
tmp=reference@misc$scPred@train
class(tmp) #list
names(tmp)

get_probabilities(reference) %>% head()

get_scpred(reference)

plot_probabilities(reference)
#scPred is built on top of the caret package and allows using a large set of prediction models (e.g. logistic regression, decision trees, bagging, neural networks, etc).

# start predction: Cell classification
query<-NormalizeData(query)
system.time(query <- scPredict(query, reference)) #12sec

#scPred now uses Harmony to align the query data onto the training low-dimensional space used as reference. Once the data is aligned, cells are classified using the pre-trained models.

DimPlot(query, group.by = "scpred_prediction", reduction = "scpred")
query <- RunUMAP(query, reduction = "scpred", dims = 1:30)
DimPlot(query, group.by = "scpred_prediction", label = TRUE, repel = TRUE)

FeaturePlot(query, c("scpred_B.cell", "scpred_CD4.T.cell", "scpred_CD8.T.cell", 
                     "scpred_cMono", "scpred_ncMono", "scpred_Plasma.cell", 
                     "scpred_cDC", "scpred_pDC"))
crossTab(query, "cell_type", "scpred_prediction")


# advanced options
## get classifier
get_classifiers(reference)
caret::plot.train(get_classifiers(reference)[["NK cell"]])
## use a different probability threshold
query <- scPredict(query, reference, recompute_alignment = FALSE, threshold = 0.9)

## Parallel training (See Caret parallel processing for more details)
library(doParallel)
cl <- makePSOCKcluster(2)
registerDoParallel(cl)
reference <- trainModel(reference, model = "mda", allowParallel = TRUE)

## extract the classifier from the Seurat obj
scpred <- get_scpred(reference)
query <- scPredict(query, scpred)