# online tutorial: http://htmlpreview.github.io/?https://github.com/aertslab/SCENIC/blob/master/inst/doc/SCENIC_Running.html

options(stringsAsFactors = FALSE)
library(Seurat)
library(SCENIC)

## read in processed wholebrain data
file="../data/wholebrain_filtered.rds";
dat=readRDS(file); #a seurat object
#df.expr0=dat@assays$RNA@counts #12616 features across 100527 samples 

# keep male cells
mydata=subset(dat,sex=='male')
mydata 
#cellInfo=mydata@meta.data

# read in cell types to be analyzed
cell.types=scan('../cell.type_common.gene/cell.types.at.least.200_both.sex.txt',what='',sep='\t')
length(cell.types) #36 cell types
#cell.types2=as.character(sapply(cell.types,function(i){stringr::str_replace_all(i, '\\/', '\\_')}))


# begin scenic
my.path='./test_scenic/';
dir.create(my.path)
setwd('./test_scenic/')

cell.type=cell.types[1] #"adult brain perineurial glial cell"
cell.type2=stringr::str_replace_all(cell.type, '\\/', '\\_')

cell=mydata[,mydata@meta.data$annotation==cell.type]
exprMat=as.matrix(cell@assays$RNA@counts) #gene by cell matrix
dim(exprMat) # 12616   735 

### Initialize SCENIC settings
org <- "dmel" # or hgnc, or mgi
dbDir <- "../cisTarget_databases" # RcisTarget databases location
myDatasetTitle <- "SCENIC on flY brain" # choose a name for your analysis
data(defaultDbNames)
dbs <- defaultDbNames[[org]]
dbs
# 5kbUpTx 
#"dm6-5kb-upstream-full-tx-11species.mc8nr.feather"
scenicOptions <- initializeScenic(org=org, dbDir=dbDir, dbs=dbs, datasetTitle=myDatasetTitle, nCores=2) 
scenicOptions@settings

# Save to use at a later time...
saveRDS(scenicOptions, file="int/scenicOptions.Rds") 

### Co-expression network
# (Adjust minimum values according to your dataset)
#minCountsPerGene, Minimum counts per gene required
#minSamples, Minimum number of samples (cells) in which the gene should be detected
genesKept <- geneFiltering(exprMat, scenicOptions=scenicOptions,
                           minCountsPerGene=3*.01*ncol(exprMat),
                           minSamples=ncol(exprMat)*.01)
#3383	genes with counts per gene > 22.05
#3383	genes detected in more than 7.35 cells
length(genesKept) #3113 #'1.1_genesKept.Rds' is saved
dim(exprMat) 

interestingGenes <- c("dsx","lncRNA:roX1","lncRNA:roX2")
interestingGenes[which(!interestingGenes %in% genesKept)]

exprMat_filtered <- exprMat[genesKept, ]
dim(exprMat_filtered) #3113  735
rm(exprMat)

### Correlation
#(This step can be run either before/after or simultaneously to GENIE3/GRNBoost)
runCorrelation(exprMat_filtered, scenicOptions) #1.2_corrMat.Rds is saved

#GENIE3
# Optional: add log (if it is not logged/normalized already)
exprMat_filtered <- log2(exprMat_filtered+1) 

# Run GENIE3 (too slow)
set.seed(123)
i=sample(1:ncol(exprMat_filtered),200,replace = F)
exprMat_filtered=exprMat_filtered[,i]
runGenie3(exprMat_filtered, scenicOptions)
#Using 309 TFs as potential regulators...
# 3113 gene in 735 cells, too slow
# 3113 gene in 200 cells, 1036~1046, 10min


### Build and score the GRN (runSCENIC_…)
scenicOptions <- readRDS("int/scenicOptions.Rds")
scenicOptions@settings$verbose <- TRUE
scenicOptions@settings$nCores <- 2
scenicOptions@settings$seed <- 123

scenicOptions@settings$dbs
# coexMethod=c("top5perTarget")
scenicOptions <- runSCENIC_1_coexNetwork2modules(scenicOptions)
scenicOptions <- runSCENIC_2_createRegulons(scenicOptions, coexMethod=c("top5perTarget")) 
# make take several mins

scenicOptions <- runSCENIC_3_scoreCells(scenicOptions, exprMat_filtered)

saveRDS(scenicOptions, file="int/scenicOptions.Rds") # To save status
########################################################################
########################################################################
########################################################################
#Optional steps:
#Binarize the network activity (regulon on/off)
aucellApp <- plotTsne_AUCellApp(scenicOptions, exprMat_log)
savedSelections <- shiny::runApp(aucellApp)

# Save the modified thresholds:
newThresholds <- savedSelections$thresholds
scenicOptions@fileNames$int["aucell_thresholds",1] <- "int/newThresholds.Rds"
saveRDS(newThresholds, file=getIntName(scenicOptions, "aucell_thresholds"))
saveRDS(scenicOptions, file="int/scenicOptions.Rds") 

#Once you have optimized the thresholds, run runSCENIC_4_aucell_binarize to binarize the AUC, and generate some extra figures and clusterings:
# scenicOptions@settings$devType="png"
scenicOptions <- runSCENIC_4_aucell_binarize(scenicOptions)

### Clustering / dimensionality reduction on the regulon activity
nPcs <- c(5) # For toy dataset
scenicOptions@settings$seed <- 123 # same seed for all of them
# Run t-SNE with different settings:
fileNames <- tsneAUC(scenicOptions, aucType="AUC", nPcs=nPcs, perpl=c(5,15,50))
fileNames <- tsneAUC(scenicOptions, aucType="AUC", nPcs=nPcs, perpl=c(5,15,50), onlyHighConf=TRUE, filePrefix="int/tSNE_oHC")
# Plot as pdf (individual files in int/):
fileNames <- paste0("int/",grep(".Rds", grep("tSNE_", list.files("int"), value=T), value=T))

#Note: The toy dataset only contains ~8 regulons; using more than 8 PCs will not provide any difference…
#and to view/compare them…

par(mfrow=c(length(nPcs), 3))
fileNames <- paste0("int/",grep(".Rds", grep("tSNE_AUC", list.files("int"), value=T, perl = T), value=T))
plotTsne_compareSettings(fileNames, scenicOptions, showLegend=FALSE, varName="CellType", cex=.5)


#The chosen t-SNE can then be saved as default to use for plots (can also be “binary”, see below):
scenicOptions@settings$defaultTsne$aucType <- "AUC"
scenicOptions@settings$defaultTsne$dims <- 5
scenicOptions@settings$defaultTsne$perpl <- 15
saveRDS(scenicOptions, file="int/scenicOptions.Rds")

### Export to loom/SCope
#DGEM (Digital gene expression matrix)
#(non-normalized counts)
exprMat <- get_dgem(open_loom(loomPath))
dgem <- exprMat
head(colnames(dgem))  #should contain the Cell ID/name

# Export:
scenicOptions@fileNames$output["loomFile",] <- "output/mouseBrain_SCENIC.loom"
export2loom(scenicOptions, exprMat)

#Loading results from a .loom file
library(SCopeLoomR)
scenicLoomPath <- getOutName(scenicOptions, "loomFile")
loom <- open_loom(scenicLoomPath)

# Read information from loom file:
regulons_incidMat <- get_regulons(loom)
regulons <- regulonsToGeneLists(regulons_incidMat)
regulonsAUC <- get_regulons_AUC(loom)
regulonsAucThresholds <- get_regulon_thresholds(loom)
embeddings <- get_embeddings(loom)


# Projection the AUC and TF expression onto t-SNEs
exprMat_log <- exprMat # Better if it is logged/normalized
aucellApp <- plotTsne_AUCellApp(scenicOptions, exprMat_log) # default t-SNE
savedSelections <- shiny::runApp(aucellApp)

print(tsneFileName(scenicOptions))
tSNE_scenic <- readRDS(tsneFileName(scenicOptions))
aucell_regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")

# Show TF expression:
par(mfrow=c(2,3))
AUCell::AUCell_plotTSNE(tSNE_scenic$Y, exprMat, aucell_regulonAUC[onlyNonDuplicatedExtended(rownames(aucell_regulonAUC))[c("Dlx5", "Sox10", "Sox9","Irf1", "Stat6")],], plots="Expression")
# Save AUC as PDF:
Cairo::CairoPDF("output/Step4_BinaryRegulonActivity_tSNE_colByAUC.pdf", width=20, height=15)
par(mfrow=c(4,6))
AUCell::AUCell_plotTSNE(tSNE_scenic$Y, cellsAUC=aucell_regulonAUC, plots="AUC")
dev.off()

library(KernSmooth)
library(RColorBrewer)
dens2d <- bkde2D(tSNE_scenic$Y, 1)$fhat
image(dens2d, col=brewer.pal(9, "YlOrBr"), axes=FALSE)
contour(dens2d, add=TRUE, nlevels=5, drawlabels=FALSE)

regulons <- loadInt(scenicOptions, "regulons")
regulons[c("Dlx5", "Irf1")]
