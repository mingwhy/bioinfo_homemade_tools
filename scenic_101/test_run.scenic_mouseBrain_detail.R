#http://htmlpreview.github.io/?https://github.com/aertslab/SCENIC/blob/master/inst/doc/SCENIC_Running.html
##############################################################
## a detailed version
library(SCENIC)

loomPath <- system.file(package="SCENIC", "examples/mouseBrain_toy.loom")
library(SCopeLoomR)
loom <- open_loom(loomPath)
exprMat <- get_dgem(loom)
cellInfo <- get_cell_annotation(loom)
close_loom(loom)

dim(exprMat) #862 x 200

head(cellInfo)
cellInfo <- data.frame(cellInfo)
cbind(table(cellInfo$CellType))
dir.create("int")
saveRDS(cellInfo, file="int/cellInfo.Rds")

# Color to assign to the variables (same format as for NMF::aheatmap)
colVars <- list(CellType=c("microglia"="forestgreen", 
                           "endothelial-mural"="darkorange", 
                           "astrocytes_ependymal"="magenta4", 
                           "oligodendrocytes"="hotpink", 
                           "interneurons"="red3", 
                           "pyramidal CA1"="skyblue", 
                           "pyramidal SS"="darkblue"))
colVars$CellType <- colVars$CellType[intersect(names(colVars$CellType), cellInfo$CellType)]
saveRDS(colVars, file="int/colVars.Rds")
plot.new(); legend(0,1, fill=colVars$CellType, legend=names(colVars$CellType))

### Initialize SCENIC settings
org <- "mgi" # or hgnc, or dmel
dbDir <- "cisTarget_databases" # RcisTarget databases location
myDatasetTitle <- "SCENIC example on Mouse brain" # choose a name for your analysis
data(defaultDbNames)
dbs <- defaultDbNames[[org]]
scenicOptions <- initializeScenic(org=org, dbDir=dbDir, dbs=dbs, datasetTitle=myDatasetTitle, nCores=2) 

# Modify if needed
scenicOptions@inputDatasetInfo$cellInfo <- "int/cellInfo.Rds"
scenicOptions@inputDatasetInfo$colVars <- "int/colVars.Rds"
# Databases:
# scenicOptions@settings$dbs <- c("mm9-5kb-mc8nr"="mm9-tss-centered-5kb-10species.mc8nr.feather")
# scenicOptions@settings$db_mcVersion <- "v8"

# Save to use at a later time...
saveRDS(scenicOptions, file="int/scenicOptions.Rds") 

### Co-expression network
# (Adjust minimum values according to your dataset)
genesKept <- geneFiltering(exprMat, scenicOptions=scenicOptions,
                           minCountsPerGene=3*.01*ncol(exprMat),
                           minSamples=ncol(exprMat)*.01)

interestingGenes <- c("Sox9", "Sox10", "Dlx5")
# any missing?
interestingGenes[which(!interestingGenes %in% genesKept)]
exprMat_filtered <- exprMat[genesKept, ]
dim(exprMat_filtered) #770 200
rm(exprMat)

### Correlation
#(This step can be run either before/after or simultaneously to GENIE3/GRNBoost)
runCorrelation(exprMat_filtered, scenicOptions)
#GENIE3

# Optional: add log (if it is not logged/normalized already)
exprMat_filtered <- log2(exprMat_filtered+1) 

# Run GENIE3
runGenie3(exprMat_filtered, scenicOptions)
#Using 8 TFs as potential regulators...

### Build and score the GRN (runSCENIC_…)
scenicOptions <- readRDS("int/scenicOptions.Rds")
scenicOptions@settings$verbose <- TRUE
scenicOptions@settings$nCores <- 10
scenicOptions@settings$seed <- 123

#Re-load the expression matrix if necessary:
loom <- open_loom(loomPath)
exprMat <- get_dgem(loom)
close_loom(loom)
# Optional: log expression (for TF expression plot, it does not affect any other calculation)
exprMat_log <- log2(exprMat+1)
dim(exprMat)

# For a very quick run: 
# coexMethod=c("top5perTarget")
scenicOptions@settings$dbs <- scenicOptions@settings$dbs["10kb"] # For toy run
# save...

scenicOptions <- runSCENIC_1_coexNetwork2modules(scenicOptions)
scenicOptions <- runSCENIC_2_createRegulons(scenicOptions, coexMethod=c("top5perTarget")) #** Only for toy run!!
scenicOptions <- runSCENIC_3_scoreCells(scenicOptions, exprMat_log)

saveRDS(scenicOptions, file="int/scenicOptions.Rds") # To save status


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
