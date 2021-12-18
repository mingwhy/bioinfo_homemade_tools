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


# begin scenic, 36 cell types, 8hr, 20core on server
my.path='./male_scenic/';
dir.create(my.path)
setwd(my.path)

for(cell.type in cell.types){
  #cell.type=cell.types[1] #"adult brain perineurial glial cell"
  cell.type2=stringr::str_replace_all(cell.type, '\\/', '\\_')
  
  dir.create(cell.type2)
  setwd(cell.type2)
  
  cell=mydata[,mydata@meta.data$annotation==cell.type]
  exprMat=as.matrix(cell@assays$RNA@counts) #gene by cell matrix
  dim(exprMat) # 12616   735 
  
  ### Initialize SCENIC settings
  org <- "dmel" # or hgnc, or mgi
  dbDir <- "../../cisTarget_databases" # RcisTarget databases location
  myDatasetTitle <- "SCENIC on flY brain" # choose a name for your analysis
  data(defaultDbNames)
  dbs <- defaultDbNames[[org]]
  dbs
  # 5kbUpTx 
  #"dm6-5kb-upstream-full-tx-11species.mc8nr.feather"
  scenicOptions <- initializeScenic(org=org, dbDir=dbDir, dbs=dbs, 
                                    datasetTitle=myDatasetTitle, nCores=20) 
  scenicOptions@settings
  
  # Save to use at a later time...
  saveRDS(scenicOptions, file="scenicOptions.Rds")
  
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
  #set.seed(123)
  #i=sample(1:ncol(exprMat_filtered),200,replace = F)
  #exprMat_filtered=exprMat_filtered[,i]
  runGenie3(exprMat_filtered, scenicOptions)
  #Using 309 TFs as potential regulators...
  # 3113 gene in 735 cells, too slow
  # 3113 gene in 200 cells, 1036~1046, 10min
  
  if(F){
    ### Build and score the GRN (runSCENIC_…)
    scenicOptions <- readRDS("scenicOptions.Rds")
    scenicOptions@settings$verbose <- TRUE
    scenicOptions@settings$nCores <- 20
    scenicOptions@settings$seed <- 123
    
    scenicOptions@settings$dbs
    # coexMethod=c("top5perTarget")
    scenicOptions <- runSCENIC_1_coexNetwork2modules(scenicOptions)
    scenicOptions <- runSCENIC_2_createRegulons(scenicOptions, coexMethod=c("top5perTarget")) 
    # make take several mins
    
    scenicOptions <- runSCENIC_3_scoreCells(scenicOptions, exprMat_filtered)
    
    saveRDS(scenicOptions, file="scenicOptions.Rds") # To save status
  }
  setwd('../')
}  
