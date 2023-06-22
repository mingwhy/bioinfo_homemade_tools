
options(stringsAsFactors = FALSE)
library(Seurat)
library(SCENIC)

setwd('./male_scenic');
files=Sys.glob('./*')
files

for(file in files){
  #setwd('./female_scenic/adult brain perineurial glial cell/')
  setwd(file);
  scenicOptions <- readRDS("scenicOptions.Rds")
  scenicOptions@settings$verbose <- TRUE
  scenicOptions@settings$nCores <- 2
  scenicOptions@settings$seed <- 123
  
  scenicOptions@settings$dbs
  # coexMethod=c("top5perTarget")
  scenicOptions <- runSCENIC_1_coexNetwork2modules(scenicOptions)
  # make take several mins
  scenicOptions <- runSCENIC_2_createRegulons(scenicOptions, coexMethod=c("top5perTarget")) 
  #scenicOptions <- runSCENIC_3_scoreCells(scenicOptions, exprMat_filtered)
  #saveRDS(scenicOptions, file="scenicOptions.Rds") # To save status
  setwd('../')
}

