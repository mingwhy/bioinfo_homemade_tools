
library(Matrix)
library(reticulate)
source("./Path2Data.R")

DataSets <- list(Path2Svensson, Path2Tian, Path2Stumpf)

Results <- vector(length = length(DataSets))
names(Results) <- c("Svensson", "Tian", "Stumpf")

for(i in 1:length(Results)){
  load(DataSets[[i]])
  CountsMatrix <- t(CountsMatrix)
  
  py_run_file("./HyperEigen.py")

  Results[i] <- as.numeric(py$large_val)
  rm(CountsMatrix)
}
