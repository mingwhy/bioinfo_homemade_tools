Sys.setenv(RETICULATE_PYTHON = "/Users/mingyang/anaconda3/envs/RCFGL/bin/python")
library(reticulate)
py_config()

library(Matrix)
source("./3clusters.R")

N <- 300
G <- 1800
Trials <- 100

Results <- matrix(nrow = 3*Trials, ncol = 2)

set.seed(1)
k <- 1
for(i in 1:Trials){
  
  for(j in 1:3){
    if(j == 1){
      CountsMatrix <- One_Cluster(N, G)
    } else if(j ==2){
      CountsMatrix <- Two_Clusters(N, G)
    } else {
      CountsMatrix <- Three_Clusters(N, G)
    }
    
    Results[k,2] <- j
    
    py_run_file("./HyperEigen.py")
    Results[k,1] <- as.numeric(py$large_val)
    
    k <- k+1
  }
}


tapply(Results[,1], Results[,2], mean)

cor.test(Results[,2], Results[,1], "two.sided", "spearman", exact = FALSE)

