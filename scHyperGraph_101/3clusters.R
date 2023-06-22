
asSparseMatrix <- function(CountsMatrix){
  
  if(all(CountsMatrix == 0)){
    CountsMatrix <- sparseMatrix(i = {}, j = {}, dims = dim(CountsMatrix))
  } else {
    
    Ind <- which(CountsMatrix != 0, arr.ind = T)
    
    CountsMatrix <- sparseMatrix(i = Ind[,1],
                                 j = Ind[,2], 
                                 x = CountsMatrix[Ind],
                                 dims = dim(CountsMatrix))
  }
  CountsMatrix
}

One_Cluster <- function(N, G){
  
  Counts <- rpois(N, lambda = rgamma(n = 1, rate = 0.3, shape = 0.6))
  for(i in 2:G){
    Counts <- cbind(Counts, rpois(N, lambda = rgamma(n = 1, rate = 0.3, shape = 0.6)))
  }
  
  asSparseMatrix(Counts)
}

Two_Clusters <- function(N, G){
  
  stopifnot(N/2 == round(N/2))
  N <- N/2
  
  Counts1 <- rpois(N, lambda = rgamma(n = 1, rate = 0.3, shape = 0.6))
  for(i in 2:(G/2)){
    Counts1 <- cbind(Counts1, rpois(N, lambda = rgamma(n = 1, rate = 0.3, shape = 0.6)))
  }
  for(i in 1:(G/2)){
    Counts1 <- cbind(Counts1, rpois(N, lambda = rgamma(n = 1, rate = 3, shape = 0.6)))
  }
  
  Counts2 <- rpois(N, lambda = rgamma(n = 1, rate = 3, shape = 0.6))
  for(i in 2:(G/2)){
    Counts2 <- cbind(Counts2, rpois(N, lambda = rgamma(n = 1, rate = 3, shape = 0.6)))
  }
  for(i in 1:(G/2)){
    Counts2 <- cbind(Counts2, rpois(N, lambda = rgamma(n = 1, rate = 0.3, shape = 0.6)))
  }
  
  Counts <- rbind(Counts1, Counts2)
  asSparseMatrix(Counts)
}


Three_Clusters <- function(N, G){

  stopifnot(N/3 == round(N/3))
  N <- N/3
  
  Counts1 <- rpois(N, lambda = rgamma(n = 1, rate = 0.3, shape = 0.6))
  for(i in 2:(G/3)){
    Counts1 <- cbind(Counts1, rpois(N, lambda = rgamma(n = 1, rate = 0.3, shape = 0.6)))
  }
  for(i in 1:(2*(G/3))){
    Counts1 <- cbind(Counts1, rpois(N, lambda = rgamma(n = 1, rate = 3, shape = 0.6)))
  }
  
  Counts2 <- rpois(N, lambda = rgamma(n = 1, rate = 3, shape = 0.6))
  for(i in 2:(G/3)){
    Counts2 <- cbind(Counts2, rpois(N, lambda = rgamma(n = 1, rate = 3, shape = 0.6)))
  }
  for(i in 1:(G/3)){
    Counts2 <- cbind(Counts2, rpois(N, lambda = rgamma(n = 1, rate = 0.3, shape = 0.6)))
  }
  for(i in 1:(G/3)){
    Counts2 <- cbind(Counts2, rpois(N, lambda = rgamma(n = 1, rate = 3, shape = 0.6)))
  }
  
  Counts3 <- rpois(N, lambda = rgamma(n = 1, rate = 3, shape = 0.6))
  for(i in 2:(2*(G/3))){
    Counts3 <- cbind(Counts3, rpois(N, lambda = rgamma(n = 1, rate = 3, shape = 0.6)))
  }
  for(i in 1:(G/3)){
    Counts3 <- cbind(Counts3, rpois(N, lambda = rgamma(n = 1, rate = 0.3, shape = 0.6)))
  }
  
  Counts <- rbind(Counts1, Counts2, Counts3)
  asSparseMatrix(Counts)
}
