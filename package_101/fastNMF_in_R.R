sessionInfo()
#R version 4.0.2 (2020-06-22)
#Platform: x86_64-apple-darwin17.0 (64-bit)
#Running under: OS X  12.1
if(F){
  devtools::install_github("zdebruine/RcppSparse")  #https://github.com/zdebruine/RcppSparse
  devtools::install_github("zdebruine/RcppML") # #https://github.com/zdebruine/RcppML
}
#https://github.com/zdebruine/RcppML
library(RcppML)
packageVersion("RcppML") # 0.5.4
?RcppML::nmf


library(Matrix)
# basic NMF
A =rsparsematrix(1000, 100, 0.1)
dim(A) #1000 100
model <- nmf(A, k = 10)
model
dim(model@w) #1000 10
model@d #10 eigenvalues
dim(model@h) #10 100

# compare rank-2 NMF to second left vector in an SVD
data(iris)
A <- as(as.matrix(iris[,1:4]), "dgCMatrix")
nmf_model <- nmf(A, 2, tol = 1e-5)
nmf_model

bipartitioning_vector <- apply(nmf_model$w, 1, diff)
second_left_svd_vector <- base::svd(A, 2)$u[,2]
abs(cor(bipartitioning_vector, second_left_svd_vector))

# compare rank-1 NMF with first singular vector in an SVD
abs(cor(nmf(A, 1)$w[,1], base::svd(A, 2)$u[,1]))

# symmetric NMF
A <- crossprod(rsparsematrix(100, 100, 0.02))
dim(A)
model <- nmf(A, 10, tol = 1e-5, maxit = 1000)
plot(model$w, t(model$h))
# see package vignette for more examples

