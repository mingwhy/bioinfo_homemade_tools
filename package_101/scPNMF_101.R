
#https://github.com/JSB-UCLA/scPNMF
#https://htmlpreview.github.io/?https://github.com/JSB-UCLA/scPNMF/blob/main/inst/docs/scPNMF2.html

library(RcppEigen)
library(ggplot2)
library(tibble)
library(irlba)
library(Rcpp)
library(parallel)
#importFrom(akmedoids,elbow_point)
#importFrom(magrittr,"%>%")
#importFrom(methods,is)
#importFrom(parallel,mclapply)
#importFrom(tibble,as_tibble)

load('scPNMF-main/data/zheng4.rda')
#data(zheng4, package = "scPNMF")
Input_zheng4 <- SingleCellExperiment::logcounts(zheng4)
dim(Input_zheng4) #2192 gene by 3994 cells

########################################
## scPNMF::PNMFfun
Rcpp::sourceCpp('scPNMF-main/src/all_func.cpp')
# return to Console panel

#res_pnmf <- scPNMF::PNMFfun(X = Input_zheng4,K = 15, method="EucDist", tol=1e-4, maxIter=1000, verboseN = TRUE)
X = Input_zheng4;K = 15; 
method="EucDist"; tol=1e-4; maxIter=1000; verboseN = TRUE;
zerotol=1e-10;maxIter=500; label=NULL;mu=1;lambda=0.01;seed=123;

if (is.null(rownames(X))) {
  stop("Gene names missing!")
}
if (is.null(colnames(X))) {
  stop("Cell names missing!")
}

#nmfmod <- NMF::nmf(X, rank)
set.seed(seed)
Init <- irlba(X, nv = K)
Winit <- Init$u
Winit <- abs(Winit)

if (method == "EucDist") {
  W <- PNMF_EucDistC(X, Winit, tol, maxIter, verboseN, zerotol) 
  W <- W/norm(W, "2")
  ld <- t(X) %*% W
}else if (method == "KL") {
  W <- PNMF_KLC(X, Winit, tol, maxIter, verboseN, zerotol) 
  W <- W/norm(W, "2")
  ld <- t(X) %*% W
}else if (method == "DPNMF") {
  if (length(label) != dim(X)[2]) {
    stop("Cluster labels must have same length as number of cells.")
  }
  cluvec <- as.factor(label)
  cluvec.num <- as.numeric(cluvec)
  cluvec.ord <- order(cluvec.num)
  Xord <- X[, cluvec.ord]
  clunum <- as.integer(table(cluvec))
  
  W <- DPNMFC(X, Winit, tol, maxIter, verboseN, zerotol, Xord, clunum, mu, lambda) 
  W <- W/norm(W, "2")
  ld <- t(X) %*% W
}

rownames(W) <- rownames(X)
rownames(ld) <- colnames(X) 
colnames(W) <- paste0("Basis", 1:K)
colnames(ld) <- paste0("Basis", 1:K)

res_pnmf=list(Weight=W, Score=ld)
#return(list(Weight=W, Score=ld))

########################################
W <- res_pnmf$Weight
S <- res_pnmf$Score


