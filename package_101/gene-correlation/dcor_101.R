
library(Pigengene) #for loading aml data
data(aml)
dim(aml)
system.time(
  dcor1 <- dcor.matrix(Data=aml[,1:10])
)
dcor1

#https://privefl.github.io/blog/a-guide-to-parallelism-in-r/
library(bigstatsr)
library(foreach)
cl <- parallel::makeCluster(4)
doParallel::registerDoParallel(cl)
mat=FBM(1000,1000)
tmp<-foreach(i = 1:999,.combine='c') %:% 
  foreach(j=(i+1):1000,.combine='c') %dopar% {
  mat[i,j]=energy::dcor(x=aml[,i],y=aml[,j])
  NULL
}
parallel::stopCluster(cl)
#mat[]


## Comparison with Pearson:
cor1 <- abs(cor(aml[,1:5]))
## With 202 samples, distance and Pearson correlations do not differ much:
dcor1-cor1 
dcor2 <- dcor.matrix(Data=aml[1:20,1:5])
cor2 <- abs(cor(aml[1:20,1:5]))
## Distance correlation is more robust if fewer samples are available:
dcor2-cor2
plot(dcor2-cor1,cor1-cor2,xlim=c(-0.5,0.5),ylim=c(-0.5,0.5))
