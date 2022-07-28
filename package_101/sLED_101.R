
#https://github.com/lingxuez/sLED
#devtools::install_github("lingxuez/sLED")
library(sLED)
n <- 50
p <- 100
set.seed(99)
X <- matrix(rnorm(n*p, mean=0, sd=1), nrow=n, ncol=p)
set.seed(42)
Y <- matrix(rnorm(n*p, mean=0, sd=1), nrow=n, ncol=p)

result <- sLED(X=X, Y=Y, npermute=50)

names(result)
result$pVal

# power of sLED
n <- 50
p <- 100
## The first population is still standard normal
set.seed(99)
X <- matrix(rnorm(n*p, mean=0, sd=1), nrow=n, ncol=p)
## For the second population, the first 10 features have different correlation structure
s <- 10
sigma.2 <- diag(p)
sigma.2[1:s, 1:s] <- sigma.2[1:s, 1:s] + 0.2
set.seed(42)
Y2 <- MASS::mvrnorm(n, mu=rep(0, p), Sigma=sigma.2)

#Now we run sLED. Note that the changes in covariance matrices occur at 10% of the features, 
#so the ideal sparsity parameter sumabs should be around (usually slightly less than) sqrt(0.1)=0.32. 
#Here, we pick sumabs=0.25, and use 100 permutations:
## for reproducibility, let's also set the seeds for permutation
result <- sLED(X=X, Y=Y2, sumabs.seq=0.25, npermute=100, seeds = c(1:100))
result$pVal
## [1] 0.04 #The p-value is near zero. 
#Further more, let's check which features are detected by sLED, 
#that is, have non-zero leverage:
which(result$leverage != 0)
#[1]  1  2  3  4  6  7  9 10 16 26 46 50 52 83

#We see that sLED correctly identifies most of the first 10 signaling features. 
#We can also run sLED across a range of sparsity parameters sumabs at once:
## here we let sumabs.seq to be a vector of 3 different sparsity parameters
result <- sLED(X=X, Y=Y2, sumabs.seq=c(0.2, 0.25, 0.3), npermute=100, seeds = c(1:100))

## we can check the 3 p-values, which are all < 0.05
result$pVal
## [1] 0.04 0.03 0.02

## let's also look at which features have non-zero leverage
detected.genes <- apply(result$leverage, 1, function(x){which(x!=0)})
names(detected.genes) <- paste0("sumabs=",result$sumabs.seq)
detected.genes
## $`sumabs=0.2`
## [1]  2  6  7  9 16 26 50 52
## 
## $`sumabs=0.25`
## [1]  1  2  3  4  6  7  9 10 16 26 46 50 52 83
## 
## $`sumabs=0.3`
## [1]  1  2  3  4  6  7  9 10 16 18 19 20 26 35 45 46 50 52 60 63 83

#Sparsity parameter
#A key tuning parameter for sLED is the sparsity parameter, sumabs. 
#This corresponds to the parameter c in equation (2.16) in our paper. 
#Larger values of sumabs correspond to denser solutions. 
#Roughly speaking, sumabs^2 provides a loose lowerbound on the proportion of features to be detected. 
#sumabs can take any value between 1/\sqrt{p} and 1. 
#In practice, a smaller value is usually preferred for better interpretability.

#Parallelization
result_multicore <- sLED(X=X, Y=Y, npermute=1000, useMC=TRUE, mc.cores=2)

