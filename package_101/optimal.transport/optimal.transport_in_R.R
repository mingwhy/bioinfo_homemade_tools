
#https://www.rdocumentation.org/packages/transport/versions/0.12-2/topics/wasserstein
library(transport)
set.seed(27)
x <- pp(matrix(runif(500),250,2))
y <- pp(matrix(runif(500),250,2))
x; #250 obs x 2 dim
y; #250 obs x 2 dim
wasserstein(x,y,p=1) #0.067
wasserstein(x,y,p=2) #0.079

#######################################################################################
#https://www.r-bloggers.com/2017/07/matching-optimal-transport-and-statistical-tests/
#https://freakonometrics.hypotheses.org/50958

#remotes::install_github("pierrejacob/winference") #for cost_matrix_Lp 
# didn't work, change to python (https://pythonot.github.io/auto_examples/plot_Intro_OT.html)
set.seed(13)
npoints <- 25
mu1 <- c(1,1)
mu2 <- c(0,2)
Sigma1 <- diag(1, 2, 2)
Sigma2 <- diag(1, 2, 2)
Sigma2[2,1] <- Sigma2[1,2] <- -0.5
Sigma1 <- 0.4 * Sigma1
Sigma2 <- 0.4 *Sigma2

library(mnormt)
X1 <- rmnorm(npoints, mean = mu1, Sigma1)
X2 <- rmnorm(npoints, mean = mu2, Sigma2)
plot(X1[,1], X1[,2], ,col="blue")
points(X2[,1], X2[,2], col = "red")


ground_p <- 2
p <- 1
w1 <- rep(1/npoints, npoints)
w2 <- rep(1/npoints, npoints)

library(transport)
library(winference)
C <- cost_matrix_Lp(t(X1), t(X2), ground_p)
a <- transport(w1, w2, costm = C^p, method = "shortsimplex")
