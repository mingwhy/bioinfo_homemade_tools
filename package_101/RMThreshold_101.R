#http://matstat.org/content_en/RMT/RMThreshold_Intro.pdf
library(RMThreshold)

# check if a matrix is well-conditioned for the RMT-based algorithm
## function: rm.matrix.validation
## eg1
set.seed(777)
random.mat <- create.rand.mat(size = 1000, distrib = "norm")$rand.matr
dim(random.mat) #1000 by 1000 symmetric matrix
isSymmetric(random.mat)
res <- rm.matrix.validation(random.mat)
str(res)

## eg2
library(igraph); library(Matrix)
g <- erdos.renyi.game(1000, 0.1)
image(as.matrix(get.adjacency(g)))
rm.matrix.validation(as.matrix(get.adjacency(g)))

## eg3
matlist = list()
for (i in 1:4) matlist[[i]] = get.adjacency(erdos.renyi.game(250, 0.1))
mat <- bdiag(matlist)
dim(mat)
image(mat)
rm.matrix.validation(as.matrix(mat))

# Finding a candidate signal-noise separating threshold for the matrix
## function: rm.get.threshold
set.seed(777)
random.mat <- create.rand.mat(size = 1000, distrib = "norm")$rand.matr
res <- rm.get.threshold(random.mat)
str(res)

## 
matlist = list()
set.seed(979)
for (i in 1:4) matlist[[i]] = get.adjacency(erdos.renyi.game(250, 0.1))
mat <- bdiag(matlist)
dim(mat) #1000 by 1000
rm.matrix.validation(as.matrix(mat))

m <- mat != 0
g <- graph.adjacency(m, mode = "undirected")
clusters(g) # 4 clusters of size 250, as expected

set.seed(979)
mat1 = add.Gaussian.noise(as.matrix(mat), mean = 0, stddev = 0.1)
rm.matrix.validation(mat1)
m1 <- mat1 != 0
g1 <- graph.adjacency(m1, mode = "undirected")
clusters(g1) # a single big cluster with 1000 nodes

res <- rm.get.threshold(mat1) # noisy matrix as input
res

# running the main algorithm on a smaller interval of thresholds
res <- rm.get.threshold(random.mat, interval = c(2.5, 3.5))
cleaned <- rm.denoise.mat(mat1, 0.6)
matr <- cleaned != 0
g <- graph.adjacency(matr, mode = "undirected")
clusters(g) # 4 clusters reconstructed !

# Applying the identified threshold to the matrix
## rm.denoise.mat
cleaned.matrix <- rm.denoise.mat(random.mat, threshold = 3.2)
cleaned.matrix <- rm.discard.zeros(cleaned.matrix)
dim(cleaned.matrix)
m3 <- cleaned.matrix != 0
g3 <- graph.adjacency(m3, mode = "undirected")
clusters(g3)

