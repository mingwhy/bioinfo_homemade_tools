
#https://meichendong.github.io/JGNsc/articles/JGNsc.html

#devtools::install_github("meichendong/JGNsc")
#BiocManager::install('matrixcalc')
library(huge)#due to 'Error in 'RunJGNsc': lapply(theta.star.t, huge.npn) : object 'huge.npn' not found'
library(JGL)#used in RunJGNsc

library(JGNsc)

set.seed(1)

# example: the first 20 genes have different structure, the rest of genes have the same structures.
scenario = "DI20" # can change this to DD/DI100/ID
nsample = 100

nivec.list.diff <- list(nivec= c(rep(2,10), rep(20,4)),nivec2 = rep(20,5))
diffblk = list( 1:10,1)
sigma.list.1 <- generateSigmaList(nivec.list = nivec.list.diff, structure = "Diff S, Identical W", diffblk = diffblk)
# Create a list of count matrices
countlist.1 <- getCountList(sigma.list = sigma.list.1, nvec = rep(nsample,2), a3=3, b3=1)

str(countlist.1)
observed.list <- list(t(countlist.1[[1]]$count), t(countlist.1[[2]]$count))

sapply(observed.list,dim)
simulation.result <- RunJGNsc(observed.list = observed.list, min.cell = 10, runNetwork = T,
                              l1.vec = seq(1,40, by=3)/100, l2.vec = seq(1,10, by=2)/100)

# check the structure of the result
str(simulation.result)
# check the AIC values
library(ggplot2)
dt <- as.data.frame(simulation.result$aic.table)
ggplot(dt, aes(x=V1, y=V3, color = V2)) + geom_point() + xlab("lambda1") + ylab("AIC") + labs(color = "lambda2")

library(pheatmap)
# ESTIMATED
partcorr.est1.trunc <- trunc_precision(simulation.result$partcorr[[1]], threshold = 0.0001)
pheatmap(partcorr.est1.trunc, cluster_rows = F, cluster_cols = F)

# SIMULATION TRUTH
partcorr.true1 <- prec2partialcorr(countlist.1[[1]]$precision)
partcorr.true1.trunc <- trunc_precision(partcorr.true1, threshold = 0.0001)
partcorr.true1[1:6,1:6]
pheatmap(partcorr.true1, cluster_rows = F, cluster_cols = F)


