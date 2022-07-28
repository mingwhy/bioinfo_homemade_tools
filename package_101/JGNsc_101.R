
#https://meichendong.github.io/JGNsc/articles/JGNsc.html

#devtools::install_github("meichendong/JGNsc")
#BiocManager::install('matrixcalc')
#BiocManager::install('huge',force=T)
#BiocManager::install('JGL',force=T)
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
length(simulation.result$partcorr)#the list of partial correlation matrices)
length(simulation.result$JGL$theta)# precision matrices 

simulation.result$partcorr[[1]][1:5,1:5]
simulation.result$JGL$theta[[1]][1:5,1:5]
# check the AIC values
library(ggplot2)
dt <- as.data.frame(simulation.result$aic.table)
ggplot(dt, aes(x=V1, y=V3, color = V2)) + geom_point() + xlab("lambda1") + ylab("AIC") + labs(color = "lambda2")

library(pheatmap)
# ESTIMATED
dim(simulation.result$partcorr[[1]])
tmp=simulation.result$partcorr[[1]];
x=tmp[upper.tri(tmp)]
summary(x)
partcorr.est1.trunc <- trunc_precision(simulation.result$partcorr[[1]], 
                                       threshold = 0.0001)
pheatmap(partcorr.est1.trunc, cluster_rows = F, cluster_cols = F)

# SIMULATION TRUTH
partcorr.true1 <- prec2partialcorr(countlist.1[[1]]$precision)
partcorr.true1.trunc <- trunc_precision(partcorr.true1, threshold = 0.0001)
partcorr.true1[1:6,1:6]
pheatmap(partcorr.true1, cluster_rows = F, cluster_cols = F)

partcorr.est1.trunc[1:6,1:6]
partcorr.true1[1:6,1:6]

########################
#https://github.com/meichendong/JGNsc/tree/main/vignettes/data
#MB.3cond <- getObservedList(mtx = as.matrix(sobj@assays$RNA@counts), group = sobj$group3, geneSet = c(toupper(meta2$GeneSymbol), "MYC", "OTX2"))
# the list of observed count matrices can be import as below
MB.3cond <- readRDS("J0001_MB3cond.rds")
length(MB.3cond)
sapply(MB.3cond,dim)
#     Intermediate Group 4 Group 3
#[1,]          775     775     775
#[2,]          952    1765    2138
MB.3cond[[1]][1:3,1:3] #gene X cell matrix

if(F){
#Run JGNsc imputation and continuization step. The current version may take hours or longer if the data dimension is high. Faster computational method is in development.
MB.3cond.continuous <- RunJGNsc(observed.list = MB.3cond, min.cell = 20, runNetwork = F)

#Tuning parameter selection by AIC, for JGL model:
res <- getJGLTuningParamResult(GauList = MB.3cond.continuous$theta.star.npn)
# transform the precision matrix to partial correlation
partcorr <- lapply(res$jgl.res, prec2partialcorr)
# you can read in the calculated partial correlation matrices for the visualization steps
}

partcorr <- readRDS("J0001_MB_partcorr.rds")
length(partcorr) #3 groups
sapply(partcorr,dim)
#     [,1] [,2] [,3]
#[1,]  653  653  653
#[2,]  653  653  653

sapply(partcorr,sum)
#Intermediate Group 4 Group 3

#GSEA
# metabolism enzyme only pathway
meta2 <- read.csv("Mammalian_Metabolic_enzyme_genes.txt") #https://github.com/meichendong/JGNsc/tree/main/vignettes/data
dim(meta2) #904 x 3

test <- Map2Pathways(partcorr.list = partcorr,
                     #conditions = unique(sobj$group3),
                     conditions=c('Group 3','Intermediate','Group 4'),
                     #GeneInterest = "MYC",
                     GeneInterest = "OTX2",
                     pathwayRef = meta2,
                     pathwayRef_geneVariable = "GeneSymbol",
                     pathwayRef_pathVariable = "SigmaMiniMap.Term",
                     threshold =0)
length(test$cond.connect)
sapply(test$cond.connect,dim)

GSEA.table <- test$GSEA.table
dim(GSEA.table) #16 x 12

#Visualize joint networks for MYC connected genes
#devtools::install_github("briatte/ggnet")
#BiocManager::install('sna',force=T)
library(ggnet);

#gene1 = "MYC"
gene1 = "OTX2"
gconnect = c(toupper(test$cond.connect[[1]]$GeneSymbol), gene1)
highlight = gconnect %in% c(gene1)
net1 = plot_onenet(partcorr[[1]][gconnect,gconnect], 
                   #gname = unique(sobj$group3)[1], 
                   gname='Group 3',
                   circlenet = T, nodecolor = c("orange","lightblue")[highlight+1], family.vec = highlight)

gconnect = c(toupper(test$cond.connect[[2]]$GeneSymbol), gene1)
highlight = gconnect %in% c(gene1)
net2 = plot_onenet(partcorr[[2]][gconnect,gconnect], 
                   #gname = unique(sobj$group3)[2], 
                   gname='Intermediate',
                   circlenet = T, nodecolor = c("orange","lightblue")[highlight+1], family.vec = highlight)

gconnect = c(toupper(test$cond.connect[[3]]$GeneSymbol), gene1)
highlight = gconnect %in% c(gene1)
net3 = plot_onenet(partcorr[[3]][gconnect,gconnect], 
                   #gname = unique(sobj$group3)[3], 
                   gname='Group 4',
                   circlenet = T, nodecolor = c("orange","lightblue")[highlight+1], family.vec = highlight)

library(cowplot)
#plot_grid(net1, net2, ncol = 2) #MYC
plot_grid(net1, net2, net3, ncol = 3) #OTX2


