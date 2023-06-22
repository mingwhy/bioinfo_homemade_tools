# locCSN use sLED to test AD gene set in single-cell co-expr network between two groups
# https://xuranw.github.io/locCSN/docs/vignettes.html#comparison-of-asd-and-control-groups
# sLED: https://github.com/lingxuez/sLED
# locCSN: https://github.com/xuranw/locCSN

# data source: 
# 1) https://github.com/xuranw/locCSN/tree/main/DataStore
# 'Velme' folder
# 2) csn_asd_loc_flat_L4.txt
# https://www.dropbox.com/sh/yhwaubuo44zs344/AAA9Di_NdwbOcD1K5oDC5pwFa?dl=0

#####################################################################
# SLED COMPARISON, Apply sLED to correlation matrix of L4 cell group.
library(sLED)
# read in gene expression and metadata files
#setwd('locCSN-main/DataStore/Velme/')
log.mc.cpm.L = read.table('locCSN-main/DataStore/Velme/Velme_log_mc_cpm_L.txt')
meta.mc.L = read.table('locCSN-main/DataStore/Velme/Velme_meta_mc_L.txt')
 
dim(log.mc.cpm.L)#942 AD genes in 1778 cells
dim(meta.mc.L) #1778 cells   4

# Let's take L4 as an example
ct.name = 'L4'
meta.mc.diag = as.numeric(meta.mc.L$diagnosis[meta.mc.L$cluster == ct.name] == 'ASD')
log.mc.L = data.matrix(log.mc.cpm.L[, meta.mc.L$cluster == ct.name])

log.mc.L[1:5, 1:5]
#          mc_L_4   mc_L_7  mc_L_10  mc_L_25  mc_L_28
#SAMD11  0.000000 0.000000 0.000000 0.000000 0.000000
#SKI     5.797950 4.036630 5.298243 0.000000 3.842033
#SLC45A1 0.000000 2.814837 0.000000 0.000000 2.269254
#RERE    6.489579 5.775307 5.702040 5.917348 5.959781
#CA6     0.000000 1.965827 0.000000 0.000000 1.894637

# rownames of expression are ASD genes
asd.genes = rownames(log.mc.L)

# input are two gene expression matrix
result.cor = sLED(X = t(log.mc.L[, meta.mc.diag == 0]), Y = t(log.mc.L[, meta.mc.diag == 1]), 
                  sumabs.seq = 0.2, npermute = 100, seeds = c(1:100), adj.beta = 0)
# 100 permutation started:
# 10 ,20 ,30 ,40 ,50 ,60 ,70 ,80 ,90 ,100 ,permutations finished.

result.cor$pVal
# [1] 0.8

####################################################
# load functions of sLED for CSNs
source('https://raw.githubusercontent.com/xuranw/locCSN/main/Rcode/sLEDmodify.R')
#For each CSN matrix, we vectorize it to a vector, then column-bind the vectors by cells. 
#The final csn flatten matrix is a gene pair * cell matrix. The flattened matrix for G genes and N cells is of size G(G−1)/2×N.
ct_name='L4'
csn.flat.temp = read.table(paste0('csn_asd_loc_flat_',ct_name, '.txt'))
dim(csn.flat.temp) # 443211    449
csn.flat.temp[1:3,1:3]
choose(942,2) # 443211
# nrow in csn.flat.temp == pairwise interaction between 942 AD genes.

csn.flat.temp = data.matrix(csn.flat.temp)
csn.t.flat = (csn.flat.temp > qnorm(0.99)) + 0 #Threshold at alpha = 0.01
dim(csn.t.flat) # 443211    449

length(meta.mc.diag) #449
dim(csn.t.flat[, meta.mc.diag == 0]) #443211    211
dim(csn.t.flat[, meta.mc.diag == 1]) #443211    238
result.csn = sLED.csn(X = csn.t.flat[, meta.mc.diag == 0], Y = csn.t.flat[, meta.mc.diag == 1], 
                      sumabs.seq = 0.2, npermute = 100, seeds = c(1:100))

result.csn$pVal
# [1] 0
names(result.csn)
dim(result.csn$leverage) #1 942

#LEVERAGE GENES AND DN GENES
# Leverage genes 
length(asd.genes) #942
lev.L4 = asd.genes[result.csn$leverage > 0]
length(lev.L4) #83 genes

# DN genes (top 90%)
num.dn = min(which(cumsum(sort(result.csn$leverage, decreasing = T)) > 0.9))  
num.dn #26

dn.L4.id = which(result.csn$leverage >= sort(result.csn$leverage, decreasing = T)[num.dn])
dn.L4.id; #26 length
dn.L4 = asd.genes[dn.L4.id]
dn.L4; #26 genes

##################################################
## sLED require input as gene expression matrix (https://github.com/lingxuez/sLED/blob/master/R/sLED.R)
## sLEDmodify.r require input as flat gene networks (https://github.com/xuranw/locCSN/blob/main/Rcode/sLEDmodify.R)
## the Differential network is recovered via 'getDiffMatrix.csn' function in sLEDmodify.r
X = csn.t.flat[, meta.mc.diag == 0];
Y = csn.t.flat[, meta.mc.diag == 1];
avg.x = rowMeans(X); avg.y = rowMeans(Y); avg.d = avg.x - avg.y;
ng = ceiling(sqrt(length(avg.d)*2));
length(avg.d) #443211 = choose(942,2)
ng #recover 942 interacting genes
D.hat = matrix(0, ng, ng); k = 0
for(i in 1:(ng-1)){
  for(j in (i+1):ng){
    k = k+1
    D.hat[i, j] <- avg.d[k]
    D.hat[j, i] <- avg.d[k]
  }
}
dim(D.hat) #942 by 942 size, differential network
