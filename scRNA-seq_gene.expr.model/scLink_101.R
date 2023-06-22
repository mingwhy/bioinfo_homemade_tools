#https://github.com/Vivianstats/scLink/blob/master/vignettes/scLink-vignette.Rmd
library(scLink)
count = readRDS(system.file("extdata", "example.rds", package = "scLink"))
genes = readRDS(system.file("extdata", "genes.rds", package = "scLink"))

dim(count); #793cell by 23341 gene
length(genes) #500 genes for network construction
count.norm = sclink_norm(count, scale.factor = 1e6, filter.genes = FALSE, gene.names = genes)
dim(count.norm)  #793cell by 500gene

#If users do not have a particular gene list for network inference, they can set `filter.genes=TRUE` to filter for the top $n$ genes with largest average expression values. For example:
#count.norm = sclink_norm(count, scale.factor = 1e6, filter.genes = TRUE, n = 500)


#After the pre-processing step, we use the function `sclink_net` 
# to calculate the robust correlation matrix and identifed sparse co-expression network of scLink.
#`expr` is the normalized count matrix output by `sclink_norm` or supplied by the users.
#`lda` is the candidate regularization parameters used in scLink's graphical model. The users can set `ncores` to take advantage of parallel computation.
seq(0.5, 0.1, -0.05)
networks = sclink_net(expr = count.norm, ncores = 1, lda = seq(0.5, 0.1, -0.05))

#`sclink_net` returns a list of results. 
#The scLink's robust correlation matrix can be retrieved from the `cor` element:
networks$cor[1:3,1:3]
robust.cor=networks$cor
x=robust.cor[upper.tri(robust.cor,diag = F)]
summary(x)

# correlation matrix, concentration matrix are different!
#The gene co-expression networks and summary statistics 
#can be retrieved from the `summary` element, 
#which is a list with the same length as `lda`: each element corresponds to one regularization parameter.
length(networks$summary)
lapply(networks$summary,function(i){ c(i$bic,i$lambda,i$nedge)} )

# 9 lda values
net1 = networks$summary[[1]]
names(net1)
# "adj"    "Sigma"  "nedge"  "bic"    "lambda"
### adjacency matrix
net1$adj[1:3,1:3]
(sum(net1$adj)-nrow(net1$adj))/2  #diag
# save as net1$nedge

### concentration matrix
net1$Sigma[1:3,1:3]
x=net1$Sigma
x1=x[upper.tri(x)]
summary(x1)
sum(x1!=0) # elements !=0 zero in the concentration matrix <-> a edge present

### BIC 
net1$bic
### number of edges
net1$nedge
### regularization parameter lambda
net1$lambda

#Since it is very difficult to infer co-expression relationships for lowly expressed genes in single-cell data, we suggest the filtering step as used in `sclink_norm` to select genes. This also reduces the computational burden. 
# However, if the users would like to infer gene networks for a large gene list (e.g., > 5000 genes), 
# we suggest that the users first use `sclink_cor` to investigate the correlation structures among these genes.
corr = sclink_cor(expr = count.norm, ncores = 1) #500 genes
# save time and increase overall accuracy.