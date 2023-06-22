
#https://hms-dbmi.github.io/scde/pagoda.html
library(scde)

data(pollen)
# remove poor cells and genes
cd <- clean.counts(pollen)
# check the final dimensions of the read count matrix
dim(cd) #11310 gene by 64 cells
str(cd)

x <- gsub("^Hi_(.*)_.*", "\\1", colnames(cd))
l2cols <- c("coral4", "olivedrab3", "skyblue2", "slateblue3")[as.integer(factor(x, levels = c("NPC", "GW16", "GW21", "GW21+3")))]
l2cols

# Find error models
#Note this step takes a considerable amount of time unless multiple cores are used. We highly recommend use of multiple cores. You can check the number of available cores available using detectCores().
if(F){knn <- knn.error.models(cd, k = ncol(cd)/4, n.cores = 2, min.count.threshold = 2, min.nonfailed = 5, max.model.plots = 10)}
data(knn)

class(knn) #data.frame
dim(knn) #64 cells x 12 attributes
colnames(knn)

# Normlizing variance
#In order to accurately quantify excess variance or overdispersion, we must normalize out expected levels of technical and intrinsic biological noise. 
varinfo <- pagoda.varnorm(knn, counts = cd, trim = 3/ncol(cd), max.adj.var = 5, n.cores = 1, plot = TRUE)
class(varinfo); #list
names(varinfo)
## list top overdispersed genes
sort(varinfo$arv, decreasing = TRUE)[1:10]

#Controlling for sequencing depth
#Below we will control for the gene coverage (estimated as a number of genes with non-zero magnitude per cell) and normalize out that aspect of cell heterogeneity:
varinfo <- pagoda.subtract.aspect(varinfo, colSums(cd[, rownames(knn)]>0))

#####################################################
#Evaluate overdispersion of pre-defined gene sets
library(org.Hs.eg.db)
# translate gene names to ids
ids <- unlist(lapply(mget(rownames(cd), org.Hs.egALIAS2EG, ifnotfound = NA), function(x) x[1]))
rids <- names(ids); names(rids) <- ids 
# convert GO lists from ids to gene names
gos.interest <- unique(c(ls(org.Hs.egGO2ALLEGS)[1:100],"GO:0022008","GO:0048699", "GO:0000280")) 
go.env <- lapply(mget(gos.interest, org.Hs.egGO2ALLEGS), function(x) as.character(na.omit(rids[x]))) 
go.env <- clean.gos(go.env) # remove GOs with too few or too many genes
go.env <- list2env(go.env) # convert to an environment

#Now, we can calculate weighted first principal component magnitudes for each GO gene set in the provided environment.
pwpca <- pagoda.pathway.wPCA(varinfo, go.env, n.components = 1, n.cores = 1)
length(pwpca) #52
names(pwpca[[1]])

#We can now evaluate the statistical significance of the observed overdispersion for each GO gene set.
df <- pagoda.top.aspects(pwpca, return.table = TRUE, plot = TRUE, z.score = 1.96)
head(df)
dim(df) #8 rows

#####################################################
# Evaluate overdispersion of 'de novo' gene sets
if(file.exists('clpca.rds')){ #make take some time
  clpca <- pagoda.gene.clusters(varinfo, trim = 7.1/ncol(varinfo$mat), n.clusters = 50, n.cores = 1, plot = TRUE);
  saveRDS(clpca,'clpca.rds')
}
clpca=readRDS('clpca.rds')
df <- pagoda.top.aspects(pwpca, clpca, return.table = TRUE, plot = TRUE, z.score = 1.96)
head(df)

#####################################################
# Visualize significant aspects of heterogeneity
# get full info on the top aspects
#tam <- pagoda.top.aspects(pwpca, clpca, n.cells = NULL, z.score = qnorm(0.01/2, lower.tail = FALSE))
tam <- pagoda.top.aspects(pwpca, n.cells = NULL, z.score = qnorm(0.01/2, lower.tail = FALSE))
str(tam)
# List of 4
# $xv: 8 GO term (new coordinates) by 64 cells
# $xvw: 8 GO term (new coordinates) by 64 cells
# $gw: 1:269, gene names
# $df, data.frame, 52GO terms of 8 attributes(name, npc, n,score,z,adj.z,sh.z,adj.sh.z)
summary(tam$xv[1,]) #negative and positive
summary(tam$xvw[1,]) #only positive
sum(tam$xvw[1,]) #I guess it's a weighing factor, as it's normlized to 1
# determine overall cell clustering
hc <- pagoda.cluster.cells(tam, varinfo) 
hc #64 cells

#Next, we will reduce redundant aspects in two steps. First we will combine pathways that are driven by the same sets of genes:
tamr <- pagoda.reduce.loading.redundancy(tam, pwpca, clpca)

#In the second step we will combine aspects that show similar patterns (i.e. separate the same sets of cells). Here we will plot the cells using the overall cell clustering determined above:
tamr2 <- pagoda.reduce.redundancy(tamr, distance.threshold = 0.9, plot = TRUE, cell.clustering = hc, labRow = NA, labCol = NA, box = TRUE, margins = c(0.5, 0.5), trim = 0)

col.cols <- rbind(groups = cutree(hc, 3))
pagoda.view.aspects(tamr2, cell.clustering = hc, box = TRUE, labCol = NA, margins = c(0.5, 20), col.cols = rbind(l2cols))


#Similar views can be obtained in the R session itself. For instance, here we'll view top 10 genes associated with the top two pathways in the neurogenesis cluster: "neurogenesis" (GO:0022008) and "generation of neurons" (GO:0048699)
pagoda.show.pathways(c("GO:0022008","GO:0048699"), varinfo, go.env, cell.clustering = hc, margins = c(1,5), show.cell.dendrogram = TRUE, showRowLabels = TRUE, showPC = TRUE)

#Adding 2D embedding
library(Rtsne);
# recalculate clustering distance .. we'll need to specify return.details=T
cell.clustering <- pagoda.cluster.cells(tam,varinfo,include.aspects=TRUE,verbose=TRUE,return.details=T)
# fix the seed to ensure reproducible results
set.seed(0); 
tSNE.pagoda <- Rtsne(cell.clustering$distance,is_distance=T,initial_dims=100,perplexity=10)
par(mfrow=c(1,1), mar = c(2.5,2.5,2.0,0.5), mgp = c(2,0.65,0), cex = 1.0);
plot(tSNE.pagoda$Y,col=adjustcolor(col.cols,alpha=0.5),cex=1,pch=19,xlab="",ylab="")

#Controlling for undesired aspects of heterogeneity
# get cell cycle signature and view the top genes
#cc.pattern <- pagoda.show.pathways(c("GO:0000280", "GO:0007067"), varinfo, go.env, show.cell.dendrogram = TRUE, cell.clustering = hc, showRowLabels = TRUE)
cc.pattern <- pagoda.show.pathways(c("GO:0000280"), varinfo, go.env, show.cell.dendrogram = TRUE, cell.clustering = hc, showRowLabels = TRUE)
# subtract the pattern
varinfo.cc <- pagoda.subtract.aspect(varinfo, cc.pattern)

#######################################################
#Creating custom pathway annotations or gene sets
#http://hms-dbmi.github.io/scde/genesets.html
