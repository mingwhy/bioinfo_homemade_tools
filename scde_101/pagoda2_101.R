# http://pklab.med.harvard.edu/peterk/p2/walkthrough.nb.html
# PBMC8K dataset download from: https://support.10xgenomics.com/single-cell-gene-expression/datasets/2.1.0/pbmc8k?

library(pagoda2)
cd <- read.10x.matrices('../Documents/single.cell_datasets/pbmc_dataset/raw_gene_bc_matrices/GRCh38/',
                        version = "V2") #version V2: genes.tsv. V3:features.tsv

str(cd)

par(mfrow=c(1,2), mar = c(3.5,3.5,2.0,0.5), mgp = c(2,0.65,0), cex = 1.0)
hist(log10(colSums(cd)+1),main='molecules per cell',col='wheat',xlab='log10[ molecules per cell]')
hist(log10(rowSums(cd)+1),main='molecules per gene',col='wheat',xlab='log10[ molecules per gene]')

#cell filter: require at least 500 molecules per cell
counts <- gene.vs.molecule.cell.filter(cd,min.cell.size=500)
str(counts)

# gene filter: omit low-expressed genes
hist(log10(rowSums(counts)+1),main='Molecules per gene',xlab='molecules (log10)',col='wheat')
abline(v=1,lty=2,col=2)
counts <- counts[rowSums(counts)>=10,]
str(counts)

# get a lean count matrix. create pagoda2 object
r <- Pagoda2$new(counts,log.scale=TRUE)
#Error in setCountMatrix(x, min.cells.per.gene = min.cells.per.gene, trim = trim,  : 
#duplicate gene names are not allowed - please reduce

#make gene names unique:
rownames(counts) <- make.unique(rownames(counts))
r <- Pagoda2$new(counts,log.scale=TRUE, n.cores=2)

# adjust the variance, normlize the extent to which genes with different 
# expression magitudes will contribute to to the downstream analysis
r$adjustVariance(plot=T,gam.k=10)

# Below we’ll use the simplest, default scenario, where we first reduce the dataset dimensions by running PCA, and then move into k-nearest neighbor graph space for clustering and visualization calculations. 
# First, the PCA reduction:
r$calculatePcaReduction(nPcs=50,n.odgenes=3e3)
# The next few steps will make kNN graph, find clusters and generate a quick largeVis embedding to visualize the subpopulations:
r$makeKnnGraph(k=40,type='PCA',center=T,distance='cosine')

r$getKnnClusters(method=infomap.community,type='PCA')
M <- 30; r$getEmbedding(type='PCA',M=M,perplexity=30,gamma=1/M,alpha=1)
#Now we can visualize the embedding using the determined clusters:
r$plotEmbedding(type='PCA',show.legend=F,mark.clusters=T,min.group.size=50,shuffle.colors=F,mark.cluster.cex=1,alpha=0.1,main='clusters (largeVis)')

#We can use the same plotEmbedding() function to show all kinds of other values. 
#For instance, let’s look at depth, or an expresson pattern of a gene:
str(r$depth)
par(mfrow=c(1,2))
r$plotEmbedding(type='PCA',#embeddingType='tSNE',
                colors=r$depth,shuffle.colors=F,mark.cluster.cex=1,alpha=0.1,main='depth')
# or look at one gene expression abundance
gene <-"LYZ"
r$plotEmbedding(type='PCA',#embeddingType='tSNE',
                colors=r$counts[,gene],shuffle.colors=F,mark.cluster.cex=1,alpha=0.1,main=gene)

#############################################################
## Pathway overdispersion analysis (a la PAGODA1)
# First, build GO->gene environment:
suppressMessages(library(org.Hs.eg.db))
# translate gene names to ids
ids <- unlist(lapply(mget(colnames(r$counts),org.Hs.egALIAS2EG,ifnotfound=NA),function(x) x[1]))
# reverse map
rids <- names(ids); names(rids) <- ids;
# list all the ids per GO category
go.env <- list2env(eapply(org.Hs.egGO2ALLEGS,function(x) as.character(na.omit(rids[x]))))
class(go.env) #"environment"

# Now run overdispersion anlaysis
r$testPathwayOverdispersion(go.env,verbose=T,correlation.distance.threshold=0.95,recalculate.pca=F,top.aspects=15)

# We’ll use hierarchical differential expression results instead:
r$getHierarchicalDiffExpressionAspects(type='PCA',clusterName='community',z.threshold=3)

#We’ll make an app with that, ordering the “differential expression aspects” explicitly (otherwise if row clustering is omitted they’ll be clustered by similarity)
app <- p2.make.pagoda1.app(r,inner.clustering=TRUE,embeddingType='tSNE',clusterType='multilevel',min.group.size=50,row.clustering=list(order=rev(1:nrow(r$misc$pathwayOD$xv))))
#Show app:
show.app(app,'pbmc',browse=T)
