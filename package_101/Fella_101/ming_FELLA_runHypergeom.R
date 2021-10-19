library(FELLA)
library(org.Dm.eg.db)
library(KEGGREST)
library(igraph)
library(magrittr)

mode='tradeoff'
#mode='insulin'
#mode='fecundity'
##################################################
## read in signaficant metabolite data
input.file=paste0(mode,'_candidate.mz.csv')
outfile1=paste0(mode,'_fella.out_diffusion.csv')
outfile2=paste0(mode,'_fella.out_hypergeom.csv')
#df=data.table::fread('tradeoff_candidate.mz.csv')
df=data.table::fread(input.file)
mz=df$metabolite

## read in assay background
library(readxl)
mz.info=read_excel('./2020-12-08_xxxx-60_Data.xlsx',sheet='Metabolite Information')
dim(mz.info) #361
sum(mz.info$`KEGG ID`=='N/A') #29 measure metabolites do not have kegg id
mz.info=mz.info[mz.info$`KEGG ID`!='N/A',]
dim(mz.info) #332
tmp=mz.info[grep('\\/',mz.info$`KEGG ID`),] #1 with multiple id
tmp #8
x=strsplit(tmp$`KEGG ID`,'\\/')
x
y=sapply(x,function(i){ i[i!='NA'][1] } )
mz.info[grep('\\/',mz.info$`KEGG ID`),]$`KEGG ID`=y
dim(mz.info)#332   6

## map sig mz to bacground
length(mz) #39
sum(mz %in% mz.info$`Current MS Compounds`) #36
mz.id=mz.info[mz.info$`Current MS Compounds` %in% mz,]
dim(mz.id) #36

#########################################################
# build the KEGG graph
#https://rdrr.io/bioc/FELLA/src/R/buildGraphFromKEGGREST.R
if(!file.exists('dem_graph.rds')){
  graph<-buildGraphFromKEGGREST(organism='dme')
  saveRDS(graph,'dem_graph.rds')
}
graph=readRDS('dem_graph.rds')#10558nodes and 34285edges

tmpdir<-'FELLA_database_hyper'
if(!dir.exists(tmpdir)){
  buildDataFromGraph(keggdata.graph=graph,databaseDir = tmpdir,
                     internalDir=FALSE,matrices='hypergeom',
                     normality='diffusion',niter = 50)
}

######################################################################################
## play with buildDataFromGraph: https://rdrr.io/bioc/FELLA/src/R/buildDataFromGraph.R
graph=readRDS('dem_graph.rds')
#buildDataFromGraph(keggdata.graph=graph,databaseDir = tmpdir,
#                     internalDir=FALSE,matrices='hypergeom',
#                     normality='diffusion',niter = 50)

keggdata.graph=graph 
E(graph)$weight <- 1/E(graph)$weight

id.pathway <- which(V(graph)$com == 1) #1,2,3,4,5 5 levels network, top1,pathway
id.compound <- which(V(graph)$com == 5)#botom 5, compound
length(id.pathway) #131 pathways
length(id.compound) #4018 compounds

vcount(graph)#10558 nodes
ecount(graph) #34285 edges

x= shortest.paths(
  graph = graph, 
  v = id.compound, 
  to = id.pathway, 
  mode = "out")
dim(x) #4081 x 131
#can one compound 'walk' to a pathway in 4 steps, as compound level5, pathway level1, 5->1, 4 steps

hypergeom.matrix <- shortest.paths(
  graph = graph, 
  v = id.compound, 
  to = id.pathway, 
  mode = "out") == 4

library(Matrix)
hypergeom.matrix <- Matrix(data = hypergeom.matrix, sparse = TRUE)
# save a sparse matrix

sum(hypergeom.matrix[,20]) #167 mz in total
colnames(hypergeom.matrix)[20] #"dme00270"

#########################################################
## runHypergeom

tmpdir<-'FELLA_database_hyper/'

fella.data <- loadKEGGdata(
  databaseDir = tmpdir,
  internalDir = FALSE,
  loadMatrix='hypergeom')
fella.data
cat(getInfo(fella.data))
fella.data@hypergeom #big matrix
colnames(fella.data@hypergeom@matrix)[20] #"dme00270"
sum(fella.data@hypergeom@matrix[,20]) #167 for "dme00270"

# map metabolites to the KEGG graph
analysis.mz <- defineCompounds(
  compounds = mz.id$`KEGG ID`,
  compoundsBackground=mz.info$`KEGG ID`, #use targeted assay as the background
  data = fella.data)

length(getInput(analysis.mz)) #33
length(getExcluded(analysis.mz)) #3
analysis.mz


analysis.mz <- runHypergeom(
  object = analysis.mz,
  data = fella.data,
  p.adjust = 'fdr')
analysis.mz
# top pathways and affected compounds
plot(analysis.mz,method='hypergeom',data=fella.data,
     main='analysis',
     plotLegend = FALSE)

# top k p.scores
tab.all <- generateResultsTable(
  method = "hypergeom",
  object = analysis.mz,
  data = fella.data)
dim(tab.all) #

myGraph <- generateResultsGraph(
  object = analysis.mz, 
  method = "hypergeom", 
  data = fella.data)
show(myGraph)

# details from the enzymes reported among the top k KEGG entries
# not applicable to method = "hypergeom" 
#tab.enzyme <- generateEnzymesTable(method = "hypergeom",object = analysis.mz,data = fella.data)

exportResults(
  format = "csv",
  #file = 'fella_out_hypergeom.csv',
  file = outfile2,
  method = "hypergeom",
  object = analysis.mz,
  data = fella.data)

#dme00480  Glutathione metabolism
#dme00270  Cysteine and methionine metabolism

###################################################################
## extract mapped mz to each sig pathway
sig.path=data.table::fread('./fella_out_hypergeom.csv')
sig.path$KEGG.id

## get KEGG database
getMatrix <- function(data, method) {
  return(slot(data, method)@matrix)
}
hypergeom.matrix=getMatrix(fella.data,'hypergeom')

object=analysis.mz
metabolites.input <- getInput(object)
metabolites.input.intersect <- intersect(
  metabolites.input, 
  rownames(hypergeom.matrix))

row_comp <- which(rownames(hypergeom.matrix) %in% metabolites.input)
length(row_comp) #33 mz hits

col_comp=which(colnames(hypergeom.matrix) %in% sig.path$KEGG.id)
out=list();
for(path in col_comp){
  #hypergeom.matrix[row_comp, path]
  id=names(which(hypergeom.matrix[row_comp, path])) #matched mz hits in this pathway
  x=mz.id[mz.id$`KEGG ID` %in% id,]
  
  kegg.id=colnames(hypergeom.matrix)[path]
  kegg.name=sig.path[sig.path$KEGG.id==kegg.id,]$KEGG.name
  x$mapped.pathway.id=kegg.id
  x$mapped.pathway.name=kegg.name
  out[[kegg.id]]=x
}

df.out=Reduce(`rbind`,out)
head(df.out)
table(df.out$mapped.pathway.id)
data.table::fwrite(df.out,'sig.mz_to_pathway.csv')

###################################################################
###################################################################
## play with FELLA: https://rdrr.io/bioc/FELLA/src/R/runHypergeom.R
object=analysis.mz
data=fella.data
p.adjust='fdr'
#https://rdrr.io/bioc/FELLA/src/R/get-.R
getMatrix <- function(data, method) {return(slot(data, method)@matrix)}
hypergeom.matrix=getMatrix(data,'hypergeom')

which(colnames(hypergeom.matrix)=='dme00270')
hypergeom.matrix[,20]
sum(hypergeom.matrix[,20]) #167 mz in total

metabolites.input <- getInput(object)
metabolites.input.intersect <- intersect(
  metabolites.input, 
  rownames(hypergeom.matrix))
length(metabolites.input.intersect) #33
# I didn't specify bacground, so it's 0
metabolites.background.intersect <- intersect(
  getBackground(object), 
  rownames(hypergeom.matrix))
length(metabolites.background.intersect) #0

#   hypergeom.matrix <- hypergeom.matrix[, colSums(hypergeom.matrix) > 0]
dim(hypergeom.matrix) #4018  131, mz by pathways
row_comp <- which(rownames(hypergeom.matrix) %in% metabolites.input)
length(row_comp) #33 mz hits

pvalues.path <- vector("double", length = dim(hypergeom.matrix)[2])
names(pvalues.path) <- colnames(hypergeom.matrix)
length(pvalues.path) #131 pathways

# Additional variables for the report
pathbackground <- pvalues.path
pathhits <- pvalues.path

# p-values calculation
colnames(hypergeom.matrix)[20] #"dme00270"
path=20
#for (path in seq_len(ncol(hypergeom.matrix))) {
  hypergeom.matrix[row_comp, path]
  names(which(hypergeom.matrix[row_comp, path])) #matched mz hits in this pathway
  
  sample_success <- sum(hypergeom.matrix[row_comp, path])
  pathhits[path] <- sample_success
  
  # sum over row or column (current)?
  total_success <- sum(hypergeom.matrix[, path])
  pathbackground[path] <- total_success
  
  total_failure <- dim(hypergeom.matrix)[1] - total_success
  sample_size <- length(row_comp)
  
  pvalues.path[path] <- stats::phyper(
    sample_success - 1, 
    total_success, 
    total_failure, 
    sample_size, 
    lower.tail = FALSE)
#}

