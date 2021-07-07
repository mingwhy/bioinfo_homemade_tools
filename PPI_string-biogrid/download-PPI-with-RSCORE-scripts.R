library(igraph)
##################################################
## store gene module identified in each cell type
## then use AUCell to score each gene set in each cell
## then you have a matrix, row=gene.set, col=cell
## depend on the gene.set.activity acorss cells to define
## basic gene modules and specific gene modules
if(F){
  source("./getPPI_String.R") #in RSCROE github package
#  PPI = 'String'; species = 7227; score_threshold=600;
  PPI = 'String'; species = 7227; score_threshold=400;
  fly_network_matrix <- getPPI_String(species = species,score_threshold = score_threshold)
  dim(fly_network_matrix)
  fly_network_matrix[10:14,10:14]
  rowSums(fly_network_matrix)
  saveRDS(fly_network_matrix,'fly_network_matrix.string_400.rds');
  #saveRDS(fly_network_matrix,'fly_network_matrix.string_600.rds');
}

fly_network_matrix=readRDS('fly_network_matrix.string_600.rds');
dim(fly_network_matrix)

if(F){
  source("./getPPI_Biogrid.R") #in RSCROE github package
  PPI = 'Biogrid'; species = 7227;
  fly_network_matrix <- getPPI_Biogrid(species = species)
  dim(fly_network_matrix)
  fly_network_matrix[10:14,10:14]
  rowSums(fly_network_matrix)
  saveRDS(fly_network_matrix,'fly_network_matrix.biogrid.rds');
  #saveRDS(fly_network_matrix,'fly_network_matrix.string_600.rds');
}

fly_network_matrix=readRDS('fly_network_matrix.biogrid.rds');
dim(fly_network_matrix)


