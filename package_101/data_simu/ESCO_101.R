
#https://github.com/JINJINT/ESCO/blob/master/vignettes/esco.Rmd

if(F){
  Sys.getenv("GITHUB_PAT")
  Sys.unsetenv("GITHUB_PAT")
  Sys.getenv("GITHUB_PAT")
  
  BiocManager::install('snow',force=T)
  
  library("devtools")
  devtools::install_github("JINJINT/ESCO")
}

library(ESCO)

sim <- escoSimulate(nGenes = 100, nCells  = 50)
sim
names(assays(sim))

# Access the counts
data = assays(sim)$TrueCounts
data[1:5, 1:5]
# Information about genes
head(rowData(sim))
# Information about cells
head(colData(sim))
# Information about paramters
str(metadata(sim)$Params)

#The main part of this object is a series of features
#by samples matrix containing the simulated true counts (accessed using `assays(sim)$TrueCounts`); 
#counts with zero-inflation noise; (accessed using `assays(sim)$counts`); 
#counts with downsample noise ((accessed using `assays(sim)$observedcounts`). 
#By default, all three kinds of data matrix will be simulated, 
#while users are allowed to simulated only one or two of them for time/space saving, via specifying the ``dropout.type`` in parameters setting.


                                                                                                                                                                                                         