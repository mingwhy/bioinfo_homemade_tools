#https://github.com/JINJINT/ESCO/blob/master/vignettes/esco.Rmd

if(F){
  Sys.getenv("GITHUB_PAT")
  Sys.unsetenv("GITHUB_PAT")
  Sys.getenv("GITHUB_PAT")
  
  BiocManager::install('snow',force=T)
  
  library("devtools")
  devtools::install_github("JINJINT/ESCO")
}

#The main part of this object is a series of features
#assays(sim)$TrueCounts: by samples matrix containing the simulated true counts; 
#assays(sim)$counts: counts with zero-inflation noise; 
#assays(sim)$observedcounts: counts with downsample noise 
#By default, all three kinds of data matrix will be simulated, 
#while users are allowed to simulated only one or two of them for time/space saving, 
#via specifying the ``dropout.type`` in parameters setting.


library(ESCO)

# in my case, there is some error when running 'makeCluster'
# I took suggestion from: https://stackoverflow.com/questions/62730783/error-in-makepsockclusternames-spec-cluster-setup-failed-3-of-3-work
# to modify my ~/.Rprofile and solved this issue.

#===== start simulation ======#
sim <- escoSimulate(nGenes = 100, nCells = 50, 
                          withcorr = TRUE,numCores=2,
                          verbose = FALSE)

#===== access the data ======#
datalist = list("simulated truth"=assays(sim)$TrueCounts,
                "zero-inflated" = assays(sim)$counts, 
                "down-sampled" = assays(sim)$observedcounts)
sapply(datalist,dim)
#====== plot the data ======#
heatdata(datalist, norm = FALSE, size = 2, ncol = 3)

#====== plot the Gene correlation ======#
# object that saved all simulation configurations
simparams = metadata(sim)$Params 

# object that particularly saved the correlation structure
rholist = slot(simparams,"corr") 
length(rholist)
dim(rholist[[1]]); #50by50

# arrange the true correlation and simulated correlation
corrgenes = rownames(rholist[[1]])
gcnlist = lapply(datalist, function(data)gcn(data, genes = corrgenes))
gcnlist = append(gcnlist, list("given truth" = rholist[[1]]), 0)
heatgcn(gcnlist, size = 3, ncol = 4)

###############################################
# from PNAS 2021-Constructing Local Cell Sepcific Networks from Single Cell Data
# supp note info
#https://github.com/xuranw/locCSN/blob/main/DataStore/ESCO/ESCO_Simulation_Example.md
corr.single = list()
corr.single[[1]] = diag(75)
#diag(corr.single[[1]]) = 1
#corr.single[[1]]
corr.single[[1]][1:15, 1:15] <- corr.single[[1]][16:30, 16:30] <- corr.single[[1]][31:45, 31:45] <- 
  corr.single[[1]][46:60, 46:60] <- corr.single[[1]][61:75, 61:75]  <- 0.9
corr.single[[1]][1:15, 16:30] <- corr.single[[1]][16:30, 31:45] <- corr.single[[1]][31:45, 46:60] <-
  corr.single[[1]][46:60, 61:75] <- corr.single[[1]][16:30, 1:15] <- corr.single[[1]][31:45, 16:30] <- 
  corr.single[[1]][46:60, 31:45] <- corr.single[[1]][61:75, 46:60] <- 0.7
corr.single[[1]][1:15, 31:45] <- corr.single[[1]][16:30, 46:60] <- corr.single[[1]][31:45, 61:75] <-
  corr.single[[1]][31:45, 1:15] <- corr.single[[1]][46:60, 16:30] <- 
  corr.single[[1]][61:75, 31:45] <- 0.5
corr.single[[1]][1:15, 46:60] <- corr.single[[1]][16:30, 61:75] <- corr.single[[1]][46:60, 1:15] <-
  corr.single[[1]][61:75, 16:30] <- 0.3
corr.single[[1]][1:15, 61:75] <- corr.single[[1]][61:75, 1:15] <- 0.1
diag(corr.single[[1]]) = 1

alpha = 0.1; lib.loc = 9;
sim.single <- escoSimulateSingle(nGenes = 100, nCells = 200,  withcorr = TRUE, corr = corr.single, verbose = TRUE)
corr.single = metadata(sim.single)$Params@corr

gene.order = colnames(corr.single[[1]])[hclust(dist(corr.single[[1]]))$order]
corr.single[[1]][gene.order[seq(1, 75, 15)], gene.order[seq(1, 75, 15)]]
gene.order = gene.order[c( 1:15, 16:30, 31:45, 61:75, 46:60)]
image(corr.single[[1]][gene.order, gene.order])

gene.order = c(gene.order, setdiff(paste0('Gene', 1:100), gene.order))

# get the data
datalist = list("simulated truth" = assays(sim.single)$TrueCounts, 
                "zero-inflated" = assays(sim.single)$counts, 
                "down-sampled" = assays(sim.single)$observedcounts)


#log.cpm = log(apply(datalist$`down-sampled`, 2, function(x){x/sum(x)})*10^6+1)
log.cpm = log(apply(datalist$`simulated truth`, 2, function(x){x/sum(x)})*10^6 + 1)
log.cpm = log.cpm[gene.order, ]
write.table(log.cpm, file = paste0(dir.save, 'simsingle_logcpm_trail.txt'))

