
########################################
raw=readRDS('genotype_old.rds')
length(raw) #the total number of chunks
sapply(raw,length) #each chunk contain 17 elements
lapply(raw,'[[',1) #the name of each chunk
sapply(raw[[1]][3][[1]],dim)

########################################
## ming's code below (don't run)
library(reticulate)
np<-import('numpy') 
(files=Sys.glob('*npy'))
for(file in files){
  #raw<-np$load('genotype_young.npy',allow_pickle = 'True')
  raw<-np$load(file,allow_pickle = 'True')
  output.file=gsub("npy",'rds',file)
  cat('read in',file,'save in',output.file,'\n')
  dim(raw)
  new.list=list();
  z=0;
  for(i in 1:length(raw)){
    x=raw[[i]]
    n=length(x)/17
    #cat(i,length(x),n,'\n')
    for(j in 1:n){
      start=(j-1)*17+1
      z=z+1;
      new.list[[z]]=x[start : (start+16)] 
    }
  }
  length(new.list)
  sapply(new.list,length)
  #lapply(1:length(new.list),function(i){print(new.list[[i]][1])} )
  #saveRDS('genotype_young.rds')   
  saveRDS(new.list,file=output.file)
}

########################################
## save npy or npz files 
#https://cran.r-project.org/web/packages/RcppCNPy/vignettes/UsingReticulate.pdf
tfile <- tempfile(fileext=".npy")
set.seed(42)
m <- matrix(sort(rnorm(6)), 3, 2)
m
np$save(tfile, m)

x <- seq(1, 10)
y <- sin(x)
np$savez("file1.npz", x, y)
np$savez("file2.npz", x=x, y=y)

