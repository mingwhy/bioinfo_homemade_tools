
## NMF (call python in R, key points: `Sys.setenv` for load modules, `.->$`, as.interger() or 'xxx' to pass parameters)
library(reticulate)
#Sys.setenv(RETICULATE_PYTHON = "/Users/mingyang/anaconda3/envs/myenv/bin/python")
Sys.setenv(RETICULATE_PYTHON ='/Users/mingyang/anaconda3/bin')
py_config()
sys=import('sys')
sys$path
nimfa=import('nimfa') #open terminal, if `import nimfa` successfully, use `import sys;sys.path` to get the env set as above
np=import('numpy')

#https://nimfa.biolab.si/
V = nimfa$examples$medulloblastoma$read(normalize='True')
lsnmf = nimfa$Lsnmf(V, seed='random_vcol', rank=as.integer(50), max_iter=as.integer(100))
lsnmf_fit = lsnmf()

cat('Iterations:' ,lsnmf_fit$n_iter,'\n')
cat('Rss:' , lsnmf_fit$fit$rss(),'\n')
cat('Evar: ', lsnmf_fit$fit$evar(),'\n')
cat('K-L divergence:' , lsnmf_fit$distance(metric='kl'),'\n')
cat('Sparseness, W: , H: ' , unlist(lsnmf_fit$fit$sparseness()),'\n')
#https://github.com/mims-harvard/nimfa/issues/54

W = lsnmf_fit$basis()
dim(W) #5839 x 50, Basis matrix
H = lsnmf_fit$coef()
dim(H) #50 x 34, Mixture matrix

target.mat=np$dot(W,H) #Target estimate
dim(target.mat) #5893   34

####################################################
#big dataset test
dim(input.mat) #2000 gene X 47898 cells
V = input.mat

#LSNMF on numpy dense matrix with quality and performance measures.
lsnmf = nimfa$Lsnmf(V, seed='random_vcol', rank=as.integer(50), max_iter=as.integer(100))
start.time=Sys.time()
lsnmf_fit = lsnmf() 
end.time=Sys.time()
end.time-start.time; # 2000 x 47898, Time difference of 27.97791 mins

cat('Iterations:' ,lsnmf_fit$n_iter,'\n')
cat('Rss:' , lsnmf_fit$fit$rss(),'\n')
cat('Evar: ', lsnmf_fit$fit$evar(),'\n')
cat('K-L divergence:' , lsnmf_fit$distance(metric='kl'),'\n')
cat('Sparseness, W: , H: ' , unlist(lsnmf_fit$fit$sparseness()),'\n')
#https://github.com/mims-harvard/nimfa/issues/54


#Standard NMF - Divergence on scipy.sparse matrix with matrix factors estimation.
V = input.mat
lsnmf = nimfa$Nmf(V, seed='random_vcol', rank=as.integer(50), max_iter=as.integer(100),
                  update='euclidean', objective='fro')
start.time=Sys.time()
lsnmf_fit = lsnmf() 
end.time=Sys.time()
end.time-start.time; # 2000 x 47898, Time difference of 49.1204 mins


