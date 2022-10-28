## NMF (call python in R, key points: `Sys.setenv` for load modules, `.->$`, as.interger() or 'xxx' to pass parameters)
library(reticulate)
#Sys.setenv(RETICULATE_PYTHON = "/Users/mingyang/anaconda3/envs/myenv/bin/python")
Sys.setenv(RETICULATE_PYTHON ='/Users/mingyang/anaconda3/bin')
py_config()
sys=import('sys')
sys$path
nimfa=import('nimfa') #open terminal, if `import nimfa` successfully, use `import sys;sys.path` to get the env set as above

#https://nimfa.biolab.si/
V = nimfa$examples$medulloblastoma$read(normalize='True')
lsnmf = nimfa$Lsnmf(V, seed='random_vcol', rank=as.integer(50), max_iter=as.integer(100))
lsnmf_fit = lsnmf()

cat('Rss:' , lsnmf_fit$fit$rss(),'\n')
cat('Evar: ', lsnmf_fit$fit$evar(),'\n')
cat('K-L divergence:' , lsnmf_fit$distance(metric='kl'),'\n')
cat('Sparseness, W: , H: ' , unlist(lsnmf_fit$fit$sparseness()),'\n')

