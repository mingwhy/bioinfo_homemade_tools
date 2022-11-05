use python inside R via reticulate

If you can run some python sucessfully in ipython terminal, and want to use the same python env inside R,
you may want to type:
```import sys
sys.path
```
inside ipython to locate the python engine.

Then go to Rstudio:
```
library(reticulate)
Sys.setenv(RETICULATE_PYTHON = "/Users/mingyang/anaconda3/envs/RCFGL/bin/python")
py_config()
pd=import('pandas')
np=import('numpy')
sc=import('scanpy') #py_install('scanpy')
sce=import('scanpy.external')
tqdm=import('tqdm')
sys=import('sys')

# https://gitlab.com/olgaibanez/decibel
sys$path=c(sys$path,'./decibel/module/') #sys.path.append('./decibel/module/')
sys$path 
dcb=import('decibel') 
dcb
dcb$enge_euclidean_dist #make sure this function works
```
to use the same python env as in the ipython terminal.
