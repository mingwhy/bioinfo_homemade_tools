Sys.setenv(RETICULATE_PYTHON = "/Users/mingyang/anaconda3/envs/RCFGL/bin/python")
library(reticulate)
py_config()
library(Matrix)

mat=matrix(c(0.1,0.9,1,0.0),ncol=2,byrow = T)
mat
CountsMatrix=Matrix(mat,sparse=TRUE)
CountsMatrix

py_run_file("./HyperEigen.py")
as.numeric(py$large_val)

##
mat=matrix(c(0.1,0.9,1,0.0),ncol=2,byrow = T)
mat
inc_mat=mat;
deg_mat= diag(apply(mat^2,1,sum))
lap_mat=solve(deg_mat) %*% inc_mat %*% t(inc_mat)
lap_mat
eigen(lap_mat)

lap_mat= t(inc_mat) %*% solve(deg_mat) %*% inc_mat 
lap_mat
eigen(lap_mat)



