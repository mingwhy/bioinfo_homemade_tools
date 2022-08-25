
library(zellkonverter)
library(SingleCellExperiment)
library(tidyverse)
library(ggExtra) #
library(scRNAseq)
library(ggplot2);theme_set(theme_classic())
library(scran)
library(ggpubr)
library(Seurat)
library(grDevices);library(RColorBrewer)
library(SingleCellExperiment)
plotcol <- brewer.pal(8,'Dark2')
library(glmpca)
library(scry)

#out=readRDS('../0619_DV_normcounts/sce.shared.rds')
#out=readRDS('../0619_DV_normcounts/sce.filtered.rds')
out=readRDS('../0823_subsample.umi/sce_minCellCounts.rds')


######################################################################
library(reticulate)
Sys.setenv(RETICULATE_PYTHON = "/Users/mingyang/anaconda3/envs/myenv/bin/python")
py_config()
sys=import('sys')
sys$path

sys$path=c(sys$path,'./') #sys.path.append('./decibel/module/')
sys$path 
gcl_lib=import('gcl_library') 
gcl_lib$jackknife

py_run_string('from pathlib import Path')
py_run_string('import sys')
py_run_string('from pathlib import Path')
py_run_string('from os import path')
# for threading:
py_run_string('import threading')
# math tools:
np=import('numpy')
py_run_string('from scipy.spatial.distance import cdist')
py_run_string('import math')
py_run_string('import random')
# gcl library import:
py_run_string('from matplotlib import pyplot as plt')


files_arr = c('old.csv', 'young.csv');
jack_knifes=as.integer(70)
num_divisions=as.integer(10)
jack_knife_percentage=0.8 

csv_mat = np$genfromtxt(files_arr[1], delimiter=',')
dim(csv_mat) #(3000, 202), 3000 gene by 202 cell
np$sum(csv_mat,axis=as.integer(1)) # rowSum
apply(csv_mat,1,sum)

np$sum(csv_mat,axis=as.integer(0)) # colSum
apply(csv_mat,2,sum)

res_old=gcl_lib$jackknife(csv_mat, jack_knifes, jack_knife_percentage, num_divisions)

csv_mat = np$genfromtxt(files_arr[2], delimiter=',')
res_young=gcl_lib$jackknife(csv_mat, jack_knifes, jack_knife_percentage, num_divisions)

df.res=data.frame(gcl=c(unlist(res_old),unlist(res_young)),age=rep(c('old','young'),each=70))
library(ggplot2)
head(df.res)
ggplot(df.res,aes(x=age,y=gcl))+geom_violin()+geom_jitter(size=0.2)+theme_classic()
