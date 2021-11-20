
## this blog is inspired by 
#https://www.nxn.se/valent/2017/11/16/droplet-scrna-seq-is-not-zero-inflated
#https://divingintogeneticsandgenomics.rbind.io/post/negative-binomial-distribution-in-scrnaseq/
#https://divingintogeneticsandgenomics.rbind.io/post/negative-bionomial-distribution-in-single-cell-rnaseq/
#https://github.com/willtownes/scrna2019/blob/master/real/svensson_2019/01_exploratory.Rmd
#https://divingintogeneticsandgenomics.rbind.io/post/modeling-single-cell-rnaseq-data-with-multinomial-distribution/
#https://www.nxn.se/valent/2018/1/30/count-depth-variation-makes-poisson-scrna-seq-data-negative-binomial

## load example dataset
library(tidyverse)
# to read in H5AD file, use https://github.com/theislab/zellkonverter
library(zellkonverter)

# example data download from:https://figshare.com/projects/Zero_inflation_in_negative_control_data/61292
svensson_data=readH5AD('svensson_chromium_control.h5ad') 
raw_counts<- svensson_data@assays@data$X
dim(raw_counts)
#[1] 24116  4000

# two datasets, each with 2000 cells
colnames(raw_counts) %>% stringr::str_extract("[0-9]+_") %>% table()
#20311_ 20312_ 
#2000   2000 

# only use the 2rd dataset '20312_'
raw_counts2<- raw_counts[, grepl(pattern = "20312_", x = colnames(raw_counts))]
dim(raw_counts2)
#[1] 24116  2000
gg<-Matrix::rowSums(raw_counts2)>0 #exclude genes that are zero across all cells
Y<-raw_counts2[gg,]
dim(Y)
#[1] 21411  2000

## look at mean ~ var relationship for each gene
# https://github.com/const-ae/sparseMatrixStats
library(sparseMatrixStats)
gene_means<- sparseMatrixStats::rowMeans2(raw_counts2)
gene_vars<- sparseMatrixStats::rowVars(raw_counts2)
length(gene_means);length(gene_vars) #24116 genes

df<- bind_cols(gene_means = gene_means, gene_vars = gene_vars)
# the quadratic polynomial mean-variance relation
df %>% ggplot(aes(x = log10(gene_means), y = log10(gene_vars))) +
  geom_point() +
  theme_classic(base_size = 14) +
  ggtitle("svensson et al 2")

# as for NB distribution (https://en.wikipedia.org/wiki/Negative_binomial_distribution)
# one way of parameterization 
# Mean: mu
# Var = mu + phi * mu^2
# plug in empirical mean and var to estimate phi
model<- lm(gene_vars ~  1* gene_means + I(gene_means^2) + 0, data =df )
x=summary(model)
x$coefficients
# 0.3725

#the same value as calculated in the paper: (Table 1)
#Valentine Svensson: Droplet scRNA-seq is not zero-inflated 
#plot fitter curve with obs data
predicted_df<- data.frame(mean = df$gene_means, var_predict = 
                            df$gene_means + 0.3725 * (df$gene_means)^2 )

df %>%  ggplot(aes(x = log10(gene_means), y = log10(gene_vars))) +
  geom_point() +
  geom_line(color = "red", data = predicted_df, aes(x = log10(gene_means), y =log10(var_predict))) + 
  theme_classic(base_size = 14) +
  ggtitle("svensson et al")


#Is single cell RNAseq data 0 inflated?
phi <- 1/0.3725
zeros_nb<- (phi/(gene_means + phi))^phi
zeros_observed<- apply(raw_counts2, 1, function(x) mean(x ==0))

data.frame(zeros_nb = zeros_nb, zeros_observed = zeros_observed, 
           gene_means = gene_means) %>%
  ggplot(aes(x =log10(gene_means), y = zeros_observed)) +
  geom_point() +
  geom_line(aes(x = log10(gene_means), y = zeros_nb), color = "red") +
  theme_classic(base_size = 14) +
  ggtitle("Svensson et al 2")

