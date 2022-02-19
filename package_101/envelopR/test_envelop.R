
library(S4Vectors)
library(Seurat)
library(ggplot2);library(gridExtra);
options(stringsAsFactors = F)


file="../single.cell_datasets/fly.brain.atlas/wholebrain_filtered_valid.rds";

dat=readRDS(file);
#dat=NormalizeData(dat);
#df.umi=dat@assays$RNA@counts #12094 features across 56192 samples within 1 assay 
df.expr=dat@assays$RNA@data #logNormal

unique(dat@meta.data$sex)
table(dat@meta.data$annotation)

## separate female and male, select cell clusters that contain>=500cells in both sexes
min.cell=200;max.cell=Inf;
i=apply(table(dat$sex,dat$annotation),2,function(i)sum(i>=min.cell & i<max.cell))
#pick.cell.clusters=names(which(i==2)) #79
pick.cell.clusters=names(which(i>=1)) #79
pick.cell.clusters #min200: sn,40. sc,38. min100: sn,54; sc,60.

dat=subset(dat,annotation %in% pick.cell.clusters)
df.expr=dat@assays$RNA@counts;
dim(df.expr); #sn:12602 50588

i.cluster='Ensheathing_glia'
#i.cluster='Astrocyte-like'
#i.cluster='G-KC'
#i.cluster='Pm1/Pm2' #no clear trend
#i.cluster='Mip' #no clear trend
#i.cluster='A/B-KC'
#i.cluster=c('Pm1/Pm2','Ensheathing_glia');
#i.cluster=pick.cell.clusters;

# select 1000 var genes
one.dat=subset(dat,annotation==i.cluster)
one.dat <- NormalizeData(one.dat, normalization.method = "LogNormalize", scale.factor = 10000)
one.dat <- FindVariableFeatures(one.dat, selection.method = "vst", nfeatures = 2000)
var.genes=one.dat@assays$RNA@var.features


mat.meta=dat@meta.data[dat$annotation %in% i.cluster,]
dim(mat.meta) #1508 cells
table(mat.meta$annotation,mat.meta$Age)

mat=df.expr[,dat$annotation %in% i.cluster]
#gene.filter=Matrix::rowSums(mat>0) >= max(15,ncol(mat)*0.005)
#mat=mat[gene.filter,]
mat=mat[var.genes,]
dim(mat); # 2000 1508, ngene>ncell, ncell>ngene


#Normalize and log-transform the data
expr.norm <- t(t(as.matrix(mat))/colSums(mat))*10000
expr.norm.log <- log(expr.norm + 1)

expr.scale <- t(scale(t(expr.norm.log))) #scale for each gene: center and scale
dim(expr.scale)  #gene by cell, every gene has 1508 'time series' data points
h=expr.scale
m=nrow(h) #ngene
t=ncol(h) #ncell
dim(h) #2000 gene by 1508 cell 

#############################################
## envelop
library(tidyverse)
library(ggExtra) #
library(mvtnorm)
library(Amelia) #
library(ggridges)
library(envelopeR)
library(modelr)
library(Matrix)
library(kableExtra)

subject_info=mat.meta
targeted=t(expr.scale)
dim(subject_info) #1508obs x 34 attributes
dim(targeted) #1508 obs x 2000 gene

# set response matrix
Y  <- targeted %>% as.matrix
dim(Y)

# set covariate matrix
names(subject_info)
table(subject_info$Age)
X <- subject_info %>% mutate(Sex = sex) %>%
  model_matrix(~ Age + Sex) %>%
  as.matrix
dim(X) #1508 x 3 covars
head(X)


Xfit  <- X[, c("Age", "Sexmale"), drop=FALSE]
Yfit  <- Y[, ]
dim(Xfit) #1508 x 2 #age, sex
dim(Yfit) #1508 2000


## Get the rank
s <- getRank(Yfit)
s #53
q <- ncol(X)
q #3 covars
X[1:3,] #intercept, age, sex


#############################
#### Envelope Fit
############################
res=readRDS('Ensheathing_glia_1508cell_2000gene_fitResult.rds');
if(F){
start.time=Sys.time()
res <- fit_envelope(Y=Yfit, X=Xfit, #X=scale(Xfit[indices, ]), 
                    #D=D, s=s,prior_counts=prior_counts,
                    s=s,
                    distn="covreg",
                    maxIters=1000, 
                    #U1=0*diag(s), U0=0*diag(r),L=0
                    Vinit="OLS", )
end.time=Sys.time()
cat(end.time-start.time,'seconds\n')
#"2022-02-13 21:42:10 PST"
#"2022-02-14 04:09:07 PST"
saveRDS(res,'Ensheathing_glia_1508cell_2000gene_fitResult.rds')
}

YVfit  <- Yfit %*% res$V
dim(YVfit) #1508 X 53 sub-space fitted features
dim(Yfit) #1508 X 2000 original features

covreg_fit=readRDS('Ensheathing_glia_1508cell_2000gene_covreg_fit.rds');
if(F){
#https://rdrr.io/cran/covreg/man/covreg.mcmc.html
covreg_fit  <- covreg::covreg.mcmc(YVfit ~ Xfit - 1,
                                   YVfit ~ Xfit,
                                   #R=5, 
                                   niter=10000,
                                   nthin=10)
saveRDS(covreg_fit,'Ensheathing_glia_1508cell_2000gene_covreg_fit.rds')
}
lapply(covreg_fit,dim)

## mean fit
mean_coefs  <- covreg_fit$B1.psamp
dim(mean_coefs) #53 subspace x 2 x 1000

## Covariance Fit
##https://rdrr.io/cran/covreg/src/R/cov.psamp.R
cov_psamp  <- covreg::cov.psamp(covreg_fit)
dim(cov_psamp) #16   53   53 1000

nrow(covreg_fit$matrix.cov) #1508 obs
nrow(unique(covreg_fit$matrix.cov)) #16 unique obs 

mean_coefs_age  <- res$V %*% mean_coefs[, 1, ]
mean_coefs_sex  <- res$V %*% mean_coefs[, 2, ]
dim(mean_coefs_age) #2000gene x 1000, 2000, original #feature
dim(mean_coefs_sex) #2000gene x 1000

head(colnames(Y)) #original gene feature names
rownames(mean_coefs_age)  <- rownames(mean_coefs_sex)  <-  colnames(Y)


## plot
dim(Xfit) #1508 x 2
X_unique <- unique(Xfit[, c("Age", "Sexmale")])
dim(Xfit) #1508 obs x 2 covars
dim(X_unique) #16  2

index_map  <- match(apply(Xfit, 1, function(x) 
  paste(x, collapse="_")),apply(X_unique, 1, function(x) paste(x, collapse="_")))
length(index_map) #1508

nms <- paste(X_unique[, 1], X_unique[, 2], sep="_")
nms <- X_unique[, 1] #age
length(nms)  #63
head(nms)

## aging_to_plot
which(X_unique[,2]==1) #male
which(X_unique[,2]==0) #female

## plot posterior distributions for min, max and quartiles of X
#dim(X_unique[which(X_unique[,2]==1),]) #8 x 2
#ix <- sort(X_unique[which(X_unique[,2]==1),1],index.return=TRUE)$ix
#ix <- sort(X_unique[,1],index.return=TRUE)$ix
#tmp=floor(c(1, 25, 50, 75, 100)*length(ix)/100)
#if(tmp[1]==0){tmp[1]=1}
#to_plot <- ix[tmp]
to_plot <- which(X_unique[,2]==1)
names(to_plot)=X_unique[to_plot, 1] 
to_plot

## values correspond to the quantiles of x
#names(to_plot)  <- c(0, 0.25, 0.5, 0.75, 1)
#to_plot  <- c(13, 21, 63, 28, 1, 52) ## MALES
#to_plot  <- c(13, 21, 62, 28, 1, 52) ## MALES
X_unique[to_plot, ] # males

#save(Yfit, Xfit, covreg_fit, cov_psamp, Vfit, res, file = paste0("targeted_covreg-", today(), ".Rdata"))

dim(res$V) #2000 metabolites x 53 subspace
dim(cov_psamp) #16   53   53 1000

cols <- rev(colorspace::sequential_hcl(length(to_plot)+1, "Viridis"))
#cols <-RColorBrewer::brewer.pal(length(to_plot), 'Dark2')
s #53

source('ming_posterior_plot.R') #for age colors
for(i in seq(1, (s-1), 2)) {
  
  combo  <- create_plots(res$V, cov_psamp,
                         n1=to_plot[1], n2=to_plot[length(to_plot)],
                         to_plot = to_plot, col_values=cols, 
                         obs_names=nms[to_plot],  #cell at selected ages
                         view=c(i, i+1), nlabeled=20,
                         labels=colnames(Y),  #gene names
                         main="Ensheathing_glia",legend.title="Age")
  combo
  
  ggsave(sprintf("aging_Plots/aging_plot-%i%i.pdf", i, i+1), combo, width=14)
}



## save regression result
tmp=map_dfr(list(Age=mean_coefs_age, Sex=mean_coefs_sex), 
        function(mat) 
          apply(mat, 1,
                function(x) {
                  frac_neg  <- mean(x < 0)
                  pval <- 2*min(frac_neg, 1-frac_neg)
                  tstat  <- mean(x) / sd(x)
                  c("P-value"=pval, "T-statistic"=tstat)
                }) %>% t %>% as_tibble(rownames="Metabolite"), .id="Type")
regression_stats <- tmp %>%
  group_by(Type) %>% 
  arrange(`P-value`, desc(abs(`T-statistic`))) %>%
  mutate(`Q-value` = `P-value`*n()/row_number()) %>%
  ungroup() 
 

regression_stats %>% filter(`Q-value` < 0.05, Type == "Sex") %>%
  kable(format="latex") %>%
  cat(., file = "Ensheathing_glia_sex.tex")

regression_stats %>% filter(`Q-value` < 0.05, Type == "Age") %>%
  kable(format="latex") %>%
  cat(., file = "Ensheathing_glia_age.tex")




