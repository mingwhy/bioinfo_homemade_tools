
library(S4Vectors)
library(Seurat)
library(ggplot2);library(gridExtra);
options(stringsAsFactors = F)

## read in data
file="~/Documents/Data_fly_FCA/fly.brain.atlas/wholebrain_filtered_valid.rds";

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

#i.cluster='Ensheathing_glia'
#i.cluster='Astrocyte-like'
#i.cluster='G-KC'
#i.cluster='Pm1/Pm2' #no clear trend
#i.cluster='Mip' #no clear trend
i.cluster='A/B-KC'
#i.cluster=c('Pm1/Pm2','Ensheathing_glia');
#i.cluster=pick.cell.clusters;

# select 2000 var genes
one.dat=subset(dat,annotation==i.cluster & dat$Age %in% c(1,9,30))
one.dat
one.dat <- NormalizeData(one.dat, normalization.method = "LogNormalize", scale.factor = 10000)
one.dat <- FindVariableFeatures(one.dat, selection.method = "vst", nfeatures = 2000)
#one.dat <- FindVariableFeatures(one.dat, selection.method = "vst", nfeatures = 3000)
var.genes=one.dat@assays$RNA@var.features

mat.meta=dat@meta.data[dat$annotation %in% i.cluster & dat$Age %in% c(1,9,30),]
dim(mat.meta) #688 cells
table(mat.meta$annotation,mat.meta$Age)

mat=df.expr[,dat$annotation %in% i.cluster & dat$Age %in% c(1,9,30)]
#gene.filter=Matrix::rowSums(mat>0) >= max(15,ncol(mat)*0.005)
#mat=mat[gene.filter,]
mat=mat[var.genes,]
dim(mat); # 2000  256, ngene>ncell, ncell>ngene
sum(colnames(mat)==rownames(mat.meta)) #256 cells

#Normalize and log-transform the data
expr.norm <- t(t(as.matrix(mat))/colSums(mat))*10000
#expr.norm.log <- log(expr.norm + 1)
t.expr.norm <- t(expr.norm)

#scale for each gene within each age group: center
mat.meta$Age[!duplicated(mat.meta$Age)]
df<-Reduce(`rbind`,lapply(mat.meta$Age[!duplicated(mat.meta$Age)],function(i){
  tmp=t.expr.norm[mat.meta$Age==i,];
  tmp=scale(tmp,center = T,scale=F)
  #summary(apply(tmp,1,mean))
  tmp
}) )
dim(df)
df<-df[rownames(mat.meta),]
sum(rownames(df)==rownames(mat.meta)) #256

expr.scale=df;
summary(apply(expr.scale,1,mean)) #cell
summary(apply(expr.scale,2,mean)) #gene

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
targeted=expr.scale
dim(subject_info) #256 cells x 34 attributes
dim(targeted) #256 obs x 2000 gene

# set response matrix
Y  <- targeted %>% as.matrix
dim(Y)

# set covariate matrix
names(subject_info)
table(subject_info$Age)
X <- subject_info %>% mutate(Sex = sex) %>%
  model_matrix(~ Age + Sex) %>%
  as.matrix
dim(X) #256 x 3 covars
head(X)


#Xfit  <- X[, c("Age", "Sexmale"), drop=FALSE]
Xfit  <- X[, c("Age"), drop=FALSE]
Yfit  <- Y[, ]
dim(Xfit) #256 x 2 #age, sex
dim(Yfit) #256 2000


## Get the rank
s <- getRank(Yfit)
s #48
q <- ncol(X)
q #3 covars
X[1:3,] #intercept, age, sex


#############################
#### Envelope Fit
############################
res=readRDS('A-B-KC_256cell_2000gene_fitResult.rds');
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
  #"2022-05-12 21:08:24 PDT"
  #"2022-05-12 21:11:38 PDT"
  saveRDS(res,'A-B-KC_256cell_2000gene_fitResult.rds')
}

names(res) #V, intercept, ...

YVfit  <- Yfit %*% res$V
dim(Yfit) #1508 cell x 200 gene, original response Y
dim(YVfit) #1508 cell x 5, project Y into reduced sub-space

## cal mahalanobis distance for projected data
dim(subject_info)
dim(YVfit)
sum(rownames(YVfit)==rownames(subject_info))
mh.out<-lapply(unique(subject_info$Age),function(age){
  x=YVfit[subject_info$Age==age,]
  d.mh=mahalanobis(x,colMeans(x),cov(x)) #https://www.geeksforgeeks.org/how-to-calculate-mahalanobis-distance-in-r/
})
df=data.frame(age=rep(unique(subject_info$Age),sapply(mh.out,length)),dist=unlist(mh.out))
head(df)
df$age=factor(df$age)
ggplot(df,aes(x=dist,group=age,col=age))+geom_density()

ggplot(df,aes(x=age,y=dist,col=age))+geom_violin()+#geom_boxplot(outlier.shape = NA)+
  scale_y_log10()+geom_jitter(size=0.1)+theme_bw()

#######
covreg_fit=readRDS('A-B-KC_256cell_2000gene_covreg_fit.rds');
if(F){
  #https://rdrr.io/cran/covreg/man/covreg.mcmc.html
  covreg_fit  <- covreg::covreg.mcmc(YVfit ~ Xfit - 1,
                                     YVfit ~ Xfit,
                                     #R=5, 
                                     niter=10000,
                                     nthin=10) #10000/10, 1000 samples extracted
  saveRDS(covreg_fit,'A-B-KC_256cell_2000gene_covreg_fit.rds')
}
lapply(covreg_fit,dim)
names(covreg_fit) #"B1.psamp"    "B2.psamp"    "A.psamp"     "matrix.mean" "matrix.cov"    

## mean fit
mean_coefs  <- covreg_fit$B1.psamp
dim(mean_coefs) #5-dim subspace x 2 est.para(age,sex) x 1000 (posterior.samples)
mean_coefs[,,1]

## Covariance Fit
##https://rdrr.io/cran/covreg/src/R/cov.psamp.R
cov_psamp  <- covreg::cov.psamp(covreg_fit)
dim(cov_psamp) # 6 5 5 1000, 6 due to unique sex_age combination, 3 age x 2 sex=6 combos
cov_psamp[,,,1]
cov_psamp[1,,,1] #one sex_age combo (predictor) <=> one posterior.sample of cov.mat
sum(diag(cov_psamp[1,,,1])) #155.7345
det(cov_psamp[1,,,1]) #4981383
sum(eigen(cov_psamp[1,,,1])$values)  #155.7345 #https://stats.stackexchange.com/questions/225434/a-measure-of-variance-from-the-covariance-matrix
combo6.1000var<-lapply(1:6,function(j){
  tmp=sapply(1:1000,function(i){
    #tmp=cov_psamp[j,,,i] #one sex_age combo (predictor) <=> one posterior.sample of cov.mat
    #sum(diag(tmp))
    det(cov_psamp[j,,,i]) 
  })
  return(tmp)
})

X_unique <- unique(Xfit[, c("Age", "Sexmale")])
dim(Xfit) #256 obs x 2 covars
dim(X_unique) #6  2
label=paste0('Age',X_unique[,1],', Sexmale',X_unique[,2])

df.group.var=data.frame(group=rep(label,sapply(combo6.1000var,length)),var=unlist(combo6.1000var))
summary(df.group.var$var)
ggplot(df.group.var,aes(x=group,y=var))+geom_jitter(size=0.01)+coord_flip()+
  theme_bw()+geom_violin(fill=NA)

##
nrow(covreg_fit$matrix.cov) #256 obs
nrow(unique(covreg_fit$matrix.cov)) #6 unique obs (sex+age)

dim(res$V) #2000 gene x 5 factors
dim(mean_coefs[, 1, ]) #5 factors X 1000 post.samples
mean_coefs_age  <- res$V %*% mean_coefs[, 1, ]
mean_coefs_sex  <- res$V %*% mean_coefs[, 2, ]
dim(mean_coefs_age) #2000gene x 1000, 2000, original #feature
dim(mean_coefs_sex) #2000gene x 1000

head(colnames(Y)) #original gene feature names
rownames(mean_coefs_age)  <- rownames(mean_coefs_sex)  <-  colnames(Y)

mean_coefs_age[1:3,1:3]

## plot
dim(Xfit) #256 x 2
X_unique <- unique(Xfit[, c("Age", "Sexmale")])
dim(Xfit) #256 obs x 2 covars
dim(X_unique) #6  2
X_unique

index_map  <- match(apply(Xfit, 1, function(x) 
  paste(x, collapse="_")),apply(X_unique, 1, function(x) paste(x, collapse="_")))
length(index_map) #256
table(index_map) #number of cell per age-sex combination

nms <- paste(X_unique[, 1], X_unique[, 2], sep="_")
nms <- X_unique[, 1] #age
length(nms)  #6
head(nms)

## aging_to_plot
which(X_unique[,2]==1) #male 2 3 4
which(X_unique[,2]==0) #female 1 5 6

## plot posterior distributions for min, max and quartiles of X
#dim(X_unique[which(X_unique[,2]==1),]) #8 x 2
#ix <- sort(X_unique[which(X_unique[,2]==1),1],index.return=TRUE)$ix
#ix <- sort(X_unique[,1],index.return=TRUE)$ix
#tmp=floor(c(1, 25, 50, 75, 100)*length(ix)/100)
#if(tmp[1]==0){tmp[1]=1}
#to_plot <- ix[tmp]
to_plot <- which(X_unique[,2]==1) #male
X_unique[to_plot, 1] #3 ages
names(to_plot)=X_unique[to_plot, 1] 
to_plot #plot 3 age groups

## values correspond to the quantiles of x
#names(to_plot)  <- c(0, 0.25, 0.5, 0.75, 1)
#to_plot  <- c(13, 21, 63, 28, 1, 52) ## MALES
#to_plot  <- c(13, 21, 62, 28, 1, 52) ## MALES
X_unique[to_plot, ] # males

#save(Yfit, Xfit, covreg_fit, cov_psamp, Vfit, res, file = paste0("targeted_covreg-", today(), ".Rdata"))

dim(res$V) #2000 gene x 5 (factors) subspace
dim(cov_psamp) # 6    5    5 1000

cols <- rev(colorspace::sequential_hcl(length(to_plot)+1, "Viridis"))
#cols <-RColorBrewer::brewer.pal(length(to_plot), 'Dark2')
s #5

#source('ming_posterior_plot.R') #for age colors
dir.create('aging_Plots')
for(i in seq(1, (s-1), 2)) {
  
  combo  <- create_plots(res$V, cov_psamp,
                         n1=to_plot[1], n2=to_plot[length(to_plot)],
                         to_plot = to_plot, col_values=cols, 
                         obs_names=nms[to_plot],  #cell at selected ages
                         view=c(i, i+1), nlabeled=20,
                         labels=colnames(Y),  #gene names
                         main="Kenyon cells",legend.title="Age")
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
dim(regression_stats) #4000 x 5
head(regression_stats)

regression_stats %>% filter(`Q-value` < 0.05, Type == "Sex") %>%
  kable(format="latex") %>%
  cat(., file = "KC_sex.tex")
head(regression_stats)

regression_stats %>% filter(`Q-value` < 0.05, Type == "Age") %>%
  kable(format="latex") %>%
  cat(., file = "KC_age.tex")
dim(regression_stats) #4000 x 5


