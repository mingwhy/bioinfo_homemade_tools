library(tidyverse)
library(ggExtra)
library(mvtnorm)
library(Amelia)
library(ggridges)
library(envelopeR)
library(modelr)
library(Matrix)
library(MetaboAnalystR)
library(kableExtra)

targeted  <- read_csv("metabolomics_data/targeted.csv")
subject_info  <- read_csv("metabolomics_data/subject_info.csv")

dim(targeted) #198 x 109 metabolites
dim(subject_info) #198 x 4

Y  <- targeted %>% as.matrix

## Remove extreme outliers (more than 4 median asbolute deviations from 0)
Y  <- apply(Y, 2, function(x) {
  mad_val <- mad(x, na.rm=TRUE)
  outliers  <- which(abs(x) > 4*mad_val)
  x[outliers]  <- NA
  x
})

## Impute missing values using
Y <- amelia(Y, m = 1, empri = 100)$imputations$imp1
dim(Y) #198 x 109

## Focus only on aging among controls
table(subject_info$Type)
control_indices  <- which(subject_info$Type %in% c("CY", "CM", "CO"))

X <- subject_info %>% 
  model_matrix(~ Age + Sex) %>%
  as.matrix
head(X)

Xfit  <- X[control_indices, c("Age", "SexM"), drop=FALSE]
Yfit  <- Y[control_indices, ]
dim(Xfit) #85 x 2
dim(Yfit) #85 x 109
head(Xfit) #only Age and SexM in the columns

Xfit[order(Xfit[, 1]), ]
xord  <- order(Xfit[, 1]) #order by age

s <- getRank(Yfit)
s #18
q <- ncol(X)
q #3

fit_obj  <- lm(Yfit ~ Xfit)
getRank(fit_obj$residuals)
# 18

indices  <- 1:nrow(Xfit)
prior_counts  <- 1

#############################
#### Envelope Fit
############################
dim(scale(Xfit[indices, ])) #85 x 2

#https://github.com/afranks86/envelopeR/blob/master/metabolomics/untargeted.R
if(!file.exists('targeted_fit.rds')){
  res <- fit_envelope( Y=Yfit, X=scale(Xfit[indices, ]), 
                              s=s, distn="covreg",
                              #D=D, prior_counts=prior_counts,
                               maxIters=1000, #U1=0*diag(s), U0=0*diag(r),L=0
                               Vinit="OLS")
  saveRDS(res,'targeted_fit.rds')
}

res=readRDS('targeted_fit.rds')

YVfit  <- Yfit %*% res$V
dim(YVfit) #85 x 18
s #18

# generater posterior.samples for covariance matrix
if(!file.exists('targeted_covreg_fit.rds')){
  covreg_fit  <- covreg::covreg.mcmc(YVfit ~ scale(Xfit[indices, ]) - 1,
                                   YVfit ~ scale(Xfit[indices, ]),
                                   R=10, niter=10000, nthin=10)
  saveRDS(covreg_fit,'targeted_covreg_fit.rds')
}

covreg_fit=readRDS('targeted_covreg_fit.rds')
names(covreg_fit) # "B1.psamp"    "B2.psamp"    "A.psamp"     "matrix.mean" "matrix.cov" 

## mean fit
mean_coefs  <- covreg_fit$B1.psamp
dim(mean_coefs) #18 2 1000, 18 factors, 2 X.covar, 1000 samples

## Covariance Fit
cov_psamp  <- covreg::cov.psamp(covreg_fit)
dim(Yfit) #85obs x 109 features
dim(res$V) #109 x 18, feature to factor
dim(mean_coefs[, 1, ]) #18 1000
head(Xfit) #column1=Age, column2=Sex

mean_coefs_age  <- res$V %*% mean_coefs[, 1, ]
mean_coefs_sex  <- res$V %*% mean_coefs[, 2, ]
rownames(mean_coefs_age)  <- rownames(mean_coefs_sex)  <-  colnames(Y) #metabolite names
dim(mean_coefs_age) #109 feature x 1000 samples
dim(mean_coefs_sex)

map_dfr(list(Age=mean_coefs_age, Sex=mean_coefs_sex), function(mat) apply(mat, 1,
                     function(x) {
                       frac_neg  <- mean(x < 0)
                       pval <- 2*min(frac_neg, 1-frac_neg)
                       tstat  <- mean(x) / sd(x)
                       c("P-value"=pval, "T-statistic"=tstat)
                     }) %>% t %>% as_tibble(rownames="Metabolite"), .id="Type") %>%
  group_by(Type) %>% 
  arrange(`P-value`, desc(abs(`T-statistic`))) %>%
  mutate(`Q-value` = `P-value`*n()/row_number()) %>%
  ungroup() ->
  regression_stats
dim(regression_stats) #218
table(regression_stats$Type) #109 age, 109 sex
head(regression_stats)

#regression_stats %>% filter(`Q-value` < 0.05, Type == "Sex") %>%
#  kable(format="latex") %>% cat(., file = "targeted_sex.tex")

#regression_stats %>% filter(`Q-value` < 0.05, Type == "Age") %>%
#  kable(format="latex") %>% cat(., file = "targeted_age.tex")

X_unique <- unique(Xfit[indices, c("Age", "SexM")])
X_unique #63 unique combinations
dim(Xfit[indices,]) #85obs x 2 covar

index_map  <- match(apply(Xfit, 1, function(x) paste(x, collapse="_")), apply(X_unique, 1, function(x) paste(x, collapse="_")))
length(index_map) #85
max(index_map) #63, each of the 85 obs map to only 1 X.covar combination

#nms <- paste(X_unique[, 1], X_unique[, 2], sep="_")
nms <- X_unique[, 1]
nms #63
head(X_unique)

head(X_unique)
x=lapply(1:dim(cov_psamp)[4],function(j){
  tr(cov_psamp[2,,,j])
})
x=unlist(x)
summary(x)
boxplot(x)

## aging_to_plot
rownames(X_unique)=as.character(1:nrow(X_unique))
sum(X_unique[,2]==1) #30 male samples
X_unique[which(X_unique[,2]==1),1]
order(X_unique[which(X_unique[,2]==1),1])
quantile(X_unique[which(X_unique[,2]==1),1])
#min 21, 25% 33, median 57, 75%72, max 86
names(which(X_unique[X_unique[,2]==1,1]==86))
to_plot  <- c(28,33, 2, 35,5) ## MALES
X_unique[to_plot, ]
#to_plot  <- c(44, 26, 43, 2, 27) ## FEMALES
#X_unique[to_plot, ]

#save(Yfit, Xfit, covreg_fit, cov_psamp, YVfit, res, file = paste0("targeted_covreg.Rdata"))

cols <- rev(colorspace::sequential_hcl(length(to_plot)+1, "Viridis"))
cols
dim(res$V)
for(i in seq(1, ncol(res$V), 2)) {
  if((i+1)>ncol(res$V)){next}
  combo  <- create_plots(res$V, cov_psamp,
                         n1=to_plot[1], n2=to_plot[length(to_plot)],
                         to_plot = to_plot, col_values=cols, nlabeled=20,
                         obs_names=nms[to_plot], view=c(i, i+1),  #view, which PC to use for rotating V matrix
                         labels=colnames(Y), main="Targeted", legend.title="Age")
  ggsave(sprintf("output_plots/aging_plot-%i%i.pdf", i, i+1), combo, width=14)
}

## gender_to_plot
to_plot  <- c(43, 32)
X_unique[to_plot, ]
nms <- X_unique[, 2]

cbind(X_unique, 1:nrow(X_unique))[order(X_unique[, 1]), ]
to_plot  <- c(49, 14, 15, 58)
to_plot  <- c(14, 49, 58, 15)
X_unique[to_plot, ]

nms <- paste(X_unique[, 1], X_unique[, 2], sep="_")

cols <- rev(colorspace::sequential_hcl(length(to_plot)+1, "Viridis"))
for(i in seq(1, s, 2)) {

  combo  <- create_plots(res$V, cov_psamp,
                         n1=to_plot[1], n2=to_plot[length(to_plot)],
                         to_plot = to_plot, col_values=cols,
                         obs_names=nms[to_plot], view=c(i, i+1), labels=colnames(Y))
  ggsave(sprintf("output_plots/sex_plot-%i%i.pdf", i, i+1), combo, width=14)
}


## Save all metabolites
#all_mets <- subject_data$Metabolite %>% unique %>% as.data.frame
#write_csv(all_mets, "all_mets.csv")

i  <- 7
rotation_result <- rotate_basis(res$V, cov_psamp, n1=to_plot[1], n2=to_plot[length(to_plot)])
V_rotated  <- rotation_result$rotV[, i:(i+1)]
rownames(V_rotated) <- colnames(Y)
sig_mets1 <- names(sort(abs(V_rotated[, 1]), decreasing=TRUE))[1:20]
sig_mets2 <- names(sort(abs(V_rotated[, 2]), decreasing=TRUE))[1:20]
sig_mets_both <- names(sort(rowSums(V_rotated[, 1:2]^2), decreasing=TRUE))[1:20]

Vall <- rowSums(rotation_result$rotV^2)
names(Vall) <- colnames(Y)
sig_mets_both <- names(sort(Vall, decreasing=TRUE)[1:10])



