# install envelopeR: devtools::install_github("afranks86/envelopeR")
# https://github.com/afranks86/envelopeR
library(tidyverse)
library(ggExtra) #
library(mvtnorm)
library(Amelia) #
library(ggridges)
library(envelopeR)
library(modelr)
library(Matrix)
#library(MetaboAnalystR) #https://github.com/xia-lab/MetaboAnalystR
library(kableExtra)

targeted  <- read_csv("data/targeted.csv")
subject_info  <- read_csv("data/subject_info.csv")
dim(targeted) #198obs x 109 metabolies
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
dim(Y) #198 obs x 109 metabolites

# set covariate matrix
names(subject_info)
table(subject_info$Type)
subject_info  <- subject_info %>% mutate(Type2 = ifelse(Type != "AD", "C", "AD"))
X <- subject_info %>% mutate(Type2 = ifelse(Type != "AD", "C", "AD")) %>%
  model_matrix(~ Age + Type2 + Sex) %>%
  as.matrix
dim(X) #198 x 4 covars

## Focus only on aging among controls
control_indices  <- which(subject_info$Type %in% c("CY", "CM", "CO"))

Xfit  <- X[control_indices, c("Age", "SexM"), drop=FALSE]
Yfit  <- Y[control_indices, ]
dim(Xfit) #85 x2 #age, sex
dim(Yfit) #85 x 109

Xfit[order(Xfit[, 1]), ]
xord  <- order(Xfit[, 1])

## Get the rank
s <- getRank(Yfit)
s #18
q <- ncol(X)
q #4 covars
X[1:3,] #intercept, age, type2c, sexM

fit_obj  <- lm(Yfit ~ Xfit)
getRank(fit_obj$residuals)
#17

indices  <- 1:nrow(Xfit)
prior_counts  <- 1

#############################
#### Envelope Fit
############################

res <- fit_envelope(Y=Yfit, X=Xfit, #X=scale(Xfit[indices, ]), 
                    #D=D, s=s,prior_counts=prior_counts,
                    s=s,
                    distn="covreg",
                    maxIters=1000, 
                    #U1=0*diag(s), U0=0*diag(r),L=0
                    Vinit="OLS", )

YVfit  <- Yfit %*% res$V
dim(YVfit) #85 X 2 sub-space fitted features
dim(Yfit) #85 X 109 original features

#https://rdrr.io/cran/covreg/man/covreg.mcmc.html
#covreg_fit  <- covreg::covreg.mcmc(YVfit ~ scale(Xfit[indices, ]) - 1,
#                                   YVfit ~ scale(Xfit[indices, ]),
#                                   #R=10, 
#                                   niter=10000, nthin=10, verb=FALSE)
covreg_fit  <- covreg::covreg.mcmc(YVfit ~ Xfit - 1,
                                   YVfit ~ Xfit,
                                   #R=5, 
                                   niter=10000,
                                   nthin=10)

lapply(covreg_fit,dim)

## mean fit
mean_coefs  <- covreg_fit$B1.psamp
dim(mean_coefs) #18 x 2 x 1000

## Covariance Fit
##https://rdrr.io/cran/covreg/src/R/cov.psamp.R
cov_psamp  <- covreg::cov.psamp(covreg_fit)
dim(cov_psamp) #63   18   18 1000

nrow(covreg_fit$matrix.cov) #85 obs
nrow(unique(covreg_fit$matrix.cov)) #63, only 63 unique obs 

mean_coefs_age  <- res$V %*% mean_coefs[, 1, ]
mean_coefs_sex  <- res$V %*% mean_coefs[, 2, ]
dim(mean_coefs_age) #109 x 1000, 109, original #feature
dim(mean_coefs_sex) #109 x 1000

head(colnames(Y)) #original metabolites feature names
rownames(mean_coefs_age)  <- rownames(mean_coefs_sex)  <-  colnames(Y)


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
  cat(., file = "targeted_sex.tex")

regression_stats %>% filter(`Q-value` < 0.05, Type == "Age") %>%
  kable(format="latex") %>%
  cat(., file = "targeted_age.tex")


## plot
#X_unique <- unique(Xfit[indices, c("Age", "SexM")])
X_unique <- unique(Xfit[, c("Age", "SexM")])
dim(Xfit) #85 obs x 2 covars
dim(X_unique) #63  2

index_map  <- match(apply(Xfit,1,function(x) 
  paste(x, collapse="_")),apply(X_unique,1,function(x) paste(x, collapse="_")))
length(index_map) #85

nms <- paste(X_unique[, 1], X_unique[, 2], sep="_")
nms <- X_unique[, 1]
length(nms)  #63

## aging_to_plot
which(X_unique[,2]==1) #male
which(X_unique[,2]==0) #female

## plot posterior distributions for min, max and quartiles of X
dim(X_unique[which(X_unique[,2]==1),]) #30 x 2
ix <- sort(X_unique[which(X_unique[,2]==1),1], index.return=TRUE)$ix
tmp=floor(c(1, 25, 50, 75, 100)*length(ix)/100)
if(tmp[1]==0){tmp[1]=1}
to_plot <- ix[tmp]

## values correspond to the quantiles of x
names(to_plot)  <- c(0, 0.25, 0.5, 0.75, 1)
#to_plot  <- c(13, 21, 63, 28, 1, 52) ## MALES
#to_plot  <- c(13, 21, 62, 28, 1, 52) ## MALES
X_unique[to_plot, ]

#save(Yfit, Xfit, covreg_fit, cov_psamp, Vfit, res, file = paste0("targeted_covreg-", today(), ".Rdata"))

dim(res$V) #109 metabolites x 18 subspace
dim(cov_psamp) #63 unique obs x 18 subspace

cols <- rev(colorspace::sequential_hcl(length(to_plot)+1, "Viridis"))

for(i in seq(1, s, 2)) {
  
  combo  <- create_plots(res$V, cov_psamp,
                         n1=to_plot[1], n2=to_plot[length(to_plot)],
                         to_plot = to_plot, col_values=cols, nlabeled=20,
                         obs_names=nms[to_plot], view=c(i, i+1), labels=colnames(Y), 
                         main="Targeted", legend.title="Age")
  ggsave(sprintf("aging_plot-%i%i.pdf", i, i+1), combo, width=14)
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
  ggsave(sprintf("sex_plot-%i%i.pdf", i, i+1), combo, width=14)
}

######################################
## Save all metabolites
all_mets <- subject_data$Metabolite %>% unique %>% as.data.frame
write_csv(all_mets, "all_mets.csv")

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

mSet <- InitDataObjects("conc", "msetora", FALSE)
mSet <- Setup.HMDBReferenceMetabolome(mSet, "all_mets.csv");
cmpd.vec <- sig_mets_both
mSet <- Setup.MapData(mSet, cmpd.vec);
mSet <- CrossReferencing(mSet, "name");
mSet <- CreateMappingResultTable(mSet)
mSet <- SetMetabolomeFilter(mSet, T);
mSet <- SetCurrentMsetLib(mSet, "smpdb_pathway", 2);
mSet <- CalculateHyperScore(mSet)
mSet <- CalculateHyperScore(mSet)

## Scores
mSet$analSet$ora.mat
mSet$name.map$map.table
mSet$analSet$ora.hits

as_tibble(mSet$dataSet$map.table)$KEGG


