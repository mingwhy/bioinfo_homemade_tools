
#https://github.com/afranks86/envelopeR
library(envelopeR)

test_data <- generate_test_data(q=2, cov_rank=2, intercept=TRUE, 
                                gamma_sd=1, error_sd=1, 
                                seed=4445)
Y <- test_data$Y
X  <- test_data$X
s  <- test_data$s
p  <- test_data$p
V  <- test_data$V
dim(Y) #100 x 10
dim(X) #100 x 2
s #2 subspace dimension
p #10 feature
dim(V) #10 x 2 (10 feature ~ 2 factors in the shared subspace)

# fit model
## `get_rank` can be used to infer the appropriate dimension of hte subspace of material variation
s_hat <- getRank(Y)
s_hat #2
envfit <- fit_envelope(Y, X, distn="covreg", s=s_hat,
                       Vinit="OLS",
                       verbose_covreg=FALSE)
print(sprintf("Subspace similarity: %f", 
              tr(envfit$V  %*% t(envfit$V)  %*% V  %*% t(V))/s))
## Sphericity test: if small pvalue, reject isotropy of immaterial subspace
Vperp <- NullC(envfit$V) 
print(sprintf("P-value of sphericity test is: %f", 
              test_sphericity(Y %*% Vperp)$p.value))

##
YV  <- Y %*% envfit$V
dim(Y) #100 x 10 feature;
dim(YV) #100 x 2 factor; 
dim(X) #100 x 2;

res  <-  covreg::covreg.mcmc(YV ~ X - 1, YV ~ X, niter=10000, nthin=10, verb=FALSE)
names(res) #"B1.psamp"    "B2.psamp"    "A.psamp"     "matrix.mean" "matrix.cov" 
cov_psamp  <- covreg::cov.psamp(res)
cov_psamp  <- cov_psamp[, , , 1:1000]

dim(YV) #100 x 2
dim(cov_psamp) #100 x 2 x 2 x 1000, 100 unique X row.combinations

## plot
type  <- "mag"

## plot posterior distributions for min, max and quartiles of X
ix <- sort(X[, 1], index.return=TRUE)$ix
obs_to_plot <- ix[c(1, 25, 50, 75, 100)]
obs_to_plot
## values correspond to the quantiles of x
names(obs_to_plot)  <- c(0, 0.25, 0.5, 0.75, 1)

# maximaize the differnece between below two to set up shared_subspace
obs_to_plot[1] #min quantile of X
obs_to_plot[length(obs_to_plot)] #max quantile of X

cols <- colorspace::sequential_hcl(5, "viridis")
post  <- create_plots(envfit$V, cov_psamp,
                      n1=obs_to_plot[1], n2=obs_to_plot[length(obs_to_plot)],
                      to_plot = obs_to_plot, col_values=cols,
                      labels=colnames(Y), plot_type="posterior", alpha=0.5)

post

# Biplot
#Contours show posterior mean covariances matrices for (\Psi_x). 
#Points represent the largest magnitue loadings on this two-dimensional subspace. We choose the subspace to maximize the difference in covariance matrices between observation (n_1) and (n_2)
rownames(envfit$V)  <- 1:p
biplot  <- create_plots(envfit$V, cov_psamp,
                        n1=obs_to_plot[1], n2=obs_to_plot[length(obs_to_plot)],
                        to_plot = obs_to_plot, col_values=cols,
                        labels=colnames(Y), plot_type="biplot")
biplot

###############################################
## create_plots subfunctions, rotate_basis
post  <- create_plots(envfit$V, cov_psamp,
                      n1=obs_to_plot[1], n2=obs_to_plot[length(obs_to_plot)],
                      to_plot = obs_to_plot, col_values=cols,
                      labels=colnames(Y), plot_type="posterior", alpha=0.5)
#rotate_basis <- function(V, samples, n1=1, n2=NULL) {
dim(cov_psamp)   #100,2,2,1000
V=envfit$V;samples=cov_psamp;
n1=obs_to_plot[1];n2=NULL;
  
  S <- ncol(V)
  
  pmPsi1 <- apply(samples[n1, , , ], 1:2, mean)
  dim(pmPsi1) #2 2
  if(is.null(n2)) {
    ## Compare subspace that explains largest source of variability
    O <- svd(pmPsi1[1:S, 1:S])$u
    dim(O)
  } else {
    ## Compare largest difference between 2 obs
    pmPsi2 <- apply(samples[n2, , , ], 1:2, mean)
    O <- svd(pmPsi1[1:S, 1:S] - pmPsi2[1:S, 1:S])$u
  }
  dim(V) # feature -> factor, 10 x 2
  Vstar <- V[, 1:S] %*% O
  dim(Vstar) # feature -> factor, 10 x 2
  list(rotV=Vstar, rotMat=O)
}


