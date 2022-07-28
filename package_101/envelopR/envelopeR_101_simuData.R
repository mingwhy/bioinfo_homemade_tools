#https://github.com/afranks86/envelopeR
#devtools::install_github("afranks86/envelopeR")

library(envelopeR)

test_data <- generate_test_data(q=2, cov_rank=2, intercept=TRUE, gamma_sd=1, error_sd=1, seed=4445)
Y <- test_data$Y
X  <- test_data$X
s  <- test_data$s
p  <- test_data$p
V  <- test_data$V

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

#We can also get cov_psamp from the list returned by fit_envelope by running cov_psamp <- cov.psamp(envfit$covariance_list$covreg_res) (default 100 iterations).
YV  <- Y %*% envfit$V
dim(YV) #100 2

#https://www.rdocumentation.org/packages/covreg/versions/1.0/topics/covreg.mcmc
res  <-  covreg::covreg.mcmc(YV ~ X - 1, YV ~ X, niter=10000, nthin=10, verb=FALSE)
cov_psamp  <- covreg::cov.psamp(res)
cov_psamp  <- cov_psamp[, , , 1:1000]
dim(cov_psamp) #100    2    2 1000

#Posterior plot
#We plot posterior samples of (\Psi_x) for minimum, maximum and quartiles of the first covariate (X_1). Non-overlapping groups is suggestive of differences in the covariances a posteriori.
type  <- "mag"

## plot posterior distributions for min, max and quartiles of X
dim(X) #100 x 2
ix <- sort(X[, 1], index.return=TRUE)$ix
obs_to_plot <- ix[c(1, 25, 50, 75, 100)]

## values correspond to the quantiles of x
names(obs_to_plot)  <- c(0, 0.25, 0.5, 0.75, 1)

cols <- colorspace::sequential_hcl(5, "viridis")

dim(envfit$V); dim(cov_psamp)
post  <- create_plots(envfit$V, cov_psamp,
                      n1=obs_to_plot[1], n2=obs_to_plot[length(obs_to_plot)],
                      to_plot = obs_to_plot, col_values=cols,
                      labels=colnames(Y), plot_type="posterior", alpha=0.5)
post

# biplot
rownames(envfit$V)  <- 1:p
biplot  <- create_plots(envfit$V, cov_psamp,
                        n1=obs_to_plot[1], n2=obs_to_plot[length(obs_to_plot)],
                        to_plot = obs_to_plot, col_values=cols,
                        labels=colnames(Y), plot_type="biplot")
biplot



