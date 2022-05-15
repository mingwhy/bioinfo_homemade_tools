options(mc.cores = parallel::detectCores())
library(tidyverse)
library(mvtnorm)
library(rstiefel)
library(glmnet)
library(covreg)
library(envelopeR)
library(lubridate)
source('envelope_functions.r') #for get_mean_cov_hat()
source('utility_functions.r') #for steinsLoss()

start.time=Sys.time()
## Number of features
p <- 25
## Number of predictors?? number of samples/observations
n <- 100
gamma_sd  <- 1

## error sd
error_sd  <- 1

##beta_sd_vec <- 0:5
beta_sd_vec <- c(3)
nreps <- 1; #nreps <- 100

## Rank of the matrix
qmax <- 2; #qmax <- 4
## dimension of v
s <- 4
#sfit <- c(2:s, 2*s, 3*s, p)
#(sfit <- c(2,s, 2*s))
(sfit <- c(10,20))

create_cov <- function(X, gammaList, s, scaler=1) {
    sig_X <- matrix(0, ncol=s, nrow=s)
    for(i in 1:s) {
        gammaX <- gammaList[[i]] %*% t(X) #(equation12,gamma_k*x_i part)
        sig_X <- sig_X + tcrossprod(gammaX) #(equation12,crossprod)
    }
    sig_X
}

file_name <- sprintf("sim_s_comparison_%s.Rdata", lubridate::today())
subspace_sim_array  <- array(dim=c(length(sfit), nreps,  qmax, length(beta_sd_vec)))
steins_loss_array  <- array(dim=c(length(sfit), nreps, qmax, length(beta_sd_vec)))
squared_error_loss_array  <- array(dim=c(length(sfit), nreps, qmax, length(beta_sd_vec)))

sink('log.txt')
for(q_cur in 1:qmax) {
  for(rep in 1:nreps) {      
    beta_count  <- 1
    
    for(beta_sd in beta_sd_vec) {
      X <- matrix(rnorm(n*q_cur, 0, 1), nrow=n, ncol=q_cur)
      gammaList <- lapply(1:s, function(i) matrix(rnorm(s*q_cur, 0, sd=gamma_sd), nrow=s, ncol=q_cur))
      sapply(gammaList,dim)
      cov_list  <- lapply(1:n, function(i){create_cov(X[i, , drop=FALSE], gammaList, s, scaler=1) 
        + error_sd^2 * diag(s)} ) #equation 12, sigma^2*I part
      sapply(cov_list,dim) #simu 100reps, each with subsapce 4x4
      
      V  <- rustiefel(p, s) #p: nfeature, s:dim of subspace
      Vnull  <- rstiefel::NullC(V)
      dim(V) #25x4
      dim(Vnull) #25x21
      
      beta_mat <- matrix(rnorm(q_cur*s, 0, 1), nrow=q_cur)
      dim(beta_mat) # predictor.combo x subspace.dim
      beta <- beta_sd * beta_mat
      
      Z <- sapply(1:n, function(i) {
        rmvnorm(1, X[i, ] %*% beta, sigma=cov_list[[i]])
      }) %>% t
      
      Y  <- Z %*% t(V)  +
        matrix(rnorm(n * (p-s), sd=error_sd), nrow=n, ncol=p-s) %*% t(Vnull)
      
      s_count <- 1
      for(s_cur in sfit) {
        start.time=Sys.time()
        
        if(s_cur == p) {
          fit <- covreg.mcmc(Y ~ X - 1, Y ~ X, R=s_cur, verb=FALSE)
          mean_cov_hat <- get_mean_cov_hat(fit, inv=TRUE)
          
          steins_loss <- sapply(1:n, function(i) {
            a <- t(V) %*% mean_cov_hat[i, , ] %*% V
            b <- cov_list[[i]]
            steinsLoss(a, solve(b))
          }) %>% mean
          
          mean_cov_hat <- get_mean_cov_hat(fit)
          se_loss <- sapply(1:n, function(i) {
            a <- t(V) %*% mean_cov_hat[i, , ] %*% V
            b <- cov_list[[i]]
            sum((a-b)^2)
          }) %>% mean
          
          post_mean <- apply(m.psamp(fit), 1:2, mean)
          Vhat <- svd(post_mean)$v[, 1:s]
        }else{
          envfit <- fit_envelope(Y, X, distn="covreg", s=s_cur,
                                 fmean = "YV ~ X - 1",
                                 fcov = "YV ~ X",
                                 Vinit="OLS")
          Vhat <- envfit$V
          YV <- Y %*% Vhat
          Vhat_perp <- rstiefel::NullC(Vhat)
          fit <- covreg.mcmc(YV ~ X - 1, YV ~ X, R=s, verb=FALSE)
          
          cov_psamp <- cov.psamp(fit)
          
          YVp <- Y %*% Vhat_perp
          sigma2_hat <- sum(YVp^2)/(n*(p-s))
          
          mean_cov_hat <- get_mean_cov_hat(fit, inv=TRUE)
          steins_loss <- sapply(1:n, function(i) {
            a <- t(V) %*% Vhat %*% mean_cov_hat[i, ,] %*% t(Vhat) %*% V  + sigma2_hat * tcrossprod(t(V) %*% Vhat_perp)
            b <- cov_list[[i]]
            steinsLoss(a, solve(b))
          }) %>% mean
          
          mean_cov_hat <- get_mean_cov_hat(fit)
          se_loss <- sapply(1:n, function(i) {
            a <- t(V) %*% Vhat %*% mean_cov_hat[i, ,] %*% t(Vhat) %*% V  + sigma2_hat * tcrossprod(t(V) %*% Vhat_perp)
            b <- cov_list[[i]]
            sum((a-b)^2)
          }) %>% mean
        }
        subspace_sim  <- tr(Vhat  %*% t(Vhat)  %*% V  %*% t(V))/s
        
        subspace_sim_array[s_count, rep, q_cur, beta_count] <- subspace_sim
        steins_loss_array[s_count, rep, q_cur, beta_count] <- steins_loss
        squared_error_loss_array[s_count, rep, q_cur, beta_count] <- se_loss
        
        print(sprintf("--------- %i, %i, %i, %i: %f, %f, %f ----------------", s_cur, rep, q_cur, beta_sd, subspace_sim, steins_loss, se_loss))
        s_count <- s_count + 1
        
        end.time=Sys.time()
        cat('#rep',rep,'s_cur',s_cur,'q_cur',q_cur,'finished\n')
        print(start.time)
        print(end.time)
      }
      save(subspace_sim_array, steins_loss_array, squared_error_loss_array, file=file_name)
    }
    
    beta_count <- beta_count + 1
  }
  
  save(subspace_sim_array, steins_loss_array, squared_error_loss_array, file=file_name)
}
sink()

#30min for 5rep x 2q_cur, sfit = 2,4,8
#which is the slowest step? q_cur (predictor in X) or s_cur (rank of subsapce)?


