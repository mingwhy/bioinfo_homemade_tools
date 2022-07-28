
#https://github.com/churchill-lab/scRATE

#devtools::install_github('churchill-lab/scRATE')
library(scRATE)
#devtools::install_github(repo = "mojaveazure/loomR", ref = "develop",force=T)
library(loomR)

#download example data: ftp://churchill-lab.jax.org/analysis/scRATE
ds <- connect('DC-like_cells.loom')
cntmat <- t(ds$matrix[,])
gsymb <- ds$row.attrs$GeneID[]
ds$close_all()

head(gsymb)
exposure <- log(colSums(cntmat))
head(exposure)

gg <- 4153
gsymb[gg] 
y <- cntmat[gg,]
y 

gexpr <- data.frame(y=y, exposure=exposure)   
#gexpr <- data.frame(y=y, exposure=exposure, celltype, sex)
gexpr

# model_fit <- fit_count_models(gexpr)

## there was some error message running above fit_count_models()
## I check source code on: https://rdrr.io/github/churchill-lab/scRATE/src/R/fit_count_models.R

## error message was generated when fitting model 3 and model 4
## which can be reproduced as below:
formula_string=NULL;
nCores=2;seed=1004; adapt_delta=0.8;model2fit=NULL

covariates <- names(gexpr)[-c(1, 2)]
if(is.null(formula_string)) {
   message(sprintf("Formulating the default additive model..."))
   formula_string <- 'y ~ 1'
   if(!identical(covariates, character(0))) {
      for (covar in covariates) {
         formula_string <- paste(formula_string, sprintf(' + (1|%s)', covar))
      }
   }
}

message(sprintf("Formula: %s", formula_string))
f12 <- as.formula(formula_string)
f34 <- as.formula(paste(formula_string, ' + offset(exposure)'))

fitting <- list()

message('Fitting data with Zero-Inflated Poisson model...')
myprior_3 <- get_prior(bf(f34, zi ~ 1),
                       family = zero_inflated_poisson(),
                       data = gexpr)
myprior_3_values <- eval(parse(text=gsub("student_t", "c", myprior_3$prior[1])))
fitting[["ZIP"]] <- brm(bf(f34, zi ~ 1),
                        family = zero_inflated_poisson(),
                        data = gexpr,
                        prior = myprior_3,
                        control = list(adapt_delta = adapt_delta),
                        cores = nCores,
                        seed = seed,
                        refresh = 500)
# error message: Error in unserialize(socklist[[n]]) : error reading from connection


#I found an online answer: https://discourse.mc-stan.org/t/persistent-error-message-error-in-unserialize-socklist-n-error-reading-from-connection/24184/2
# suggested using: options(brms.backend = "cmdstanr")
# I need to install 'cmdstanr' first

############################################
# install cmdstanr
#https://mc-stan.org/cmdstanr/articles/cmdstanr.html
#https://githubmemory.com/repo/stan-dev/cmdstanr/issues/552
#remotes::install_github("stan-dev/cmdstanr")
# test cmdstanr
library(cmdstanr)
library(posterior)
library(bayesplot)
color_scheme_set("brightblue")

check_cmdstan_toolchain()
install_cmdstan(cores = 2)
#* Finished installing CmdStan to /Users/ming/.cmdstan/cmdstan-2.28.1
#CmdStan path set to: /Users/ming/.cmdstan/cmdstan-2.28.1
# check path
cmdstan_path()
cmdstan_version() #2.28.1

file <- file.path(cmdstan_path(), "examples", "bernoulli", "bernoulli.stan")
mod <- cmdstan_model(file)
mod$print()
mod$exe_file()
data_list <- list(N = 10, y = c(0,1,0,0,0,0,0,0,0,1))

fit <- mod$sample(
   data = data_list, 
   seed = 123, 
   chains = 4, 
   parallel_chains = 4,
   refresh = 500
)
fit$summary()

###########################################################################
## set up brms.backend and test subfunciton in `fit_count_models`
## on my mac laptop, after installing 'cmdstanr'. directly 'library(scRATE)` 
## led Rstudio to quite automatically. 
## My solution is: options(xxx'cmdstanr') first, then load scRATE library

options(brms.backend = "cmdstanr")
library(scRATE)

message('Fitting data with Zero-Inflated Poisson model...')
myprior_3 <- get_prior(bf(f34, zi ~ 1),
                       family = zero_inflated_poisson(),
                       data = gexpr)
myprior_3_values <- eval(parse(text=gsub("student_t", "c", myprior_3$prior[1])))
fitting[["ZIP"]] <- brm(bf(f34, zi ~ 1),
                        family = zero_inflated_poisson(),
                        data = gexpr,
                        prior = myprior_3,
                        control = list(adapt_delta = adapt_delta),
                        cores = nCores,
                        seed = seed,
                        refresh = 500)

#########################################
## go back to original function in scRATE
model_fit <- fit_count_models(gexpr)
model_fit
model_fit$NB
model_fit$NB$fitted.values

elpd_loo <- compare_count_models(model_fit)
elpd_loo   

select_model(elpd_loo, margin=1)

