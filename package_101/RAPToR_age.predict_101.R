
#RAPToR (Real Age Prediction from Transcriptome staging on Reference) 
#https://github.com/LBMC/RAPToR
################################
## install 
# CRAN packages
#install.packages(c("ica", "mgcv", "parallel", "data.table", "pryr", "beeswarm", "Rdpack", "R.rsp"))
library(ica);library(mgcv)
#install.packages(c("pryr", "R.rsp"))
library(pryr)
library(R.rsp)
library(limma)
library(GEOquery)

# to use already built ref
#devtools::install_github("LBMC/wormRef")
library(wormRef)
list_refs("wormRef")

#devtools::install_github("LBMC/RAPToR", build_vignettes = T)
library(RAPToR)
vignette("RAPToR")

################################
## setup data: dsaeschimann2017
data_folder <- "~/Downloads/"

requireNamespace("wormRef", quietly = T)
requireNamespace("utils", quietly = T)
requireNamespace("GEOquery", quietly = T) # May need to be installed with bioconductor
requireNamespace("Biobase", quietly = T) # same

raw2tpm <- function(rawcounts, genelengths){
  if(nrow(rawcounts) != length(genelengths))
    stop("genelengths must match nrow(rawcounts).")
  x <- rawcounts/genelengths
  return(t( t(x) * 1e6 / colSums(x) ))
}

fpkm2tpm <- function(fpkm){
  return(exp(log(fpkm) - log(colSums(fpkm)) + log(1e6)))
}

geo_dsaeschimann2017 <- "GSE80157"

g_url_dsaeschimann2017 <- GEOquery::getGEOSuppFiles(geo_dsaeschimann2017, makeDirectory = FALSE, fetch_files = FALSE)
g_file_dsaeschimann2017 <- paste0(data_folder, "dsaeschimann2017.txt.gz")
utils::download.file(url = as.character(g_url_dsaeschimann2017$url[2]), destfile = g_file_dsaeschimann2017)

X_dsaeschimann2017 <- read.table(gzfile(g_file_dsaeschimann2017), h=T, sep = '\t', stringsAsFactors = F, row.names = 1)

# convert to tpm & wb_id
X_dsaeschimann2017 <- X_dsaeschimann2017[rownames(X_dsaeschimann2017)%in%wormRef::Cel_genes$wb_id,]
X_dsaeschimann2017 <- raw2tpm(rawcounts = X_dsaeschimann2017, 
                              genelengths = wormRef::Cel_genes$transcript_length[match(rownames(X_dsaeschimann2017),
                                                                                       wormRef::Cel_genes$wb_id)])

# pheno data
P_dsaeschimann2017 <- Biobase::pData(GEOquery::getGEO(geo_dsaeschimann2017, getGPL = F)[[1]])
P_dsaeschimann2017[,10:34] <- NULL
P_dsaeschimann2017[, 3:8] <- NULL

colnames(P_dsaeschimann2017)[4] <- "strain"
P_dsaeschimann2017$strain <- factor(P_dsaeschimann2017$strain)
P_dsaeschimann2017$title <- make.names(P_dsaeschimann2017$title)

colnames(X_dsaeschimann2017) <- gsub('RNASeq_riboM_', '', colnames(X_dsaeschimann2017), fixed = T)
P_dsaeschimann2017$title <- gsub('RNASeq_riboM_', '', P_dsaeschimann2017$title, fixed = T)

# get age 
P_dsaeschimann2017$age <- as.numeric(sub('(\\d+)\\shours', '\\1', P_dsaeschimann2017$`time in development:ch1`))


X_dsaeschimann2017 <- X_dsaeschimann2017[, P_dsaeschimann2017$title]

dsaeschimann2017 <- list(g = X_dsaeschimann2017, p = P_dsaeschimann2017)
save(dsaeschimann2017, file = paste0(data_folder, "dsaeschimann2017.RData"), compress = "xz")

# cleanup
file.remove(g_file_dsaeschimann2017)
rm(geo_dsaeschimann2017, g_url_dsaeschimann2017, g_file_dsaeschimann2017, X_dsaeschimann2017, P_dsaeschimann2017)
################################
## basic usage
# load data
names(dsaeschimann2017)
dim(dsaeschimann2017$g) #19595    43
dim(dsaeschimann2017$p) #43 6

dsaeschimann2017$g[1:5,1:4]
head(dsaeschimann2017$p, n = 5)

# quantile-normalize
dsaeschimann2017$g <- limma::normalizeBetweenArrays(dsaeschimann2017$g, method = "quantile")
dsaeschimann2017$g <- log1p(dsaeschimann2017$g) # log1p(x) = log(x + 1)

# choose and load a reference
plot_refs("wormRef")
r_larv <- prepare_refdata("Cel_larval", "wormRef", n.inter = 600)
dim( r_larv$interpGE) # 18718   600
length(r_larv$time.series) #600

# age estimation
ae_dsaeschimann2017 <- ae(samp = dsaeschimann2017$g,                         # input gene expression matrix
                        refdata = r_larv$interpGE, # reference gene expression matrix
                        ref.time_series = r_larv$time.series) 
# result
plot(ae_dsaeschimann2017, groups = dsaeschimann2017$p$strain, show.boot_estimates = T)
names(ae_dsaeschimann2017)
head(ae_dsaeschimann2017$age.estimates)

lm_dsaeschimann2017 <- lm(ae_dsaeschimann2017$age.estimates[,1] ~ dsaeschimann2017$p$age)
summary(lm_dsaeschimann2017)

# understanding the output
summary(ae_dsaeschimann2017)
head(ae_dsaeschimann2017$age.estimates)
#age.estimate, the global estimate for the sample (whole gene set).
#lb, ub, the lower and upper bounds of the bootstrapped age estimates’ confidence interval (Median Absolute Deviation).
#cor.score, the correlation score between the sample and reference at the age estimate.
par(mfrow=c(2,2))
plot_cor.ae(ae_dsaeschimann2017, subset = c(1,4,9,11))

################################
## metabolome
# build reference
X=readRDS('~/Documents/aging_clock/mz.filter.combat.rds')
dim(X) #feature by sample
apply(X,1,mean) #Data must not be gene-centered, as this destroys the relationship between gene levels within a sample.
X=as.data.frame(X)

inp=readRDS('~/Documents/aging_clock/raw_mz.rds')
meta=as.data.frame(inp$meta)
meta$age=meta$age_at_doc
dog_ref<- ge_im(X = X, 
              p = meta, 
              formula = "X ~ s(age, bs = 'ts') ", nc = 32)
dog_ref

# find #pc and specify model type
dim(X)
dog_pca <- stats::prcomp(t(X), center = TRUE, scale = FALSE, rank = 166)
nc <- sum(summary(dog_pca)$importance[3,] < .99) + 1
nc #93

smooth_methods <- c("tp", "ts", "cr", "ps")
flist <- as.list(paste0("X ~ s(age, bs = \'", smooth_methods, "\') ")) #, k=", k, ", fx=TRUE
flist

m_cv <- ge_imCV(X = X, p = meta, formula_list = flist,
                    cv.n = 20, nc = nc, nb.cores = 4)

plot(m_cv, names = paste0("bs = ", smooth_methods), outline = F,
     swarmargs = list(cex = .8))

# predict from the model 
# setup new data
n.inter <- 100
ndat <- data.frame(age = seq(min(meta$age), max(meta$age),  l = n.inter))
               
# predict
pred_m_comp <- predict(dog_ref, ndat, as.c = TRUE) # in component space
pred_m_ge <- predict(dog_ref, ndat)
dim(pred_m_comp) #100 x 32
dim(pred_m_ge) #166 x 100

# make a 'reference object' 
r_dog <- list(interpGE = pred_m_ge, 
              time.series = ndat$age)

# age prediction via ae() function
ae_test_dog <- ae(X, r_dog$interpGE, r_dog$time.series)
names(ae_test_dog)
head(ae_test_dog$age.estimates)

plot(meta$age,ae_test_dog$age.estimates[,1],pch=16)
abline(a=0,b=1,lty=5)

