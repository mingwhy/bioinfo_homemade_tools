## parameter choice suggestion: https://support.bioconductor.org/p/p134550/
#https://bioconductor.org/packages/release/bioc/vignettes/zinbwave/inst/doc/intro.html

#BiocManager::install("zinbwave")
library(zinbwave)
library(scRNAseq) #remotes::install_github("LTLA/scRNAseq")
library(matrixStats)
library(magrittr)
library(ggplot2)
library(biomaRt)

# Register BiocParallel Serial Execution
BiocParallel::register(BiocParallel::SerialParam())

fluidigm <- ReprocessedFluidigmData(assays = "tophat_counts")
fluidigm
#class: SingleCellExperiment 

table(colData(fluidigm)$Coverage_Type)

# filter gene
filter <- rowSums(assay(fluidigm)>5)>5
table(filter)
fluidigm <- fluidigm[filter,]

# log1p transform data 
assay(fluidigm) %>% log1p %>% rowVars -> vars
names(vars) <- rownames(fluidigm)
vars <- sort(vars, decreasing = TRUE)
head(vars)

# use the top/highly 100 expressed genes
fluidigm <- fluidigm[names(vars)[1:100],]

assayNames(fluidigm)[1] <- "counts"

# fit model with k=2, unknown sample-level, or cell-level covariates
fluidigm_zinb <- zinbwave(fluidigm, K = 2, epsilon=1000)

W <- reducedDim(fluidigm_zinb)
data.frame(W, bio=colData(fluidigm)$Biological_Condition,
           coverage=colData(fluidigm)$Coverage_Type) %>%
  ggplot(aes(W1, W2, colour=bio, shape=coverage)) + geom_point() + 
  scale_color_brewer(type = "qual", palette = "Set1") + theme_classic()

#fit model with k=2, and add cell-level covariates
#Sample-level covariates
#The column Coverage_Type in the colData of fluidigm contains the coverage information.
tmp=colData(fluidigm)
colnames(tmp)
table(tmp$Coverage_Type)

fluidigm_cov <- zinbwave(fluidigm, K=2, X="~Coverage_Type", epsilon=1000)
W <- reducedDim(fluidigm_cov)
dim(W)
data.frame(W, bio=colData(fluidigm)$Biological_Condition,
           coverage=colData(fluidigm)$Coverage_Type) %>%
  ggplot(aes(W1, W2, colour=bio, shape=coverage)) + geom_point() + 
  scale_color_brewer(type = "qual", palette = "Set1") + theme_classic()


#fit model with k=2, and add gene-level covariates
#Gene-level covariates
mart <- useMart("ensembl")
mart <- useDataset("hsapiens_gene_ensembl", mart = mart)
bm <- getBM(attributes=c('hgnc_symbol', 'start_position',
                         'end_position', 'percentage_gene_gc_content'),
            filters = 'hgnc_symbol',
            values = rownames(fluidigm),
            mart = mart)

bm.bp=bm;
fluidigm.bp=fluidigm;

dim(fluidigm);dim(bm)
sum(rownames(fluidigm) %in% bm$hgnc_symbol)
sum(bm$hgnc_symbol %in% rownames(fluidigm))
length(unique(bm$hgnc_symbol ))
bm[duplicated(bm$hgnc_symbol),]
bm[bm$hgnc_symbol=='LDHA',]

bm=bm[!duplicated(bm$hgnc_symbol),]
fluidigm=fluidigm[rownames(fluidigm) %in% bm$hgnc_symbol,]
dim(fluidigm)
dim(bm)

bm$length <- bm$end_position - bm$start_position
len <- tapply(bm$length, bm$hgnc_symbol, mean)
len <- len[rownames(fluidigm)]
gcc <- tapply(bm$percentage_gene_gc_content, bm$hgnc_symbol, mean)
gcc <- gcc[rownames(fluidigm)]
rowData(fluidigm) <- data.frame(gccontent = gcc, length = len)
fluidigm
head(rowData(fluidigm))

fluidigm_gcc <- zinbwave(fluidigm, K=2, V="~gccontent + log(length)", epsilon=1000)

## instead of using a wrapper function, use zinbFit to get more model fitting information
#https://github.com/drisso/zinbwave/blob/master/R/zinb_fit.R
dim(fluidigm) #99gene x 130 cell
zinb <- zinbFit(fluidigm, K=2, V="~gccontent + log(length)", epsilon=1000,verbose=T)
zinb <- zinbFit(fluidigm, K=2,  X="~Coverage_Type", epsilon=1000,verbose=T)

#https://bioconductor.org/packages/devel/bioc/vignettes/zinbwave/inst/doc/intro.html
#Here, we also specify observationalWeights = TRUE to compute observational weights, useful for differential expression (see next section).
fluidigm_zinb <- zinbwave(fluidigm, fitted_model = zinb, K = 2, epsilon=1000,
                          observationalWeights = TRUE)

#The zinbwave package can be used to compute observational weights to “unlock” bulk RNA-seq tools for single-cell applications, as illustrated in (Van den Berge et al. 2018).
#zinbwave optionally computes the observational weights when specifying observationalWeights = TRUE as in the code chuck above. 
#See the man page of zinbwave. The weights are stored in an assay named weights and can be accessed with the following call.

# use these weights as posterior probability for 'measurement confidence'
# https://genomebiology.biomedcentral.com/articles/10.1186/s13059-018-1406-4#Sec2
#... From the ZINB-WaVE density of Eq. 1, one can read- ily derive the posterior probability that a count yij was generated from the NB count component: ...
weights <- assay(fluidigm_zinb, "weights")
dim(fluidigm); #99 130
dim(weights); #99 130 
summary(as.numeric(weights))
#    Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 0.001206 0.798709 1.000000 0.788453 1.000000 1.000000 
apply(weights,2,function(i){sum(i<0.99)}) #if require posterior probability >=0.99, confidence cutoff
hist(apply(weights,2,function(i){sum(i<0.99)}))

assayNames(Y)  
dataY <- assay(Y, "counts")
dim(dataY) # 99 130

################################################
## look at the output fitting result
head(zinb@X) #130 cell, design matrix for cell
dim(zinb@V) #99gene x 3 attributes
head(zinb@V) #covarite for each gene
dim(zinb@O_mu) #130cell by 99 gene
dim(zinb@O_pi)
zinb@O_mu[1:3,1:3]
zinb@O_pi[1:3,1:3]

dim(zinb@X) #130 x 2
dim(zinb@beta_mu) #2 x 99
dim(zinb@V) #99x1
dim(zinb@gamma_mu) #1x130
dim(zinb@W) #130x2 (k=2)
dim(zinb@alpha_mu) #2x99
dim(zinb@O_mu)

# based on paper method: https://www.nature.com/articles/s41467-017-02554-5#Sec12
log.u.ij=zinb@X %*% zinb@beta_mu + 
  t(zinb@V %*% zinb@gamma_mu) +
  zinb@W %*% zinb@alpha_mu +
  zinb@O_mu
dim(log.u.ij) #130 x 99
u.ij=t(exp(log.u.ij)) #99x130
u.ij[1:3,1:3]



logit.pi.ij= zinb@X %*% zinb@beta_pi + 
  t(zinb@V %*% zinb@gamma_pi) +
  zinb@W %*% zinb@alpha_pi +
  zinb@O_pi
pi.ij=exp(logit.pi.ij)/(1+exp(logit.pi.ij))
pi.ij=t(pi.ij)
dim(pi.ij) #99 by 130
pi.ij[1:3,1:3]
max(pi.ij)
min(pi.ij)
## use zinbwave to estiamte pi and calculate posterior.prob of non-detection or drop-out rate
## can you simulate data based on estimated values to confirm model fitting

############################
dim(zinb@beta_pi) # 2 x 99
dim(zinb@gamma_pi) #1x130
dim(zinb@alpha_pi) #2x99
dim(zinb@O_pi) #130 x 99

length(zinb@zeta) #99 genes


zinb@epsilon_min_logit


###########################
# 5, t-SNE representation
set.seed(93024)

library(Rtsne)
W <- reducedDim(fluidigm_cov)
tsne_data <- Rtsne(W, pca = FALSE, perplexity=10, max_iter=5000)

data.frame(Dim1=tsne_data$Y[,1], Dim2=tsne_data$Y[,2], 
           bio=colData(fluidigm)$Biological_Condition,
           coverage=colData(fluidigm)$Coverage_Type) %>%
  ggplot(aes(Dim1, Dim2, colour=bio, shape=coverage)) + geom_point() + 
  scale_color_brewer(type = "qual", palette = "Set1") + theme_classic()

#6 Normalized values and deviance residuals
fluidigm_norm <- zinbwave(fluidigm, K=2, epsilon=1000, normalizedValues=TRUE,
                          residuals = TRUE)
fluidigm_norm

#7 The zinbFit function
zinb <- zinbFit(fluidigm, K=2, epsilon=1000)
fluidigm_zinb <- zinbwave(fluidigm, fitted_model = zinb, K = 2, epsilon=1000,
                          observationalWeights = TRUE)

#8 Differential Expression
#Differential expression with DESeq2
library(DESeq2)

dds <- DESeqDataSet(fluidigm_zinb, design = ~ Biological_Condition)

dds <- DESeq(dds, sfType="poscounts", useT=TRUE, minmu=1e-6)
res <- lfcShrink(dds, contrast=c("Biological_Condition", "NPC", "GW16"),
                 type = "normal")
head(res)


weights <- assay(fluidigm_zinb, "weights")

