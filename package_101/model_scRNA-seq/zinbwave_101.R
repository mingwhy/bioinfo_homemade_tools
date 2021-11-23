#https://bioconductor.org/packages/release/bioc/vignettes/zinbwave/inst/doc/intro.html

#BiocManager::install("zinbwave")

library(zinbwave)
library(scRNAseq)
library(matrixStats)
library(magrittr)
library(ggplot2)
library(biomaRt)

# Register BiocParallel Serial Execution
BiocParallel::register(BiocParallel::SerialParam())

fluidigm <- ReprocessedFluidigmData(assays = "tophat_counts")
fluidigm

table(colData(fluidigm)$Coverage_Type)

filter <- rowSums(assay(fluidigm)>5)>5
table(filter)

fluidigm <- fluidigm[filter,]

assay(fluidigm) %>% log1p %>% rowVars -> vars
names(vars) <- rownames(fluidigm)
vars <- sort(vars, decreasing = TRUE)
head(vars)

fluidigm <- fluidigm[names(vars)[1:100],]

assayNames(fluidigm)[1] <- "counts"

fluidigm_zinb <- zinbwave(fluidigm, K = 2, epsilon=1000)

W <- reducedDim(fluidigm_zinb)

data.frame(W, bio=colData(fluidigm)$Biological_Condition,
           coverage=colData(fluidigm)$Coverage_Type) %>%
  ggplot(aes(W1, W2, colour=bio, shape=coverage)) + geom_point() + 
  scale_color_brewer(type = "qual", palette = "Set1") + theme_classic()
#Adding covariates

#Sample-level covariates
fluidigm_cov <- zinbwave(fluidigm, K=2, X="~Coverage_Type", epsilon=1000)
W <- reducedDim(fluidigm_cov)
dim(W)

data.frame(W, bio=colData(fluidigm)$Biological_Condition,
           coverage=colData(fluidigm)$Coverage_Type) %>%
  ggplot(aes(W1, W2, colour=bio, shape=coverage)) + geom_point() + 
  scale_color_brewer(type = "qual", palette = "Set1") + theme_classic()

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
# is there model likelihood value?
# blup, how to compare two model, with BIC or likelihood est?

#https://github.com/drisso/zinbwave/blob/master/R/zinb_fit.R
dim(fluidigm) #99gene x 130 cell
zinb <- zinbFit(fluidigm, K=2, V="~gccontent + log(length)", epsilon=1000,verbose=T)

class(fluidigm)
Y=fluidigm
K=2; V="~gccontent + log(length)"; epsilon=1000;verbose=T
commondispersion=TRUE;zeroinflation = TRUE
nb.repeat.initialize=2; maxiter.optimize=25
stop.epsilon.optimize=.0001;
BPPARAM=BiocParallel::bpparam()

assayNames(Y)  
dataY <- assay(Y, "counts")
dim(dataY) # 99 130

if(!is.matrix(V)) {
  tryCatch({
    f <- as.formula(V)
    V <- model.matrix(f, data=rowData(Y))
  },
  error = function(e) {
    stop("V must be a matrix or a formula with variables in rowData(Y)")
  })
}
f
dim(V) #99 x 3

# Apply zinbFit on the assay of SummarizedExperiment
res <- zinbFit(dataY, V, K, commondispersion, zeroinflation,
               verbose, nb.repeat.initialize, maxiter.optimize,
               stop.epsilon.optimize, BPPARAM)

################################################
head(zinb@X) #130 cell, design matrix for cell
dim(zinb@V) #99gene x 3
head(zinb@V) #covarite for each gene
dim(zinb@O_mu) #130cell by 99 gene
dim(zinb@O_pi)
zinb@O_mu[1:3,1:3]
zinb@O_pi[1:3,1:3]

dim(zinb@W) #130 x 2
zinb@beta_mu
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