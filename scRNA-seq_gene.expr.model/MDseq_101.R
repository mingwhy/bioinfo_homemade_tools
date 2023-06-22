
library(MDSeq)
data(sampleData)

# raw count table 
dat <- sample.exprs
dim(dat)

dat[1:6, 1:6]


# covariates
X <- sample.pheno[,c("X1","X2")]
head(X)
        

# group information
group <- sample.pheno$group
summary(factor(group))
 


# Filtering and Normalization

#We could filter out lowly expressed genes that provide little information for gene expression analysis. The function *filter.counts()* is provided in the MDSeq to filter lowly expressed genes. By default, it will remove genes whose mean counts per million (cpm) is less than 0.05. cpm is calculated using the edgeR package [@robinson2010edger]. In this example, we applied a more stringent threshold with mean cpm greater than 0.1. As a result, 21992 high-quality genes were obtained from the original 49722 genes. We highly recommend users to remove lowly expressed genes for robust analysis of gene expression mean and variability. Different criteria or other methods can also be applied. For example, one could filter out more genes by setting a higher cutoff point of cpm; or one could include more genes by setting a smaller cutoff point. 

#Since the MDSeq evaluates relative differences between treatments, results may be effected by technical biases due to compositional differences among libraries, and technical biasness should be taken into account by normalization. In this example, we computed the TMM normalization factors [@robinson2010edger].  The normalization factors can be incorporated, in the MDSeq, as either normalized counts or offsets in the mean-dispersion GLMs.  The MDSeq provides the function *normalize.counts()* which can calculate normalized counts based on different normalization methods, such as "upperquartile" [@bullard2010evaluation] and "RLE" [@anders2010differential]. User could also apply more complex normalization approaches such as the "cqn", which considers adjustment for both gene length and GC content. This optional method depends on the implementation in the R package "cqn" [@hansen2012removing]. Offsets in the GLMs can also be applied in the function *MDSeq()* by incorporating normalization factors with the *offsets* parameter.  In the analysis of large-scale RNA-seq data, the approaches of using normalized counts and offsets in the GLMs are largely equivalent, as effects of discretizing normalized counts are minute. However, in computing the mean-dispersion model, the dispersion of normalized counts are assumed to be overdispersed with $\phi>1$ when normalized counts are used, whereas the dispersion of raw counts are assumed to be overdispersed with $\phi_{raw}>1$ when offsets in the GLMs are applied.  Further details are provided in Supplementary Methods of Ran and Daye (2017).

# lowly expressed genes were filtered by mean cpm value across all samples
dat.filtered <- filter.counts(dat, mean.cpm.cutoff = 0.1)  
dim(dat.filtered)


# Using normalized counts 
dat.normalized <- normalize.counts(dat.filtered, group=group, method="TMM")
#Using TMM normaliztion.
#Calculating normalization(scaling) factors with the TMM method.
#Adjusting counts to rescaling factors.

# including normalization factor as an offset in the mean-dispersion GLMs
cnf <- calcNormFactors(dat.filtered, method="TMM") 
libsize <- colSums(dat.filtered)              #normalization factor
rellibsize <- libsize/exp(mean(log(libsize))) #relative library size
nf <- cnf * rellibsize                        #normalization factor including library size


# Design matrix and constrast

#The MDSeq allows users to specify the design matrix using various contrast settings. Function *get.model.matrix()* allows one to generate the design matrix and form the contrast. For example, in the analysis of the sample dataset, the goal is to make a comparison between two treatments or groups, which is the most common and simple situation in biological research. It is necessary to form a contrast between the case and control groups and assign it to each sample within its group in the design matrix. First of all, treatment/group indicator needs to be in factor formats with predefined levels and labels. User could use the default contrast matrix or provide any proper contrast matrix. The default considers the contrast setting to make comparison such as $"-1*control + 1*case"$. It can also be extended to multi-group comparisons, greater than or equal to 3 groups. User can compare any of the treatment groups using different contrast settings. For example, a three-groups comparison can use a sum to zero contrast (default in MDSeq) as follows.

get.model.matrix(factor(rep(1:3,each=4)))

get.model.matrix(factor(rep(1:3,each=4)), contrast.type = "contr.treatment")


# treatment group assignment
group <- factor(sample.pheno$group, labels = c("Control", "Case"))
table(group)

# make design matrix with proper contrast setting
groups <- get.model.matrix(group)
groups


# Checking and removing outliers

#Outlier detection is an important step in RNA-seq studies to allow robust model estimation and inference. The MDSeq provides a computationally efficient procedure to detect outliers that are influential for statistical inference on a given set of parameters of interest [@MDSeq2016]. The function *remove.outlier()* will remove outliers and provide a summary of outlier detection.  A count matrix of RNA-seq counts will be outputted, in which outliers are set to NA.  A summary indicates the status (label=0 for successful detection) and numbers of outliers found at each gene.  Outliers will be ignored in ensuing analyses.

#After outlier detection, outliers within a gene were replaced by *NA*.  The mean-dispersion GLM utilizes the observed Fisher information to compute the variance-covariance matrix.  In some rare scenarios, the observed Fisher information may be non-invertible.  An error is flagged in these scenarios, even though an estimate may still be available.  We recommend removing those genes with nonzero error flags from downstream analysis. Using the first 1,000 genes, we removed 2 genes and obtained 998 genes. Outliers were replaced by *NA* in the output.


# Check outliers using parallel process with 4 threads
# For the sake of simplicity, we demonstrate outlier detection 
# by using the first 1,000 genes.

dat.checked <- remove.outlier(dat.normalized[1:1000, ], X=X, U=X, 
                              contrast = groups, mc.cores = 4)

# status of outlier checking
table(dat.checked$outliers$status)
  
# frequency distribtuion of outliers
table(dat.checked$outliers$num.outliers)


# remove genes with status flag other than 0
counts <- dat.checked$count[dat.checked$outliers$status==0,]
dim(counts)


# Differential mean and dispersion analysis

#Differential mean and dispersion analysis can be performed simultaneously using the main function *MDSeq()*, which allows parallel processing with multiple threads. For experiments with many genes, one can take advantage of parallelized computation. For instance, an example data with 200 samples and around 20,000 genes will cost less than 30 minutes with 4 threads on a lab desktop. MDSeq requires at least several inputs, such as expression count matrix, optional covariates of the mean(X) and dispersion(U), design matrix with proper contrast setting. By default, this function will simultaneously test zero-inflation using a likelihood ratio test. One can also assess the goodness-of-fit of zero-inflation using function *gof.ZI()*. MDSeq will return an "ZIMD" class object. In this example, estimations are stored in *fit$Dat*, which contains not only the estimations of GLM coefficients but also the results for testing zero-inflation. For example, one can find "s" the estimated proportion of excess zeros, "ZI.pval" p-value of likelihood ratio test, and "ZI" zero-inflation indicator. 

# using parallel process with 4 threads
fit <- MDSeq(counts, X=X, U=X, contrast = groups, mc.cores = 4)

# simultaneously test zero-inflation
colnames(fit$Dat)
head(fit$Dat[,c(1:8)])
head(fit$Dat[c("s","ZI.pval", "ZI")], 20)

## Testing with a given |log2|-fold-change

#It is often of interest to identify genes with differential changes beyond a given threshold level that could more readily allow for experimental replication and biological interpretations. Traditional methods of statistical inference often focus on evaluating the compliant hypothesis that any change in differential expressions $Ha: |log2FC| \neq 0$ may occur. This would often result in the selection of a large proportion of genes that are only mildly differentially expressed and cannot be replicated experimentally, especially in RNA-seq studies with moderate to large numbers of samples. Therefore, we provided inequality hypothesis tests $Ha: |log2FC| > \tau$ for both differential mean and dispersion that evaluate whether absolute log-fold changes are above a given threshold level $\tau$. We applied a rigorous development based on one-sided hypothesis tests within restricted parameter spaces and union-intersection principle [@MDSeq2016]. The function *extract.ZIMD()* is provided in the MDSeq that incorporates inequality hypothesis tests. One can specify any |log2|-fold-change threshold "log2FC.threshold", which will be applied on both mean and dispersion testing. In addition. this function offers the comparison between any possible single pair. One can test any comparison, such as "2 vs. 1" or "3 vs. 2" in a three-levels treatment, for example. To make a comparison  "2 vs. 1", one needs to specify "A='2'" and "B='1'" in the "compare" argument.


# given log2-fold-change threshold = 1
result <- extract.ZIMD(fit, compare = list(A="Case", B="Control"), log2FC.threshold = 1)
head(result)

