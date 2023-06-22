#https://lbozhilova.github.io/COGENT/tutorial/tutorial.html
# install
if(F){
Sys.getenv("GITHUB_PAT")
Sys.unsetenv("GITHUB_PAT")
Sys.getenv("GITHUB_PAT")
devtools::install_github("lbozhilova/COGENT")
}

library(COGENT);library(ggplot2);
library(ggthemes);library(parallel);
set.seed(101019)

load("tutorialData.RData")
head(yeastData[,1:6])
dim(yeastData) #200gene x 27 samples
class(yeastData) #data.frame

buildPearson <- function(df, quant=0.90){
  check <- checkExpressionDF(df)
  A <- cor(t(df[,colnames(df)!="Name"]), use="pairwise.complete.obs", method="pearson")
  threshold <- quantile(A[upper.tri(A)], quant, na.rm=TRUE)
  A <- 1*(A>=threshold); diag(A) <- 0
  colnames(A) <- rownames(A) <- df$Name
  return(A)
}
pearsonA <- buildPearson(yeastData)
dim(pearsonA) #200 x 200
pearsonA[1:3,1:3]

buildKendall<- function(df, quant=0.90){
  check <- checkExpressionDF(df)
  A <- cor(t(df[,colnames(df)!="Name"]), use="pairwise.complete.obs", method="kendall")
  threshold <- quantile(A[upper.tri(A)], quant, na.rm=TRUE)
  A <- 1*(A>=threshold); diag(A) <- 0
  colnames(A) <- rownames(A) <- df$Name
  return(A)
}
kendallA <- buildKendall(yeastData)

PKcomparison <- getEdgeSimilarity(list(pearsonA, kendallA), align=FALSE)
PKcomparison$nodeCount
PKcomparison$globalSimilarity 
PKcomparison$localSimilarity

hist(PKcomparison$localSimilarity,
     main="Local similarity between Pearson and Kendall networks",
     xlab="Similarity",
     breaks=seq(0, 1, .05), col="cornflowerblue")

#COGENT evaluates network consistency by repeatedly splitting the set of samples of gene expression data in two, 
#and building a separate network from each sample subset. The more similar the resulting networks are, the more consistent the network construction method is judged to be. 
calculateDegrees <- function(A){
  deg <- rowSums(A)
  names(deg) <- colnames(A)
  return(deg)
}

x <- cogentSingle(df=yeastData, netwkFun=buildPearson, propShared=0.50)
x$globalSimilarity


x <- cogentSingle(df=yeastData, netwkFun=buildKendall, propShared=0.50,
                  nodeFun=calculateDegrees, nodeModes="all",
                  use="pairwise.complete.obs", method="pearson",
                  k.or.p=0.10,
                  scale=TRUE)
c("globalSimilarity"=x$globalSimilarity, "corDegrees"=x$corSimilarity)

stabilityPearson <- cogentLinear(df=yeastData, netwkFun=buildPearson, propShared=0.50,
                                 repCount=100, nodeFun=calculateDegrees, nodeModes="all",
                                 use="pairwise.complete.obs", method="pearson", 
                                 k.or.p=0.10, scale=TRUE)
head(stabilityPearson)
dim(stabilityPearson) #100 x 6

stabilityKendall <- cogentParallel(df=yeastData, netwkFun=buildKendall, propShared=0.50, 
                                   repCount=100, threadCount=2, nodeFun=calculateDegrees, 
                                   nodeModes="all", use="pairwise.complete.obs", 
                                   method="pearson", k.or.p=0.10, scale=TRUE)
head(stabilityKendall)

## Choosing a threshold
getThresholdStability <- function(th){
  dfSplit <- splitExpressionData(yeastData, propShared=0)
  A <- lapply(dfSplit, function(df) buildPearson(df, th))
  return(getEdgeSimilarityCorrected(A, type="expected")) 
}

aggregateThresholdStability <- function(th, repCount=100){
  thS <- replicate(repCount, getThresholdStability(th), simplify=FALSE)
  thS <- do.call("rbind", thS); thS <- apply(thS, 2, unlist)
  return(as.data.frame(cbind(thS, threshold=th)))
}
thresholds <- seq(0.5, 0.99, 0.01)
thresholdComparisonDF <- mclapply(thresholds, aggregateThresholdStability, mc.cores=6)
thresholdComparisonDF <- do.call("rbind", thresholdComparisonDF)

thresholdComparisonDF <- subset(thresholdComparisonDF, 
                                !is.na(thresholdComparisonDF$correctedSimilarity))
ggplot(thresholdComparisonDF, aes(x=1-threshold, y=correctedSimilarity)) +
  geom_smooth() +
  theme_economist_white() +
  ggtitle("Threshold choice for Pearson correlation networks") +
  scale_y_continuous("Density-adjusted consistency") +
  scale_x_continuous("Network density", breaks=seq(0, 0.5, .05))


