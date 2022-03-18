
#https://uc-r.github.io/random_forests
library(rsample)      # data splitting 
library(randomForest) # basic implementation
library(ranger)       # a faster implementation of randomForest
library(caret)        # an aggregator package for performing many machine learning models
library(h2o)          # an extremely fast java-based platform

# Create training (70%) and test (30%) sets for the AmesHousing::make_ames() data.
# Use set.seed for reproducibility

set.seed(123)
ames_split <- initial_split(AmesHousing::make_ames(), prop = .7)
ames_train <- training(ames_split)
ames_test  <- testing(ames_split)

dim(ames_train) #2051 81
sort(colnames(ames_train))
ncol(ames_train)/3 #26 features randomly subset every time

# what to predict: a continuous variable
summary(ames_train$Sale_Price)
set.seed(123)
# default RF model
m1 <- randomForest(
  formula = Sale_Price ~ .,
  data    = ames_train
)
m1
plot(m1)
which.min(m1$mse) #number of trees with lowest MSE
sqrt(min(m1$mse))


OOB_RMSE <- vector(mode = "numeric", length = 100)

for(i in seq_along(OOB_RMSE)) {

  optimal_ranger <- ranger(
    formula         = Sale_Price ~ ., 
    data            = ames_train, 
    num.trees       = 500,
    mtry            = 24,
    min.node.size   = 5,
    sample.fraction = .8,
    importance      = 'impurity'
  )
  
  OOB_RMSE[i] <- sqrt(optimal_ranger$prediction.error)
}

hist(OOB_RMSE, breaks = 20)


optimal_ranger$variable.importance %>% 
  tidy() %>%
  dplyr::arrange(desc(x)) %>%
  dplyr::top_n(25) %>%
  ggplot(aes(reorder(names, x), x)) +
  geom_col() +
  coord_flip() +
  ggtitle("Top 25 important variables")
###############################################
## https://stats.stackexchange.com/questions/56895/do-the-predictions-of-a-random-forest-model-have-a-prediction-interval
# predict a continuous variable
set.seed(1)

x1 <- rep(0:1, each=500)
x2 <- rep(0:1, each=250, length=1000)

y <- 10 + 5*x1 + 10*x2 - 3*x1*x2 + rnorm(1000)

fit1 <- lm(y ~ x1 * x2)

newdat <- expand.grid(x1=0:1, x2=0:1)

(pred.lm.ci <- predict(fit1, newdat, interval='confidence'))
(pred.lm.pi <- predict(fit1, newdat, interval='prediction'))

library(randomForest)
fit2 <- randomForest(y ~ x1 + x2, ntree=1001)

pred.rf <- predict(fit2, newdat, predict.all=TRUE)

pred.rf.int <- apply(pred.rf$individual, 1, function(x) {
  c(mean(x) + c(-1, 1) * sd(x), 
    quantile(x, c(0.025, 0.975)))
})

t(pred.rf.int)

# Create random forest
# For classification
library(randomForest)
iris.rf <- randomForest(Species ~ ., 
                        data = iris, 
                        importance = TRUE,
                        proximity = TRUE)

# Print classification model
print(iris.rf)

#########################################
## read in gene expression data and sample meta information 
expr.mat=data.table::fread('../2019_paper_reproduce.result/gene_by_sample_log2TPM.txt')
expr.mat=as.matrix(expr.mat,rownames=1)

sample.meta=data.table::fread('../2019_paper_reproduce.result/sample.meta_sex.label.txt')

dim(expr.mat) #8934, 54
dim(sample.meta) #54, 4
sum(colnames(expr.mat)==sample.meta$GSM.id) #54


# pick female and male samples
expr.mat=as.data.frame(expr.mat)
sample.meta=as.data.frame(sample.meta)
expr.mat=expr.mat[,sample.meta$cluster!='PB']
sample.meta=sample.meta[sample.meta$cluster!='PB',]

dim(expr.mat) #gene by sample
input.dat=as.data.frame(t(expr.mat))
input.dat[1:3,1:3]
input.dat$sex=sample.meta$cluster
input.dat$sex=factor(input.dat$sex)

rf2019 <- randomForest(sex ~ ., 
                        data = input.dat, 
                        importance = TRUE,
                        proximity = TRUE)
rf2019
model=rf2019

# feature importance in random forest
#https://www.statistik.uni-dortmund.de/useR-2008/slides/Strobl+Zeileis.pdf
# Gini importance
model$importance
importance(model,type=2)

# permutation importance
#mean decrease in classification accuracy after permuting Xj over all trees
importance(model,type=1)

# prediction confidence
# continous case: Do the predictions of a Random Forest model have a prediction interval?
# https://stats.stackexchange.com/questions/56895/do-the-predictions-of-a-random-forest-model-have-a-prediction-interval
# multiclass case: Using randomForest package in R, how to get probabilities from classification model?
# https://stackoverflow.com/questions/25715502/using-randomforest-package-in-r-how-to-get-probabilities-from-classification-mo
# https://stackoverflow.com/questions/28641298/what-do-xtest-and-ytest-do-in-the-randomforest-algorithm-in-r
# https://uc-r.github.io/random_forests

library(rsample)      # data splitting 
library(randomForest) # basic implementation
library(ranger)       # a faster implementation of randomForest
library(caret)        # an aggregator package for performing many machine learning models
library(h2o)          # an extremely fast java-based platform

# create training and validation data 
set.seed(123)
valid_split <- initial_split(input.dat, .8)
ames_train <- training(valid_split)
ames_test  <- testing(valid_split)

# training data
ames_train_v2 <- analysis(valid_split)

# validation data
ames_valid <- assessment(valid_split)
x_test <- ames_valid[setdiff(names(ames_valid), "sex" )]
y_test <- ames_valid$sex

model <- randomForest(formula=sex ~ ., 
                      data = ames_train_v2, 
                      xtest=x_test,
                      ytest=y_test,
                      keep.forest=TRUE,
                      importance = TRUE,
                      proximity = TRUE)
prob=predict(model,input.dat,type="prob")
#model=randomForest(x,y,xtest=x,ytest=y,keep.forest=TRUE). 
#prob=predict(model,x,type="prob")
dim(prob) #49 x 2
sum(prob[,1]>0.5) #29 assigned to female
sum(prob[,2]>0.5) #20 assigned to male


