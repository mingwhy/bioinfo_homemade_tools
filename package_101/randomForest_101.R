
###############################################
#UC Business Analytics R Programming Guide
#https://uc-r.github.io/random_forests
library(rsample)      # data splitting 
library(randomForest) # basic implementation
library(ranger)       # a faster implementation of randomForest
library(caret)        # an aggregator package for performing many machine learning models
library(h2o)          # an extremely fast java-based platform

library(tidyverse)
# Create training (70%) and test (30%) sets for the AmesHousing::make_ames() data.
# Use set.seed for reproducibility

set.seed(123)
dat=AmesHousing::make_ames()
dim(dat) #2930 obs x 81 feature
nrow(dat) * 0.7 #2051 for trainning
ames_split <- rsample::initial_split(AmesHousing::make_ames(), prop = .7)
ames_split
#<Analysis/Assess/Total>
#  <2051/879/2930>

ames_train <- training(ames_split)
ames_test  <- testing(ames_split)

dim(ames_train) #2051 81
ncol(ames_train)/3 #~27 features randomly subset every time

# what to predict: a continuous variable
summary(ames_train$Sale_Price)
set.seed(123)
# default RF model
m1 <- randomForest(
  formula = Sale_Price ~ .,
  data    = ames_train
)
m1
plot(m1) #OOB error
which.min(m1$mse) #number of trees with lowest MSE
# 500
sqrt(min(m1$mse))
# 25924.37

# randomForest also allows us to use a validation set to measure predictive accuracy if we did not want to use the OOB samples. 
# create training and validation data 
nrow(ames_train)*0.8 #1640 for training
set.seed(123)
valid_split <- rsample::initial_split(ames_train, .8)
valid_split
#<1640/411/2051>

# training data
ames_train_v2 <- rsample::analysis(valid_split)
dim(ames_train_v2) #1640   81
ames_train_v2

# validation data
ames_valid <- assessment(valid_split)
dim(ames_valid) #411 for validation
x_test <- ames_valid[setdiff(names(ames_valid), "Sale_Price")] #as 'Sale_Price' to be predicted
y_test <- ames_valid$Sale_Price
dim(x_test) #411 obs x 80 predictors
length(y_test) #411 y to be predicted

rf_oob_comp <- randomForest(
  formula = Sale_Price ~ .,
  data    = ames_train_v2,
  xtest   = x_test,
  ytest   = y_test
)

# extract OOB & validation errors
oob <- sqrt(rf_oob_comp$mse)
validation <- sqrt(rf_oob_comp$test$mse)
length(oob) #500, one with one ntree value

# compare error rates 

tibble::tibble(
  `Out of Bag Error` = oob,
  `Test error` = validation,
  ntrees = 1:rf_oob_comp$ntree
) %>%  #gather, wide to long
  gather(Metric, RMSE, -ntrees) %>%
  ggplot(aes(ntrees, RMSE, color = Metric)) +
  geom_line() +
  scale_y_continuous(labels = scales::dollar) +
  xlab("Number of trees")

## Tuning, mtry, #var used in each bootstrap
# names of features
features <- setdiff(names(ames_train), "Sale_Price")
features #80 predictors
#tuneRf will start at a value of mtry that you supply and increase by a certain step factor until the OOB error stops improving be a specified amount. For example, the below starts with mtry = 5 and increases by a factor of 1.5 until the OOB error stops improving by 1%. 
set.seed(123)
m2 <- tuneRF(
  x          = ames_train[features],
  y          = ames_train$Sale_Price,
  ntreeTry   = 500,
  mtryStart  = 5,
  stepFactor = 1.5,
  improve    = 0.01,
  trace      = FALSE      # to not show real-time progress 
)
plot(m2) #optimal mtry ~ 15

#Full grid search with ranger, which is faster than randomForest
# randomForest speed
system.time(
  ames_randomForest <- randomForest(
    formula = Sale_Price ~ ., 
    data    = ames_train, 
    ntree   = 500,
    mtry    = floor(length(features) / 3)
  )
)
#user  system elapsed 
#27.788   0.051  27.828 
# ranger speed
system.time(
  ames_ranger <- ranger(
    formula   = Sale_Price ~ ., 
    data      = ames_train, 
    num.trees = 500,
    mtry      = floor(length(features) / 3)
  )
)
#user  system elapsed 
#5.404   0.023   0.301 

# after parameter tuning is done, Lets repeat this model to get a better expectation of our error rate. 
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

#Furthermore, you may have noticed we set importance = 'impurity' in the above modeling, which allows us to assess variable importance. 
# Variable importance is measured by recording the decrease in MSE each time a variable is used as a node split in a tree. 
optimal_ranger$variable.importance %>% 
  tidy() %>%
  dplyr::arrange(desc(x)) %>%
  dplyr::top_n(25) %>%
  ggplot(aes(reorder(names, x), x)) +
  geom_col() +
  coord_flip() +
  ggtitle("Top 25 important variables")

## predicting
# randomForest
pred_randomForest <- predict(ames_randomForest, ames_test)
head(pred_randomForest)

# ranger
pred_ranger <- predict(ames_ranger, ames_test)
head(pred_ranger$predictions)

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

###############################################
# Create random forest
# For classification
library(randomForest)
iris.rf <- randomForest(Species ~ ., 
                        data = iris, 
                        importance = TRUE,
                        proximity = TRUE)

# Print classification model
print(iris.rf)

########################################################################################
## read in gene expression data and sample meta information 
expr.mat=data.table::fread('../external_data/2019_paper_reproduce.result/gene_by_sample_log2TPM.txt')
expr.mat=as.matrix(expr.mat,rownames=1)

sample.meta=data.table::fread('../external_data/2019_paper_reproduce.result/sample.meta_sex.label.txt')

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

## if you are using ranger for prediction
# https://stackoverflow.com/questions/55654644/predicted-probabilities-in-r-ranger-package
library("ranger")
#You need to train a "probabilistic classifier"-type ranger object:
iris.ranger = ranger(Species ~ ., data = iris, probability = TRUE)
#This object computes a matrix (n_samples, n_classes) when used in the predict.ranger function:
probabilities = predict(iris.ranger, data = iris)$predictions


##############################################
# training and testing, prediction confidence

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
dim(ames_train) #13713  6074
dim(ames_test) #3429 6074

system.time(
  rf_ranger<- ranger(
    formula=sex ~ ., 
    data = ames_train, 
    num.trees = 500,
    importance = 'impurity',
    mtry = floor(length(features) / 3),
    probability = TRUE
  )
)

rf_ranger
rf_ranger$variable.importance %>% 
  tidy() %>%
  dplyr::arrange(desc(x)) %>%
  dplyr::top_n(25) %>%
  ggplot(aes(reorder(names, x), x)) +
  geom_col() +
  coord_flip() +
  ggtitle("Top 25 important variables")

pred_ranger = predict(rf_ranger, ames_test)
head(pred_ranger$predictions) #probabilities
predict.labels=ifelse(pred_ranger$predictions[,1]>0.5,'embryoFemale','embryoMale')
table(predict.labels,ames_test$sex)
#embryoFemale embryoMale
#embryoFemale         2040        131
#embryoMale             53       1205

library(caret)

y <- as.factor(ames_test$sex) # factor of positive / negative cases
predictions <- as.factor(predict.labels) # factor of predictions
# comfusion matrix
cm <- confusionMatrix(predictions, reference =y)
cm$byClass

(precision <- posPredValue(predictions, y))
(recall <- sensitivity(predictions, y))
(F1 <- (2 * precision * recall) / (precision + recall))

#https://stackoverflow.com/questions/8499361/easy-way-of-counting-precision-recall-and-f1-score-in-r
library (ROCR);
y=ames_test$sex=='embryoFemale'# logical array of positive / negative cases
predictions <- pred_ranger$predictions[,1] # array of predictions， 1nd <=> embryoFemale

pred <- prediction(predictions, y);

# Recall-Precision curve             
RP.perf <- performance(pred, "prec", "rec");
plot (RP.perf);

# ROC curve
ROC.perf <- performance(pred, "tpr", "fpr");
plot (ROC.perf);

# ROC area under the curve
auc.tmp <- performance(pred,"auc");
(auc <- as.numeric(auc.tmp@y.values))

# F1-score performance(pred,"f") gives a vector of F1-scores 
performance(pred,"f")

