#https://uc-r.github.io/random_forests

library(rsample)      # data splitting 
library(randomForest) # basic implementation
library(ranger)       # a faster implementation of randomForest
library(caret)        # an aggregator package for performing many machine learning models
library(h2o)          # an extremely fast java-based platform

# Create training (70%) and test (30%) sets for the AmesHousing::make_ames() data.
# Use set.seed for reproducibility

set.seed(123)
inp=AmesHousing::make_ames()
ames_split <- initial_split(AmesHousing::make_ames(), prop = .7)
ames_train <- training(ames_split)
ames_test  <- testing(ames_split)
dim(ames_test) 
dim(ames_train)

# for reproduciblity
set.seed(123)
hist(ames_train$Sale_Price) #response var, a numeric variable

# default RF model
m1 <- randomForest(
  formula = Sale_Price ~ .,
  data    = ames_train
)
m1
plot(m1)

