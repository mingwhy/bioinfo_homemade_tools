##########################################
# svm with caret (https://rpubs.com/uky994/593668)
library(tidyverse)
library(caret)

# Load the data
data("PimaIndiansDiabetes2", package = "mlbench")
pima.data <- na.omit(PimaIndiansDiabetes2)
# Inspect the data
sample_n(pima.data, 3)

# Set up Repeated k-fold Cross Validation
train_control <- trainControl(method="repeatedcv", number=5, repeats=1)
# Fit the model (SVM linear classifier)
svm1 <- train(diabetes ~., data = pima.data, method = "svmLinear", trControl = train_control,  preProcess = c("center","scale"))
#View the model
svm1


##########################################
#Example for svm feature selection in R
#https://stackoverflow.com/questions/17529537/example-for-svm-feature-selection-in-r
#http://topepo.github.io/caret/recursive-feature-elimination.html
#http://topepo.github.io/caret/recursive-feature-elimination.html#rfe

data(BloodBrain, package="caret")
logBBB
x <- scale(bbbDescr[,-nearZeroVar(bbbDescr)])
x <- x[, -findCorrelation(cor(x), .8)]
x <- as.data.frame(x)
dim(x); #208 samples, 71 features
length(logBBB) #208 responses
svmProfile <- rfe(x, logBBB,
                  sizes = c(2, 5, 10, 20),
                  rfeControl = rfeControl(functions = caretFuncs,
                                          method='cv',
                                          number = 5,repeats=10),
                  ## pass options to train()
                  method = "svmRadial")
svmProfile
names(svmProfile)
svmProfile$results
head(svmProfile$variables)
table(svmProfile$variables$Variables)
head(svmProfile$variables)

tmp=svmProfile$variables
tmp[tmp$Variables==2,]

head(sort(table(tmp$var),decreasing=T),10) #most frequently selected variables
svmProfile
#The top 5 variables (out of 71):
#fpsa3, tcsa, prx, tcpa, most_positive_charge

#At each iteration of feature selection, the Si top ranked predictors are retained, the model is refit and performance is assessed. The value of Si with the best performance is determined and the top Si predictors are used to fit the final model. 
svmProfile$obsLevels
