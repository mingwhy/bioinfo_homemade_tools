# supporting vector machine
# https://www.geeksforgeeks.org/classifying-data-using-support-vector-machinessvms-in-r/
# linear
set.seed(100)
x <- matrix(rnorm(40),20,2)
y <- rep(c(-1,1),c(10,10))
x[y == 1,] = x[y == 1,] + 1
plot(x, col = y + 3, pch = 19)

library(e1071)
data = data.frame(x, y = as.factor(y))
head(data)

data.svm = svm(y ~ ., data = data, kernel = "linear", cost = 10, scale = FALSE)
print(data.svm)

plot(data.svm, data)

#test on real data
rownames(sub.df)
library(e1071)
data = data.frame(t(sub.df), y = substr(names(sub.df),0,1))
head(data)
data=data[1:22,] #remove 2 NA
data$y=as.numeric(factor(data$y)) #1=female, 2=male
data.svm = svm(y ~ ., data = data, kernel = "linear", cost = 10, scale = FALSE)
print(data.svm)

plot(data.svm, data)

library(caTools)
dataset=data
set.seed(123)
split = sample.split(dataset$y, SplitRatio = 0.7)

training_set = subset(dataset, split == TRUE)
test_set = subset(dataset, split == FALSE)

# Feature Scaling for train and test separately
training_set[-5] = scale(training_set[-5])
test_set[-5] = scale(test_set[-5])

classifier = svm(formula = y ~ .,
                 data = training_set,
                 type = 'C-classification',
                 kernel = 'linear')
# Predicting the Test set results
y_pred = predict(classifier, newdata = test_set[-5])
table(test_set$y,y_pred)

## extend to single cell data
## read in embryo data
x=readRDS('../2019_paper_reproduce.result/prepared_embryo.data.rds')
test=as.matrix(x$umi.mat)
test.meta=x$sample.meta
class(test)

gene.names=rownames(test);
pick.genes=gene.names[grep('Sxl|msl|roX|mle|mof',gene.names)]
pick.genes=gene.names[grep('Sxl|msl_2|roX',gene.names)]
test=test[pick.genes,]
dim(test)
test_set=t(test)
test_set = scale(test_set)
colnames(test_set)
colnames(test_set)=c('msl.2','roX1','roX2','Sxl')

# Predicting the Test set results
y_pred = predict(classifier, newdata = test_set)
table(y_pred)
#y_pred
#1   2 
#803 494 
table(y_pred,test.meta$sex)
