#http://pnp.mathematik.uni-stuttgart.de/isa/steinwart/software/R/demo.html
#install.packages("liquidSVM", repos="http://pnp.mathematik.uni-stuttgart.de/isa/steinwart/software/R")

library(liquidSVM)
#Least Squares Regression
reg <- liquidData('reg-1d')
model <- svm(Y~., reg$train)
result <- test(model, reg$test)
errors(result)
plot(reg$train$X1, reg$train$Y, ylim=c(-.2,.8), ylab='Y', xlab='X1', axes=T, pch='.', cex=2.5)
curve(predict(model, x), add=T, col='red', lwd=2)
model <- svm(Y~., reg)
errors(model$last_result)[1]

#Multi-Class Classification
#https://github.com/liquidSVM/liquidSVM
banana <- liquidData('banana-mc')
banana
model <- mcSVM(Y~., banana$train,display=1, threads=2)
result <- test(model, banana$test)
errors(result)

mycol <- c('red', 'blue', 'cyan', 'green')
plot(banana$train$X1, banana$train$X2, pch=20, col=mycol[banana$train$Y], ylab='', xlab='', axes=F, lwd=0.25)
x <- seq(-1,1,.01)
z <- matrix(predict(model,expand.grid(x,x)),length(x))
contour(x,x,z, add=T, levels=1:4, col=1, lwd=3)
errors(test(model, banana$test))

#https://rdrr.io/cran/liquidSVM/man/liquidSVM-package.html
#https://www.rdocumentation.org/packages/liquidSVM/versions/1.2.4/topics/mcSVM
set.seed(123)
## Multiclass classification
modelIris <- liquidSVM::mcSVM(Species ~ ., iris,predict.prob = T)
modelIris

y <- predict(modelIris, iris, probability=TRUE)
names(y)
head(y)

pred=apply(y,1,function(x) names(which.max(x)))
table(pred,iris$Species)

modelIris$train_data
