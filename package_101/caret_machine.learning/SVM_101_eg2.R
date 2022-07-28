#http://topepo.github.io/caret/train-models-by-tag.html#Support_Vector_Machines
#https://uc-r.github.io/svm
# set pseudorandom number generator
set.seed(10)

# Attach Packages
library(tidyverse)    # data manipulation and visualization
library(kernlab)      # SVM methodology
library(e1071)        # SVM methodology
#library(ISLR)         # contains example data set "Khan"
library(RColorBrewer) # customized coloring of plots

# Construct sample data set - completely separated
x <- matrix(rnorm(20*2), ncol = 2)
y <- c(rep(-1,10), rep(1,10))
x[y==1,] <- x[y==1,] + 3/2
dat <- data.frame(x=x, y=as.factor(y))

# Plot data
ggplot(data = dat, aes(x = x.2, y = x.1, color = y, shape = y)) + 
  geom_point(size = 2) +
  scale_color_manual(values=c("#000000", "#FF0000")) +
  theme(legend.position = "none")

# Fit Support Vector Machine model to data set
svmfit <- svm(y~., data = dat, kernel = "linear", scale = FALSE)
# Plot Results
plot(svmfit, dat)

# fit model and produce plot
kernfit <- ksvm(x, y, type = "C-svc", kernel = 'vanilladot')
plot(kernfit, data = x)


