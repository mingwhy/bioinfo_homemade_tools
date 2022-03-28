#https://docs.h2o.ai/h2o/latest-stable/h2o-docs/downloading.html
if(F){
# Downloading & Installing H2O
1. Download and Run from the Command Line
 you may need to install java runtime before hand

1) Click the Download H2O button on the http://h2o-release.s3.amazonaws.com/h2o/latest_stable.html page. This downloads a zip file that contains everything you need to get started.
2) From your terminal, unzip and start H2O as in the example below.
# cd ~/Downloads
# unzip h2o-3.30.0.6.zip
# cd h2o-3.30.0.6
# java -jar h2o.jar
3) Point your browser to http://localhost:54321 to open up the H2O Flow web GUI.

2. Install in R
Perform the following steps in R to install H2O. Copy and paste these commands one line at a time.

1) The following two commands remove any previously installed H2O packages for R.
if ("package:h2o" %in% search()) { detach("package:h2o", unload=TRUE) }
if ("h2o" %in% rownames(installed.packages())) { remove.packages("h2o") }

2) Next, download packages that H2O depends on.
pkgs <- c("RCurl","jsonlite")
for (pkg in pkgs) {
  if (! (pkg %in% rownames(installed.packages()))) { install.packages(pkg) }
}

3) Download and install the H2O package for R.
install.packages("h2o", type="source", repos=(c("http://h2o-release.s3.amazonaws.com/h2o/latest_stable_R")))

4) Optionally initialize H2O and run a demo to see H2O at work.
library(h2o)
localH2O = h2o.init()
demo(h2o.kmeans)

}

## test: https://docs.h2o.ai/h2o/latest-stable/h2o-docs/data-science/svm.html

library(h2o)
h2o.init()

# Import the splice dataset into H2O:
splice <- h2o.importFile("https://s3.amazonaws.com/h2o-public-test-data/smalldata/splice/splice.svm")

# Build and train the model:
svm_model <- h2o.psvm(gamma = 0.01,
                      rank_ratio = 0.1,
                      y = "C1",
                      training_frame = splice,
                      disable_training_metrics = FALSE)

# Eval performance:
perf <- h2o.performance(svm_model)


## used case: https://www.r-bloggers.com/2018/07/dalex-and-h2o-machine-learning-model-interpretability-and-feature-explanation/
# together with R package `DALEX`
# load required packages
library(rsample)
library(dplyr)
library(purrr)
library(ggplot2)
library(h2o)
#install.packages("DALEX")
library(DALEX) 

# initialize h2o session
h2o.no_progress()
h2o.init()

# classification data
#https://rsample.tidymodels.org/reference/attrition.html
#https://uc-r.github.io/dalex
#install.packages('modeldata')
library(modeldata)
data(attrition)
head(attrition)
dim(attrition) #1470 x 31
df <- attrition %>% 
  mutate_if(is.ordered, factor, ordered = FALSE) %>%
  mutate(Attrition = recode(Attrition, "Yes" = "1", "No" = "0") %>% factor(levels = c("1", "0")))
head(df)
dim(df) #1470 x 31

# convert to h2o object
df.h2o <- as.h2o(df)

# create train, validation, and test splits
set.seed(123)
splits <- h2o.splitFrame(df.h2o, ratios = c(.7, .15), destination_frames = c("train","valid","test"))
names(splits) <- c("train","valid","test")
dim(splits$train)#1026   31
dim(splits$valid) # 233  31
dim(splits$test)#211  31

# variable names for resonse & features
y <- "Attrition"
x <- setdiff(names(df), y) 
x

# elastic net model 
glm <- h2o.glm(
  x = x, 
  y = y, 
  training_frame = splits$train,
  validation_frame = splits$valid,
  family = "binomial",
  seed = 123
)

# random forest model
rf <- h2o.randomForest(
  x = x, 
  y = y,
  training_frame = splits$train,
  validation_frame = splits$valid,
  ntrees = 1000,
  stopping_metric = "AUC",    
  stopping_rounds = 10,         
  stopping_tolerance = 0.005,
  seed = 123
)

# gradient boosting machine model
gbm <-  h2o.gbm(
  x = x, 
  y = y,
  training_frame = splits$train,
  validation_frame = splits$valid,
  ntrees = 1000,
  stopping_metric = "AUC",    
  stopping_rounds = 10,         
  stopping_tolerance = 0.005,
  seed = 123
)
h2o.varimp(gbm) #https://docs.h2o.ai/h2o/latest-stable/h2o-docs/variable-importance.html

# model performance
h2o.auc(glm, valid = TRUE)
#0.7870935
h2o.auc(rf, valid = TRUE)
# 0.7681021
h2o.auc(gbm, valid = TRUE)
#0.7468242

## DALEX procedures
# convert feature data to non-h2o objects
x_valid <- as.data.frame(splits$valid)[, x]

# make response variable numeric binary vector
y_valid <- as.vector(as.numeric(as.character(splits$valid$Attrition)))
head(y_valid)

# create custom predict function
pred <- function(model, newdata)  {
  results <- as.data.frame(h2o.predict(model, as.h2o(newdata)))
  return(results[[3L]])
}

pred(rf, x_valid) %>% head()


# elastic net explainer
explainer_glm <- explain(
  model = glm,
  data = x_valid,
  y = y_valid,
  predict_function = pred,
  label = "h2o glm"
)

# random forest explainer
explainer_rf <- explain(
  model = rf,
  data = x_valid,
  y = y_valid,
  predict_function = pred,
  label = "h2o rf"
)

# GBM explainer
explainer_gbm <- explain(
  model = gbm,
  data = x_valid,
  y = y_valid,
  predict_function = pred,
  label = "h2o gbm"
)

# example of explainer object
class(explainer_glm)

# residual diagnostics
# compute predictions & residuals
resids_glm <- model_performance(explainer_glm)
resids_rf  <- model_performance(explainer_rf)
resids_gbm <- model_performance(explainer_gbm)

# assess quantiles for residuals
resids_glm
resids_rf
resids_gbm

# create comparison plot of residuals for each model
p1 <- plot(resids_glm, resids_rf, resids_gbm)
p2 <- plot(resids_glm, resids_rf, resids_gbm, geom = "boxplot")

gridExtra::grid.arrange(p1, p2, nrow = 1)

#  Variable importance
## DALEX uses a model agnostic variable importance measure computed via permutation.
# compute permutation-based variable importance
# open http://localhost:54321/ in the browser before hand
dim(explainer_glm$data) #233  30
dim(explainer_rf$data) #233 30
vip_glm <- variable_importance(explainer_glm,  loss_function = loss_root_mean_square) 
#vip_glm <- variable_importance(explainer_glm, n_sample = 233, loss_function = loss_root_mean_square) 
#vip_glm <- variable_importance(explainer_glm, n_sample = -1, loss_function = loss_root_mean_square) 
#vip_rf  <- variable_importance(explainer_rf, n_sample = -1, loss_function = loss_root_mean_square)
#vip_gbm <- variable_importance(explainer_gbm, n_sample = -1, loss_function = loss_root_mean_square)
vip_rf  <- variable_importance(explainer_rf, loss_function = loss_root_mean_square)
vip_gbm <- variable_importance(explainer_gbm, loss_function = loss_root_mean_square)


plot(vip_glm, vip_rf, vip_gbm, max_vars = 10)


