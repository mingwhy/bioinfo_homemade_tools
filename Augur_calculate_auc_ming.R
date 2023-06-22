
calculate_auc_ming = function(input=train_dat, # expression matrix for training the machine-learning model 
                         input_test=test_dat, ## expression matrix for testing the machine-learning model 
                         meta = meta,
                         meta_test = meta_test,
                         label_col = "label",
                         cell_type_col = "cell_type",
                         n_subsamples = 50,
                         subsample_prop = 0.8, #subsample_size = 20,
                         folds = 3,
                         min_cells = 20,
                         var_quantile = 0.5,
                         feature_perc = 0.5,
                         n_threads = 4,
                         show_progress = T,
                         classifier='rf', #classifier = c("rf", "lr"),
                         rf_engine='randomForest',
                         # random forest parameters
                         rf_params = list(trees = 100,
                                          mtry = 2,
                                          min_n = NULL,
                                          importance = 'accuracy'),
                         # logistic regression parameters
                         lr_params = list(mixture = 1, penalty = 'auto')
) {  

  # extract cell types and label from metadata
  expr = input
  meta %<>% droplevels()
  cell_types = meta[[cell_type_col]]
  labels = meta[[label_col]]
  labels_test = meta_test[[label_col]]
  
  subsample_size = min(as.integer(subsample_prop * table(cell_types,labels)))
  if(subsample_size<min_cells){ stop('not enough cells in the training dataset\n') }
  
  # make sure `label` is not a rowname in `input` (feature in RF)
  if ("label" %in% rownames(expr)) {
    warning("row `label` exists in input; changing ...")
    to_fix = which(rownames(expr) == "label")
    rownames(expr)[to_fix] = paste0("label", seq_along(rownames(expr)[to_fix]))
  }

  # remove missing values
  missing = is.na(expr)
  if (any(missing)) {
    stop("matrix contains ", sum(missing), "missing values")
  }

  mode = "classification"
  # check whether we are working with multiclass data
  multiclass = n_distinct(labels) > 2

  # check if showing progress or not
  if (show_progress == T) {
    apply_fun = pbmclapply
  } else {
    apply_fun = mclapply
  }

    
  # iterate over cell type clusters
  res <- apply_fun(unique(cell_types),
    mc.cores = n_threads, function(cell_type) {
      # skip this cell type if there aren't enough cells
      y = labels[cell_types == cell_type]
      if (mode == 'classification') {
        if (min(table(y)) < min_cells) {
          warning("skipping cell type ", cell_type,
                  ": minimum number of cells (", min(table(y)),
                  ") is less than ", min_cells)
          return(list())
        }
      }
      # subset the entire expression matrix to this cell type
      X = expr[, cell_types == cell_type]

      # select features by variance
      min_features_for_selection = 1000
      if (nrow(X) >= min_features_for_selection) {
        X %<>% select_variance(var_quantile, filter_negative_residuals = F)
      }

      # set up subsamples and results bin
      tmp_results = data.frame()
      tmp_results_on_new = data.frame()

      n_iter = ifelse(n_subsamples < 1, 1, n_subsamples)
      
      for (subsample_idx in seq_len(n_iter)) {
        # seed RNG for reproducibility
        set.seed(subsample_idx)

        # optionally, skip the subsampling process
        if (n_subsamples < 1) {
          # randomly select features
          if (nrow(X) >= min_features_for_selection &
              feature_perc < 1) {
            X0 = select_random(X, feature_perc)
          } else {
            X0 = X
          }
          X0 %<>%
            t() %>%
            as.matrix() %>%
            as.data.frame() %>%
            # fix any non-unique columns
            repair_names() %>%
            # add labels
            mutate(label = y)
        } else {
          if (mode == 'regression') {
            subsample_idxs = data.frame(label = y,
                                        position = seq_along(y)) %>%
              do(sample_n(., subsample_size)) %>%
              pull(position)
          } else {
            subsample_idxs = data.frame(label = y,
                                        position = seq_along(y)) %>%
              group_by(label) %>%
              do(sample_n(., subsample_size)) %>%
              pull(position)
          }
          y0 = y[subsample_idxs]
          # randomly select features
          if (nrow(X) >= min_features_for_selection &
              feature_perc < 1) {
            X0 = select_random(X, feature_perc)
          } else {
            X0 = X
          }
          # coerce to data frame
          X0 %<>%
            magrittr::extract(, subsample_idxs) %>%
            t() %>%
            #magrittr::extract(, colVars(.) > 0) %>%
            as.matrix() %>%
            as.data.frame() %>%
            # fix any non-unique columns
            repair_names() %>%
            # convert labels back to a factor
            mutate(label = y0) #sample by (feature+label) columns
        }
        X0_test = test_dat %>% t() %>%
            #magrittr::extract(, colVars(.) > 0) %>%
            as.matrix() %>%
            as.data.frame() %>%
            mutate(label = meta_test$label) #sample by (feature+label) columns
        X0_test$label=gsub('control','normal',X0_test$label)
        X0_test$label=gsub('rapa','starved',X0_test$label)        
        X0_test$label=factor(X0_test$label,levels=c('normal','starved'))

        # set up model
        if (classifier == "rf") {
          importance = T
          clf = rand_forest(trees = !!rf_params$trees,
                            mtry = !!rf_params$mtry,
                            min_n = !!rf_params$min_n,
                            mode = mode) %>%
            set_engine(rf_engine, seed = 1, importance = T, localImp = T)
        } 

        # fit models in cross-validation
        if (mode == "classification") {
          cv = vfold_cv(X0, v = folds, strata = 'label')
        } else {
          cv = vfold_cv(X0, v = folds)
        }
        
        #cv$splits
        #dim(X0);dim(cv$splits[[1]]$data)
        #length(cv$splits[[1]]$in_id)        
        #dim(assessment(folded$splits[[1]]))

        withCallingHandlers({
          folded = cv %>%
            mutate(
              recipes = splits %>%
                map(~ prepper(., recipe = recipe(.$data, label ~ .))),
              test_data = splits %>% map(analysis),
              fits = map2(
                recipes,
                test_data, #use trained data to fit a RF model
                ~ fit(
                  clf,
                  label ~ .,
                  data = bake(object = .x, new_data = .y)
                )
              )
            )
        }, warning = function(w) {
          if (grepl("dangerous ground", conditionMessage(w)))
            invokeRestart("muffleWarning")
        })   
        #folded$recipes   
        #folded$fits #already fitted model

        # predict on the left-out data
        retrieve_class_preds = function(split, recipe, model) {
          test = bake(recipe, assessment(split))
          tbl = tibble(
            true = test$label,
            pred = predict(model, test)$.pred_class,
            prob = predict(model, test, type = 'prob')) %>%
            # convert prob from nested df to columns
            cbind(.$prob) %>%
            select(-prob)
          return(tbl)
        }

        retrieve_class_preds_on_new_data = function(split, recipe, model) {
          #model=folded$fits[[1]]
          #split=folded$splits[[1]]
          #recipe=folded$recipes[[1]]
          #test = bake(recipe, assessment(split))
          #ncell=nrow(assessment(folded$splits[[1]]))
          test = bake(recipe, X0_test) #test on all rapa cells
          tbl = tibble(
            true = test$label,
            pred = predict(model, test)$.pred_class,
            prob = predict(model, test, type = 'prob')) %>%
            # convert prob from nested df to columns
            cbind(.$prob) %>%
            select(-prob)
          return(tbl)
        }

        predictions = folded %>%
          mutate(pred = list(splits,recipes,fits))
        #predictions$pred[[1]] #three splits
        #predictions$pred[[2]] #three recipes
        #predictions$pred[[2]] #three fitted RF model
        predictions = predictions %>% #do actual prediction
          mutate(pred = pmap(pred, retrieve_class_preds))        
        #dim(predictions$test_data[[1]]) #train data point
        #dim(predictions$pred[[1]]) #test data point

        predictions_on_new = folded %>%
          mutate(pred = list(splits,recipes,fits))
        predictions_on_new = predictions_on_new %>%
          mutate(pred = pmap(pred, retrieve_class_preds_on_new_data))


        # evalulate the predictions
        if (mode == 'regression') {
          multi_metric = metric_set(ccc, huber_loss_pseudo, huber_loss,
            mae, mape, mase, rpd, rpiq, rsq_trad, rsq, smape, rmse)
        } else {
          multi_metric = metric_set(accuracy, precision, recall, sens,
                                    spec, npv, ppv, roc_auc)
        }

        prob_select = 3
        if (mode == "classification") {
          estimator = ifelse(multiclass, "macro", "binary")
          if (multiclass)
            prob_select = seq(3, 3 + n_distinct(labels) - 1)
          metric_fun = function(x)
            multi_metric(x,
                         truth = true,
                         estimate = pred,
                         prob_select,
                         estimator = estimator)
        } else {
          metric_fun = function(x)
            multi_metric(x,
                         truth = true,
                         estimate = pred,
                         prob_select)
        }
        
        eval = predictions %>%
          mutate(
            metrics = pred %>%
              map(metric_fun)
          ) %>%
          extract2('metrics')

        eval_on_new = predictions_on_new %>%
          mutate(
            metrics = pred %>%
              map(metric_fun)
          ) %>%
          extract2('metrics')

        # clean up the results
        names(eval)=as.character(1:length(eval))
        result = eval %>%
          map2_df(., names(.), ~ mutate(.x, fold = .y)) %>%
          set_colnames(gsub("\\.", "", colnames(.))) %>%
          mutate(cell_type = cell_type,
                 subsample_idx = subsample_idx)
        
        names(eval_on_new)=as.character(1:length(eval_on_new))
        result_on_new = eval_on_new %>%
          map2_df(., names(.), ~ mutate(.x, fold = .y)) %>%
          set_colnames(gsub("\\.", "", colnames(.))) %>%
          mutate(cell_type = cell_type,
                 subsample_idx = subsample_idx)

        # rearrange columns
        result %<>%
          dplyr::select(cell_type, subsample_idx, fold, metric, estimator,
                        estimate)
        result_on_new %<>%
          dplyr::select(cell_type, subsample_idx, fold, metric, estimator,
                        estimate)
        # add to results
        tmp_results %<>% bind_rows(result)
        tmp_results_on_new %<>% bind_rows(result_on_new)
      }
      res=list('results' = tmp_results, 'results_on_new' = tmp_results_on_new)
      return(res)
    }
  )
  #names(res[[1]])
  
  # ignore warnings from yardstick
  if (any(map_lgl(res, ~ "warning" %in% class(.)))) {
    res = res$value
  }

  # make sure at least one cell type worked
  if (all(lengths(res) == 0))
    stop("no cell type had at least ", min_cells, " cells in all conditions")

  # summarise AUCs (or CCCs) per cell type
  if (mode == "classification") {
    AUCs = res %>%
      map("results") %>%
      bind_rows() %>%
      filter(metric == "roc_auc") %>%
      group_by(cell_type, subsample_idx) %>%
      summarise(estimate = mean(estimate)) %>%
      ungroup() %>%
      group_by(cell_type) %>%
      summarise(auc = mean(estimate)) %>%
      ungroup() %>%
      arrange(desc(auc))
    
    AUCs_on_new = res %>%
      map("results_on_new") %>%
      bind_rows() %>%
      filter(metric == "roc_auc") %>%
      group_by(cell_type, subsample_idx) %>%
      summarise(estimate = mean(estimate)) %>%
      ungroup() %>%
      group_by(cell_type) %>%
      summarise(auc = mean(estimate)) %>%
      ungroup() %>%
      arrange(desc(auc))
  } 
  
  results = res %>%
    map("results") %>%
    bind_rows()
  
  results_on_new = res %>%
    map("results_on_new") %>%
    bind_rows()
  
  # create an object for return
  params = list(
    n_subsamples = n_subsamples,
    subsample_size = subsample_size,
    folds = folds,
    min_cells = min_cells,
    #var_quantile = var_quantile,
    #feature_perc = feature_perc,
    n_threads = n_threads,
    classifier = classifier
  )
  if (classifier == "rf")
    params$rf_params = rf_params
  if (classifier == "lr")
    params$lr_params = lr_params
  obj = list(
    #X = expr,
    #y = labels,
    #cell_types = cell_types,
    parameters = params,
    results = results,
    results_on_new = results_on_new
    #feature_importance = feature_importances
  )
  if (mode == "classification") {
    obj$AUC = AUCs
    obj$AUC_on_new = AUCs_on_new
  } else if (mode == "regression") {
    obj$CCC = CCCs
  }

  return(obj)
}


