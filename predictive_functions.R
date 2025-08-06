

library(tidyverse)
library(caret)
library(pROC)
library(smotefamily)
library(UBL)

library(wesanderson)

mycols <- colorRampPalette(wes_palette("Royal1"))(5)

options(ggplot2.discrete.fill = mycols)
options(ggplot2.discrete.colour = mycols)



## GGPLOT THEMES

#options(ggplot2.discrete.fill = c("#1E90FF", "#00AA5A", "#F75A3F", "#8E7BFF","#00D1EE", "#00E6BB", "#F9F871", "#F45680", "#A5ABBD", "#B60E50"))
#options(ggplot2.discrete.colour = c("#1E90FF","#00AA5A", "#F75A3F",  "#8E7BFF","#00D1EE", "#00E6BB", "#F9F871", "#F45680", "#A5ABBD", "#B60E50"))

#options(ggplot2.discrete.fill = c("#A1C6EA","#FD8B2F", "#00AA5A", "#8E7BFF","#00D1EE", "#00E6BB", "#F9F871", "#F45680", "#A5ABBD", "#B60E50"))
#options(ggplot2.discrete.colour = c("#A1C6EA","#FD8B2F","#00AA5A",   "#8E7BFF","#00D1EE", "#00E6BB", "#F9F871", "#F45680", "#A5ABBD", "#B60E50"))


randomforest_params = list(ntree = 500, 
                           mtry = 1, 
                           nodesize = 1, 
                           balance_weights = FALSE)

xgboost_params = list(max_depth = 2, 
                      learning_rate = 1, 
                      nrounds = 2,
                      nthread = 2, 
                      objective = "binary:logistic", 
                      balance_weights=FALSE)

catboost_params <- list(
  iterations = 50,              # Number of boosting rounds
  learning_rate = 0.02,           # Step size shrinkage
  depth = 2,                     # Depth of the trees
  loss_function = "Logloss",     # Binary classification
  eval_metric = "AUC",           # Evaluation metric
  random_seed = 123,            
  use_best_model = TRUE,         # Stop early if no improvemelnt
  od_type = "Iter",              # Overfitting detector type
  od_wait = 20,                  # Rounds to wait before stopping
  verbose = FALSE,               
  thread_count = 1,              
  balance_weights = TRUE,  
  bootstrap_type = "Bayesian",
  l2_leaf_reg = 1,
  subsample = 1,  # only if bootstrap type ="Bernouilli"
  grow_policy = "SymmetricTree"
)

smote_params=list(K=5, dup_size="balance")

makeKmeans <- function(datasc, levs, varnames, SEED=123, folds=c()){
  library(stats)
  library(caret)
  set.seed(SEED)
  train_df <- datasc %>% dplyr::select(-class, -sample)  %>% dplyr::select(all_of(varnames))
  mod_kmeans <- kmeans(train_df, centers=length(levs), iter.max = 100, nstart=100)
  predict_kmeans <-levels(datasc$class)[mod_kmeans$cluster] %>% factor(levels=levs)
  confmat_kmeans <- confusionMatrix(predict_kmeans, datasc$class, positive = levs[2])
  return(list(confmat_no_l1o=confmat_kmeans, mod=mod_kmeans, predicted=predict_kmeans))
}

makeKmeans_l1o <- function(datasc, levs, varnames, SEED=123, folds=c(), 
                           do_smote=FALSE, 
                           smote_params=list(K=5, dup_size="balance")){
  library(stats)
  library(caret)
  set.seed(SEED)
  train_df_all <- datasc %>% dplyr::select(-class, -sample)  %>% dplyr::select(all_of(varnames))
  
  if(length(folds)==0){
    folds <- 1:nrow(datasc)
  }
  
  preds <- c()
  all_dists <- list()
  for(i in folds){
    train_df <- train_df_all[-i, ]
    test_df <- train_df_all[i, ]
    # Separar clases
    train_labels <- datasc$class[-i]
    test_labels <- datasc$class[i]
    
    if(do_smote){
      form <- as.formula(paste0("class ~ ", paste(varnames, sep="+", collapse= "+")))
      smote_df <- datasc[-i, ] %>% select(class, all_of(varnames))
      smoteData <- SmoteClassif(form, smote_df, 
                                C.perc = smote_params$dup_size, 
                                k = smote_params$K, repl = FALSE,
                                dist = "Euclidean", p = 2)
      train_df <- smoteData %>% select(-class)
      train_labels <- factor(smoteData$class, levels=levs)
    }else{
      smoteData = NULL
    }
    
    mod_kmeans <- kmeans(train_df, centers=length(levs), iter.max = 100, nstart=100)
    
    assign_class <- table(mod_kmeans$cluster, train_labels)
    cents <- mod_kmeans$centers
    rownames(cents) <- paste0("C_", rownames(cents))
    train_dists <- dist(rbind(cents, train_df %>% as.matrix )) %>% 
      as.matrix %>% 
      as.data.frame %>% 
      rownames_to_column("sample") %>% 
      dplyr::select("sample", all_of(rownames(cents))) %>% 
      filter(! sample %in% rownames(cents)) %>% 
      dplyr::mutate(class = train_labels) %>% 
      group_by(class) %>% 
      dplyr::summarise(across(all_of(rownames(cents)), .fns = mean)) %>% 
      rowwise() %>%
      dplyr::mutate(min_dist = min(c_across(all_of(rownames(cents))))) %>%
      ungroup() %>% 
      dplyr::arrange(min_dist)

    free_groups <- rownames(cents)
    for(ir in 1:nrow(train_dists)){
      newgr <- free_groups[which.min(train_dists[ir, free_groups] %>% as.vector %>% unlist)]
      train_dists$cluster[ir] <- newgr
      free_groups <- free_groups[free_groups != newgr]
    }
  
    rownames(cents) <- train_dists$class[match(rownames(cents), train_dists$cluster)]
    
    dists <- dist(rbind(cents, test_df)) %>% as.matrix
    assign_class <- rownames(cents)
    
    if(length(i)> 1){
      pred_i <- apply(dists[rownames(test_df),  assign_class], MAR=1, \(x) assign_class[which.min(x)])
      all_dists[[i]] <- dists[rownames(test_df),  assign_class]
    }else{
      pred_i <- assign_class[which.min(dists[rownames(test_df),  assign_class])]
      all_dists[[i]] <- dists[rownames(test_df),  assign_class]
    }
    preds  <- c(preds, pred_i )
  }
  preds <- factor(preds, levels=levels(datasc$class))
  confmat_kmeans <- confusionMatrix(preds, datasc$class, positive = levs[2])
  if(length(levs) == 2){
    all_dists <- all_dists %>% bind_rows
    probs <- 1 - (all_dists %>% pull(!!sym(levs[2])))/rowSums(all_dists)
    roc1 <- roc(response=as.numeric(datasc$class)-1, predictor=probs)
  }else{
    probs <- all_dists %>% bind_rows %>% 
      dplyr::mutate_all(\(x) 1 - x/rowSums(.))
    roc1 <- multiclass.roc(response=datasc$class, predictor=as.matrix(probs))
  }
  
  mod_kmeans_all <- kmeans(train_df_all, centers=length(levs), iter.max = 100, nstart=100)
  predict_kmeans_nol1o <-levels(datasc$class)[mod_kmeans_all$cluster] %>% factor(levels=levs)
  confmat_kmeans_nol1o <- confusionMatrix(predict_kmeans_nol1o, datasc$class, positive = levs[2])
  
  return(list(confmat = confmat_kmeans, 
              confmat_no_l1o=confmat_kmeans_nol1o, 
              preds=preds,
              pred_probs=probs,
              preds_no_l1o=predict_kmeans_nol1o, 
              mod=mod_kmeans_all,
              roc_obj_no_l1o=NULL,
              roc_auc_no_l1o=NULL,
              roc_obj=roc1,
              roc_auc=as.numeric(roc1$auc))
         )
}

makeKnn_l1o <- function(datasc, levs, varnames, 
                        different_ks=c(1, 3, 5, 7, 9, 11, 13), 
                        folds=c(), 
                        do_smote=FALSE, 
                        smote_params=list(K=5, dup_size="balance")){
  library(class)
  #library(gmodels)
  
  results <- list()
  datasc <- datasc %>% dplyr::mutate(class=factor(class))
  train_df_all <- datasc %>% 
    dplyr::select(-class, -sample) %>% 
    dplyr::select(all_of(varnames))
  
  if(length(folds)==0){
    folds <- 1:nrow(datasc)
  }
  for(k in  different_ks){
    kname = paste("K=", as.character(k), sep="", collapse="")
    results[[kname]] <- list()
    preds <- c()
    pred_probs <- list()
    for(i in folds){
      train_df <- train_df_all[-i, ]
      test_df <- train_df_all[i, ]
      
      # Separar clases
      train_labels <- datasc$class[-i]
      test_labels <- datasc$class[i]
      
      if(do_smote){
        form <- as.formula(paste0("class ~ ", paste(varnames, sep="+", collapse= "+")))
        smote_df <- datasc[-i, ] %>% select(class, all_of(varnames))
        smoteData <- SmoteClassif(form, smote_df, 
                                  C.perc = smote_params$dup_size, 
                                  k = smote_params$K, repl = FALSE,
                                  dist = "Euclidean", p = 2)
        train_df <- smoteData %>% select(-class)
        train_labels <- factor(smoteData$class, levels=levs)
      }else{
        smoteData = NULL
      }
      kires <- knn(train_df, test_df, train_labels, k = k, prob = T)
      preds  <- c(preds, kires)
      pred_probs[[i]] <- kires
      
    }
    if(length(levs == 2)){
      prob_vec <- map_vec(pred_probs, \(x){ifelse(x==levs[1], 1-attr(x, "prob"), attr(x, "prob"))})
      roc1 <- roc(response=as.numeric(datasc$class)-1, predictor=prob_vec)
      results[[kname]][["roc_obj"]] <- roc1
      results[[kname]][["roc_auc"]] <- as.numeric(roc1$auc)
    }else{
      prob_vec <- map_vec(pred_probs, \(x){ifelse(x==levs[1], 1-attr(x, "prob"), attr(x, "prob"))})
      roc1 <- multiclass.roc(response=datasc$class, predictor=prob_vec)
      results[[kname]][["roc_obj"]] <- roc1
      results[[kname]][["roc_auc"]] <- as.numeric(roc1$auc)
    }
    
    results[[kname]][["preds"]] <- levs[preds] %>% factor(levels=levs)
    results[[kname]][["confmat"]] <- confusionMatrix(results[[kname]][["preds"]], 
                                                     factor(datasc$class), 
                                                     positive = levs[2])
    results[[kname]][["pred_probs"]] <- prob_vec
     
    
    
    preds2 <- knn(train_df_all, train_df_all, datasc$class, k = k, prob = T)
    results[[kname]][["preds_no_l1o"]] <- levs[preds2] %>% factor(levels=levs)
    results[[kname]][["confmat_no_l1o"]] <- confusionMatrix(results[[kname]][["preds_no_l1o"]], 
                                                            factor(datasc$class), 
                                                            positive = levs[2])
  }
  
  return(results)
}

makeKnn <- function(datasc, levs, nvars, different_ks=c(1, 3, 5, 7, 9, 11, 13)){
  library(class)
  #library(gmodels)
  
  test_pred <- list()
  conf_matrices_knn <- list()
  train_df <- datasc %>% dplyr::select(-class, -sample) %>% dplyr::select(all_of(varnames))
  train_labels <- datasc$class %>% factor(levels=levs)
  
  for(k in  different_ks){
    kname = paste("K=", as.character(k), sep="", collapse="")
    test_pred[[kname]] <- knn(train_df, train_df, train_labels, k = k, prob = T)
    conf_matrices_knn[[kname]] <- confusionMatrix(test_pred[[kname]], 
                                                  train_labels, 
                                                  positive = levs[2])
  }
  
  return(list(confmats=conf_matrices_knn, mods=test_pred))
}

makeNaiveBayes_l1o <- function(datasc, levs, varnames, 
                               SEED=123, folds=c(), 
                               do_smote=FALSE, 
                               smote_params=list(K=5, dup_size="balance")){
  library(e1071)
  set.seed(SEED)
  
  predict_bayes1 <- factor()
  predict_bayes1_probs <- list()
  df <- datasc %>% dplyr::select(-class, -sample) %>% dplyr::select(all_of(varnames))
  if(length(folds)==0){
    folds <- 1:nrow(datasc)
  }
  for(i in folds){
    # Separar datos
    train_df <- df[-i, ]
    test_df <- df[i, ]
    
    # Separar clases
    train_labels <- datasc$class[-i]
    test_labels <- datasc$class[i]
    
    if(do_smote){
      form <- as.formula(paste0("class ~ ", paste(varnames, sep="+", collapse= "+")))
      smote_df <- datasc[-i, ] %>% select(class, all_of(varnames))
      smoteData <- SmoteClassif(form, smote_df, 
                                C.perc = smote_params$dup_size, 
                                k = smote_params$K, repl = FALSE,
                                dist = "Euclidean", p = 2)
      train_df <- smoteData %>% select(-class)
      train_labels <- factor(smoteData$class, levels=levs)
    }else{
      smoteData = NULL
    }
    
    mod_bayes2 <- naiveBayes(train_df, train_labels, laplace = 0)
    predict_bayes1 <- c(predict_bayes1, predict(mod_bayes2, test_df))
    predict_bayes1_probs[[i]] <- predict(mod_bayes2, test_df, type="raw")
    
  }
  confusionMatrix_bayes1 <- confusionMatrix(predict_bayes1, datasc$class, 
                                            positive = levs[2])
  if(length(levs)==2){
    probs_vector <- map(predict_bayes1_probs, \(xx) xx %>% as.data.frame) %>% 
      bind_rows %>% pull(!!sym(levs[2]))
    roc1 <- roc(response=as.numeric(datasc$class)-1, predictor=probs_vector)
  }else{
    probs_vector <- map(predict_bayes1_probs, \(xx) xx %>% as.data.frame) %>% 
      bind_rows 
    #probs_vector2 <- map(predict1_probs, \(xx) attr(xx, "probabilities") %>% as.data.frame) %>% 
    #   map2(datasc$class, \(x, nn) x[as.character(nn)]) %>% unlist
    roc1 <- multiclass.roc(response=datasc$class, predictor=probs_vector)
  }
  modwithall <- naiveBayes(df, datasc$class, laplace = 0)
  predict_bayes2 <- predict(modwithall, df)
  confMatrix_bayes2_nol1o <- confusionMatrix(predict_bayes2, datasc$class, 
                                             positive = levs[2])
  return(list(confmat=confusionMatrix_bayes1, 
              confmat_no_l1o=confMatrix_bayes2_nol1o,
              preds=predict_bayes1, 
              pred_probs=probs_vector,
              preds_no_l1o=predict_bayes2, 
              mod=modwithall,
              roc_obj_no_l1o=NULL,
              roc_auc_no_l1o=NULL,
              roc_obj=roc1,
              roc_auc=as.numeric(roc1$auc))
         )
}

make_classifTree_l1o <- function(datasc, levs, varnames, 
                                 folds=c(), 
                                 balance_weights = TRUE,
                                 do_smote=FALSE, 
                                 smote_params=list(K=5, dup_size="balance")){
  library(C50)
  predict_tree1 <- factor()
  predict_probs <- list()
  df <- datasc %>% dplyr::select(-class, -sample)  %>% dplyr::select(all_of(varnames))
  if(balance_weights & !do_smote){
    warning("Using weights in C5.0 not implemented")
    classweights <- table(datasc$class)
    sweights <- max(classweights)/classweights[datasc$class]
  }else{
    sweights <- rep(1, nrow(datasc))
  }
  if(length(folds)==0){
    folds <- 1:nrow(datasc)
  }
  for(i in folds){
    # Separar datos
    train_df <- df[-i, ]
    test_df <- df[i, ]
    #train_weighs <- sweights[-i]
    
    # Separar clases
    train_labels <- datasc$class[-i]
    test_labels <- datasc$class[i]
    
    if(do_smote){
      form <- as.formula(paste0("class ~ ", paste(varnames, sep="+", collapse= "+")))
      smote_df <- datasc[-i, ] %>% select(class, all_of(varnames))
      smoteData <- SmoteClassif(form, smote_df, 
                                C.perc = smote_params$dup_size, 
                                k = smote_params$K, repl = FALSE,
                                dist = "Euclidean", p = 2)
      train_df <- smoteData %>% select(-class)
      train_labels <- factor(smoteData$class, levels=levs)
    }else{
      smoteData = NULL
    }
    
    mod_tree1 <- C5.0(train_df, train_labels, trials = 20) # , weights=train_weighs # makes it crash!
    predict_tree1 <- c(predict_tree1, predict(mod_tree1, test_df))
    predict_probs[[i]] <- predict(mod_tree1, test_df, type = "prob")
  }
  
  confmat_tree1 <- confusionMatrix(predict_tree1, datasc$class, positive = levs[2])
  if(length(levs)==2){
    probs_vector <- map(predict_probs, \(xx) xx %>% as.data.frame) %>% 
      bind_rows %>% pull(!!sym(levs[2]))
    roc1 <- roc(response=as.numeric(datasc$class)-1, predictor=probs_vector)
  }else{
    probs_vector <- map(predict_probs, \(xx) xx %>% as.data.frame) %>% 
      bind_rows #%>% pull(!!sym(levs[2]))
    roc1 <- multiclass.roc(response=datasc$class, predictor=probs_vector)
  }
  mod_all <- C5.0(df, datasc$class, trials = 20) # , weights=sweighs
  predict_tree2 <- predict(mod_all, df)
  confmat_tree2 <- confusionMatrix(predict_tree2, datasc$class, positive = levs[2])
  return(list(confmat=confmat_tree1, 
              mod=mod_all, 
              preds=predict_tree1, 
              pred_probs=probs_vector,
              preds_no_l1o=predict_tree2,
              confmat_no_l1o=confmat_tree2,
              roc_obj_no_l1o=NULL,
              roc_auc_no_l1o=NULL,
              roc_obj=roc1,
              roc_auc=as.numeric(roc1$auc)
  ))
}

make_randomForest_l1o <- function(datasc, levs, varnames, 
                                  folds=folds(), 
                                  randomforest_params = randomforest_params,
                                  do_smote=FALSE, 
                                  smote_params=list(K=5, dup_size="balance")
                                  ){
  library(randomForest)
  df <- datasc %>% dplyr::select(-class, -sample)  %>% dplyr::select(all_of(varnames))
  
  if(randomforest_params$balance_weights & !do_smote){
    classweights <- table(datasc$class)
    sweights <- max(classweights)/classweights[datasc$class]
  }else{
    sweights <- rep(1, nrow(datasc))
  }
  
  
  predict_tree1 <- factor()
  predict_tree1_probs <- list()
  
  
  if(length(folds)==0){
    folds <- 1:nrow(datasc)
  }
  for(i in folds){
    # Separar datos
    train_df <- df[-i, ]
    test_df <- df[i, ]
    train_weighs <- sweights[-i]
  
    # Separar clases
    train_labels <- datasc$class[-i]
    test_labels <- datasc$class[i]
    
    if(do_smote){
      form <- as.formula(paste0("class ~ ", paste(varnames, sep="+", collapse= "+")))
      smote_df <- datasc[-i, ] %>% select(class, all_of(varnames))
      smoteData <- SmoteClassif(form, smote_df, 
                                C.perc = smote_params$dup_size, 
                                k = smote_params$K, repl = FALSE,
                                dist = "Euclidean", p = 2)
      train_df <- smoteData %>% select(-class)
      train_labels <- factor(smoteData$class, levels=levs)
      train_weighs <- rep(1, nrow(train_df))
    }else{
      smoteData = NULL
    }
    mod_tree1 <- randomForest(x=train_df, y=train_labels, levels=levs, 
                              weights = train_weighs,
                              ntree = randomforest_params$ntree, 
                              mtry = randomforest_params$mtry, 
                              nodesize = randomforest_params$nodesize)
    predict_tree1 <- c(predict_tree1, predict(mod_tree1, test_df))
    predict_tree1_probs[[i]] <- predict(mod_tree1, test_df, type = "prob")
    
  }
  confmat_tree1 <- confusionMatrix(predict_tree1, datasc$class, positive = levs[2])
  if(length(levs)==2){
    probs_vector <- map(predict_tree1_probs, \(xx) xx %>% as.data.frame) %>% 
      bind_rows %>% pull(!!sym(levs[2]))
    roc1 <- roc(response=as.numeric(datasc$class)-1, predictor=probs_vector)
  }else{
    probs_vector <- predict_tree1_probs %>% map(as.data.frame) %>% bind_rows
    roc1 <- multiclass.roc(response=datasc$class, predictor=probs_vector)
    roc_auc <- as.numeric(roc1$auc)
  }
  
  mod_tree1 <- randomForest(x=df, y=datasc$class, levels=levs, 
                            weights = sweights,
                            ntree = randomforest_params$ntree, 
                            mtry = randomforest_params$mtry, 
                            nodesize = randomforest_params$nodesize)
  predict_tree2 <- predict(mod_tree1, df)
  confmat_tree2 <- confusionMatrix(predict_tree2, datasc$class, positive = levs[2])
  return(list(confmat=confmat_tree1, 
              confmat_no_l1o=confmat_tree2,
              mod=mod_tree1, 
              preds=predict_tree1, 
              pred_probs =probs_vector,
              preds_no_l1o=predict_tree2,
              roc_obj_no_l1o=NULL,
              roc_auc_no_l1o=NULL,
              roc_obj=roc1,
              roc_auc=as.numeric(roc1$auc), 
              smoteData=smoteData, 
              params =  randomforest_params))
}


make_xgboost_l1o <- function(datasc, levs, varnames, 
                                  folds=folds(), 
                             xgboost_params = xgboost_params,
                                  do_smote=FALSE, 
                                  smote_params=list(K=5, dup_size="balance")
){
  library(xgboost)
  df <- datasc %>% dplyr::select(-class, -sample)  %>% dplyr::select(all_of(varnames))
  if(length(folds)==0){
    folds <- 1:nrow(datasc)
  }
  predict_probs = numeric(0)
  
  if(xgboost_params$balance_weights){
    
    #class_weights <- class_weights/min(class_weights)
    #weights <- class_weights[datasc$class]
    
    class_weights <- 1/table(datasc$class)
    if(length(levs)==2){
      posweight <- class_weights[levs[1]]/class_weights[levs != levs[1]]
    } else{
      posweight <- as.vector(class_weights)
      names(posweight) <- names(class_weights)
    }
    
  }else{
    class_weights <- NULL
    if(length(levs)==2){
      posweight <- 1
    }else{
      posweight <- rep(1, length(levs))
      names(posweight) <- levs
    }
  } 
  
  for(i in folds){
    # Separar datos
    train_df <- df[-i, ]
    test_df <- df[i, ]
    
    # Separar clases
    train_labels <- datasc$class[-i]
    test_labels <- datasc$class[i]
    #if(xgboost_params$balance_weights){
    #  train_weights <- weights[-i]
    #}else{
    #  train_weights <- NULL
    #}
    
    if(do_smote){
      form <- as.formula(paste0("class ~ ", paste(varnames, sep="+", collapse= "+")))
      smote_df <- datasc[-i, ] %>% select(class, all_of(varnames))
      smoteData <- SmoteClassif(form, smote_df, 
                                C.perc = smote_params$dup_size, 
                                k = smote_params$K, repl = FALSE,
                                dist = "Euclidean", p = 2)
      train_df <- smoteData %>% select(-class)
      train_labels <- factor(smoteData$class, levels=levs)
      #if(xgboost_params$balance_weights){
      #  train_weights <- class_weights[train_labels]
      #}
    }else{
      smoteData = NULL
    }
    mod_tree1 <- xgboost(x = train_df, y = train_labels,
                         #weights=train_weights,
                         scale_pos_weight = posweight,
                         max_depth = xgboost_params$max_depth, 
                         learning_rate = xgboost_params$learning_rate,
                         nrounds = xgboost_params$nrounds,
                         min_child_weight = xgboost_params$min_child_weight,
                         subsample = xgboost_params$subsample,
                         colsample_bytree = xgboost_params$colsample_bytree,
                         gamma = xgboost_params$gamma,
                         reg_lambda=xgboost_params$reg_lambda,
                         reg_alpha= xgboost_params$reg_alpha,
                         nthread = xgboost_params$nthread, 
                         objective = xgboost_params$objective)
    if(length(levs) == 2){
      predict_probs <- c(predict_probs, predict(mod_tree1, test_df))
    }else{
      predict_probs <- rbind(predict_probs, predict(mod_tree1, test_df))
    }
    
    
  }
  if(length(levs) == 2){
    predict_tree1 <- factor(levs[as.integer(round(predict_probs))+1], levels=levs)
    roc1 <- roc(response=as.numeric(datasc$class)-1, predictor=predict_probs)
  }else{
    predict_tree1 <- factor(colnames(predict_probs)[apply(predict_probs, MAR=1, which.max)], levels=levs)
    roc1 <- multiclass.roc(response=datasc$class, predictor=predict_probs)
  }
  confmat_tree1 <- confusionMatrix(predict_tree1, datasc$class, positive = levs[2])
  
  mod_tree1 <-  xgboost(x = df, y = datasc$class,
                        #weights=weights,
                        scale_pos_weight = posweight,
                        max_depth = xgboost_params$max_depth, 
                        learning_rate = xgboost_params$learning_rate,
                        nrounds = xgboost_params$nrounds,
                        min_child_weight = xgboost_params$min_child_weight,
                        subsample = xgboost_params$subsample,
                        colsample_bytree = xgboost_params$colsample_bytree,
                        gamma = xgboost_params$gamma,
                        reg_lambda=xgboost_params$reg_lambda,
                        reg_alpha= xgboost_params$reg_alpha,
                        nthread = xgboost_params$nthread, 
                        objective = xgboost_params$objective)
  
  predict_tree2 <- predict(mod_tree1, df)
  if(length(levs) == 2){
    predict_tree2 <- factor(levs[as.integer(round(predict_tree2))+1], levels=levs)
  }else{
    predict_tree2 <- factor(levs[apply(predict_tree2, MAR=1, which.max)], levels=levs)
  }
  confmat_tree2 <- confusionMatrix(predict_tree2, datasc$class, positive = levs[2])
  
  return(list(confmat=confmat_tree1, 
              confmat_no_l1o=confmat_tree2,
              mod=mod_tree1, 
              preds=predict_tree1, 
              pred_probs = predict_probs, 
              preds_no_l1o=predict_tree2,
              roc_obj_no_l1o=NULL,
              roc_auc_no_l1o=NULL,
              roc_obj=roc1,
              roc_auc=as.numeric(roc1$auc),
              xgboost_params = xgboost_params,
              smoteData=smoteData))
}


make_catboost_l1o <- function(datasc, levs, varnames, 
                             folds=folds(), 
                             catboost_params = catboost_params,
                             do_smote=FALSE, 
                             smote_params=list(K=5, dup_size="balance")
){
  library(catboost)
  datasc$class <- factor(datasc$class, levels=levs)
  df <- datasc %>% dplyr::select(-class, -sample)  %>% dplyr::select(all_of(varnames))
  if(length(folds)==0){
    folds <- 1:nrow(datasc)
  }
  predict_probs = numeric(0)
  
  if(catboost_params$balance_weights & ! do_smote){
    class_weights <- table(datasc$class)
    class_weights_vec <- max(class_weights)/class_weights %>% as.vector
    names(class_weights_vec) <- levs
  } else {
    class_weights_vec <- rep(1, length(levs))
  }
  
  for(i in folds){
    # Separar datos
    train_df <- df[-i, ]
    test_df <- df[i, ]
    
    # Separar clases
    train_labels <- datasc$class[-i]
    test_labels <- datasc$class[i]
    #if(xgboost_params$balance_weights){
    #  train_weights <- weights[-i]
    #}else{
    #  train_weights <- NULL
    #}
    
    if(do_smote){
      form <- as.formula(paste0("class ~ ", paste(varnames, sep="+", collapse= "+")))
      smote_df <- datasc[-i, ] %>% select(class, all_of(varnames))
      #smoteData <- SMOTE(train_df, train_labels, K=smote_params$K, dup_size = smote_params$dup_size)
      smoteData <- SmoteClassif(form, smote_df, 
                                C.perc = smote_params$dup_size, 
                                k = smote_params$K, repl = FALSE,
                   dist = "Euclidean", p = 2)
      
      train_df <- smoteData %>% select(-class)
      train_labels <- factor(smoteData$class, levels=levs)
      #if(xgboost_params$balance_weights){
      #  train_weights <- class_weights[train_labels]
      #}
    }else{
      smoteData = NULL
    }
    if(length(levs) == 2){
      train_pool <- catboost.load_pool(data = train_df, label = as.integer(train_labels == levs[2]))
    }else{
      train_pool <- catboost.load_pool(data = train_df, label = as.integer(train_labels)-1)
    }
    test_pool <- catboost.load_pool(data = test_df)
    
    if(catboost_params$bootstrap_type == "Bernoulli"){
      model <- catboost.train(learn_pool = train_pool, params = list(
        depth = catboost_params$depth,
        learning_rate = catboost_params$learning_rate,
        iterations = catboost_params$iterations,
        loss_function = catboost_params$loss_function,
        eval_metric = catboost_params$eval_metric,
        bootstrap_type = catboost_params$bootstrap_type,
        l2_leaf_reg = catboost_params$l2_leaf_reg,
        subsample = catboost_params$subsample,
        grow_policy = catboost_params$grow_policy,
        auto_class_weights= catboost_params$auto_class_weights,
        thread_count = catboost_params$thread_count,
        #class_weights = class_weights_vec,
        logging_level = "Silent"
      ))
    }else{
      model <- catboost.train(learn_pool = train_pool, params = list(
        depth = catboost_params$depth,
        learning_rate = catboost_params$learning_rate,
        iterations = catboost_params$iterations,
        loss_function = catboost_params$loss_function,
        eval_metric = catboost_params$eval_metric,
        bootstrap_type = catboost_params$bootstrap_type,
        l2_leaf_reg = catboost_params$l2_leaf_reg,
        #subsample = catboost_params$subsample,
        grow_policy = catboost_params$grow_policy,
        auto_class_weights= catboost_params$auto_class_weights,
        thread_count = catboost_params$thread_count,
        #class_weights = class_weights_vec,
        logging_level = "Silent"
      ))
    }
    pred_prob <- catboost.predict(model, test_pool, prediction_type = "Probability")
    
    if(length(levs) == 2){
       predict_probs <- c(predict_probs, pred_prob)
    }else{
      predict_probs <- rbind(predict_probs, pred_prob)
    }
    
  }
  if(length(levs) == 2){
    predict_tree1 <- factor(levs[as.integer(round(predict_probs))+1], levels=levs)
    roc1 <- roc(response=as.numeric(datasc$class)-1, predictor=predict_probs)
  }else{
    predict_tree1 <- factor(levs[apply(predict_probs, MAR=1, which.max)], levels=levs)
    colnames(predict_probs) <- as.character(levels(datasc$class))
    roc1 <- multiclass.roc(response=datasc$class, predictor=predict_probs)
  }
  confmat_tree1 <- confusionMatrix(predict_tree1, datasc$class, positive = levs[2])
  
  train_pool <- catboost.load_pool(data = df, label = as.integer(datasc$class == levs[2]))
  test_pool <- catboost.load_pool(data = df)
  
  if(catboost_params$bootstrap_type == "Bernoulli"){
   model2 <- catboost.train(learn_pool = train_pool, params = list(
     depth = catboost_params$depth,
     learning_rate = catboost_params$learning_rate,
     iterations = catboost_params$iterations,
     loss_function = catboost_params$loss_function,
     eval_metric = catboost_params$eval_metric,
     bootstrap_type = catboost_params$bootstrap_type,
     l2_leaf_reg = catboost_params$l2_leaf_reg,
     subsample = catboost_params$subsample,
     grow_policy = catboost_params$grow_policy,
     auto_class_weights= catboost_params$auto_class_weights,
     thread_count = catboost_params$thread_count,
     #class_weights = class_weights_list,
     logging_level = "Silent"
  ))
  }else{
    model2 <- catboost.train(learn_pool = train_pool, params = list(
      depth = catboost_params$depth,
      learning_rate = catboost_params$learning_rate,
      iterations = catboost_params$iterations,
      loss_function = catboost_params$loss_function,
      eval_metric = catboost_params$eval_metric,
      bootstrap_type = catboost_params$bootstrap_type,
      l2_leaf_reg = catboost_params$l2_leaf_reg,
      #subsample = catboost_params$subsample,
      grow_policy = catboost_params$grow_policy,
      auto_class_weights= catboost_params$auto_class_weights,
      thread_count = catboost_params$thread_count,
      #class_weights = class_weights_list,
      logging_level = "Silent"
    ))
  }
  predict_tree2 <- catboost.predict(model2, test_pool, prediction_type = "Probability")
  if(length(levs) == 2){
    predict_tree2 <- factor(levs[as.integer(round(predict_tree2))+1], levels=levs)
  }else{
    predict_tree2 <- factor(levs[apply(predict_tree2, MAR=1, which.max)], levels=levs)
  }
  confmat_tree2 <- confusionMatrix(predict_tree2, datasc$class, positive = levs[2])
  
  return(list(confmat=confmat_tree1, 
              confmat_no_l1o=confmat_tree2,
              mod=model2, 
              preds=predict_tree1, 
              pred_probs = predict_probs, 
              preds_no_l1o=predict_tree2,
              roc_obj_no_l1o=NULL,
              roc_auc_no_l1o=NULL,
              roc_obj=roc1,
              roc_auc=as.numeric(roc1$auc),
              xgboost_params = catboost_params,
              smoteData=smoteData))
}


make_svm_l1o <- function(datasc, levs, varnames, kernel="linear", SEED=123, folds=c(), 
                         do_smote=FALSE, 
                         smote_params=list(K=5, dup_size=2), 
                         balance_classes=TRUE){
  library(e1071)
  datasc$class <- factor(datasc$class)
  df <- datasc %>% dplyr::select(-class, -sample)  %>% dplyr::select(all_of(varnames))
  predict1 <- factor(levels=levs)
  predict1_probs <- list()
  
  if(balance_classes & ! do_smote){
    class_weight <- "inverse"
  }else{
    class_weight <- rep(1, length(levs))
    names(class_weight) <- levs
  }
  
  if(length(folds)==0){
    folds <- 1:nrow(datasc)
  }
  set.seed(SEED)
  for(i in folds){
    # Separar datos
    train_df <- df[-i, ]
    test_df <- df[i, ]
    
    # Separar clases
    train_labels <- datasc$class[-i]
    test_labels <- datasc$class[i]
    
    if(do_smote){
      form <- as.formula(paste0("class ~ ", paste(varnames, sep="+", collapse= "+")))
      smote_df <- datasc[-i, ] %>% select(class, all_of(varnames))
      smoteData <- SmoteClassif(form, smote_df, 
                                C.perc = smote_params$dup_size, 
                                k = smote_params$K, repl = FALSE,
                                dist = "Euclidean", p = 2)
      train_df <- smoteData %>% dplyr::select(-class)
      train_labels <- factor(smoteData$class, levels=levs)
    }else{
      smoteData = NULL
    }
    
    mod <- e1071::svm(x = train_df, y = train_labels, scale=TRUE, kernel=kernel, 
                      class.weights = class_weight,
                      probability = TRUE)
    predict1 <- c(predict1, predict(mod, test_df))
    predict1_probs[[i]] <- predict(mod, test_df, probability = TRUE)
    
  }
  
  confmat1 <- confusionMatrix(predict1, datasc$class, positive = levs[2])
  if(length(levs)==2){
    probs_vector <- map(predict1_probs, \(xx) attr(xx, "probabilities") %>% as.data.frame) %>% 
      bind_rows %>% pull(!!sym(levs[2]))
    roc1 <- roc(response=as.numeric(datasc$class)-1, predictor=probs_vector)
  }else{
    probs_vector <- map(predict1_probs, \(xx) attr(xx, "probabilities") %>% as.data.frame) %>% 
      bind_rows #%>% pull(!!sym(levs[2]))
    roc1 <- multiclass.roc(response=datasc$class, predictor=probs_vector)
    roc_auc <- as.numeric(roc1$auc)
  }
  
  mod_all <- e1071::svm(x = df, y = datasc$class, scale=TRUE, kernel=kernel, 
                        probability = TRUE, class.weights = class_weight)
  predict2 <- predict(mod_all, df)
  predict2_probs <- predict(mod_all, df, probability = TRUE)
  confmat2 <- confusionMatrix(predict2, datasc$class, positive = levs[2])
  
  mod_all_noscale <- e1071::svm(x = df[, varnames[1:2]], y = datasc$class, scale=FALSE, 
                         kernel=kernel, class.weights = "inverse")
  predict2_noscale <- predict(mod_all_noscale, df[, varnames[1:2]])
  confmat2_noscale <- confusionMatrix(predict2_noscale, datasc$class, positive = levs[2])
  
  return(list(confmat=confmat1, 
              confmat_no_l1o=confmat2,
              mod=mod_all, 
              preds=predict1, 
              pred_probs = probs_vector,
              pred_probs_obj = predict1_probs,
              preds_no_l1o=predict2,
              mod_noscale=mod_all_noscale, 
              preds_noscale=predict2_noscale, 
              confmat_noscale=confmat2_noscale,
              roc_obj_no_l1o=NULL,
              roc_auc_no_l1o=NULL,
              roc_obj=roc1,
              roc_auc=as.numeric(roc1$auc))
         )
}

make_glm_l1o <- function(datasc, levs, varnames, folds= c(), 
                         do_smote=FALSE, 
                         smote_params=list(K=5, dup_size=2)){
  predict_glm1 <- c()
  df <- datasc %>% dplyr::select(-sample) %>% dplyr::select(class, all_of(varnames))
  formula <- paste0("class ~ ", paste(varnames, sep="+", collapse="+")) %>% as.formula()
  
  if(length(folds)==0){
    folds <- 1:nrow(datasc)
  }
  
  for(i in folds){
    # Separar datos
    train_df <- df[-i, ]
    test_df <- df[i, ]
    
    if(do_smote){
      form <- as.formula(paste0("class ~ ", paste(varnames, sep="+", collapse= "+")))
      smote_df <- datasc[-i, ] %>% select(class, all_of(varnames))
      smoteData <- SmoteClassif(form, smote_df, 
                                C.perc = smote_params$dup_size, 
                                k = smote_params$K, repl = FALSE,
                                dist = "Euclidean", p = 2)
      train_df <- smoteData %>% dplyr::mutate(class=factor(class, levels=levs))
    }else{
      smoteData = NULL
    }
    
    mod_glm <- glm(formula, data=train_df, family = binomial)
    predict_glm1 <- c(predict_glm1, predict(mod_glm, test_df, type = "response"))
    
  }
  predict1 <- ifelse(predict_glm1 > 0.5, levs[2], levs[1]) %>% factor(levels=levs)
  confmat1 <- confusionMatrix(predict1, datasc$class, positive = levs[2])
  roc1 <- roc(response=as.numeric(datasc$class)-1, predictor=predict_glm1)
  
  mod_all <- glm(formula, data=datasc, family = binomial)
  predict2_probs <- predict(mod_all, df, type = "response")
  predict2 <- ifelse(predict2_probs > 0.5, levs[2], levs[1]) %>% factor(levels=levs)
  confmat2 <- confusionMatrix(predict2, datasc$class, positive = levs[2])
  roc2 <- roc(response=as.numeric(datasc$class)-1, predictor=predict2_probs)
  return(list(confmat=confmat1, 
              confmat_no_l1o=confmat2,
              mod=mod_all, 
              preds=predict1, 
              pred_probs = predict_glm1,
              preds_no_l1o=predict2, 
              pred_probs_no_l1o=predict2_probs,
              roc_obj_no_l1o=roc2,
              roc_auc_no_l1o=as.numeric(roc2$auc),
              roc_obj=roc1,
              roc_auc=as.numeric(roc1$auc)
              ))
}


make_glm_l1o_multiclass <- function(datasc, levs, varnames, folds=c(), 
                                    do_smote=FALSE, 
                                    smote_params=list(K=5, dup_size=2)){
  library(nnet)
  predict_glm1 <- c()
  df <- datasc %>% dplyr::select(-sample) %>% dplyr::select(class, all_of(varnames))
  formula <- paste0("class ~ ", paste(varnames, sep="+", collapse="+")) %>% as.formula()
  
  if(length(folds)==0){
    folds <- 1:nrow(datasc)
  }
  
  for(i in folds){
    # Separar datos
    train_df <- df[-i, ]
    test_df <- df[i, ]
    
    if(do_smote){
      form <- as.formula(paste0("class ~ ", paste(varnames, sep="+", collapse= "+")))
      smote_df <- datasc[-i, ] %>% select(class, all_of(varnames))
      smoteData <- SmoteClassif(form, smote_df, 
                                C.perc = smote_params$dup_size, 
                                k = smote_params$K, repl = FALSE,
                                dist = "Euclidean", p = 2)
      train_df <- smoteData %>% dplyr::mutate(class=factor(class, levels=levs))
    }else{
      smoteData = NULL
    }
    
    mod_glm <- multinom(formula, data=train_df)
    newpred <- predict(mod_glm, test_df, type="prob")
    predict_glm1 <- rbind(predict_glm1, newpred)
    
  }
  
  predict1<- colnames(predict_glm1)[apply(predict_glm1, MAR=1, which.max)] %>% factor
  confmat1 <- confusionMatrix(predict1, factor(datasc$class))
  
  assertthat::assert_that(all(colnames(predict_glm1) == levels(datasc$class)))
  roc1 <- multiclass.roc(response=datasc$class, predictor=predict_glm1)
  roc_auc <- as.numeric(roc1$auc)
  
  mod_all <- multinom(formula, data=datasc, family = binomial)
  predict2 <- predict(mod_all, df, type="probs")
  classes2 <- colnames(predict2)[apply(predict2, MAR=1, \(x)which(x==max(x)))] %>% factor
  confmat2 <- confusionMatrix(classes2, factor(datasc$class))
  
  roc_obj_fullmod <- apply(predict2, MAR=2, \(x) multiclass.roc(datasc$class, x))
  roc_auc_fullmod <- sapply(roc_obj_fullmod, \(x)x$auc) %>% mean
  return(list(confmat=confmat1, 
              confmat_no_l1o=confmat2,
              mod=mod_all, 
              preds=predict1, 
              preds_no_l1o=classes2,
              roc_obj_no_l1o=roc_obj_fullmod,
              roc_auc_no_l1o=roc_auc_fullmod,
              roc_obj=roc1,
              roc_auc=roc_auc
  ))
}
get_signif_components <- function(datasc, levs){
  df <- datasc %>% dplyr::select(-sample) 
  varnames <- names(df)[names(df)!="class"]
  res <- data.frame()
  ps <- c()
  for(i in varnames){
    df$aux <- df[, i]
    mod_glm <- glm(class ~ aux, data=df, family = binomial)
    ps <- c(ps, summary(mod_glm)$coefficients[2, 4])
    predict_glm1 <- predict(mod_glm, df)
    predict_glm1 <- ifelse(predict_glm1 > 0.5, levs[2], levs[1]) %>% factor(levels=levs)
    confmat <- confusionMatrix(predict_glm1, datasc$class, positive = levs[2])
    auxdf <- data.frame(var=i, 
                        pval=summary(mod_glm)$coefficients[2, 4],
                        Accuracy=confmat$overall["Accuracy"],
                        Sensitivity = confmat$byClass["Sensitivity"],
                        Specificity = confmat$byClass["Specificity"],
                        PPV = confmat$byClass["Pos Pred Value"],
                        NPV = confmat$byClass["Neg Pred Value"]
    )
    res <- rbind(res, auxdf)
  }
  return(res %>% dplyr::arrange(pval))
}

get_signif_components_multiclass <- function(datasc, levs, plim=0.05){
  df <- datasc %>% dplyr::select(-sample) %>% dplyr::filter(!is.na(class))
  varnames <- names(df)[names(df)!="class"]
  res <- data.frame()
  ps <- c()
  for(i in varnames){
    df$aux <- df[, i]
    mod <- mod <- lm(aux ~ class, data=df)
    modsum <- summary(mod)
    
    any_sig <- any(modsum$coefficients[2:length(levs), 4] < plim)
    which_sig = which(modsum$coefficients[2:length(levs), 4] < 0.05) %>% names %>% paste(collapse="_")
    auxdf <- data.frame(var=i, 
                        any_sig = any_sig,
                        which_sig = which_sig
    ) %>% cbind(broom::glance(mod))
    res <- rbind(res, auxdf) 
  }
  return(res %>% dplyr::rename(pval = p.value) %>% dplyr::arrange(pval))
}

getTableFromConfmatrices <- function(modlist){
  res <- map(modlist, \(mod){
    data.frame(
      Accuracy_l1out=if(is.null(mod$confmat)) NA else mod$confmat$overall["Accuracy"],
      Kappa_l1out=if(is.null(mod$confmat)) NA else mod$confmat$overall["Kappa"],
      Sensitivity_l1out = if(is.null(mod$confmat)) NA else mod$confmat$byClass["Sensitivity"],
      Specificity_l1out = if(is.null(mod$confmat)) NA else mod$confmat$byClass["Specificity"],
      PPV_l1out = if(is.null(mod$confmat)) NA else mod$confmat$byClass["Pos Pred Value"],
      NPV_l1out = if(is.null(mod$confmat)) NA else mod$confmat$byClass["Neg Pred Value"],
      Precision_l1out = if(is.null(mod$confmat)) NA else mod$confmat$byClass["Precision"],
      Recall_l1out = if(is.null(mod$confmat)) NA else mod$confmat$byClass["Recall"],
      BalancedAccuracy_l1out = if(is.null(mod$confmat)) NA else mod$confmat$byClass["Balanced Accuracy"],
      AUC_l1out = if(is.null(mod$roc_auc)) NA else mod$roc_auc,
      Accuracy=if(is.null(mod$confmat_no_l1o)) NA else mod$confmat_no_l1o$overall["Accuracy"],
      Kappa=if(is.null(mod$confmat_no_l1o)) NA else mod$confmat_no_l1o$overall["Kappa"],
      Sensitivity = if(is.null(mod$confmat_no_l1o)) NA else mod$confmat_no_l1o$byClass["Sensitivity"],
      Specificity = if(is.null(mod$confmat_no_l1o)) NA else mod$confmat_no_l1o$byClass["Specificity"],
      PPV = if(is.null(mod$confmat_no_l1o)) NA else mod$confmat_no_l1o$byClass["Pos Pred Value"],
      NPV = if(is.null(mod$confmat_no_l1o)) NA else mod$confmat_no_l1o$byClass["Neg Pred Value"],
      Precision = if(is.null(mod$confmat_no_l1o)) NA else mod$confmat_no_l1o$byClass["Precision"],
      Recall = if(is.null(mod$confmat_no_l1o)) NA else mod$confmat_no_l1o$byClass["Recall"],
      BalancedAccuracy = if(is.null(mod$confmat_no_l1o)) NA else mod$confmat_no_l1o$byClass["Balanced Accuracy"]
    ) 
  }) %>% bind_rows() %>% 
    dplyr::mutate(model = names(modlist)) %>% 
    dplyr::select(model, everything()) %>% 
    arrange(desc(Accuracy_l1out))
  rownames(res) <- NULL
  return(res)
}

getTableFromConfmatrices_multiclass <- function(modlist){
  res <- map(modlist, \(mod){
    data.frame(
      Accuracy_l1out=if(is.null(mod$confmat)) NA else mod$confmat$overall["Accuracy"] %>% mean,
      Kappa_l1out=if(is.null(mod$confmat)) NA else mod$confmat$overall["Kappa"] %>% mean,
      Sensitivity_l1out = if(is.null(mod$confmat)) NA else mod$confmat$byClass[, "Sensitivity"] %>% mean,
      Specificity_l1out = if(is.null(mod$confmat)) NA else mod$confmat$byClass[, "Specificity"] %>% mean,
      PPV_l1out = if(is.null(mod$confmat)) NA else mod$confmat$byClass[, "Pos Pred Value"] %>% mean,
      NPV_l1out = if(is.null(mod$confmat)) NA else mod$confmat$byClass[, "Neg Pred Value"] %>% mean,
      Precision_l1out = if(is.null(mod$confmat)) NA else mod$confmat$byClass[, "Precision"] %>% mean,
      Recall_l1out = if(is.null(mod$confmat)) NA else mod$confmat$byClass[, "Recall"] %>% mean,
      BalancedAccuracy_l1out = if(is.null(mod$confmat)) NA else mod$confmat$byClass[, "Balanced Accuracy"] %>% mean,
      AUC_l1out = if(is.null(mod$roc_auc)) NA else mod$roc_auc,
      Accuracy=if(is.null(mod$confmat_no_l1o)) NA else mod$confmat_no_l1o$overall["Accuracy"] %>% mean,
      Kappa=if(is.null(mod$confmat_no_l1o)) NA else mod$confmat_no_l1o$overall["Kappa"] %>% mean,
      Sensitivity = if(is.null(mod$confmat_no_l1o)) NA else mod$confmat_no_l1o$byClass[, "Sensitivity"] %>% mean,
      Specificity = if(is.null(mod$confmat_no_l1o)) NA else mod$confmat_no_l1o$byClass[, "Specificity"] %>% mean,
      PPV = if(is.null(mod$confmat_no_l1o)) NA else mod$confmat_no_l1o$byClass[, "Pos Pred Value"] %>% mean,
      NPV = if(is.null(mod$confmat_no_l1o)) NA else mod$confmat_no_l1o$byClass[, "Neg Pred Value"] %>% mean,
      Precision = if(is.null(mod$confmat_no_l1o)) NA else mod$confmat_no_l1o$byClass[, "Precision"] %>% mean,
      Recall = if(is.null(mod$confmat_no_l1o)) NA else mod$confmat_no_l1o$byClass[, "Recall"] %>% mean,
      BalancedAccuracy = if(is.null(mod$confmat_no_l1o)) NA else mod$confmat$byClass[, "Balanced Accuracy"] %>% mean
    ) 
  }) %>% bind_rows() %>% 
    dplyr::mutate(model = names(modlist)) %>% 
    dplyr::select(model, everything()) %>% 
    arrange(desc(Accuracy_l1out))
  rownames(res) <- NULL
  return(res)
}

make_ensemble_votes <- function(datasc, levs, modlist, model_res, param="Kappa_l1out", 
                          min_val=0, prop=TRUE, 
                          only_1_knn=FALSE){ # 0.65
  remove_knn <- model_res %>% 
    dplyr::arrange(desc(!!sym(param))) %>% 
    dplyr::filter(grepl("KNN", model)) %>% pull(model)
  remove_knn <- remove_knn[2:length(remove_knn)]
  m2use <- model_res %>% 
    dplyr::filter(!!sym(param) >= min_val) %>% 
    dplyr::filter(! (model %in% remove_knn & only_1_knn)) %>%  
    pull(model)
  preddf <- map(m2use, \(x) tibble( !!x := modlist[[x]]$preds))  %>% bind_cols()
  #names(preddf) <- m2use
  if(prop){
    ponderfac <- model_res[match(m2use, model_res$model), param]
    ponderfac <- (ponderfac - min(ponderfac))/(max(ponderfac) - min(ponderfac)) + 0.1
  }else{
    ponderfac <- rep(1, length(m2use))
  }
  preds <- c()
  votes <- list()
  for(i in 1:nrow(preddf)){
    classcore <- map_vec(levs, \(ll) sum(ponderfac[preddf[i, ] == ll]))
    names(classcore) <- levs
    l1 <- length(preds)
    preds <- c(preds, levs[which.max(classcore)] )
    l2 <- length(preds)
    cat(i, ": L1=", l1, ", L2=", l2, ifelse(l1==l2, " --WARNING--", ""),  "\n")
    votes[[i]] <- classcore
  }
  preds <- factor(preds, levels=levs)
  confmat1 <- confusionMatrix(preds, datasc$class, positive = levs[2])
  
  if(length(levs) == 2){
    probs <- map_vec(votes, \(x)x[levs[2]]/sum(x) )
    roc1 <- roc(response=as.numeric(datasc$class)-1, predictor=probs)
  }else{
    probs <- purrr::map(votes, .f = \(x) x/sum(x)) %>% 
      bind_rows %>% as.matrix #%>% pull(!!sym(levs[2]))
    roc1 <- multiclass.roc(response=datasc$class, predictor=probs)
  }
  return(list(confmat=confmat1, 
              confmat_no_l1o=NULL,
              mod=NULL, 
              preds=preds, 
              pred_df=preddf,
              preds_no_l1o=NULL,
              roc_obj_no_l1o=NULL,
              roc_auc_no_l1o=NULL,
              roc_obj=roc1,
              roc_auc=as.numeric(roc1$auc)
  ))
}

make_ensemble_probs <- function(datasc, levs, modlist, model_res, param="BalancedAccuracy_l1out", #"Kappa_l1out", 
                                min_val=0, prop=TRUE, 
                                only_1_knn=FALSE){ # 0.65
  remove_knn <- model_res %>% 
    dplyr::arrange(desc(!!sym(param))) %>% 
    dplyr::filter(grepl("KNN", model)) %>% pull(model)
  remove_knn <- remove_knn[2:length(remove_knn)]
  m2use <- model_res %>% 
    dplyr::filter(!!sym(param) >= min_val) %>% 
    dplyr::filter(! (model %in% remove_knn & only_1_knn)) %>%  
    pull(model)
  preddf <- map(m2use, \(x) modlist[[x]]$pred_probs)  %>% bind_cols()
  names(preddf) <- m2use
  if(prop){
    ponderfac <- model_res[match(m2use, model_res$model), param]
    ponderfac <- (ponderfac - min(ponderfac))/(max(ponderfac) - min(ponderfac)) + 0.1
    
  }else{
    ponderfac <- rep(1, length(m2use))
  }
  preds <- c()
  avg_probs <- c()
  for(i in 1:nrow(preddf)){
    classcore <-sum(preddf[i, ]*ponderfac)/sum(ponderfac)
    avg_probs <- c(avg_probs, classcore)
    preds <- c(preds, levs[as.integer(round(classcore))+1] )
  }
  preds <- factor(preds, levels=levs)
  confmat1 <- confusionMatrix(preds, datasc$class, positive = levs[2])
  
  roc1 <- roc(response=as.numeric(datasc$class)-1, predictor=avg_probs)
  
  return(list(confmat=confmat1, 
              confmat_no_l1o=NULL,
              mod=NULL, 
              preds=preds, 
              pred_probs=avg_probs, 
              pred_df=preddf,
              preds_no_l1o=NULL,
              roc_obj_no_l1o=NULL,
              roc_auc_no_l1o=NULL,
              roc_obj=roc1,
              roc_auc=as.numeric(roc1$auc)
  ))
}

makeAllModels <- function(datasc, plim=0.01, opt, name="Condition", nfolds=0, 
                          xgboost_params = xgboost_params,
                          catboost_params = catboost_params, 
                          randomforest_params = randomforest_params,
                          do_smote=FALSE,
                          smote_params=smote_params, 
                          ensemble_param = "BalancedAccuracy_l1out",
                          ensemble_minval = 0, 
                          ensemble_1knn = FALSE, 
                          do_ensemble_probs=TRUE){
  levs <- datasc %>% pull(class) %>% as.factor %>% levels
  # Select features
  if(length(levs)==2){
    compsig <- get_signif_components(datasc, levs)
  }else{
    compsig <- get_signif_components_multiclass(datasc, levs)
  }
  
  tryCatch({readr::write_tsv(compsig, file=paste0(opt$out, "significant_PCAcomponents_", name,".tsv"))},
           error = function(x){print("ERROR writting sig Components")})
  
  varnames <- c(compsig$var[compsig$pval <= plim])
  if(length(varnames) < 2){
    varnames <- compsig %>% arrange(pval) %>% head(2) %>% pull(var)
  }
  
  if(nfolds == 0){
    folds <- c() ## leave 1 out
  }else{
    folds <- createFolds(datasc$class, k = nfolds, list = TRUE, returnTrain = FALSE)
  }
  
  
  if(length(levs)==2){
    res_glms <- make_glm_l1o(datasc, levs, varnames, folds = folds, do_smote = do_smote, smote_params = smote_params)
  }else{
    res_glms <- make_glm_l1o_multiclass(datasc, levs, varnames, folds = folds, do_smote = do_smote, smote_params = smote_params)
  }
  cat("-- GLM finished\n")
  res_svm_lin <- make_svm_l1o(datasc, levs, varnames, kernel="linear", folds = folds, do_smote = do_smote, smote_params = smote_params, balance_classes = TRUE)
  res_svm_rad <- make_svm_l1o(datasc, levs, varnames, kernel="radial", folds = folds, do_smote = do_smote, smote_params = smote_params, balance_classes = FALSE)
  cat("-- SVMs finished\n")
  res_randfor <- make_randomForest_l1o(datasc, levs, varnames, folds = folds, do_smote = do_smote, smote_params = smote_params, randomforest_params = randomforest_params)
  cat("-- RandomForest finished\n")
  res_tree <- make_classifTree_l1o(datasc, levs, varnames, folds = folds, do_smote = do_smote, smote_params = smote_params, balance_weights = TRUE)
  cat("-- C5.0 Tree finished\n")
  res_naivebayes <- makeNaiveBayes_l1o(datasc, levs, varnames, SEED=SEED, folds = folds, do_smote = do_smote, smote_params = smote_params)
  cat("-- NaiveBayes finished\n")
  res_knn_l1o <- makeKnn_l1o(datasc, levs, varnames, different_ks=seq(3,11, by=2), folds = folds, do_smote = do_smote, smote_params = smote_params)
  #res_knn_no_l1o <- makeKnn(datasc, levs, varnames, different_ks=seq(1,13, by=2))
  cat("-- KNN finished\n")
  res_kmeans_l1o <- makeKmeans_l1o(datasc, levs, varnames, SEED=SEED, folds = folds, do_smote = do_smote, smote_params = smote_params)
  cat("-- K-Means finished\n")
  res_xgboost <- make_xgboost_l1o(datasc, levs, varnames, xgboost_params = xgboost_params, folds = folds, do_smote = do_smote, smote_params = smote_params)
  cat("-- XGBoost finished\n")
  res_catboost <- make_catboost_l1o(datasc, levs, varnames, catboost_params = catboost_params, folds = folds, do_smote = do_smote, smote_params = smote_params)
  cat("-- CatBoost finished\n")
  
  modlist <- list("logistic_regression" = res_glms, 
                  "SVM-linear"=res_svm_lin, 
                  "SVM-radial"=res_svm_rad,
                  "RandomForest"=res_randfor,
                  "C5.0 Tree"=res_tree,
                  "NaiveBayes"=res_naivebayes,
                  "XGBoost"=res_xgboost,
                  "CatBoost"=res_catboost,
                  "KMeans"=res_kmeans_l1o)
  if(length(levs)>2)names(modlist)[[1]] <- "Multinom"
  for(k in names(res_knn_l1o)) modlist[[paste0("KNN-", k)]] <- res_knn_l1o[[k]]
  save(modlist, file=paste0(opt$out, "all_models_", name, ".RData"))
  
  if(length(levs)==2){
    model_res <- getTableFromConfmatrices(modlist)
  }else{
    model_res <- getTableFromConfmatrices_multiclass(modlist)
  }
  
  modlist$Ensemble <- make_ensemble_votes(datasc, levs, modlist, model_res, 
                                          param = ensemble_param, 
                                          min_val = ensemble_minval, 
                                          only_1_knn = ensemble_1knn)
  cat("-- Ensemble finished\n")
  if(do_ensemble_probs){
    modlist$Ensemble2 <- make_ensemble_probs(datasc, levs, modlist, model_res, param = ensemble_param, min_val = ensemble_minval, only_1_knn = ensemble_1knn)
  }
  if(length(levs)==2){
    model_res <- getTableFromConfmatrices(modlist)
  }else{
    model_res <- getTableFromConfmatrices_multiclass(modlist)
  }
  
  write_tsv(model_res, file=paste0(opt$out,"summary_all_models", name, ".tsv"))
  return(list(models=modlist, modummary=model_res, component_pvals=compsig, varnames=varnames))
}


plotSVM<-function(modelo_svm, datasc, varnames, opt, name){
  #Sacado de: https://rpubs.com/Joaquin_AR/267926
  datos <- datasc[,varnames[1:2]] 
  names(datos) <- paste("X", 1:ncol(datos), sep="")
  datos$y <- datasc$class
  rangos <- datos %>% dplyr::select(matches("^X[0-9]+")) %>% map(range)
  new_xs <- map(rangos, \(x)seq(from = x[1], to = x[2], length = 75))
  
  # Interpolación de puntos
  nuevos_puntos <- expand.grid(new_xs)
  
  # Predicción según el modelo
  predicciones <- predict(object = modelo_svm, newdata = nuevos_puntos)
  
  # Se almacenan los puntos predichos para dar color a las regiones
  color_regiones <- data.frame(nuevos_puntos, y = predicciones)
  
  # Para extraer la ecuación del hiperplano y del margen es necesario aplicar 
  # algebra lineal.
  beta <- drop(t(modelo_svm$coefs) %*% as.matrix(datos[, c("X1", "X2")])[modelo_svm$index,])
  beta0 <- modelo_svm$rho
  
  
  g1 <- ggplot() +
    # Representación de las 2 regiones empleando los puntos y coloreándolos
    # según la clase predicha por el modelo
    geom_point(data = color_regiones, aes(x = X1, y = X2, color = as.factor(y)),
               size = 0.2, alpha=0.5) +
    # Se añaden las observaciones
    geom_point(data = datos, aes(x = X1, y = X2, color = as.factor(y)),
               size = 2) +
    # Se identifican aquellas observaciones que son vectores soporte del modelo
    geom_point(data = datos[modelo_svm$index, ],
               aes(x = X1, y = X2, color = as.factor(y)),
               shape = 21, colour = "black",
               size = 2) +
    #scale_color_lancet()+
    #scale_fill_lancet() +
    theme_bw() #+theme(legend.position = "none")
    
    if(modelo_svm$kernel == 0 & length(unique(datasc$class))==2){
    # Se añaden las rectas del hiperplano y los márgenes
    g1<-g1 + geom_abline(intercept = beta0/beta[2], slope = -beta[1]/beta[2]) +
            geom_abline(intercept = (beta0 - 1)/beta[2], slope = -beta[1]/beta[2],
                linetype = "dashed") +    
            geom_abline(intercept = (beta0 + 1)/beta[2], slope = -beta[1]/beta[2],
                linetype = "dashed") 
    }
  
  ggsave(filename = paste0(opt$out, name, "_SVMplot.pdf"), plot = g1, width = 8, height = 5)
  return(g1)
}

callDoAllModelsFromALLPCAs <- function(all_pcas, name, metadata, vars2pca=c("Condition"), 
                                       variable_plim=0.01, 
                                       meta_vars = c(),
                                       nfolds = 0,
                                       xgboost_params = xgboost_params,
                                       catboost_params = catboost_params, 
                                       randomforest_params = randomforest_params,
                                       do_smote=FALSE,
                                       smote_params=list(K=5, dup_size="balance")){
  datasc <- all_pcas[[1]]$pca$x %>% 
    as.data.frame %>% 
    rownames_to_column("sample") %>% 
    dplyr::mutate(class=unlist(metadata[match(sample, metadata$sampleID), vars2pca[1]])) %>% 
    dplyr::filter(!is.na(class)) %>% 
    dplyr::mutate(class=factor(class))
  if(length(meta_vars) > 0){
    meta_filt <- metadata %>% dplyr::select(sampleID, all_of(meta_vars))
    byy <- join_by(sample == sampleID)
    datasc <- datasc %>% inner_join(meta_filt, by=byy)
    
  }
  allmodssumm <- makeAllModels(datasc, plim=variable_plim, opt, name= name, nfolds = nfolds, 
                               xgboost_params = xgboost_params,
                               catboost_params = catboost_params, 
                               randomforest_params = randomforest_params,
                               do_smote = do_smote, smote_params = smote_params, 
                               do_ensemble_probs = FALSE)
  
  modelo_svm <- allmodssumm$models$`SVM-linear`$mod_noscale
  allmodssumm$plot_svm_rad <-plotSVM(modelo_svm, datasc, allmodssumm$varnames, opt, paste0(name, "_linear"))
  
  modelo_svm <- allmodssumm$models$`SVM-radial`$mod_noscale
  allmodssumm$plot_svm_rad <- plotSVM(modelo_svm, datasc, allmodssumm$varnames, opt, paste0(name, "_radial"))
  
  return(allmodssumm)
}


callDoAllModelsFromALLPCAsOriginalVars <- function(all_pcas, PCs, modelo_svm, vstdf, 
                                                   name, vars2pca=c("Condition"), metadata, 
                                                   daares, topns = c(5, 10, 20),
                                                   variable_plim=0.01, 
                                                   meta_vars = c() ,
                                                   nfolds = 0,
                                                   xgboost_params = xgboost_params,
                                                   catboost_params = catboost_params, 
                                                   randomforest_params = randomforest_params,
                                                   do_smote=FALSE,
                                                   smote_params=list(K=5, dup_size="balance")){
  
  pcts <- summary(all_pcas[[1]]$pca)$importance[2, PCs]
  pcslope <- pcts[1]/pcts[2]
  beta <- drop(t(modelo_svm$coefs) %*% all_pcas$Condition$pca$x[modelo_svm$index,PCs])
  bslope <- -beta[1]/beta[2]
  
  rotvals <- all_pcas[[1]]$pca$rotation %>% as.data.frame %>% dplyr::select(all_of(PCs)) %>% 
    rownames_to_column("taxon") %>% 
    dplyr::mutate(
      score1 = abs(all_pcas[[1]]$pca$rotation[, PCs[1]]),
      score2 = abs(all_pcas[[1]]$pca$rotation[, PCs[2]]), 
      score3 = abs(all_pcas[[1]]$pca$rotation[, PCs[1]]) + abs(all_pcas[[1]]$pca$rotation[, PCs[2]]),
      score4 = pcslope*abs(all_pcas[[1]]$pca$rotation[, PCs[1]]) + abs(all_pcas[[1]]$pca$rotation[, PCs[2]]),
      score5 = bslope*abs(all_pcas[[1]]$pca$rotation[, PCs[1]]) + abs(all_pcas[[1]]$pca$rotation[, PCs[2]]),
      daa_pval  = -10*log10(daares$padj[match(taxon, daares$taxon)])
    )
  names(rotvals)[!names(rotvals) %in% c("taxon", PCs)] <- c(paste0(PCs[1], " score"), 
                                                            paste0(PCs[2], " score"), 
                                                            paste0(PCs, collapse="+"),
                                                            paste0(PCs[1], " and ", PCs[2], " combined 2"),
                                                            paste0(PCs[1], " and ", PCs[2], " combined "),
                                                            "DESeq pval"
  )
  modresults <- list()
  for(score in names(rotvals)[!names(rotvals) %in% c("taxon", PCs)]){
    for(topn in topns){
      print(paste0("Fitting models to ", score, '_', topn))
      toptaxa <- rotvals[order(rotvals[, score], decreasing = T), "taxon"][1:topn]
      df2pred <- vstdf %>% dplyr::filter(gene %in% toptaxa) %>% column_to_rownames("gene") %>% as.matrix %>% t %>% 
        as.data.frame() %>% rownames_to_column("sample") %>% 
        dplyr::mutate(class=unlist(metadata[match(sample, metadata$sampleID), vars2pca[1]]))
      names(df2pred) <- gsub("[\\.\\-\\[\\]()]", "", names(df2pred), perl=T)
      modresults[[paste0(score, ' top ', as.character(topn))]] <- makeAllModels(df2pred, plim=1, opt, name= paste0(name, "_modsIndBacs_", score, "_top", topn), 
                                                                                nfolds = nfolds, 
                                                                                xgboost_params = xgboost_params,
                                                                                catboost_params = catboost_params, 
                                                                                randomforest_params = randomforest_params,
                                                                                do_smote = do_smote, smote_params = smote_params)
      modresults[[paste0(score, ' top ', as.character(topn))]]$taxa <- toptaxa
      
    }
  }##make models
  print("Merging models")
  modall_table <- map(names(modresults), \(x){
    res <- modresults[[x]]$modummary %>% 
      dplyr::mutate(sel_method=x,
                    varsused=paste0(modresults[[x]]$taxa, collapse="|"))
    
  }) %>% bind_rows()
  
  return(list(fullresults=modresults, allmodsum=modall_table))
}


makeLinePlotComparingPhobjs <- function(all_model_results, opt, models_name1="padj_taxa_res", models_name2="praw_taxa_res"){
  
  all_sig_tables <- names(all_model_results) %>% map(\(name){
    a <- all_model_results[[name]][[models_name1]]$modummary %>% 
      dplyr::mutate(taxa_group="padj")
    b <- all_model_results[[name]][[models_name2]]$modummary %>% 
      dplyr::mutate(taxa_group="praw")
    rbind(a, b) %>% dplyr::mutate(input_data=name)
  }) %>% bind_rows() %>% dplyr::arrange(desc(Accuracy_l1out)) %>% 
    dplyr::select(input_data, model, taxa_group, everything())
  
  write_tsv(all_sig_tables, file = paste0(opt$out, "/all_model_summaries.tsv"))
  
  (g1 <- ggplot(all_sig_tables, aes(x=model, 
                                    y=Accuracy_l1out, 
                                    col=input_data,
                                    fill=input_data, 
                                    group=input_data))+
      facet_grid(taxa_group~.)+
      geom_point()+
      geom_line()+
      scale_color_cosmic()+
      scale_fill_cosmic() +
      ggpubr::theme_pubr() +
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
  )
  ggsave(paste0(opt$out, "/all_model_accuracy.pdf"), g1, 
         width = 8, height = 8)
  return(g1)
}

makeLinePlotComparingSamePhobjModels<- function(phname, all_model_results, opt,
                                                w=8, h=12, get_pcnames_from="padj_taxa_res", plot_extra=FALSE, 
                                                plot_indiv = TRUE, plot_normal_with_smote=FALSE,
                                                filter_out = c("Ensemble2"),
                                                order_by_measure = "Accuracy_l1out", from_smote=FALSE, name="", 
                                                sel_method_name="PCA DESeq"){
  outdir <- paste0(opt$out, "/", phname, "_", name, "/")
  if(!dir.exists(outdir)) dir.create(outdir)
  
  indivname <- "padj_taxa_res_indiv"
  lindaname <- "padj_taxa_res05linda"
  if(from_smote){
    #get_pcnames_from <- paste0(get_pcnames_from, "_SMOTE")
    indivname <- paste0(indivname, "_SMOTE")
    lindaname <- paste0(lindaname, "_SMOTE")
  }
  
  resph <- all_model_results[[phname]]
  pcnames <- resph[[get_pcnames_from]]$varnames
  tabs <- resph[[get_pcnames_from]]$modummary %>% dplyr::mutate(sel_method = sel_method_name, varsused = paste(pcnames, collapse="|"))
  if(lindaname %in% names(resph) & plot_extra){
    linda_pcnames <- resph[[lindaname]]$varnames
    aux <- resph[[lindaname]]$modummary %>% dplyr::mutate(sel_method = "PCA LinDA", varsused = paste(linda_pcnames, collapse="|"))
    tabs <- rbind(
      tabs, 
      aux
    )
  }
  if(indivname %in% names(resph) & plot_indiv){
    tabs <- rbind(
      tabs,
      resph[[indivname]]$allmodsum
    )
  }
  if(plot_normal_with_smote){
    get_pcnames_from_smote <- paste0(get_pcnames_from, "_SMOTE")
    aux <- resph[[get_pcnames_from_smote]]$modummary %>% dplyr::mutate(sel_method = paste0(sel_method_name, " SMOTE"), varsused = paste(pcnames, collapse="|"))
    tabs <- rbind(
      tabs, 
      aux
    )
  }
  
  
  write_tsv(tabs, file = paste0(outdir, phname, "_modelSummariesWithIndividualSpecies.tsv"))
  tabs2plot <- tabs %>% 
    dplyr::filter(!grepl("\\+", sel_method)) %>% 
    dplyr::filter(!grepl("combined 2", sel_method)) %>% 
    dplyr::filter(!grepl("PC11 score", sel_method)) %>% 
    dplyr::filter(!model %in% filter_out) %>% 
    dplyr::mutate(model=gsub("logistic_regression", "Logistic Regr.", model),
                sel_method = gsub("  ", " ", sel_method),
                sel_method = gsub("DESeq", "DESeq2", sel_method)
    )
  #modorder <- tabs2plot %>% group_by(model) %>% dplyr::summarise(maxacc = max(!!sym(order_by_measure))) %>% 
  #  dplyr::arrange(desc(maxacc)) %>% pull(model)
  
  levorder <- tabs2plot %>% group_by(sel_method) %>% dplyr::summarise(maxacc = max(!!sym(order_by_measure))) %>% 
    dplyr::arrange(desc(maxacc)) %>% pull(sel_method)
  
  modorder <- tabs2plot %>% filter(sel_method == levorder[1]) %>% 
    dplyr::arrange(desc(!!sym(order_by_measure))) %>% pull(model)

  tabs2plot <- tabs2plot %>% dplyr::mutate(model = factor(model, levels=modorder),
                                         sel_method = factor(sel_method, levels=levorder)) 
  (g1 <- ggplot(tabs2plot, aes(x=model, 
                             y=Accuracy_l1out, 
                             col=sel_method,
                             fill=sel_method, 
                             group=sel_method
                             #shape=sel_method,
                             #linetype=sel_method)
  ))+
    #facet_grid(sel_method~.)+
    geom_col(alpha=1, width = 0.8, position="dodge")+
    #geom_point(size=2)+
    #geom_line(alpha=1)+
    scale_color_cosmic()+
    scale_fill_cosmic() +
    ggpubr::theme_pubr() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
  )
  ggsave(paste0(outdir, "/all_model_accuracy_1.pdf"), g1, 
       width = w, height = w*0.5)
  tabs2plot2 <- tabs2plot %>%  dplyr::filter(!grepl("top (5|10)", sel_method, perl=T)) 
  g2 <- ggplot(tabs2plot2, aes(x=model, 
                             y=Accuracy_l1out, 
                             col=sel_method,
                             fill=sel_method, 
                             group=sel_method
                             #shape=sel_method,
                             #linetype=sel_method)
  ))+
  #facet_grid(sel_method~.)+
  #geom_col(alpha=1, width = 0.8, position="dodge")+
  geom_point(size=2)+
  geom_line(alpha=1)+
  #scale_color_cosmic()+
  #scale_fill_cosmic() +
  ggpubr::theme_pubr() +
  ylab("Accuracy") + 
  xlab("")+ 
  theme(legend.position="right")+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

  ggsave(paste0(outdir, "/", phname, "all_model_accuracy_2.pdf"), g2, 
       width = 12, height = 8)
  g3 <- ggplot(tabs2plot2, aes(x=model, 
                             y=Sensitivity_l1out, 
                             col=sel_method,
                             fill=sel_method, 
                             group=sel_method
                             #shape=sel_method,
                             #linetype=sel_method)
  ))+
  #facet_grid(sel_method~.)+
  #geom_col(alpha=1, width = 0.8, position="dodge")+
  geom_point(size=2)+
  geom_line(alpha=1)+
  #scale_color_cosmic()+
  #scale_fill_cosmic() +
  ggpubr::theme_pubr() +
  ylab("Sensitivity")+ 
  xlab("")+ 
  theme(legend.position="right")+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

  ggsave(paste0(outdir, "/", phname, "all_model_sensitivity.pdf"), g3, 
       width = 12, height = 8)

  g4 <- ggplot(tabs2plot2, aes(x=model, 
                             y=Specificity_l1out, 
                             col=sel_method,
                             fill=sel_method, 
                             group=sel_method
                             #shape=sel_method,
                             #linetype=sel_method)
  ))+
  #facet_grid(sel_method~.)+
  #geom_col(alpha=1, width = 0.8, position="dodge")+
  geom_point(size=2)+
  geom_line(alpha=1)+
  #scale_color_cosmic()+
  #scale_fill_cosmic() +
  ggpubr::theme_pubr() +
  ylab("Specificity")+ 
  xlab("") + 
  theme(legend.position="right") + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

  ggsave(paste0(outdir, "/", phname, "all_model_Specificity.pdf"), g4, 
       width = 12, height = 8)

  g5 <- ggplot(tabs2plot2, aes(x=model, 
                             y=Kappa_l1out, 
                             col=sel_method,
                             fill=sel_method, 
                             group=sel_method
                             #shape=sel_method,
                             #linetype=sel_method)
  ))+
  #facet_grid(sel_method~.)+
  #geom_col(alpha=1, width = 0.8, position="dodge")+
  geom_point(size=2)+
  geom_line(alpha=1)+
  #scale_color_cosmic()+
  #scale_fill_cosmic() +
  ggpubr::theme_pubr() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  ylab("Kappa")+ 
  xlab("") + 
  scale_y_continuous(n.breaks = 6)+
  theme(legend.position="right") + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

  ggsave(paste0(outdir, "/", phname, "all_model_Kappa.pdf"), g5, 
       width = 6, height = 3)

  
  g6 <- ggplot(tabs2plot2, aes(x=model, 
                               y=BalancedAccuracy_l1out, 
                               col=sel_method,
                               fill=sel_method, 
                               group=sel_method
                               #shape=sel_method,
                               #linetype=sel_method)
  ))+
    #facet_grid(sel_method~.)+
    #geom_col(alpha=1, width = 0.8, position="dodge")+
    geom_point(size=2)+
    geom_line(alpha=1)+
    #scale_color_cosmic()+
    #scale_fill_cosmic() +
    ggpubr::theme_pubr() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
    ylab("Balanced Accuracy")+ 
    xlab("") + 
    scale_y_continuous(n.breaks = 6)+
    theme(legend.position="right") + 
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
  
  ggsave(paste0(outdir, "/", phname, "all_model_BalancedAccuracy.pdf"), g6, 
         width = 6, height = 3)
  
  
  g7 <- ggplot(tabs2plot2, aes(x=model, 
                               y=AUC_l1out, 
                               col=sel_method,
                               fill=sel_method, 
                               group=sel_method
                               #shape=sel_method,
                               #linetype=sel_method)
  ))+
    #facet_grid(sel_method~.)+
    #geom_col(alpha=1, width = 0.8, position="dodge")+
    geom_point(size=2)+
    geom_line(alpha=1)+
    #scale_color_cosmic()+
    #scale_fill_cosmic() +
    ggpubr::theme_pubr() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
    ylab("AUC")+ 
    xlab("") + 
    scale_y_continuous(n.breaks = 6)+
    theme(legend.position="right") + 
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
  
  ggsave(paste0(outdir, "/", phname, "all_model_AUC.pdf"), g7, 
         width = 6, height = 3)
  
  cw <- cowplot::plot_grid(plotlist=list(g2, g3, g4), ncol = 1)
  pdf(paste0(outdir, "/", phname, ifelse(plot_extra, "_LinDA_", "") , "_all_model_combined.pdf"), width = w, height = h)
  print(cw)
  dev.off()

  cw <- cowplot::plot_grid(plotlist=list(g2, g5), ncol = 1)
  pdf(paste0(outdir, "/", phname, "_all_model_combined2.pdf"), width = w, height = w)
  print(cw)
  dev.off()

  cw <- cowplot::plot_grid(plotlist=list(g2, g5, g3, g4), ncol = 1)
  pdf(paste0(outdir, "/", phname,ifelse(plot_extra, "_LinDA_", "") , "_all_model_combined3.pdf"), width = w, height = w*1.7)
  print(cw)
  dev.off()
  
  cw <- cowplot::plot_grid(plotlist=list(g7, g5, g3, g4), ncol = 1)
  pdf(paste0(outdir, "/", phname,ifelse(plot_extra, "_LinDA_", "") , "_all_model_combined5.pdf"), width = w, height = w*1.7)
  print(cw)
  dev.off()
  
  cw <- cowplot::plot_grid(plotlist=list(g6, g5, g3, g4), ncol = 1)
  pdf(paste0(outdir, "/", phname,ifelse(plot_extra, "_LinDA_", "") , "_all_model_combined6.pdf"), width = w, height = w*1.7)
  print(cw)
  dev.off()
  
  cw <- cowplot::plot_grid(plotlist=list(g5 + theme(legend.position = "none"), 
                                         g6 + theme(legend.position = "none"), 
                                         g3 + theme(legend.position = "none"), 
                                         g4 + theme(legend.position = "none"), 
                                         g7 + theme(legend.position = "none")), ncol = 2)
  pdf(paste0(outdir, "/", phname,ifelse(plot_extra, "_LinDA_", "") , "_all_model_combined4.pdf"), width = w*1.7, height = w*1.5)
  print(cw)
  dev.off()

}

makeLinePlotComparingSamePhobjModels_Cov<- function(phname, condnames, 
                                                    all_model_results, 
                                                    name, opt, w=8, h=12){
  outdir <- paste0(opt$out, phname)
  if(!file.exists(outdir)) dir.create(outdir)
  all_sig_tables <- map(condnames, \(name){
    all_model_results[[phname]][[name]]$modummary %>% 
      dplyr::mutate(taxa_group=name, input_data=name)
  }) %>% bind_rows() %>% dplyr::arrange(desc(Accuracy_l1out)) %>% 
    dplyr::select(input_data, model, taxa_group, everything())
  
  model_levels <- all_sig_tables %>% group_by(model) %>% 
    dplyr::summarise(acc=max(Accuracy_l1out)) %>% 
    dplyr::arrange(desc(acc)) %>% pull(model)
  all_sig_tables <- all_sig_tables %>% dplyr::mutate(model = factor(model, levels=model_levels))
  write_tsv(all_sig_tables, file = paste0(outdir, "/all_model_summaries.tsv"))
  
  (g1 <- ggplot(all_sig_tables, aes(x=model, 
                                    y=Accuracy_l1out, 
                                    col=input_data,
                                    fill=input_data, 
                                    group=input_data))+
      #facet_grid(taxa_group~.)+
      geom_point()+
      geom_line()+
      scale_color_cosmic()+
      scale_fill_cosmic() +
      ggpubr::theme_pubr() +
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
  )
  ggsave(paste0(outdir, "/all_model_accuracy_", name, ".pdf"), g1, 
         width = w, height = h)
  return(g1)
  
}


getPCnamesFromAllresults <- function(phname, all_model_results, 
                                     get_pcnames_from="padj_taxa_res", 
                                     pca_name="padj_taxa_pcas",
                                     varname="Condition"){
  PCs <- all_model_results[[phname]][[get_pcnames_from]]$varnames
  sumPCA <- summary(all_model_results[[phname]][[pca_name]][[varname]]$pca)$importance
  PC_names <- paste(PCs, ' (', round(100*sumPCA[2, PCs], 1),'%)', sep="")
  PCs_newnames <- PCs
  names(PCs_newnames) <- PC_names
  return(PCs_newnames)
}

makePCsBoxplot <- function(phname, all_model_results, opt,
                           get_pcnames_from="padj_taxa_res", 
                           pca_name="padj_taxa_pcas", 
                           varname="Condition", w=4, h=6){
  outdir <- paste0(opt$out, "/", phname)
  if(!dir.exists(outdir)) dir.create(outdir)

  PCs <- all_model_results[[phname]][[get_pcnames_from]]$varnames 
  pc_order<- gsub("PC", "", PCs) %>% as.numeric %>% order
  PCs <- PCs[pc_order]
  PCs_newnames <- getPCnamesFromAllresults(phname, all_model_results, get_pcnames_from, pca_name, varname)
  PCs_newnames <- PCs_newnames[pc_order]
  metadata <- all_model_results[[phname]]$metadata
  dfsamples <- all_model_results$remove_tanda2[[pca_name]][[varname]]$pca$x %>% 
    as.data.frame() %>% 
    dplyr::select(all_of(PCs)) %>% 
    rownames_to_column("sampleID") %>% 
    dplyr::mutate(Condition = metadata[[varname]][match(sampleID, metadata$sampleID)]) %>% 
    tidyr::gather("PC", "score", -sampleID, -Condition) %>% 
    group_by(PC) 
  
  # Multiplicar por -1 si el componente es menor en deprimidos, para plotear más fácil
  
    pcfactors <- dfsamples %>% group_by(PC, Condition) %>% dplyr::summarise(media = mean(score)) %>% 
      tidyr::spread(Condition, media) 
    if(all(as.character(unique(dfsamples$Condition)) %in% c("Control", "Depression"))){
      pcfactors <- pcfactors %>% dplyr::mutate(factor = ifelse(Depression < Control, -1, 1))
    }else{
      pcfactors <- pcfactors %>% dplyr::mutate(factor = 1)
    }
  
  dfsamples <- dfsamples %>% 
    dplyr::mutate(score = score*pcfactors$factor[match(PC, pcfactors$PC)]) %>% 
    dplyr::mutate(PC = factor(PC, levels = PCs),
           PC = fct_recode(PC, !!!PCs_newnames),
           Condition = gsub("Depression", "Depr.", Condition))
  
  comp <- combn(unique(dfsamples$Condition), 2, simplify = F)
  signif_levels=c("***"=0.001, "**"=0.01, "*"=0.05, "ns"=1.1)
  
  gpcbox <- ggplot(dfsamples, aes(x=Condition, y=score)) +
    facet_grid(~PC)+
    geom_violin(aes(fill=Condition)) +
    geom_boxplot(width=0.2)+
    ggsignif::stat_signif(test="t.test", na.rm=T, comparisons = comp, 
                          step_increase=0.03,
                          tip_length = 0.01,
                          map_signif_level=signif_levels,
                          vjust=0.4,
                          color = "black"
    )+
    mytheme +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
  
  ggsave(filename = paste0(outdir, "/", get_pcnames_from, "_boxplot_signif_PCs.pdf"), gpcbox, width = w, height = h)
  write_tsv(dfsamples, file = paste0(outdir, "/", get_pcnames_from, "_boxplot_signif_PCs_data.tsv"))
  write_tsv(pcfactors, file = paste0(outdir, "/", get_pcnames_from, "_boxplot_signif_PCs_PCFactors.tsv"))
  return(list(plot=gpcbox, tab=dfsamples, pcfactors=pcfactors))
}

makePCBarplot <- function(phname, all_model_results, pcBoxplots, daa_all, opt,
                          get_pcnames_from="padj_taxa_res", 
                          pca_name="padj_taxa_pcas", 
                          varname="Condition", w=8, h=14){
  outdir <- paste0(opt$out, "/", phname)
  if(!dir.exists(outdir)) dir.create(outdir)
  predictions <- all_model_results[[phname]][[get_pcnames_from]]$models$`KNN-K=5`$preds_no_l1o
  PCs <- all_model_results[[phname]][[get_pcnames_from]]$varnames
  pc_order<- gsub("PC", "", PCs) %>% as.numeric %>% order
  PCs <- PCs[pc_order]
  PCs_newnames <- getPCnamesFromAllresults(phname, all_model_results, get_pcnames_from, pca_name, varname)
  PCs_newnames <- PCs_newnames[pc_order]
  daatab <- daa_all[[phname]]$resdf
  pcfactors <- pcBoxplots[[phname]]$pcfactors %>% column_to_rownames("PC") %>% dplyr::select(factor) 
  
  df <- all_model_results$remove_tanda2[[pca_name]][[varname]]$pca$rotation %>% 
    as.data.frame() %>% 
    dplyr::select(all_of(PCs)) %>% 
    rownames_to_column("taxon") %>% 
    dplyr::mutate(across(all_of(PCs), \(x)x*pcfactors[cur_column(), 1])) ## Multiplicar por factor para que coincida con LFC
  
  dfmerged <- merge(df, daatab, by="taxon", all.x=T, all.y=F) %>% 
    dplyr::arrange(desc(!!sym(PCs[1]))) %>% 
    dplyr::mutate(taxon = gsub("_", " ", taxon),
                  taxon = gsub("[\\[\\]]", "", taxon),
                  taxon = factor(taxon, levels = taxon))
  
  newlevnames <- c( PCs_newnames, "LFC"="log2FoldChangeShrink","-10log(adj. p)"="padj")
  dflong <- dfmerged %>% 
    dplyr::select(all_of(c("taxon", "padj", PCs, "log2FoldChangeShrink"))) %>% 
    dplyr::mutate(padj = -10*log10(padj)) %>% 
    tidyr::gather(key="variable", "value", -taxon) %>% 
    dplyr::mutate(color = ifelse(variable == "padj",C_NS, ifelse(value < 0, C_CTRL, C_CASE)),
                  variable = fct_recode(variable, !!!newlevnames),
                  variable = factor(variable, levels = names(newlevnames))
    )
  
  gbars <- ggplot(dflong, aes(y=value, x=taxon, fill=color)) +
    facet_wrap(~ variable, nrow=1, scales = "free_x")+
    geom_col()+
    coord_flip() +
    theme_classic() +
    scale_fill_manual(values = c(C_CTRL, C_CTRL_LINK2, C_CASE))+
    theme(axis.text.y = element_text(size = 8, 
                                     colour = "black", angle = 0, 
                                     face = "italic"))+
    theme(axis.text.x = element_text(size = 10, 
                                     colour = "black", angle = 0, 
                                     face = "plain"))+
    theme(strip.text.x = element_text(size = 14, 
                                      colour = "black", angle = 0, face = "plain")) +
    thin_barplot_lines +
    theme(legend.position="none")
  ggsave(filename = paste0(outdir, "/", get_pcnames_from, "_barplots_PCs_and_LFC.pdf"), gbars, width = w, height = h)
  write_tsv(dflong,  paste0(outdir, "/", get_pcnames_from, "_barplots_PCs_and_LFC_data.tsv"))
  return(gbars)
}

plotPrediction<-function(phname, mod2plot, all_model_results, opt,
                         get_pcnames_from="padj_taxa_res", 
                         pca_name="padj_taxa_pcas", 
                         varname="Condition", pred_mode="l1o", w=6, h=4){
  outdir <- paste0(opt$out, "/", phname)
  if(!dir.exists(outdir)) dir.create(outdir)
  if(pred_mode=="l1o"){
    predictions <- all_model_results[[phname]][[get_pcnames_from]]$models[[mod2plot]]$preds
  }else{
    predictions <- all_model_results[[phname]][[get_pcnames_from]]$models[[mod2plot]]$preds_no_l1o
  }
  snames <- all_model_results[[phname]][[pca_name]][[varname]]$pca$x %>% rownames
  PCs <- all_model_results[[phname]][[get_pcnames_from]]$varnames 
  pc_order <- gsub("PC", "", PCs) %>% as.numeric %>% order
  PCs <- PCs[pc_order]
  PCs_newnames <- getPCnamesFromAllresults(phname, all_model_results, get_pcnames_from, pca_name, varname) %>% names %>% sort
  PCs_newnames <- PCs_newnames[pc_order]
  
  df <- pcBoxplots[[phname]]$tab %>% 
    dplyr::mutate(Predicted = predictions[match(sampleID, snames)]) %>% 
    spread(key=PC, value=score)
  if(all(unique(df$Condition) %in% c("Control", "Depression", "Depr."))){
    df <- df %>% dplyr::mutate(Condition = fct_recode(Condition, "Control"="Control", "Depression"="Depr."),
                  Good = ifelse(Condition==Predicted, TRUE, FALSE),
                  size2 = ifelse(Good, 0, 1)) 
  }else{
    df <- df %>% dplyr::mutate(Predicted = gsub("Depression", "D", as.character(Predicted)),
                               Predicted = gsub("Control", "C", as.character(Predicted)),
                               Good = ifelse(Condition==Predicted, TRUE, FALSE),
                               size2 = ifelse(Good, 0, 1)) 
  }
  confmat <- caret::confusionMatrix(factor(df$Condition), factor(df$Predicted))
  
  gm <- ggplot(df, aes(x=!!sym(PCs_newnames[1]), y =!!sym(PCs_newnames[2]), col=Condition))+
    geom_point(size=2)+
    geom_point(col="black", alpha=df$size2, size=0.6)+
    mytheme +
    ggtitle(paste0(mod2plot, ", Acc=", as.character(round(confmat$overall["Accuracy"], 3))))
  ggsave(filename = paste0(outdir, "/", mod2plot, "_predictionPlot.pdf"), width = w, height = h)
  gm <- gm + theme(legend.position = 'none')
  return(gm)
}

plotAllModelPredictions <- function(phname, all_model_results, opt,
                                    get_pcnames_from="padj_taxa_res", 
                                    pca_name="padj_taxa_pcas", 
                                    varname="Condition", w=16, h=16,
                                    pred_mode="l1o"){
  outdir <- paste0(opt$out, "/", phname)
  if(!file.exists(outdir)) dir.create(outdir)
  modelnames <- names(all_model_results[[phname]][[get_pcnames_from]]$models)
  modnames_order <- map_vec(modelnames, \(x)all_model_results[[phname]][[get_pcnames_from]]$models[[x]]$confmat$overall["Accuracy"]) %>% order(decreasing = T)
  modelnames <- modelnames[modnames_order]
  modplots <- map(modelnames, 
                  \(modpred){
                    tryCatch(plotPrediction(phname, modpred, all_model_results, opt, 
                                            get_pcnames_from, pca_name, varname, pred_mode), error=\(x)list())
                  })
  names(modplots) <- names(all_model_results[[phname]][[get_pcnames_from]]$models)
  modplots2 <- list()
  for(mn in names(modplots)){
    m <- modplots[[mn]]
    if(class(m)[1] == "list" & length(m)==0) next
    modplots2[[mn]] <- m
  }
  if(length(modplots2) > 1){
  cw <- cowplot::plot_grid(plotlist = modplots2)
  pdf(paste0(outdir, "/all_model_predictions.pdf"), width = w, height = h)
  print(cw)
  dev.off()
  }else{
    cat("FAILED prediction plots for", phname)
  }
  return(list(cw=cw, plots=modplots))
}

