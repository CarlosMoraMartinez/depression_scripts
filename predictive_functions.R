library(tidyverse)
library(caret)
library(pROC)

makeKmeans <- function(datasc, levs, varnames, SEED=123){
  library(stats)
  library(caret)
  set.seed(SEED)
  train_df <- datasc %>% dplyr::select(-class, -sample)  %>% dplyr::select(all_of(varnames))
  mod_kmeans <- kmeans(train_df, centers=length(levs), iter.max = 100, nstart=100)
  predict_kmeans <-levels(datasc$class)[mod_kmeans$cluster] %>% factor(levels=levs)
  confmat_kmeans <- confusionMatrix(predict_kmeans, datasc$class, positive = levs[2])
  return(list(confmat_no_l1o=confmat_kmeans, mod=mod_kmeans, predicted=predict_kmeans))
}

makeKnn_l1o <- function(datasc, levs, varnames, different_ks=c(1, 3, 5, 7, 9, 11, 13)){
  library(class)
  #library(gmodels)
  
  results <- list()
  datasc <- datasc %>% dplyr::mutate(class=factor(class))
  train_df_all <- datasc %>% 
    dplyr::select(-class, -sample) %>% 
    dplyr::select(all_of(varnames))
  for(k in  different_ks){
    kname = paste("K=", as.character(k), sep="", collapse="")
    results[[kname]] <- list()
    preds <- c()
    for(i in 1:nrow(datasc)){
      train_df <- train_df_all[-i, ]
      test_df <- train_df_all[i, ]
      # Separar clases
      train_labels <- datasc$class[-i]
      test_labels <- datasc$class[i]
      
      preds  <- c(preds, knn(train_df, test_df, train_labels, k = k, prob = T))
      
    }
    
    results[[kname]][["preds"]] <- levs[preds] %>% factor(levels=levs)
    results[[kname]][["confmat"]] <- confusionMatrix(results[[kname]][["preds"]], 
                                                     factor(datasc$class), 
                                                     positive = levs[2])
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

makeNaiveBayes_l1o <- function(datasc, levs, varnames, SEED=123){
  library(e1071)
  set.seed(SEED)
  
  predict_bayes1 <- factor()
  df <- datasc %>% dplyr::select(-class, -sample) %>% dplyr::select(all_of(varnames))
  for(i in 1:nrow(datasc)){
    # Separar datos
    train_df <- df[-i, ]
    test_df <- df[i, ]
    
    # Separar clases
    train_labels <- datasc$class[-i]
    test_labels <- datasc$class[i]
    mod_bayes2 <- naiveBayes(train_df, train_labels, laplace = 0)
    predict_bayes1 <- c(predict_bayes1, predict(mod_bayes2, test_df))
    
  }
  confusionMatrix_bayes1 <- confusionMatrix(predict_bayes1, datasc$class, 
                                            positive = levs[2])
  modwithall <- naiveBayes(df, datasc$class, laplace = 0)
  predict_bayes2 <- predict(modwithall, df)
  confMatrix_bayes2_nol1o <- confusionMatrix(predict_bayes2, datasc$class, 
                                             positive = levs[2])
  return(list(confmat=confusionMatrix_bayes1, confmat_no_l1o=confMatrix_bayes2_nol1o,
              preds=predict_bayes1, preds_no_l1o=predict_bayes2, mod=modwithall))
}

make_classifTree_l1o <- function(datasc, levs, varnames){
  library(C50)
  predict_tree1 <- factor()
  df <- datasc %>% dplyr::select(-class, -sample)  %>% dplyr::select(all_of(varnames))
  for(i in 1:nrow(datasc)){
    # Separar datos
    train_df <- df[-i, ]
    test_df <- df[i, ]
    
    # Separar clases
    train_labels <- datasc$class[-i]
    test_labels <- datasc$class[i]
    mod_tree1 <- C5.0(train_df, train_labels, trials = 20)
    predict_tree1 <- c(predict_tree1, predict(mod_tree1, test_df))
    
  }
  
  confmat_tree1 <- confusionMatrix(predict_tree1, datasc$class, positive = levs[2])
  mod_all <- C5.0(df, datasc$class, trials = 20)
  predict_tree2 <- predict(mod_all, df)
  confmat_tree2 <- confusionMatrix(predict_tree2, datasc$class, positive = levs[2])
  return(list(confmat=confmat_tree1, 
              mod=mod_all, 
              preds=predict_tree1, 
              preds_no_l1o=predict_tree2,
              confmat_no_l1o=confmat_tree2,
              roc_obj_no_l1o=NULL,
              roc_auc_no_l1o=NULL,
              roc_obj=NULL,
              roc_auc=NULL
              ))
}

make_randomForest_l1o <- function(datasc, levs, varnames){
  library(randomForest)
  df <- datasc %>% dplyr::select(-class, -sample)  %>% dplyr::select(all_of(varnames))
  predict_tree1 <- factor()
  for(i in 1:nrow(datasc)){
    # Separar datos
    train_df <- df[-i, ]
    test_df <- df[i, ]
    # Separar clases
    train_labels <- datasc$class[-i]
    test_labels <- datasc$class[i]
    mod_tree1 <- randomForest(x=train_df, y=train_labels)
    predict_tree1 <- c(predict_tree1, predict(mod_tree1, test_df))
    
  }
  confmat_tree1 <- confusionMatrix(predict_tree1, datasc$class, positive = levs[2])
  mod_tree1 <- randomForest(x=df, y=datasc$class)
  predict_tree2 <- predict(mod_tree1, df)
  confmat_tree2 <- confusionMatrix(predict_tree2, datasc$class, positive = levs[2])
  return(list(confmat=confmat_tree1, 
              confmat_no_l1o=confmat_tree2,
              mod=mod_tree1, 
              preds=predict_tree1, 
              preds_no_l1o=predict_tree2,
              roc_obj_no_l1o=NULL,
              roc_auc_no_l1o=NULL,
              roc_obj=NULL,
              roc_auc=NULL))
}


make_svm_l1o <- function(datasc, levs, varnames, kernel="linear", SEED=123){
  library(e1071)
  datasc$class <- factor(datasc$class)
  df <- datasc %>% dplyr::select(-class, -sample)  %>% dplyr::select(all_of(varnames))
  predict1 <- factor(levels = levs)
  set.seed(SEED)
  for(i in 1:nrow(datasc)){
    # Separar datos
    train_df <- df[-i, ]
    test_df <- df[i, ]
    
    # Separar clases
    train_labels <- datasc$class[-i]
    test_labels <- datasc$class[i]
    mod <- svm(x = train_df, y = train_labels, scale=TRUE, kernel=kernel)
    predict1 <- c(predict1, predict(mod, test_df))
    
  }
  confmat1 <- confusionMatrix(predict1, datasc$class, positive = levs[2])
  mod_all <- svm(x = df, y = datasc$class, scale=TRUE, kernel=kernel)
  predict2 <- predict(mod_all, df)
  confmat2 <- confusionMatrix(predict2, datasc$class, positive = levs[2])
  
  mod_all_noscale <- svm(x = df[, varnames[1:2]], y = datasc$class, scale=FALSE, kernel=kernel)
  predict2_noscale <- predict(mod_all_noscale, df[, varnames[1:2]])
  confmat2_noscale <- confusionMatrix(predict2_noscale, datasc$class, positive = levs[2])
  
  
  predict2 <- predict(mod_all, df, type="probs")
  confmat2 <- confusionMatrix(predict2, factor(datasc$class))
  
  return(list(confmat=confmat1, 
              confmat_no_l1o=confmat2,
              mod=mod_all, preds=predict1, 
              preds_no_l1o=predict2,
              mod_noscale=mod_all_noscale, 
              preds_noscale=predict2_noscale, 
              confmat_noscale=confmat2_noscale,
              roc_obj_no_l1o=NULL,
              roc_auc_no_l1o=NULL,
              roc_obj=NULL,
              roc_auc=NULL))
}

make_glm_l1o <- function(datasc, levs, varnames){
  predict_glm1 <- c()
  df <- datasc %>% dplyr::select(-sample) 
  formula <- paste0("class ~ ", paste(varnames, sep="+", collapse="+")) %>% as.formula()
  for(i in 1:nrow(datasc)){
    # Separar datos
    train_df <- df[-i, ]
    test_df <- df[i, ]
    
    mod_glm <- glm(formula, data=train_df, family = binomial)
    predict_glm1 <- c(predict_glm1, predict(mod_glm, test_df))
    
  }
  predict1 <- ifelse(predict_glm1 > 0.5, levs[2], levs[1]) %>% factor(levels=levs)
  confmat1 <- confusionMatrix(predict1, datasc$class, positive = levs[2])
  mod_all <- glm(formula, data=datasc, family = binomial)
  predict2 <- predict(mod_all, df)
  predict2 <- ifelse(predict2 > 0.5, levs[2], levs[1]) %>% factor(levels=levs)
  confmat2 <- confusionMatrix(predict2, datasc$class, positive = levs[2])
  return(list(confmat=confmat1, confmat_no_l1o=confmat2,
              mod=mod_all, preds=predict1, preds_no_l1o=predict2))
}


make_glm_l1o_multiclass <- function(datasc, levs, varnames){
  library(nnet)
  predict_glm1 <- c()
  df <- datasc %>% dplyr::select(-sample) 
  formula <- paste0("class ~ ", paste(varnames, sep="+", collapse="+")) %>% as.formula()
  for(i in 1:nrow(datasc)){
    # Separar datos
    train_df <- df[-i, ]
    test_df <- df[i, ]
    
    mod_glm <- multinom(formula, data=train_df)
    newpred <- predict(mod_glm, test_df, type="prob")
    predict_glm1 <- rbind(predict_glm1, newpred)
    
  }
  
  predict1<- colnames(predict_glm1)[apply(predict_glm1, MAR=1, \(x)which(x==max(x)))] %>% factor
  confmat1 <- confusionMatrix(predict1, factor(datasc$class))
  roc_obj <- apply(predict_glm1, MAR=2, \(x) multiclass.roc(datasc$class, x))
  roc_auc <- sapply(roc_obj, \(x)x$auc) %>% mean
  
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
              roc_obj=roc_obj,
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
  return(res %>% arrange(pval))
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
      Accuracy=if(is.null(mod$confmat_no_l1o)) NA else mod$confmat_no_l1o$overall["Accuracy"],
      Kappa=if(is.null(mod$confmat_no_l1o)) NA else mod$confmat_no_l1o$overall["Kappa"],
      Sensitivity = if(is.null(mod$confmat_no_l1o)) NA else mod$confmat_no_l1o$byClass["Sensitivity"],
      Specificity = if(is.null(mod$confmat_no_l1o)) NA else mod$confmat_no_l1o$byClass["Specificity"],
      PPV = if(is.null(mod$confmat_no_l1o)) NA else mod$confmat_no_l1o$byClass["Pos Pred Value"],
      NPV = if(is.null(mod$confmat_no_l1o)) NA else mod$confmat_no_l1o$byClass["Neg Pred Value"],
      Precision = if(is.null(mod$confmat)) NA else mod$confmat_no_l1o$byClass["Precision"],
      Recall = if(is.null(mod$confmat)) NA else mod$confmat_no_l1o$byClass["Recall"]
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
      Accuracy=if(is.null(mod$confmat_no_l1o)) NA else mod$confmat_no_l1o$overall["Accuracy"] %>% mean,
      Kappa=if(is.null(mod$confmat_no_l1o)) NA else mod$confmat_no_l1o$overall["Kappa"] %>% mean,
      Sensitivity = if(is.null(mod$confmat_no_l1o)) NA else mod$confmat_no_l1o$byClass[, "Sensitivity"] %>% mean,
      Specificity = if(is.null(mod$confmat_no_l1o)) NA else mod$confmat_no_l1o$byClass[, "Specificity"] %>% mean,
      PPV = if(is.null(mod$confmat_no_l1o)) NA else mod$confmat_no_l1o$byClass[, "Pos Pred Value"] %>% mean,
      NPV = if(is.null(mod$confmat_no_l1o)) NA else mod$confmat_no_l1o$byClass[, "Neg Pred Value"] %>% mean,
      Precision = if(is.null(mod$confmat)) NA else mod$confmat_no_l1o$byClass[, "Precision"] %>% mean,
      Recall = if(is.null(mod$confmat)) NA else mod$confmat_no_l1o$byClass[, "Recall"] %>% mean
    ) 
  }) %>% bind_rows() %>% 
    dplyr::mutate(model = names(modlist)) %>% 
    dplyr::select(model, everything()) %>% 
    arrange(desc(Accuracy_l1out))
  rownames(res) <- NULL
  return(res)
}
makeAllModels <- function(datasc, plim=0.01, opt, name="Condition"){
  levs <- datasc %>% pull(class) %>% as.factor %>% levels
  # Select features
  if(length(levs)==2){
    compsig <- get_signif_components(datasc, levs)
  }else{
    compsig <- get_signif_components_multiclass(datasc, levs)
  }
  write_tsv(compsig, file=paste0(opt$out, "significant_PCAcomponents_", name,".tsv"))
  varnames <- c(compsig$var[compsig$pval <= plim])
  if(length(varnames) < 2){
    varnames <- compsig %>% arrange(pval) %>% head(2) %>% pull(var)
  }
  
  if(length(levs)==2){
    res_glms <- make_glm_l1o(datasc, levs, varnames)
  }else{
    res_glms <- make_glm_l1o_multiclass(datasc, levs, varnames)
  }
  res_svm_lin <- make_svm_l1o(datasc, levs, varnames, kernel="linear")
  res_svm_rad <- make_svm_l1o(datasc, levs, varnames, kernel="radial")
  res_randfor <- make_randomForest_l1o(datasc, levs, varnames)
  res_tree <- make_classifTree_l1o(datasc, levs, varnames)
  res_naivebayes <- makeNaiveBayes_l1o(datasc, levs, varnames, SEED=SEED)
  res_knn_l1o <- makeKnn_l1o(datasc, levs, varnames, different_ks=seq(1,13, by=2))
  #res_knn_no_l1o <- makeKnn(datasc, levs, varnames, different_ks=seq(1,13, by=2))
  res_kmeans <- makeKmeans(datasc, levs, varnames, SEED=SEED)
  
  modlist <- list("logistic_regression" = res_glms, 
                  "SVM-linear"=res_svm_lin, 
                  "SVM-radial"=res_svm_rad,
                  "RandomForest"=res_randfor,
                  "Tree" = res_tree,
                  "NaiveBayes"=res_naivebayes,
                  "KMeans"=res_kmeans)
  if(length(levs)>2)names(modlist)[[1]] <- "Multinom"
  for(k in names(res_knn_l1o)) modlist[[paste0("KNN-", k)]] <- res_knn_l1o[[k]]
  save(modlist, file=paste0(opt$out, "all_models_", name, ".RData"))
  
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

callDoAllModelsFromALLPCAs <- function(all_pcas, name, metadata, vars2pca=c("Condition")){
  datasc <- all_pcas[[1]]$pca$x %>% 
    as.data.frame %>% 
    rownames_to_column("sample") %>% 
    dplyr::mutate(class=unlist(metadata[match(sample, metadata$sampleID), vars2pca[1]])) %>% 
    dplyr::filter(!is.na(class)) %>% 
    dplyr::mutate(class=factor(class))
  allmodssumm <- makeAllModels(datasc, plim=0.01, opt, name= name)
  
  modelo_svm <- allmodssumm$models$`SVM-linear`$mod_noscale
  allmodssumm$plot_svm_rad <-plotSVM(modelo_svm, datasc, allmodssumm$varnames, opt, paste0(name, "_linear"))
  
  modelo_svm <- allmodssumm$models$`SVM-radial`$mod_noscale
  allmodssumm$plot_svm_rad <- plotSVM(modelo_svm, datasc, allmodssumm$varnames, opt, paste0(name, "_radial"))
  
  return(allmodssumm)
}


callDoAllModelsFromALLPCAsOriginalVars <- function(all_pcas, PCs, modelo_svm, vstdf, 
                                                name, vars2pca=c("Condition"), metadata, daares, topns = c(5, 10, 20)){
  
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
      modresults[[paste0(score, ' top ', as.character(topn))]] <- makeAllModels(df2pred, plim=1, opt, name= paste0(name, "_modsIndBacs_", score, "_top", topn))
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
  }) %>% bind_rows() %>% arrange(desc(Accuracy_l1out)) %>% 
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
                                                w=8, h=12, get_pcnames_from="padj_taxa_res"){
  outdir <- paste0(opt$out, phname)
  
  resph <- all_model_results[[phname]]
  pcnames <- resph[[get_pcnames_from]]$varnames
  tabs <- resph[[get_pcnames_from]]$modummary %>% dplyr::mutate(sel_method = "PCA", varsused = paste(pcnames, collapse="|"))
  if("padj_taxa_res_indiv" %in% names(resph)){
    tabs <- rbind(
      tabs,
      resph$padj_taxa_res_indiv$allmodsum
    )
  }
  write_tsv(tabs, file = paste0(outdir, phname, "_modelSummariesWithIndividualSpecies.tsv"))
  tabs2plot <- tabs %>% 
    dplyr::filter(!grepl("\\+", sel_method)) %>% 
    dplyr::filter(!grepl("combined 2", sel_method)) %>% 
    dplyr::filter(!grepl("PC11 score", sel_method)) %>% 
    dplyr::mutate(model=gsub("logistic_regression", "Logistic Regr.", model),
                sel_method = gsub("  ", " ", sel_method)
    )
  modorder <- tabs2plot %>% group_by(model) %>% dplyr::summarise(maxacc = max(Accuracy_l1out)) %>% 
    arrange(desc(maxacc)) %>% pull(model)
  levorder <- tabs2plot %>% group_by(sel_method) %>% dplyr::summarise(maxacc = max(Accuracy_l1out)) %>% 
    arrange(desc(maxacc)) %>% pull(sel_method)

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

  cw <- cowplot::plot_grid(plotlist=list(g2, g3, g4), ncol = 1)
  pdf(paste0(outdir, "/", phname, "_all_model_combined.pdf"), width = w, height = h)
  print(cw)
  dev.off()

  cw <- cowplot::plot_grid(plotlist=list(g2, g5), ncol = 1)
  pdf(paste0(outdir, "/", phname, "_all_model_combined2.pdf"), width = w, height = w)
  print(cw)
  dev.off()

  cw <- cowplot::plot_grid(plotlist=list(g2, g5, g3, g4), ncol = 1)
  pdf(paste0(outdir, "/", phname, "_all_model_combined3.pdf"), width = w, height = w*1.7)
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
  }) %>% bind_rows() %>% arrange(desc(Accuracy_l1out)) %>% 
    dplyr::select(input_data, model, taxa_group, everything())
  
  model_levels <- all_sig_tables %>% group_by(model) %>% 
    dplyr::summarise(acc=max(Accuracy_l1out)) %>% 
    arrange(desc(acc)) %>% pull(model)
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

