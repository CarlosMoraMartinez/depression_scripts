library(tidyverse)
library(pROC)
# library(doParallel)
# library(foreach)
# library(dplyr)

outdir <- paste0("/home/carlos/Escritorio/202311_DEPRESION/ReviewJune2025/Results_rstudio/results2/PredictDAA_10Fold/remove_tanda2/tuneCatBoost_4groups/")
if(!dir.exists(outdir)) dir.create(outdir)

#write_tsv(datasc, file = paste0(outdir, "remove_tanda2PCA_DiffTaxaPadjCovOrSingle05_datasc_input.tsv"))

datasc <- read_tsv(paste0(outdir, "/remove_tanda2PCA_DiffTaxaPadjCovOrSingle05_datasc_input.tsv"))
varnames <- c("PC3", "PC2", "PC1", "PC8", "PC26", "PC67", "PC103", "PC6", "PC50", "PC37", "PC11")
levs <- datasc %>% pull(class) %>% as.factor %>% levels
folds <- c()
datasc$class <- factor(datasc$class, levels=levs)

#all_pcas <- all_pcalists$PCA_DiffTaxaPadj


catboost_default <- list(
  iterations = 50,              # Number of boosting rounds
  learning_rate = 0.02,           # Step size shrinkage
  depth = 2,                     # Depth of the trees
  loss_function = "MultiClass",
  eval_metric = "MultiClass",
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
  grow_policy = "SymmetricTree",
  auto_class_weights = "Balanced"
)

params_grid <- expand.grid(
  iterations = c(50, 100, 200),              # Number of boosting rounds
  learning_rate = c(0.01, 0.02, 0.05, 0.1),           # Step size shrinkage
  depth = 2:4,                     # Depth of the trees
  loss_function = "MultiClass",
  eval_metric = "MultiClass",
  random_seed = 123,            
  use_best_model = TRUE,         # Stop early if no improvemelnt
  od_type = "Iter",              # Overfitting detector type
  od_wait = 20,                  # Rounds to wait before stopping
  verbose = FALSE,               
  thread_count = 4,              
  balance_weights = TRUE,  
  bootstrap_type = c("Bayesian","Bernoulli"),
  l2_leaf_reg = c(1, 3, 5),
  subsample = c(0.6, 0.8),  # only if bootstrap type ="Bernouilli"
  grow_policy = 	c("SymmetricTree", "Depthwise"),
  auto_class_weights = "Balanced"
) %>% 
  dplyr::filter(!(bootstrap_type == "Bayesian" & subsample != 0.6)) %>% 
  dplyr::mutate_if(is.factor, as.character)


#num_cores <- parallel::detectCores() 
#cl <- makeCluster(num_cores)
#registerDoParallel(cl)

results <- list()

res_catboost <- make_catboost_l1o(datasc, levs, varnames, 
                                catboost_params = catboost_default, 
                                folds = folds, do_smote = FALSE)

results[["0"]] <- res_catboost #baseline
cat("\tAccuracy=", round(res_catboost$confmat$overall[1], 3), ", Kappa=", round(res_catboost$confmat$overall[2], 3), ", AUC=", round(res_catboost$roc_auc, 2), "\n")

#results <- foreach(i = 1:10, .combine = list, .packages = c("xgboost", "dplyr", "caret", "tidyverse" "pRoc")) %dopar% {
 for(i in 1:nrow(params_grid)){  # 
   
  cat(i, " of ", nrow(params_grid), ": ", round(100*i/nrow(params_grid), 3), "%")
  catboost_params <- as.list(params_grid[i, ])
  
  res_catboost <- make_catboost_l1o(datasc, levs, varnames, 
                                catboost_params = catboost_params, 
                                folds = folds, do_smote = FALSE)

  results[[as.character(i)]] <- res_catboost

  cat("\tAccuracy=", round(res_catboost$confmat$overall[1], 3), ", Kappa=", round(res_catboost$confmat$overall[2], 3), ", AUC=", round(res_catboost$roc_auc, 2), "\n")
}
#  stopCluster(cl)


save(results, file = paste0(outdir, "tuneCatBoost_4classes_DiffTaxaPadjCovOrSingle05_1", ".RData"))
tunn <-getTableFromConfmatrices(results)
write_tsv(tunn, file=paste0(outdir, "tuneCatBoost_4classes_WithScale_pos_weight.tsv"))

tunn2 <- cbind(rbind(data.frame(catboost_default), params_grid), tunn %>% arrange(as.numeric(model))) %>% 
  arrange(desc(AUC_l1out))
write_tsv(tunn2, file=paste0(outdir, "tuneCatBoost_4classes_WithScale_pos_weight_full.tsv"))


tunn2 %>% arrange(desc(Kappa_l1out))

best_params <- tunn2  %>% arrange(desc(Kappa_l1out)) %>% dplyr::select(all_of(names(params_grid))) %>% head(1) %>% as.list()
best_params2 <- tunn2  %>% arrange(desc(AUC_l1out)) %>% dplyr::select(all_of(names(params_grid))) %>% head(1) %>% as.list()

tunn2  %>% filter(model %in% c("464")) %>%  #view() #, "538", "375"
  dplyr::select(all_of(names(params_grid))) %>% as.list
