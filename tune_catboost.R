
library(doParallel)
library(foreach)
library(dplyr)

outdir <- paste0("/home/carlos/Escritorio/202311_DEPRESION/ReviewJune2025/Results_rstudio/results2/PredictDAA_10Fold/remove_tanda2/tuneCatboost/")
if(!dir.exists(outdir)) dir.create(outdir)


params_grid <- expand.grid(
  learning_rate = c(0.01, 0.05, 0.1, 0.3),
  iterations = c(100, 500),
  depth = c(2, 4, 8),
  loss_function = "Logloss",     # Binary classification
  eval_metric = "AUC",           # Evaluation metric
  random_seed = 123,             # For reproducibility
  use_best_model = TRUE,         # Stop early if no improvement
  od_type = "Iter",              # Overfitting detector type
  od_wait = 20,                  # Rounds to wait before stopping
  verbose = FALSE,               # Suppress training output
  thread_count = 1,              # Number of threads
  balance_weights = TRUE         # Class balancing (e.g., c(1, 3) for binary imbalance)
) %>%
  dplyr::mutate_if(is.factor, as.character)


num_cores <- parallel::detectCores() 
#cl <- makeCluster(num_cores)
#registerDoParallel(cl)

results <- list()
#results <- foreach(i = 1:nrow(params_grid), .combine = list, .packages = c("xgboost", "dplyr", "caret", "tidyverse")) %dopar% {
 for(i in 1:nrow(params_grid)){  # 1:nrow(params_grid)
   
  cat(i, " of ", nrow(params_grid), ": ", round(100*i/nrow(params_grid), 3), "%")
   catboost_params <- as.list(params_grid[i, ])
  
  res_catboot <- make_catboost_l1o(datasc, levs, varnames, 
                                folds = folds, 
                                catboost_params = catboost_params, 
                                do_smote = FALSE)

  results[[as.character(i)]] <- res_catboot
  cat("\tAccuracy=", round(res_catboot$confmat$overall[1], 3), ", Kappa=", round(res_catboot$confmat$overall[2], 3), "\n")
}
#  stopCluster(cl)

save(results, file = paste0(outdir, "tuneCatboost_2classes_1.RData"))
tunn <-getTableFromConfmatrices(results)
write_tsv(tunn, file=paste0(outdir, "tuneCatboost_2classes_WithScale_pos_weight.tsv"))

##### iter 2



params_grid <- expand.grid(
  learning_rate = c(0.01, 0.02),
  iterations = c(50, 100, 20),
  depth = c(2, 3),
  loss_function = "Logloss",     # Binary classification
  eval_metric = "AUC",           # Evaluation metric
  random_seed = 123,             # For reproducibility
  use_best_model = TRUE,         # Stop early if no improvement
  od_type = "Iter",              # Overfitting detector type
  od_wait = 20,                  # Rounds to wait before stopping
  verbose = FALSE,               # Suppress training output
  thread_count = 1,              # Number of threads
  balance_weights = TRUE         # Class balancing (e.g., c(1, 3) for binary imbalance)
) %>%
  dplyr::mutate_if(is.factor, as.character)

results <- list()

for(i in 1:nrow(params_grid)){  # 1:nrow(params_grid)
  
  cat(i, " of ", nrow(params_grid), ": ", round(100*i/nrow(params_grid), 3), "%")
  catboost_params <- as.list(params_grid[i, ])
  
  res_catboot <- make_catboost_l1o(datasc, levs, varnames, 
                                   folds = folds, 
                                   catboost_params = catboost_params, 
                                   do_smote = FALSE)
  
  results[[as.character(i)]] <- res_catboot
  cat("\tAccuracy=", round(res_catboot$confmat$overall[1], 3), ", Kappa=", round(res_catboot$confmat$overall[2], 3), "\n")
}

save(results, file = paste0(outdir, "tuneCatboost_2classes_2.RData"))
tunn2 <-getTableFromConfmatrices(results)
write_tsv(tunn2, file=paste0(outdir, "tuneCatboost_2classes_WithScale_pos_weight_2.tsv"))

params_grid[2, ]



##### iter 3

params_grid <- expand.grid(
  learning_rate = c(0.01, 0.015, 0.02, 0.03, 0.04, 0.05),
  iterations = c(30, 50, 80, 100),
  depth = c(2),
  loss_function = "Logloss",     # Binary classification
  eval_metric = "AUC",           # Evaluation metric
  random_seed = 123,             # For reproducibility
  use_best_model = TRUE,         # Stop early if no improvement
  od_type = "Iter",              # Overfitting detector type
  od_wait = 20,                  # Rounds to wait before stopping
  verbose = FALSE,               # Suppress training output
  thread_count = 1,              # Number of threads
  balance_weights = TRUE         # Class balancing (e.g., c(1, 3) for binary imbalance)
) %>%
  dplyr::mutate_if(is.factor, as.character)

results <- list()

for(i in 1:nrow(params_grid)){  # 1:nrow(params_grid)
  
  cat(i, " of ", nrow(params_grid), ": ", round(100*i/nrow(params_grid), 3), "%")
  catboost_params <- as.list(params_grid[i, ])
  
  res_catboot <- make_catboost_l1o(datasc, levs, varnames, 
                                   folds = folds, 
                                   catboost_params = catboost_params, 
                                   do_smote = TRUE)
  
  results[[as.character(i)]] <- res_catboot
  cat("\tAccuracy=", round(res_catboot$confmat$overall[1], 3), ", Kappa=", round(res_catboot$confmat$overall[2], 3), "\n")
}

save(results, file = paste0(outdir, "tuneCatboost_2classes_3.RData"))
tunn3 <-getTableFromConfmatrices(results)
write_tsv(tunn3, file=paste0(outdir, "tuneCatboost_2classes_WithScale_pos_weight_3.tsv"))

params_grid[9, ]


##### iter 4, all vars

params_grid <- expand.grid(
  learning_rate = c(0.01, 0.015, 0.02, 0.05),
  iterations = c( 50, 100, 200, 500),
  depth = c(2, 4, 8),
  loss_function = "Logloss",     # Binary classification
  eval_metric = "AUC",           # Evaluation metric
  random_seed = 123,             # For reproducibility
  use_best_model = TRUE,         # Stop early if no improvement
  od_type = "Iter",              # Overfitting detector type
  od_wait = 20,                  # Rounds to wait before stopping
  verbose = FALSE,               # Suppress training output
  thread_count = 1,              # Number of threads
  balance_weights = TRUE         # Class balancing (e.g., c(1, 3) for binary imbalance)
) %>%
  dplyr::mutate_if(is.factor, as.character)

results <- list()

for(i in 1:nrow(params_grid)){  # 1:nrow(params_grid)
  
  cat(i, " of ", nrow(params_grid), ": ", round(100*i/nrow(params_grid), 3), "%")
  catboost_params <- as.list(params_grid[i, ])
  
  res_catboot <- make_catboost_l1o(datasc, levs, names(datasc)[grep("PC", names(datasc))], 
                                   folds = folds, 
                                   catboost_params = catboost_params, 
                                   do_smote = TRUE)
  
  results[[as.character(i)]] <- res_catboot
  cat("\tAccuracy=", round(res_catboot$confmat$overall[1], 3), ", Kappa=", round(res_catboot$confmat$overall[2], 3), "\n")
}

save(results, file = paste0(outdir, "tuneCatboost_2classes_4_allVars.RData"))
tunn4 <-getTableFromConfmatrices(results)
write_tsv(tunn4, file=paste0(outdir, "tuneCatboost_2classes_WithScale_pos_weight_4_allVars.tsv"))

params_grid[10, ]
