
library(doParallel)
library(foreach)
library(tidyverse)

outdir <- paste0("/home/carlos/Escritorio/202311_DEPRESION/ReviewJune2025/Results_rstudio/results2/PredictDAA_10Fold/remove_tanda2/tuneXboost_4groups/")
if(!dir.exists(outdir)) dir.create(outdir)

write_tsv(datasc, file = paste0(outdir, "remove_tanda2PCA_DiffTaxaPadjCovOrSingle05_datasc_input.tsv"))

datasc <- read_tsv(paste0(outdir, "remove_tanda2PCA_DiffTaxaPadjCovOrSingle05_datasc_input.tsv"))
varnames <- c("PC3", "PC2", "PC1", "PC8", "PC26", "PC67", "PC103", "PC6", "PC50", "PC37", "PC11")
levs <- datasc %>% pull(class) %>% as.factor %>% levels
folds <- c()
#all_pcas <- all_pcalists$PCA_DiffTaxaPadj


params_grid <- expand.grid(
  learning_rate = c(0.05, 0.1, 0.3),
  max_depth = c(2, 3, 4),
  nrounds = c(30, 50, 100),
  min_child_weight = c(1, 3),
  subsample = c(0.6, 0.8, 1),
  colsample_bytree = c(0.6, 0.8, 1),
  reg_lambda = c(1, 3, 10),
  reg_alpha = c(0, 1),
  nthread = 1,
  num_class = 4,
  objective = "multi:softprob", 
  balance_weights = TRUE
) %>% 
  dplyr::mutate(objective = as.character(objective))


#num_cores <- parallel::detectCores() 
#cl <- makeCluster(num_cores)
#registerDoParallel(cl)

results <- list()

res_xgboost <- make_xgboost_l1o(datasc, levs, varnames, 
                                xgboost_params = xgboost_params_mult, 
                                folds = folds, do_smote = FALSE)

results[["0"]] <- res_xgboost #baseline


#results <- foreach(i = 1:10, .combine = list, .packages = c("xgboost", "dplyr", "caret", "tidyverse" "pRoc")) %dopar% {
 for(i in 1:nrow(params_grid)){  # 
   
  cat(i, " of ", nrow(params_grid), ": ", round(100*i/nrow(params_grid), 3), "%")
  xgboost_params <- as.list(params_grid[i, ])
  
  res_xgboost <- make_xgboost_l1o(datasc, levs, varnames, 
                                xgboost_params = xgboost_params, 
                                folds = folds, do_smote = FALSE)

  results[[as.character(i)]] <- res_xgboost
 #if(i %% 50 == 0 ){
 #  #names(results) <- (i-50+1):i
 #  save(results, file = paste0(outdir, "tuneXboost_2classes_1_", as.character(as.integer(i/50 )), ".RData"))
 #  tunn <-getTableFromConfmatrices(results)
 #  write_tsv(tunn, file=paste0(outdir, "tuneXboost_2classes_1_", as.character(as.integer(i/50 )), ".tsv"))
 #  results <- list()
 #} 
  
  cat("\tAccuracy=", round(res_xgboost$confmat$overall[1], 3), ", Kappa=", round(res_xgboost$confmat$overall[2], 3), ", AUC=", round(res_xgboost$roc_auc, 2), "\n")
}
#  stopCluster(cl)


save(results, file = paste0(outdir, "tuneXboost_4classes_DiffTaxaPadjCovOrSingle05_1", ".RData"))
tunn <-getTableFromConfmatrices(results)
write_tsv(tunn, file=paste0(outdir, "tuneXboost_4classes_WithScale_pos_weight.tsv"))
tunn2 <- cbind(rbind(data.frame(xgboost_params_mult), params_grid), tunn %>% arrange(as.numeric(model))) %>% 
  arrange(desc(AUC_l1out))
write_tsv(tunn2, file=paste0(outdir, "tuneXboost_4classes_WithScale_pos_weight_full.tsv"))


tunn2 %>% arrange(desc(Kappa_l1out))

best_params <- tunn2 %>% dplyr::select(all_of(names(params_grid))) %>% head(1) %>% as.list()

