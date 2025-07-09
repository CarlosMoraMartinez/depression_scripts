library(tidyverse)
library(pROC)
# library(doParallel)
# library(foreach)
# library(dplyr)

outdir <- paste0("/home/carlos/Escritorio/202311_DEPRESION/ReviewJune2025/Results_rstudio/results2/PredictDAA_10Fold/remove_tanda2/tuneRandomForest_4groups/")
if(!dir.exists(outdir)) dir.create(outdir)

#write_tsv(datasc, file = paste0(outdir, "remove_tanda2PCA_DiffTaxaPadjCovOrSingle05_datasc_input.tsv"))

datasc <- read_tsv(paste0(outdir, "/remove_tanda2PCA_DiffTaxaPadjCovOrSingle05_datasc_input.tsv"))
varnames <- c("PC3", "PC2", "PC1", "PC8", "PC26", "PC67", "PC103", "PC6", "PC50", "PC37", "PC11")
levs <- datasc %>% pull(class) %>% as.factor %>% levels
folds <- c()
datasc$class <- factor(datasc$class, levels=levs)

#all_pcas <- all_pcalists$PCA_DiffTaxaPadj

randomforest_params_mult = list(ntree = 500, 
                                mtry = 4, 
                                nodesize = 5, 
                                balance_weights = TRUE)


params_grid <- expand.grid(
  ntree = 500, 
  mtry = c(1:10), 
  nodesize = c(1, 3, 5, 7, 10), 
  balance_weights = TRUE
) %>% 
  dplyr::mutate_if(is.factor, as.character)


#num_cores <- parallel::detectCores() 
#cl <- makeCluster(num_cores)
#registerDoParallel(cl)

results <- list()

res_rf <- make_randomForest_l1o(datasc, levs, varnames, 
                                randomforest_params = randomforest_params_mult, 
                                folds = folds, do_smote = FALSE)

results[["0"]] <- res_rf #baseline
cat("\tAccuracy=", round(res_rf$confmat$overall[1], 3), ", Kappa=", round(res_rf$confmat$overall[2], 3), ", AUC=", round(res_rf$roc_auc, 2), "\n")

#results <- foreach(i = 1:10, .combine = list, .packages = c("xgboost", "dplyr", "caret", "tidyverse" "pRoc")) %dopar% {
 for(i in 1:nrow(params_grid)){  # 
   
  cat(i, " of ", nrow(params_grid), ": ", round(100*i/nrow(params_grid), 3), "%")
   rf_params <- as.list(params_grid[i, ])
  
  res_rf <- make_randomForest_l1o(datasc, levs, varnames, 
                                  randomforest_params = rf_params, 
                                folds = folds, do_smote = FALSE)

  results[[as.character(i)]] <- res_rf

  cat("\tAccuracy=", round(res_rf$confmat$overall[1], 3), ", Kappa=", round(res_rf$confmat$overall[2], 3), ", AUC=", round(res_rf$roc_auc, 2), "\n")
}
#  stopCluster(cl)


save(results, file = paste0(outdir, "tuneRandomForest_4classes_DiffTaxaPadjCovOrSingle05_1", ".RData"))
tunn <-getTableFromConfmatrices(results)
write_tsv(tunn, file=paste0(outdir, "tuneRandomForest_4classes_WithScale_pos_weight.tsv"))

tunn2 <- cbind(rbind(data.frame(randomforest_params_mult), params_grid), tunn %>% arrange(as.numeric(model))) %>% 
  arrange(desc(AUC_l1out))
write_tsv(tunn2, file=paste0(outdir, "tuneRandomForest_4classes_WithScale_pos_weight_full.tsv"))


tunn2 %>% arrange(desc(Kappa_l1out))

best_params <- tunn2  %>% arrange(desc(Kappa_l1out)) %>% dplyr::select(all_of(names(params_grid))) %>% head(1) %>% as.list()
best_params2 <- tunn2  %>% arrange(desc(AUC_l1out)) %>% dplyr::select(all_of(names(params_grid))) %>% head(1) %>% as.list()

tunn2  %>% filter(model %in% c("464")) %>%  #view() #, "538", "375"
  dplyr::select(all_of(names(params_grid))) %>% as.list


tunn2  %>% arrange(desc(Kappa_l1out)) %>% head
tunn2  %>% arrange(desc(AUC_l1out)) %>% head
