
library(doParallel)
library(foreach)
library(dplyr)

outdir <- paste0("/home/carlos/Escritorio/202311_DEPRESION/ReviewJune2025/Results_rstudio/results2/PredictDAA_10Fold/remove_tanda2/tuneXboost/")
if(!dir.exists(outdir)) dir.create(outdir)

params_grid <- expand.grid(
  learning_rate = c(0.01, 0.1, 0.3),
  max_depth = c(2, 4, 8),
  nrounds = seq(50, 800, by=100),
  min_child_weight = c(1, 5, 10),
  subsample = c(0.6, 0.8, 1),
  colsample_bytree = c(0.6, 0.8, 1),
  reg_lambda = seq(0, 6, by=3),
  reg_alpha = seq(0, 6, by=3),
  nthread = 1,
  objective = "binary:logistic"
) %>% dplyr::filter(
  (learning_rate < 0.05 & nrounds > 200) |
    (learning_rate > 0.05 & learning_rate <= 0.2 & nrounds >= 100 &  nrounds <= 500) |  
    (learning_rate > 0.2 & nrounds <= 150) 
) %>% 
  dplyr::mutate(objective = as.character(objective))


num_cores <- parallel::detectCores() 
#cl <- makeCluster(num_cores)
#registerDoParallel(cl)

results <- list()
#results <- foreach(i = 1:nrow(params_grid), .combine = list, .packages = c("xgboost", "dplyr", "caret", "tidyverse")) %dopar% {
 for(i in ros2get_reserva){  # 1:nrow(params_grid)
   
  cat(i, " of ", nrow(params_grid), ": ", round(100*i/nrow(params_grid), 3), "%")
  xgboost_params <- as.list(params_grid[i, ])
  xgboost_params$balance_weights=TRUE
  
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
  
  cat("\tAccuracy=", round(res_xgboost$confmat$overall[1], 3), ", Kappa=", round(res_xgboost$confmat$overall[2], 3), "\n")
}
#  stopCluster(cl)

names(results) <- (as.integer(i/50 )*50 + 1):((as.integer(i/50 )*50 ) + i%%50)
save(results, file = paste0(outdir, "tuneXboost_2classes_1_", as.character(as.integer(i/50 )+1), ".RData"))
tunn <-getTableFromConfmatrices(results)
write_tsv(tunn, file=paste0(outdir, "tuneXboost_2classes_WithScale_pos_weight.tsv"))


# read

ff <- list.files(outdir, pattern=".tsv", full.names = TRUE)

alldf <- map(ff, \(x) {
  aux <- read_tsv(x) %>% 
    filter(!is.na(Kappa_l1out))
  nn <- strsplit(basename(x), "_")[[1]]
  nn <- gsub(".tsv", "", nn[length(nn)])
  aux$set <- nn
  return(aux)
  }) %>% bind_rows()

alldf %>% nrow
params_grid %>% nrow

best100 <- alldf %>% arrange(desc(Kappa_l1out)) %>% head(3)
sets2get <- best100$set %>% unique %>% as.numeric
rows2get <- map(sets2get, \(x){
  a <- (x-1)*50 + 1
  b <- a+49
  return(a:b)
  }) %>% unlist %>% 
  sort

params_good <- params_grid[rows2get, ]

map(params_good, unique)


best <- tunn %>% arrange(desc(Kappa_l1out)) %>% filter(Kappa_l1out >= 0.725)

rows2get <- best$model %>% unique %>% as.numeric
params_good <- params_grid[rows2get, ]

map(params_good, unique)

#### based on that, repeat the process


outdir <- paste0("/home/carlos/Escritorio/202311_DEPRESION/ReviewJune2025/Results_rstudio/results2/PredictDAA_10Fold/remove_tanda2/tuneXboost_round2/")
if(!dir.exists(outdir)) dir.create(outdir)


params_grid <- expand.grid(
  learning_rate = c(0.01, 0.1, 0.3),
  max_depth = c(2, 3, 4),
  nrounds = seq(50, 450, by=50),
  min_child_weight = c(1, 2),
  subsample = c(0.5, 0.6),
  colsample_bytree = c(0.8,1),
  reg_lambda = seq(2,4, by=1),
  reg_alpha = 0,
  nthread = 1,
  objective = "binary:logistic"
) %>% dplyr::filter(
  (learning_rate < 0.05 & nrounds > 200) |
    (learning_rate > 0.05 & learning_rate <= 0.2 & nrounds >= 100 &  nrounds <= 500) |  
    (learning_rate > 0.2 & nrounds <= 150) 
) %>% 
  dplyr::mutate(objective = as.character(objective))

params_grid %>% nrow
num_cores <- parallel::detectCores() 
#cl <- makeCluster(num_cores)
#registerDoParallel(cl)

results <- list()
#results <- foreach(i = 1:nrow(params_grid), .combine = list, .packages = c("xgboost", "dplyr", "caret", "tidyverse")) %dopar% {
for(i in 1:nrow(params_grid)){  # 1:nrow(params_grid)
  
  cat(i, " of ", nrow(params_grid), ": ", round(100*i/nrow(params_grid), 3), "%")
  xgboost_params <- as.list(params_grid[i, ])
  xgboost_params$balance_weights=TRUE
  
  res_xgboost <- make_xgboost_l1o(datasc, levs, varnames, 
                                  xgboost_params = xgboost_params, 
                                  folds = folds, do_smote = FALSE)
  
  results[[as.character(i)]] <- res_xgboost
  if(i %% 50 == 0 ){
    #names(results) <- (i-50+1):i
    save(results, file = paste0(outdir, "tuneXboost_2classes_1_", as.character(as.integer(i/50 )), ".RData"))
    tunn <-getTableFromConfmatrices(results)
    write_tsv(tunn, file=paste0(outdir, "tuneXboost_2classes_1_", as.character(as.integer(i/50 )), ".tsv"))
    results <- list()
  } 
  
  cat("\tAccuracy=", round(res_xgboost$confmat$overall[1], 3), ", Kappa=", round(res_xgboost$confmat$overall[2], 3), "\n")
}

save(results, file = paste0(outdir, "tuneXboost_2classes_1_", as.character(as.integer(i/50 )+1), ".RData"))
tunn <-getTableFromConfmatrices(results)
write_tsv(tunn, file=paste0(outdir, "tuneXboost_2classes_1_", as.character(as.integer(i/50 )+1), ".tsv"))

ff <- list.files(outdir, pattern=".tsv", full.names = TRUE)

alldf <- map(ff, \(x) {
  aux <- read_tsv(x) 
  nn <- strsplit(basename(x), "_")[[1]]
  nn <- gsub(".tsv", "", nn[length(nn)])
  aux$set <- nn
  return(aux)
}) %>% bind_rows()

alldf %>% nrow
params_grid %>% nrow

allmerged <- cbind(alldf, params_grid[as.numeric(alldf$model) , ]) %>% 
  arrange(desc(Kappa_l1out))
write_tsv(allmerged, file=paste0(outdir, "tuneXboost_2classes_1_All", ".tsv"))

allmerged %>% head(1)



#model Accuracy_l1out Kappa_l1out Sensitivity_l1out Specificity_l1out PPV_l1out NPV_l1out Precision_l1out Recall_l1out  Accuracy     Kappa Sensitivity Specificity
#1   673      0.8761905   0.7467063         0.8837209         0.8709677  0.826087 0.9152542        0.826087    0.8837209 0.9238095 0.8435754   0.9302326   0.9193548
#PPV  NPV Precision    Recall set learning_rate max_depth nrounds min_child_weight subsample colsample_bytree reg_lambda reg_alpha nthread       objective
#1 0.8888889 0.95 0.8888889 0.9302326  14           0.3         2      50                1       0.6                1          3         0       1 binary:logistic

xgboost_params <- as.list(params_grid[673, ])

res_xgboost1 <- make_xgboost_l1o(datasc, levs, varnames, 
                                xgboost_params = xgboost_params, 
                                folds = folds, do_smote = FALSE)

res_xgboost2 <- make_xgboost_l1o(datasc, levs, names(datasc)[grep("^PC", names(datasc), perl=T)], 
                                 xgboost_params = xgboost_params, 
                                 folds = folds, do_smote = FALSE)


# test with all vars

outdir <- paste0("/home/carlos/Escritorio/202311_DEPRESION/ReviewJune2025/Results_rstudio/results2/PredictDAA_10Fold/remove_tanda2/tuneXboost_round2_allPCs/")
if(!dir.exists(outdir)) dir.create(outdir)

allvarnames <- names(datasc)[grep("^PC", names(datasc), perl=T)]

params_grid <- expand.grid(
  learning_rate = c(0.01, 0.1, 0.3),
  max_depth = c(2, 4, 6),
  nrounds = seq(50, 450, by=50),
  min_child_weight = c(1, 2),
  subsample = c(0.5, 0.6),
  colsample_bytree = c(0.8,1),
  reg_lambda = seq(2,4, by=1),
  reg_alpha = 0,
  nthread = 1,
  objective = "binary:logistic"
) %>% dplyr::filter(
  (learning_rate < 0.05 & nrounds > 200) |
    (learning_rate > 0.05 & learning_rate <= 0.2 & nrounds >= 100 &  nrounds <= 500) |  
    (learning_rate > 0.2 & nrounds <= 150) 
) %>% 
  dplyr::mutate(objective = as.character(objective))

params_grid %>% nrow
num_cores <- parallel::detectCores() 
#cl <- makeCluster(num_cores)
#registerDoParallel(cl)

results <- list()
#results <- foreach(i = 1:nrow(params_grid), .combine = list, .packages = c("xgboost", "dplyr", "caret", "tidyverse")) %dopar% {
for(i in 1:nrow(params_grid)){  # 1:nrow(params_grid)
  
  cat(i, " of ", nrow(params_grid), ": ", round(100*i/nrow(params_grid), 3), "%")
  xgboost_params <- as.list(params_grid[i, ])
  xgboost_params$balance_weights=TRUE
  
  res_xgboost <- make_xgboost_l1o(datasc, levs, allvarnames, 
                                  xgboost_params = xgboost_params, 
                                  folds = folds, do_smote = FALSE)
  
  results[[as.character(i)]] <- res_xgboost
  if(i %% 50 == 0 ){
    #names(results) <- (i-50+1):i
    save(results, file = paste0(outdir, "tuneXboost_2classes_1_", as.character(as.integer(i/50 )), ".RData"))
    tunn <-getTableFromConfmatrices(results)
    write_tsv(tunn, file=paste0(outdir, "tuneXboost_2classes_1_", as.character(as.integer(i/50 )), ".tsv"))
    results <- list()
  } 
  
  cat("\tAccuracy=", round(res_xgboost$confmat$overall[1], 3), ", Kappa=", round(res_xgboost$confmat$overall[2], 3), "\n")
}

save(results, file = paste0(outdir, "tuneXboost_2classes_1_", as.character(as.integer(i/50 )+1), ".RData"))
tunn <-getTableFromConfmatrices(results)
write_tsv(tunn, file=paste0(outdir, "tuneXboost_2classes_1_", as.character(as.integer(i/50 )+1), ".tsv"))

ff <- list.files(outdir, pattern=".tsv", full.names = TRUE)

alldf <- map(ff, \(x) {
  aux <- read_tsv(x) 
  nn <- strsplit(basename(x), "_")[[1]]
  nn <- gsub(".tsv", "", nn[length(nn)])
  aux$set <- nn
  return(aux)
}) %>% bind_rows()

alldf %>% nrow
params_grid %>% nrow

allmerged <- cbind(alldf, params_grid[as.numeric(alldf$model) , ]) %>% 
  arrange(desc(Kappa_l1out))
write_tsv(allmerged, file=paste0(outdir, "tuneXboost_2classes_1_All", ".tsv"))

allmerged %>% head(1)

