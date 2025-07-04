

res_all <- list()
nn <- "BASE__"
cat(nn, "\n")


res_all[[paste0("randomForest", nn)]] <-  make_randomForest_l1o(datasc, levs, varnames, folds=c(), do_smote=FALSE,
                                                                smote_params=list(K=K, dup_size=dup_size))
res_all[[paste0("GLM", nn)]] <-  make_randomForest_l1o(datasc, levs, varnames, 
                                                       folds=c(), do_smote=FALSE,
                                                       smote_params=list(K=K, dup_size=dup_size))
res_all[[paste0("SVMlin", nn)]] <-  make_svm_l1o(datasc, levs, varnames, , kernel="linear",
                                                 folds=c(), do_smote=FALSE,
                                                 smote_params=list(K=K, dup_size=dup_size))
res_all[[paste0("SVMrad", nn)]] <-  make_svm_l1o(datasc, levs, varnames, , kernel="radial",
                                                 folds=c(), do_smote=FALSE,
                                                 smote_params=list(K=K, dup_size=dup_size))
res_all[[paste0("C50tree", nn)]] <-  make_classifTree_l1o(datasc, levs, varnames, 
                                                          folds=c(), do_smote=FALSE,
                                                          smote_params=list(K=K, dup_size=dup_size))
res_all[[paste0("NaiveBayes", nn)]] <-  makeNaiveBayes_l1o(datasc, levs, varnames,  SEED=234324,
                                                           folds=c(), do_smote=FALSE,
                                                           smote_params=list(K=K, dup_size=dup_size))
res_all[[paste0("KNN5", nn)]] <-  makeKnn_l1o(datasc, levs, varnames, different_ks=5,
                                              folds=c(), do_smote=FALSE,
                                              smote_params=list(K=K, dup_size=dup_size))$`K=5`

res_all[[paste0("KMEANS", nn)]] <-  makeKmeans_l1o(datasc, levs, varnames, SEED=234324, 
                                                   folds=c(), do_smote=FALSE,
                                                   smote_params=list(K=K, dup_size=dup_size))

for(K in 3:7) {
  for(dup_size in 1:3){
    nn <- paste0("K", K, "_", "DS", dup_size)
    cat(nn, "\n")
    res_all[[paste0("randomForest__", nn)]] <-  make_randomForest_l1o(datasc, levs, varnames, folds=c(), do_smote=TRUE,
                                            smote_params=list(K=K, dup_size=dup_size))
    res_all[[paste0("GLM__", nn)]] <-  make_randomForest_l1o(datasc, levs, varnames, 
                                                           folds=c(), do_smote=TRUE,
                                                                    smote_params=list(K=K, dup_size=dup_size))
    res_all[[paste0("SVMlin__", nn)]] <-  make_svm_l1o(datasc, levs, varnames, , kernel="linear",
                                                              folds=c(), do_smote=TRUE,
                                                                    smote_params=list(K=K, dup_size=dup_size))
    res_all[[paste0("SVMrad__", nn)]] <-  make_svm_l1o(datasc, levs, varnames, , kernel="radial",
                                                              folds=c(), do_smote=TRUE,
                                                                    smote_params=list(K=K, dup_size=dup_size))
    res_all[[paste0("C50tree__", nn)]] <-  make_classifTree_l1o(datasc, levs, varnames, 
                                                               folds=c(), do_smote=TRUE,
                                                                    smote_params=list(K=K, dup_size=dup_size))
    res_all[[paste0("NaiveBayes__", nn)]] <-  makeNaiveBayes_l1o(datasc, levs, varnames,  SEED=234324,
                                                                  folds=c(), do_smote=TRUE,
                                                               smote_params=list(K=K, dup_size=dup_size))
    res_all[[paste0("KNN5__", nn)]] <-  makeKnn_l1o(datasc, levs, varnames, different_ks=5,
                                                            folds=c(), do_smote=TRUE,
                                                               smote_params=list(K=K, dup_size=dup_size))$`K=5`
    res_all[[paste0("KMEANS__", nn)]] <-  makeKmeans_l1o(datasc, levs, varnames, SEED=234324, 
                                                       folds=c(), do_smote=TRUE,
                                                            smote_params=list(K=K, dup_size=dup_size))
  
  }}

RF_table <-getTableFromConfmatrices(res_all)

RF_table2 <- RF_table %>% 
  filter(!grepl("BASE__", model)) %>% 
  tidyr::separate(model, into=c("model", "smote_params"), sep="__")

RF_tablesum <- RF_table2 %>% 
  group_by(smote_params) %>% 
  dplyr::summarise(mean_acc=mean(Accuracy_l1out), mean_kappa = mean(Kappa_l1out), 
                   sd_acc=sd(Accuracy_l1out), sd_kappa = sd(Kappa_l1out), 
                   max_acc = max(Accuracy_l1out), max_kappa = max(Kappa_l1out), 
                   min_acc = min(Accuracy_l1out), min_kappa = min(Kappa_l1out)) %>% 
  dplyr::arrange(desc(mean_kappa))

# select K5_DS1 
RF_table2 %>% filter(smote_params == "K5_DS1")

