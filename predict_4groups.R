#####################################
## Predict DEPR + OBESITY


randomforest_params_mult = list(ntree = 500, 
                           mtry = 4, 
                           nodesize = 5, 
                           balance_weights = TRUE)

xgboost_params_mult =  list(learning_rate=0.3,
                       max_depth=2,
                       nrounds =30,
                       min_child_weight=1, 
                       subsample =1,
                       colsample_bytree =0.6,
                       reg_lambda =1,
                       reg_alpha =0,
                       nthread=1,
                       objective= "multi:softprob",
                       num_class = 4, 
                       balance_weights = TRUE
)
catboost_params_mult <- list(
  iterations = 100,
  learning_rate = 0.05,
  depth = 4,
  loss_function = "MultiClass",
  eval_metric = "MultiClass",
  random_seed = 123,   
  use_best_model = TRUE,
  od_type = "Iter",
  od_wait = 20,
  verbose = FALSE,
  thread_count = 4,
  balance_weights = TRUE, #not used anymore
  bootstrap_type = "Bernoulli",
  l2_leaf_reg = 3,
  subsample = 0.6,  # only if bootstrap type ="Bernouilli"
  grow_policy = "Depthwise",
  auto_class_weights = "Balanced"
)

source(opt$predictive_functions)
load("/home/carlos/Escritorio/202311_DEPRESION/202311_DEPRESION/results_rstudio_10/DeSEQ2/DESEQ2_all.RData")
load("/home/carlos/Escritorio/202311_DEPRESION/ReviewJune2025/Results_rstudio/results2/phyloseq_original/phyloseq_all_list.RData")
load("/home/carlos/Escritorio/202311_DEPRESION/ReviewJune2025/Results_rstudio/results2/DAA_linda/lindalist.RData")
load("/home/carlos/Escritorio/202311_DEPRESION/202311_DEPRESION/results_rstudio_10/DESeq2_ControlVarsAlone/DESEQ2_controlVarsAlone_all.RData")
load("/home/carlos/Escritorio/202311_DEPRESION/202311_DEPRESION/results_rstudio_10/DESeq2_ControlVars/DESEQ2_controlVars_all.RData")

editObesity <- function(vec){
  return(ifelse(is.na(vec), NA, ifelse(vec == "Normopeso", "normal weight", "overweight")))
}
combineClasses <- function(phobj, v1, v2, newname, sep=":"){
  v1vec <- sample_data(phobj)[[v1]] %>% gsub("Depression", "D", .) %>% gsub("Control", "C", .)
  sample_data(phobj)[[newname]] <- paste(v1vec, sample_data(phobj)[[v2]], sep=sep)
  sample_data(phobj)[[newname]][is.na(sample_data(phobj)[[v1]]) | is.na(sample_data(phobj)[[v2]])] <- NA
  
  return(phobj)
}


interestvar <- "Condition"
var2add <- "obesidad"
newname <- "Depr_and_Ob"
vars2pca <- c(interestvar, var2add, newname)
phseq_to_use <- c("remove_tanda2")

ONAME <- "PredictDAA_multiclass"
opt$out <- paste0(opt$out, ONAME)
if(!dir.exists(opt$out)) dir.create(opt$out)
opt <- restaurar(opt)

all_model_multi_results <- list()
i <- phseq_to_use
#for(i in phseq_to_use){
  cat("Doing Predictive models for: ", i, "\n")
  all_model_multi_results[[i]] <- list()
  phobj <- all_phyloseq[[i]]
  
  sample_data(phobj)[[var2add]] <- editObesity(sample_data(phobj)[[var2add]])
  phobj <- combineClasses(phobj, interestvar, var2add, newname)
  samples <- sample_data(phobj)$sampleID[! is.na(sample_data(phobj)[, var2add])]
  phobj_filt <- phyloseq::prune_samples(samples, phobj)
  this_metadata <- sample_data(phobj_filt) %>% data.frame
  outdir <- paste0(opt$out, "PredictDAA_multiclass/", i, "/")
  opt$reserva <- opt$out
  opt$out <- outdir
  if(!dir.exists(opt$out)) dir.create(opt$out)
  
  ## Get data for PCA 
  df2pca <- if(is.null(daa_all[[i]]$vst_counts_df)){ daa_all[[i]]$norm_counts_df}else{ daa_all[[i]]$vst_counts_df }
  
  ## Get taxa 
  taxa_padj <- daa_all[[i]]$resdf %>% filter_taxa_padj
  taxa_praw <- daa_all[[i]]$resdf %>% filter_taxa_praw
  taxa_cov_padj <- daa_all_corrected_only[[i]]$BMI_log %>% filter_taxa_padj
  taxa_cov_praw <- daa_all_corrected_only[[i]]$BMI_log %>% filter_taxa_praw
  taxa_padj01 <- daa_all[[i]]$resdf %>% filter_taxa_padj(plim=0.01)
  taxa_padj001 <- daa_all[[i]]$resdf %>% filter_taxa_padj(plim=0.001)
  taxa_cov_padj01 <- daa_all_corrected_only[[i]]$BMI_log %>% filter_taxa_padj(plim=0.01)
  taxa_cov_padj001 <- daa_all_corrected_only[[i]]$BMI_log %>% filter_taxa_padj(plim=0.001)
  
  taxa_padj_corr <- daa_all_corrected[[i]]$BMI_log %>% filter_taxa_padj 
  taxa_padj_corr_01 <- daa_all_corrected[[i]]$BMI_log %>% filter_taxa_padj(plim=0.01)
  taxa_padj_corr_001 <- daa_all_corrected[[i]]$BMI_log %>% filter_taxa_padj(plim=0.001)
  #imcfname <- paste0(opt$reserva, "DESeq2_ControlVars/DeSEQ2/", i, "_BMI_log/",i, "_BMI_log_BMI_log_DAAshrinkNormal.tsv" )
  imcfname <- paste0("/home/carlos/Escritorio/202311_DEPRESION/202311_DEPRESION/results_rstudio_10/DESeq2_ControlVars/DeSEQ2/", i, "_BMI_log/",i, "_BMI_log_BMI_log_DAAshrinkNormal.tsv" )
  
  imctab <- read_tsv(imcfname)
  taxa_padj_cov_corr <- imctab %>% filter_taxa_padj
  taxa_padj_cov_corr_01 <- imctab %>% filter_taxa_padj(plim = 0.01)
  taxa_padj_cov_corr_001 <- imctab %>% filter_taxa_padj(plim = 0.001)
  
  taxa_padj_05 <- unique(c(taxa_padj, taxa_cov_padj))
  taxa_praw_05 <- unique(c(taxa_praw, taxa_cov_praw))
  taxa_padj_01 <- unique(c(taxa_padj01, taxa_cov_padj01))
  taxa_padj_001 <- unique(c(taxa_cov_padj001, taxa_padj001))
  taxa_padj_cov_05 <- unique(c(taxa_padj_corr, taxa_padj_cov_corr))
  taxa_padj_cov_01 <- unique(c(taxa_padj_corr_01, taxa_padj_cov_corr_01))
  taxa_padj_cov_001 <- unique(c(taxa_padj_corr_001, taxa_padj_cov_corr_001))
  
  taxa_padj_covorsingle_05 <- unique(c(taxa_padj_cov_05, taxa_padj_05))
  taxa_padj_covorsingle_01 <- unique(c(taxa_padj_01, taxa_padj_cov_01))
  taxa_padj_covorsingle_001 <- unique(c(taxa_padj_001, taxa_padj_cov_001))
  
  ## LinDA
  taxa_linda_padj_05 <- lindalist$firstContrast$resdf %>% filter_taxa_padj(plim=0.05, fc=0)
  taxa_linda_cov_padj_1 <- lindalist$contrastlist2$BMI_alone$resdf %>% filter_taxa_padj(plim=0.1, fc=0)
  taxa_linda_praw_05 <- lindalist$firstContrast$resdf %>% filter_taxa_praw(plim=0.05, fc=0)
  taxa_linda_cov_praw_05 <- lindalist$contrastlist2$BMI_alone$resdf %>% filter_taxa_praw(plim=0.05, fc=0)
  taxa_linda_padj_corr_05 <- lindalist$contrastlist2$`D_vs_C_adj.Age*Sex+BMI`$resdf %>% filter_taxa_padj(plim=0.05, fc=0)
  taxa_linda_cov_padj_corr_1 <- lindalist$contrastlist2$`BMI_adj.Age*Sex+Depr`$resdf %>% filter_taxa_padj(plim=0.1, fc=0)
  taxa_linda_praw_corr_05 <- lindalist$contrastlist2$`D_vs_C_adj.Age*Sex+BMI`$resdf %>% filter_taxa_praw(plim=0.01, fc=0)
  taxa_linda_cov_praw_corr_05 <- lindalist$contrastlist2$`BMI_adj.Age*Sex+Depr`$resdf %>% filter_taxa_praw(plim=0.01, fc=0)
  
  taxa_linda_padj = unique(taxa_linda_padj_05, taxa_linda_cov_padj_1)
  taxa_linda_praw = unique(taxa_linda_praw_05, taxa_linda_cov_praw_05)
  taxa_linda_padj_corr = unique(taxa_linda_padj_corr_05, taxa_linda_cov_padj_corr_1)
  taxa_linda_praw_corr = unique(taxa_linda_praw_corr_05, taxa_linda_cov_praw_corr_05)
  
  taxa_linda_padj_covorsingle <- unique(taxa_linda_padj, taxa_linda_padj_corr)
  taxa_linda_praw_covorsingle <- unique(taxa_linda_praw, taxa_linda_praw_corr)
  
  ## Make PCAs
  
  taxa_list <- list(
    "DiffTaxaPadjBase" = taxa_padj,
    "DiffTaxaPraw" = taxa_praw_05,
    "DiffTaxaPadj" = taxa_padj_05,
    "DiffTaxaPadj01" = taxa_padj_01,
    "DiffTaxaPadj001" = taxa_padj_001,
    "DiffTaxaPadjCov05" = taxa_padj_cov_05,
    "DiffTaxaPadjCov01" = taxa_padj_cov_01,
    "DiffTaxaPadjCov001" = taxa_padj_cov_001,
    "DiffTaxaPadjCovOrSingle05" = taxa_padj_covorsingle_05,
    "DiffTaxaPadjCovOrSingle01" = taxa_padj_covorsingle_01,
    "DiffTaxaPadjCovOrSingle001" = taxa_padj_covorsingle_001, 
    "LinDADiffTaxaPraw" = taxa_linda_praw,
    "LinDADiffTaxaPadj" = taxa_linda_padj,
    "LinDADiffTaxaPrawCov" = taxa_linda_praw_corr,
    "LinDADiffTaxaPadjCov" = taxa_linda_padj_corr,
    "LinDADiffTaxaPrawCovOrSingle" = taxa_linda_praw_covorsingle,
    "LinDADiffTaxaPadjCovOrSingle" = taxa_linda_padj_covorsingle
  )
  
  pcafilename <- paste0(opt$out, i, "_all_pcalists.RData")
  if(!file.exists(pcafilename)){
    all_pcalists <- map(names(taxa_list), 
                        \(x)makeAllPCAs(phobj_filt, df2pca, taxa_list[[x]], vars2pca, opt, x))
    names(all_pcalists) <- paste0("PCA_", names(taxa_list))
    save(all_pcalists, file = pcafilename)
  }else{
    load(file = pcafilename)
  }
  
  all_models_this <- map2(all_pcalists, names(taxa_list),
                          \(x, y)callDoAllModelsFromALLPCAs(x, name=paste0(i, y), 
                                                            metadata=this_metadata, 
                                                            vars2pca=c("Depr_and_Ob"), 
                                                            variable_plim=0.01,
                                                            xgboost_params= xgboost_params_mult,
                                                            catboost_params = catboost_params_mult,
                                                            randomforest_params = randomforest_params_mult,
                                                            do_smote=FALSE))
  names(all_models_this) <- names(taxa_list)
  save(all_models_this, file =paste0(outdir, "all_models_variableplim0.01.RData"))
  all_models_this2 <- map2(all_pcalists, names(taxa_list),
                          \(x, y)callDoAllModelsFromALLPCAs(x, name=paste0(i, y), 
                                                            metadata=this_metadata, 
                                                            vars2pca=c("Depr_and_Ob"), 
                                                            variable_plim=0.05,
                                                            xgboost_params= xgboost_params_mult,
                                                            catboost_params = catboost_params_mult,
                                                            randomforest_params = randomforest_params_mult,
                                                            do_smote=FALSE))
  names(all_models_this2) <- names(taxa_list)
  save(all_models_this2, file=paste0(outdir, "all_models_variableplim0.05.RData"))
  all_models_this3 <- map2(all_pcalists, names(taxa_list),
                           \(x, y)callDoAllModelsFromALLPCAs(x, name=paste0(i, y), 
                                                             metadata=this_metadata, 
                                                             vars2pca=c("Depr_and_Ob"), 
                                                             variable_plim=0.1,
                                                             xgboost_params= xgboost_params_mult,
                                                             catboost_params = catboost_params_mult,
                                                             randomforest_params = randomforest_params_mult,
                                                             do_smote=FALSE))
  names(all_models_this3) <- names(taxa_list)
  save(all_models_this3, file=paste0(outdir, "all_models_variableplim0.1.RData"))
  all_models_this3_smote <- map2(all_pcalists, names(taxa_list),
                           \(x, y)callDoAllModelsFromALLPCAs(x, name=paste0(i, y), 
                                                             metadata=this_metadata, 
                                                             vars2pca=c("Depr_and_Ob"), 
                                                             variable_plim=0.1,
                                                             xgboost_params= xgboost_params_mult,
                                                             catboost_params = catboost_params_mult,
                                                             randomforest_params = randomforest_params_mult,
                                                             do_smote=TRUE))
  names(all_models_this3_smote) <- names(taxa_list)
  save(all_models_this3_smote, file=paste0(outdir, "all_models_variableplim0.1_SMOTE.RData"))
  all_models_this4 <- map2(all_pcalists, names(taxa_list),
                           \(x, y)callDoAllModelsFromALLPCAs(x, name=paste0(i, y), 
                                                             metadata=this_metadata, 
                                                             vars2pca=c("Depr_and_Ob"), 
                                                             variable_plim=1,
                                                             xgboost_params= xgboost_params_mult,
                                                             catboost_params = catboost_params_mult,
                                                             randomforest_params = randomforest_params_mult,
                                                             do_smote=FALSE))
  names(all_models_this4) <- names(taxa_list)
  save(all_models_this4, file = paste0(outdir, "all_models_variableplim1.0.RData"))
  
  #see quickly which conditions are best
  #all_models_this %>% map(\(x) x$modummary %>% select(model, Accuracy_l1out, Kappa_l1out, BalancedAccuracy_l1out, AUC_l1out) %>% head(2))
  #all_models_this2 %>% map(\(x) x$modummary %>% select(model, Accuracy_l1out, Kappa_l1out, BalancedAccuracy_l1out, AUC_l1out) %>% head(2))
  all_models_this3 %>% map(\(x) x$modummary %>% select(model, Accuracy_l1out, Kappa_l1out, BalancedAccuracy_l1out, AUC_l1out) %>% head(2))
  all_models_this3_smote %>% map(\(x) x$modummary %>% select(model, Accuracy_l1out, Kappa_l1out, BalancedAccuracy_l1out, AUC_l1out) %>% head(2))
  #all_models_this4 %>% map(\(x) x$modummary %>% select(model, Accuracy_l1out, Kappa_l1out, BalancedAccuracy_l1out, AUC_l1out) %>% head(2))
  # all_models_this3 is the best (PCs with p<0.1)
  
  names(taxa_list) <- paste0("taxa_", names(taxa_list))
  save(taxa_list, file = paste0(outdir, "taxa_list.0.RData"))
  
  all_model_multi_results[[paste0(i, "_1")]] <- c(all_models_this, all_pcalists, taxa_list)
  all_model_multi_results[[paste0(i, "_1")]]$metadata <- sample_data(phobj_filt) %>% data.frame
  
  all_model_multi_results[[paste0(i, "_2")]] <- c(all_models_this2, all_pcalists, taxa_list)
  all_model_multi_results[[paste0(i, "_2")]]$metadata <- sample_data(phobj_filt) %>% data.frame
  
  all_model_multi_results[[paste0(i, "_3")]] <- c(all_models_this3, all_pcalists, taxa_list)
  all_model_multi_results[[paste0(i, "_3")]]$metadata <- sample_data(phobj_filt) %>% data.frame
  
  all_model_multi_results[[paste0(i, "_3smote")]] <- c(all_models_this3_smote, all_pcalists, taxa_list)
  all_model_multi_results[[paste0(i, "_3smote")]]$metadata <- sample_data(phobj_filt) %>% data.frame
  
  all_model_multi_results[[paste0(i, "_4")]] <- c(all_models_this4, all_pcalists, taxa_list)
  all_model_multi_results[[paste0(i, "_4")]]$metadata <- sample_data(phobj_filt) %>% data.frame
  
  opt <- restaurar(opt)
#}
opt <- restaurar(opt)
multi_model_fname <- paste0(opt$out, ONAME, "/all_model_multi_results.RData")
save(all_model_multi_results, file=multi_model_fname)
#load(multi_model_fname)

#for(k in names(all_model_multi_results[[i]]$modresults)) {cat("\n\n", k);all_model_multi_results[[i]]$modresults[[k]]$modummary %>% head(2) %>% print}

#Integrate
opt$out <- paste0(opt$out, "PredictDAA_multiclass_plots2/")
if(!dir.exists(opt$out)) dir.create(opt$out)
makeLinePlotComparingPhobjs( all_model_multi_results, opt, models_name1 = "DiffTaxaPadjBase", models_name2 = "DiffTaxaPadjCov05")
## Compare with Bacteria in componets
condnames <- names(all_model_multi_results$remove_tanda2_3)[c(1,6,7,9,10)]
phname <- "remove_tanda2_3"
makeLinePlotComparingSamePhobjModels_Cov(phname, condnames, all_model_multi_results, "TaxaGroups_trim", opt, w=8, h=5)
condnames <- names(all_model_multi_results$remove_tanda2_3)[c(9)]
phname <- "remove_tanda2_3"
makeLinePlotComparingSamePhobjModels_Cov(phname, condnames, all_model_multi_results, "TaxaGroups_trim2", opt, w=8, h=5)

# add SMOTE to same list
for(nn in names(all_model_multi_results$remove_tanda2_3)){
  newname <- paste0(nn, "_SMOTE")
  all_model_multi_results$remove_tanda2_3[[newname]] <- all_model_multi_results$remove_tanda2_3smote[[nn]]
}

condnames <- names(all_model_multi_results$remove_tanda2_3)[c(1,3,6,7,9,10, 13:17)]
for(nn in condnames){
#walk(names(all_model_multi_results), \(x)makeLinePlotComparingSamePhobjModels(x, all_model_multi_results, 
#                                                                        opt, w=6, h=8, get_pcnames_from = nn, 
#                                                                        order_by_measure = "BalancedAccuracy_l1out",
#                                                                        name=paste0(nn,"_sortBAcc")))
#walk(names(all_model_multi_results), \(x)makeLinePlotComparingSamePhobjModels(x, all_model_multi_results, 
#                                                                              opt, w=6, h=8, get_pcnames_from = nn, 
#                                                                              order_by_measure = "AUC_l1out",
#                                                                              name=paste0(nn,"_sortAUC")))
#walk(names(all_model_multi_results), \(x)makeLinePlotComparingSamePhobjModels(x, all_model_multi_results, 
#                                                                              opt, w=6, h=8, get_pcnames_from = nn, 
#                                                                              order_by_measure = "Kappa_l1out",
#                                                                              name=paste0(nn,"_sortKappa")))
  
  cat(nn, "\n")
  sel_method_name <- ifelse(grepl("LinDA", nn), "PCA LinDA","PCA DESeq")
  makeLinePlotComparingSamePhobjModels("remove_tanda2_1", all_model_multi_results, 
                                       opt, w=6, h=8, get_pcnames_from = nn, 
                                       order_by_measure = "BalancedAccuracy_l1out",
                                       sel_method_name = sel_method_name,
                                       name=paste0(nn,"_sortBAcc"))
  
  makeLinePlotComparingSamePhobjModels("remove_tanda2_1", all_model_multi_results, 
                                       opt, w=6, h=8, get_pcnames_from = nn, 
                                       order_by_measure = "AUC_l1out",
                                       sel_method_name = sel_method_name,
                                       name=paste0(nn,"_sortAUC"))
  
  makeLinePlotComparingSamePhobjModels("remove_tanda2_1", all_model_multi_results, 
                                       opt, w=6, h=8, get_pcnames_from = nn, 
                                       order_by_measure = "Kappa_l1out",
                                       sel_method_name = sel_method_name,
                                       name=paste0(nn,"_sortKappa"))
  
  makeLinePlotComparingSamePhobjModels("remove_tanda2_2", all_model_multi_results, 
                                       opt, w=6, h=8, get_pcnames_from = nn, 
                                       order_by_measure = "BalancedAccuracy_l1out",
                                       sel_method_name = sel_method_name,
                                       name=paste0(nn,"_sortBAcc"))
  
  makeLinePlotComparingSamePhobjModels("remove_tanda2_2", all_model_multi_results, 
                                       opt, w=6, h=8, get_pcnames_from = nn, 
                                       sel_method_name = sel_method_name,
                                       order_by_measure = "AUC_l1out",
                                       name=paste0(nn,"_sortAUC"))
  
  makeLinePlotComparingSamePhobjModels("remove_tanda2_2", all_model_multi_results, 
                                       opt, w=6, h=8, get_pcnames_from = nn, 
                                       order_by_measure = "Kappa_l1out",
                                       sel_method_name = sel_method_name,
                                       name=paste0(nn,"_sortKappa"))
  
makeLinePlotComparingSamePhobjModels("remove_tanda2_3", all_model_multi_results, 
                                         opt, w=6, h=8, get_pcnames_from = nn, 
                                         order_by_measure = "BalancedAccuracy_l1out",
                                         sel_method_name = sel_method_name,
                                         name=paste0(nn,"_sortBAcc"))
  
makeLinePlotComparingSamePhobjModels("remove_tanda2_3", all_model_multi_results, 
                                         opt, w=6, h=8, get_pcnames_from = nn, 
                                         sel_method_name = sel_method_name,
                                         order_by_measure = "AUC_l1out",
                                         name=paste0(nn,"_sortAUC"))

makeLinePlotComparingSamePhobjModels("remove_tanda2_3", all_model_multi_results, 
                                        opt, w=6, h=8, get_pcnames_from = nn, 
                                        order_by_measure = "Kappa_l1out",
                                     sel_method_name = sel_method_name,
                                        name=paste0(nn,"_sortKappa"))

makeLinePlotComparingSamePhobjModels("remove_tanda2_3", all_model_multi_results, 
                                     opt, w=6, h=8, get_pcnames_from = nn, 
                                     plot_normal_with_smote=TRUE,
                                     sel_method_name = sel_method_name,
                                     order_by_measure = "BalancedAccuracy_l1out",
                                     name=paste0(nn,"_sortBAccCompSMOTE"))

makeLinePlotComparingSamePhobjModels("remove_tanda2_3", all_model_multi_results, 
                                     opt, w=6, h=8, get_pcnames_from = nn, 
                                     plot_normal_with_smote=TRUE,
                                     order_by_measure = "AUC_l1out",
                                     sel_method_name = sel_method_name,
                                     name=paste0(nn,"_sortAUCCompSMOTE"))

makeLinePlotComparingSamePhobjModels("remove_tanda2_3", all_model_multi_results, 
                                     opt, w=6, h=8, get_pcnames_from = nn, 
                                     plot_normal_with_smote=TRUE,
                                     order_by_measure = "Kappa_l1out",
                                     sel_method_name = sel_method_name,
                                     name=paste0(nn,"_sortKappaCompSMOTE"))
}

outdir <- paste0(opt$out, "/all_performances/")
if(!dir.exists(outdir)) dir.create(outdir)

for(n1 in names(all_model_multi_results)){
  ns2 <- names(all_model_multi_results[[n1]])
  ns2 <- ns2[grepl("^Diff", ns2, perl=T) | grepl("^LinDA", ns2, perl=T)]
  for(n2 in ns2){
    mmsum <- all_model_multi_results[[n1]][[n2]]$modummary %>% 
      arrange(desc(BalancedAccuracy_l1out))
    mname <- paste0(outdir, n1, "_", n2, "_modsummary.tsv")
    write_tsv(mmsum, file = mname)
  }
  
  
}

### Fins aci
#####################################
walk(names(all_model_multi_results), \(x)makeLinePlotComparingSamePhobjModels(x, all_model_multi_results, 
                                                                              opt, w=6, h=8, get_pcnames_from = "DiffTaxaPadjCov05", 
                                                                              name = "DiffTaxaPadjCov05"))
## Plot boxplot PCs
pca2use <- "DiffTaxaPadjCov05"
pcBoxplots <- map(names(all_model_multi_results), \(x){
  makePCsBoxplot(x, all_model_multi_results, opt, 
                 get_pcnames_from = pca2use,
                 pca_name =  paste0("PCA_", pca2use), 
                 varname = "Depr_and_Ob",
                 w=12, h=8)
})
names(pcBoxplots) <- names(all_model_multi_results)

## Plot barplot PCs and LFC
pcBarplots <- map(names(all_model_multi_results), makePCBarplot, all_model_multi_results, pcBoxplots, daa_all, opt,
                  get_pcnames_from = pca2use,
                  pca_name =  paste0("PCA_", pca2use), 
                  varname = "Depr_and_Ob",
                  w=12, h=20)
names(pcBarplots) <- names(all_model_multi_results)

# Plot KNN (best model)

phname <- "remove_tanda2"
predplots <- map(names(all_model_multi_results), plotAllModelPredictions, all_model_multi_results, opt,
                 get_pcnames_from = pca2use,
                 pca_name =  paste0("PCA_", pca2use), 
                 varname = "Depr_and_Ob",
                 pred_mode="l1o")


opt <- restaurar(opt)

#Contingency table
library(ggmosaic)
df <- sample_data(phobj) %>% data.frame %>% 
  dplyr::filter(!is.na(IMC)) %>% 
  dplyr::mutate(overweight = ifelse(IMC>25, "overweight", "normal weight")) %>% 
  dplyr::mutate(Depr_ow = paste(Condition, overweight, sep=":"))
var1 <- "Condition"
var2 <- "overweight"

makeContingencyPlot(df, var1, var2, opt$out, "ObesityVsDepression", w=6, h=5)


######################### Use data from obesity

opt$out <- paste0(opt$out, "PredictDAA_multiclass_fromOb")
if(!dir.exists(opt$out)) dir.create(opt$out)
opt <- restaurar(opt)

all_model_multi_results <- list()
for(i in phseq_to_use){
  cat("Doing Predictive models for: ", i, "\n")
  all_model_multi_results[[i]] <- list()
  phobj <- all_phyloseq[[i]]
  
  sample_data(phobj)[[var2add]] <- editObesity(sample_data(phobj)[[var2add]])
  phobj <- combineClasses(phobj, interestvar, var2add, newname)
  samples <- sample_data(phobj)$sampleID[! is.na(sample_data(phobj)[, var2add])]
  phobj_filt <- phyloseq::prune_samples(samples, phobj)
  
  outdir <- paste0(opt$out, "PredictDAA_multiclass_fromOb/", i, "/")
  opt$reserva <- opt$out
  opt$out <- outdir
  if(!dir.exists(opt$out)) dir.create(opt$out)
  
  ## Get data for PCA 
  df2pca <- if(is.null(daa_all[[i]]$vst_counts_df)){ daa_all[[i]]$norm_counts_df}else{ daa_all[[i]]$vst_counts_df }
  
  ## Get taxa 
  taxa_padj <- daa_all[[i]]$resdf %>% filter_taxa_padj
  taxa_praw <- daa_all[[i]]$resdf %>% filter_taxa_praw
  taxa_cov_padj <- daa_all_corrected_only[[i]]$ob_o_sobrepeso %>% filter_taxa_padj
  taxa_cov_praw <- daa_all_corrected_only[[i]]$ob_o_sobrepeso %>% filter_taxa_praw
  taxa_padj01 <- daa_all[[i]]$resdf %>% filter_taxa_padj(plim=0.01)
  taxa_padj001 <- daa_all[[i]]$resdf %>% filter_taxa_padj(plim=0.001)
  taxa_cov_padj01 <- daa_all_corrected_only[[i]]$ob_o_sobrepeso %>% filter_taxa_padj(plim=0.01)
  taxa_cov_padj001 <- daa_all_corrected_only[[i]]$ob_o_sobrepeso %>% filter_taxa_padj(plim=0.001)
  
  taxa_padj_corr <- daa_all_corrected[[i]]$ob_o_sobrepeso %>% filter_taxa_padj 
  taxa_padj_corr_01 <- daa_all_corrected[[i]]$ob_o_sobrepeso %>% filter_taxa_padj(plim=0.01)
  taxa_padj_corr_001 <- daa_all_corrected[[i]]$ob_o_sobrepeso %>% filter_taxa_padj(plim=0.001)
  imcfname <- paste0(opt$reserva, "DESeq2_ControlVars/DeSEQ2/", i, "_ob_o_sobrepeso/",i, "_ob_o_sobrepeso_ob_o_sobrepeso_normal.weight_vs_overweight_DAAshrinkNormal.tsv" )
  imctab <- read_tsv(imcfname)
  taxa_padj_cov_corr <- imctab %>% filter_taxa_padj
  taxa_padj_cov_corr_01 <- imctab %>% filter_taxa_padj(plim = 0.01)
  taxa_padj_cov_corr_001 <- imctab %>% filter_taxa_padj(plim = 0.001)
  
  taxa_padj_05 <- unique(c(taxa_padj, taxa_cov_padj))
  taxa_praw_05 <- unique(c(taxa_praw, taxa_cov_praw))
  taxa_padj_01 <- unique(c(taxa_padj01, taxa_cov_padj01))
  taxa_padj_001 <- unique(c(taxa_cov_padj001, taxa_padj001))
  taxa_padj_cov_05 <- unique(c(taxa_padj_corr, taxa_padj_cov_corr))
  taxa_padj_cov_01 <- unique(c(taxa_padj_corr_01, taxa_padj_cov_corr_01))
  taxa_padj_cov_001 <- unique(c(taxa_padj_corr_001, taxa_padj_cov_corr_001))
  
  taxa_padj_covorsingle_05 <- unique(c(taxa_padj_cov_05, taxa_padj_05))
  taxa_padj_covorsingle_01 <- unique(c(taxa_padj_01, taxa_padj_cov_01))
  taxa_padj_covorsingle_001 <- unique(c(taxa_padj_001, taxa_padj_cov_001))
  ## Make PCAs
  
  taxa_list <- list(
    "DiffTaxaPadjBase" = taxa_padj,
    "DiffTaxaPraw" = taxa_praw_05,
    "DiffTaxaPadj" = taxa_padj_05,
    "DiffTaxaPadj01" = taxa_padj_01,
    "DiffTaxaPadj001" = taxa_padj_001,
    "DiffTaxaPadjCov05" = taxa_padj_cov_05,
    "DiffTaxaPadjCov01" = taxa_padj_cov_01,
    "DiffTaxaPadjCov001" = taxa_padj_cov_001,
    "DiffTaxaPadjCovOrSingle05" = taxa_padj_covorsingle_05,
    "DiffTaxaPadjCovOrSingle01" = taxa_padj_covorsingle_01,
    "DiffTaxaPadjCovOrSingle001" = taxa_padj_covorsingle_001
  )
  all_pcalists <- map(names(taxa_list), 
                      \(x)makeAllPCAs(phobj_filt, df2pca, taxa_list[[x]], vars2pca, opt, x))
  names(all_pcalists) <- paste0("PCA_", names(taxa_list))
  
  all_models_this <- map2(all_pcalists, names(taxa_list),
                          \(x, y)callDoAllModelsFromALLPCAs(x, name=paste0(i, y), 
                                                            metadata=this_metadata, 
                                                            vars2pca=c("Depr_and_Ob")))
  names(all_models_this) <- names(taxa_list)
  names(taxa_list) <- paste0("taxa_", names(taxa_list))
  
  
  all_model_multi_results[[i]] <- c(all_models_this, all_pcalists, taxa_list)
  all_model_multi_results[[i]]$metadata <- sample_data(phobj_filt) %>% data.frame
  
  opt <- restaurar(opt)
}
opt <- restaurar(opt)
save(all_model_multi_results, file=paste0(opt$out, "PredictDAA_multiclass_fromOb/all_model_multi_results.RData"))
#load(file=paste0(opt$out, "PredictDAA/all_model_multi_results.RData"))

#for(k in names(all_model_results[[i]]$modresults)) {cat("\n\n", k);all_model_results[[i]]$modresults[[k]]$modummary %>% head(2) %>% print}

#Integrate
opt$out <- paste0(opt$out, "PredictDAA_multiclass_fromOb/")
makeLinePlotComparingPhobjs( all_model_multi_results, opt, models_name1 = "DiffTaxaPadjBase", models_name2 = "DiffTaxaPadjCov05")
## Compare with Bacteria in componets
condnames <- names(all_model_multi_results$remove_tanda2)[c(1,6,7,9,10)]
makeLinePlotComparingSamePhobjModels_Cov(phname, condnames, all_model_multi_results, "TaxaGroups_trim", opt, w=8, h=5)

walk(names(all_model_multi_results), \(x)makeLinePlotComparingSamePhobjModels(x, all_model_multi_results, 
                                                                        opt, w=6, h=8, get_pcnames_from = "DiffTaxaPadjCov05"))
## Plot boxplot PCs
pca2use <- "DiffTaxaPadjCov05"
pcBoxplots <- map(names(all_model_multi_results), \(x){
  makePCsBoxplot(x, all_model_multi_results, opt, 
                 get_pcnames_from = pca2use,
                 pca_name =  paste0("PCA_", pca2use), 
                 varname = "Depr_and_Ob",
                 w=12, h=8)
})
names(pcBoxplots) <- names(all_model_multi_results)

## Plot barplot PCs and LFC
pcBarplots <- map(names(all_model_multi_results), makePCBarplot, all_model_multi_results, pcBoxplots, daa_all, opt,
                  get_pcnames_from = pca2use,
                  pca_name =  paste0("PCA_", pca2use), 
                  varname = "Depr_and_Ob",
                  w=12, h=20)
names(pcBarplots) <- names(all_model_multi_results)

# Plot KNN (best model)

phname <- "remove_tanda2"
predplots <- map(names(all_model_multi_results), plotAllModelPredictions, all_model_multi_results, opt,
                 get_pcnames_from = pca2use,
                 pca_name =  paste0("PCA_", pca2use), 
                 varname = "Depr_and_Ob",
                 pred_mode="l1o")


opt <- restaurar(opt)

#Contingency table
library(ggmosaic)
df <- sample_data(phobj) %>% data.frame %>% 
  dplyr::filter(!is.na(IMC)) %>% 
  dplyr::mutate(overweight = ifelse(IMC>25, "overweight", "normal weight")) %>% 
  dplyr::mutate(Depr_ow = paste(Condition, overweight, sep=":"))
var1 <- "Condition"
var2 <- "overweight"

makeContingencyPlot(df, var1, var2, opt$out, "ObesityVsDepression", w=6, h=5)
