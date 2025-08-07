# Predict

mycols <- grDevices::colorRampPalette( wesanderson::wes_palette("Royal1"))(5)

options(ggplot2.discrete.fill = mycols)
options(ggplot2.discrete.colour = mycols)

load("/home/carlos/Escritorio/202311_DEPRESION/202311_DEPRESION/results_rstudio_10/DeSEQ2/DESEQ2_all.RData")
load("/home/carlos/Escritorio/202311_DEPRESION/ReviewJune2025/Results_rstudio/results2/phyloseq_original/phyloseq_all_list.RData")
load("/home/carlos/Escritorio/202311_DEPRESION/ReviewJune2025/Results_rstudio/results2/DAA_linda/lindalist.RData")

#opt$out <- "/home/carmoma/Desktop/202311_DEPRESION/results_rstudio_v2_1/"

OUTDNAME <- "PredictDAA_L1O_3"
vars2pca <- c("Condition", "Sexo", "Edad")
phseq_to_use <- c("remove_tanda2")

# These are defined in the G4Micro package
NFOLDS = 0
randomforest_params = randomforest_params_default
xgboost_params =  xgboost_params_default
catboost_params <- catboost_params_default
smote_params = smote_params_default


opt <- restaurar(opt)
opt$out <- paste0(opt$out, OUTDNAME)
if(!dir.exists(opt$out)) dir.create(opt$out)
opt <- restaurar(opt)

all_model_results <- list()
#for(i in phseq_to_use){
i <- phseq_to_use
  cat("Doing Predictive models for: ", i, "\n")
  all_model_results[[i]] <- list()
  phobj <- all_phyloseq[[i]]
  outdir <- paste0(opt$out, OUTDNAME, "/", i, "/")
  opt$reserva <- opt$out
  opt$out <- outdir
  if(!dir.exists(opt$out)) dir.create(opt$out)

  taxa_padj <- daa_all[[i]]$resdf %>% dplyr::filter(padj <= opt$pval & abs(log2FoldChangeShrink) >= log2(opt$fc) ) %>%
    pull(taxon)
  taxa_praw <- daa_all[[i]]$resdf %>% dplyr::filter(pvalue <= opt$pval & abs(log2FoldChangeShrink) >= log2(opt$fc) ) %>%
    pull(taxon)

  taxa_linda_padj <- lindalist$firstContrast$resdf %>% dplyr::filter(padj <= opt$pval & abs(log2FoldChangeShrink) >= log2(opt$fc) ) %>%
    pull(taxon)

  df2pca <- if(is.null(daa_all[[i]]$vst_counts_df)){ daa_all[[i]]$norm_counts_df}else{ daa_all[[i]]$vst_counts_df }
  all_pcas_adj <- makeAllPCAs(phobj, df2pca, taxa_padj, vars2pca, opt, "DiffTaxaPadj")
  all_pcas_praw <- makeAllPCAs(phobj, df2pca, taxa_praw, vars2pca, opt, "DiffTaxaPraw")
  all_pcas_linda_padj <- makeAllPCAs(phobj, df2pca, taxa_linda_padj, vars2pca, opt, "DiffTaxaLindaPadj")

  taxa_padj01 <- daa_all[[i]]$resdf %>% dplyr::filter(padj <= 0.01 & abs(log2FoldChangeShrink) >= log2(opt$fc) ) %>%
    pull(taxon)
  taxa_padj001 <- daa_all[[i]]$resdf %>% dplyr::filter(padj <= 0.001 & abs(log2FoldChangeShrink) >= log2(opt$fc) ) %>%
    pull(taxon)
  taxa_padj01_linda <- lindalist$firstContrast$resdf %>% dplyr::filter(padj <= 0.01 & abs(log2FoldChangeShrink) >= log2(opt$fc) ) %>%
    pull(taxon)
  all_pcas_adj01 <- makeAllPCAs(phobj, df2pca, taxa_padj01, vars2pca, opt, "DiffTaxaPadj01")
  all_pcas_adj001 <- makeAllPCAs(phobj, df2pca, taxa_padj001, vars2pca, opt, "DiffTaxaPadj001")
  all_pcas_adj01_linda <- makeAllPCAs(phobj, df2pca, taxa_padj01_linda, vars2pca, opt, "DiffTaxaLindaPadj01")

  this_metadata <- sample_data(phobj) %>% data.frame
  all_model_results[[i]][["padj_taxa_taxa"]] <- taxa_padj
  all_model_results[[i]][["praw_taxa_taxa"]] <- taxa_praw
  all_model_results[[i]][["padjlinda_taxa_taxa"]] <- taxa_linda_padj
  all_model_results[[i]][["padj_taxa_pcas"]] <- all_pcas_adj
  all_model_results[[i]][["praw_taxa_pcas"]] <- all_pcas_praw
  all_model_results[[i]][["padjlinda_taxa_pcas"]] <- all_pcas_linda_padj
  all_model_results[[i]][["padj_taxa_res"]] <- callDoAllModelsFromALLPCAs(all_pcas_adj, name=paste0(i, "ConditionPadj"),
                                                                          metadata=this_metadata, vars2pca=c("Condition"), nfolds = NFOLDS,
                                                                          xgboost_params = xgboost_params, catboost_params = catboost_params,
                                                                          randomforest_params = randomforest_params)

  all_model_results[[i]][["padj_taxa_res_SMOTE"]] <- callDoAllModelsFromALLPCAs(all_pcas_adj, name=paste0(i, "ConditionPadjSMOTE"),
                                                                          metadata=this_metadata, vars2pca=c("Condition"), nfolds = NFOLDS,
                                                                          xgboost_params = xgboost_params, , catboost_params = catboost_params,
                                                                          randomforest_params = randomforest_params,
                                                                          do_smote = TRUE, smote_params = smote_params)

  all_model_results[[i]][["praw_taxa_res"]] <- callDoAllModelsFromALLPCAs(all_pcas_praw, name=paste0(i, "ConditionPraw"),
                                                                          metadata=this_metadata, vars2pca=c("Condition"), nfolds = NFOLDS,
                                                                          xgboost_params = xgboost_params, catboost_params = catboost_params,
                                                                          randomforest_params = randomforest_params)

  all_model_results[[i]][["praw_taxa_res_SMOTE"]] <- callDoAllModelsFromALLPCAs(all_pcas_praw, name=paste0(i, "ConditionPrawSMOTE"),
                                                                          metadata=this_metadata, vars2pca=c("Condition"), nfolds = NFOLDS,
                                                                          xgboost_params = xgboost_params, catboost_params = catboost_params,
                                                                          randomforest_params = randomforest_params,
                                                                          do_smote = TRUE, smote_params = smote_params)

  all_model_results[[i]][["padj_taxa_res01"]] <- callDoAllModelsFromALLPCAs(all_pcas_adj01, name=paste0(i, "ConditionPadj01"),
                                                                            metadata=this_metadata, vars2pca=c("Condition"), nfolds = NFOLDS,
                                                                            xgboost_params = xgboost_params, catboost_params = catboost_params,
                                                                            randomforest_params = randomforest_params)

  all_model_results[[i]][["padj_taxa_res001"]] <- callDoAllModelsFromALLPCAs(all_pcas_adj001, name=paste0(i, "ConditionPadj001"),
                                                                             metadata=this_metadata, vars2pca=c("Condition"), nfolds = NFOLDS,
                                                                             xgboost_params = xgboost_params, catboost_params = catboost_params,
                                                                             randomforest_params = randomforest_params)

  all_model_results[[i]][["padj_taxa_resAll"]] <- callDoAllModelsFromALLPCAs(all_pcas_adj, name=paste0(i, "ConditionPadjAllPCs"),
                                                                          metadata=this_metadata, vars2pca=c("Condition"), nfolds = NFOLDS,
                                                                          variable_plim = 1,
                                                                          xgboost_params = xgboost_params, catboost_params = catboost_params,
                                                                          randomforest_params = randomforest_params)


  all_model_results[[i]][["padj_taxa_res05linda"]] <- callDoAllModelsFromALLPCAs(all_pcas_linda_padj, name=paste0(i, "ConditionPadjLinda05"),
                                                                                 metadata=this_metadata, vars2pca=c("Condition"), nfolds = NFOLDS,
                                                                                 xgboost_params = xgboost_params, catboost_params = catboost_params,
                                                                                 randomforest_params = randomforest_params)

  all_model_results[[i]][["padj_taxa_res01linda"]] <- callDoAllModelsFromALLPCAs(all_pcas_adj01_linda, name=paste0(i, "ConditionPadjLinda01"),
                                                                                 metadata=this_metadata, vars2pca=c("Condition"), nfolds = NFOLDS,
                                                                                 xgboost_params = xgboost_params, catboost_params = catboost_params,
                                                                                 randomforest_params = randomforest_params)

  all_model_results[[i]][["padj_taxa_res05linda_SMOTE"]] <- callDoAllModelsFromALLPCAs(all_pcas_linda_padj, name=paste0(i, "ConditionPadjLinda05SMOTE"),
                                                                                 metadata=this_metadata, vars2pca=c("Condition"), nfolds = NFOLDS,
                                                                                 xgboost_params = xgboost_params, catboost_params = catboost_params,
                                                                                 randomforest_params = randomforest_params,
                                                                                 do_smote = TRUE, smote_params = smote_params)

  all_model_results[[i]][["padj_taxa_res01linda_SMOTE"]] <- callDoAllModelsFromALLPCAs(all_pcas_adj01_linda, name=paste0(i, "ConditionPadjLinda01SMOTE"),
                                                                                 metadata=this_metadata, vars2pca=c("Condition"), nfolds = NFOLDS,
                                                                                 xgboost_params = xgboost_params, catboost_params = catboost_params,
                                                                                 randomforest_params = randomforest_params,
                                                                                 do_smote = TRUE, smote_params = smote_params)

  PCs <- all_model_results[[i]][["padj_taxa_res"]]$varnames
  modelo_svm <- all_model_results[[i]][["padj_taxa_res"]]$models$`SVM-linear`$mod_noscale
  all_model_results[[i]][["padj_taxa_res_indiv"]] <- callDoAllModelsFromALLPCAsOriginalVars(all_pcas_adj, PCs, modelo_svm,
                                                                                            df2pca,
                                                                                            paste0(i, "_ConditionPadjIndiv"),
                                                                                            vars2pca=c("Condition"), this_metadata,
                                                                                            daa_all[[i]]$resdf,
                                                                                            topns = c(5, 10, 20), nfolds = NFOLDS,
                                                                                            xgboost_params = xgboost_params, catboost_params = catboost_params,
                                                                                            randomforest_params = randomforest_params)

  all_model_results[[i]][["padj_taxa_res_indiv_SMOTE"]] <- callDoAllModelsFromALLPCAsOriginalVars(all_pcas_adj, PCs, modelo_svm,
                                                                                            df2pca,
                                                                                            paste0(i, "_ConditionPadjIndivSMOTE"),
                                                                                            vars2pca=c("Condition"), this_metadata,
                                                                                            daa_all[[i]]$resdf,
                                                                                            topns = c(5, 10, 20), nfolds = NFOLDS,
                                                                                            xgboost_params = xgboost_params, catboost_params = catboost_params,
                                                                                            randomforest_params = randomforest_params,
                                                                                            do_smote = TRUE, smote_params = smote_params)

  all_model_results[[i]]$metadata <- sample_data(phobj) %>% data.frame
  opt <- restaurar(opt)
#}
opt <- restaurar(opt)
save(all_model_results, file=paste0(opt$out, OUTDNAME, "/all_model_results_L1O.RData"))
#load(file=paste0(opt$out, "PredictDAA/all_model_results.RData"))

#Integrate
opt$out <- paste0(opt$out, OUTDNAME, "/")
makeLinePlotComparingPhobjs(all_model_results, opt)


## Compare with Bacteria in componets

for(orderby in c("Accuracy_l1out", "Kappa_l1out", "BalancedAccuracy_l1out", "AUC_l1out")){

  cat(orderby, "\n\n")
  walk(names(all_model_results), makeLinePlotComparingSamePhobjModels,
       all_model_results, opt, plot_extra=FALSE, name=paste0("orderBy_", orderby, "_normal"),
       order_by_measure = orderby,
       sel_method_name="PCA")

  walk(names(all_model_results), makeLinePlotComparingSamePhobjModels,
       get_pcnames_from = "padj_taxa_res_SMOTE",
       all_model_results, opt, plot_extra=FALSE, name=paste0("orderBy_", orderby, "_SMOTEonly"),
       from_smote=TRUE,
       order_by_measure = orderby,
       sel_method_name="PCA")

  walk(names(all_model_results), makeLinePlotComparingSamePhobjModels,
       all_model_results, opt, plot_extra=FALSE, name=paste0("orderBy_", orderby, "_SMOTEComp"),
       plot_normal_with_smote=TRUE, from_smote=FALSE, plot_indiv = FALSE,
       order_by_measure = orderby,
       sel_method_name="PCA")

  walk(names(all_model_results), makeLinePlotComparingSamePhobjModels,
       get_pcnames_from="padj_taxa_res05linda",
       all_model_results, opt, plot_extra=FALSE, name= paste0("orderBy_", orderby, "_LinDA05_SMOTEComp"),
       plot_normal_with_smote=TRUE, from_smote=FALSE, plot_indiv = FALSE,
       order_by_measure = orderby,
       sel_method_name= "PCA LinDA")

  walk(names(all_model_results), makeLinePlotComparingSamePhobjModels,
       get_pcnames_from="padj_taxa_res01linda",
       all_model_results, opt, plot_extra=FALSE, name=paste0("orderBy_", orderby, "_LinDA01_SMOTEComp"),
       plot_normal_with_smote=TRUE, from_smote=FALSE, plot_indiv = FALSE,
       order_by_measure = orderby,
       sel_method_name= "PCA LinDA")
}



outdir <- paste0(opt$out, "/all_performances/")
if(!dir.exists(outdir)) dir.create(outdir)


ns2 <- names(all_model_results$remove_tanda2)
ns2 <- ns2[grepl("_res", ns2, perl=T)]
for(n2 in ns2){

  if("allmodsum" %in% names(all_model_results$remove_tanda2[[n2]])){
  mmsum <- all_model_results$remove_tanda2[[n2]]$allmodsum %>%
    arrange(desc(BalancedAccuracy_l1out))
  }else{
    mmsum <- all_model_results$remove_tanda2[[n2]]$modummary %>%
      arrange(desc(BalancedAccuracy_l1out))
  }
  mname <- paste0(outdir, "remove_tanda2_", n2, "_modsummary.tsv")
  write_tsv(mmsum, file = mname)
}

