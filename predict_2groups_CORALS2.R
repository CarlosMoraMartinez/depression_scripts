# Predict  
source(opt$predictive_functions)


#opt$out <- "/home/carmoma/Desktop/202311_DEPRESION/results_rstudio_v2_1/"
var2predict <- "status_c2"
vars2pca <- c("status_c2", "Category_T0", "hospital", "Sex", "edad_00_meses")
phseq_to_use <- names(daa_all)[c(2,7,9)]
opt <- restaurar(opt)                  
opt$out <- paste0(opt$out, "PredictDAA_onlyGain")
if(!dir.exists(opt$out)) dir.create(opt$out)
opt <- restaurar(opt)

all_model_results <- list()

food_variables<- names(s_meta)[4:29]
extra_variables <- c("edad_00_meses_logscale", "af_extraesc_m_00_logscale")

NFOLDS <- 10

for(i in phseq_to_use){
  cat("Doing Predictive models for: ", i, "\n")
  all_model_results[[i]] <- list()
  phobj <- all_phyloseq[[i]]
  outdir <- paste0(opt$out, "PredictDAA_onlyGain/", i, "/")
  opt$reserva <- opt$out
  opt$out <- outdir
  if(!dir.exists(opt$out)) dir.create(opt$out)
  
  taxa_padj <- daa_all[[i]]$all_contrasts$status_c2_Normal_vs_Excessive.gain$resdf %>% 
    dplyr::filter(padj <= opt$pval & abs(log2FoldChangeShrink) >= log2(opt$fc) ) %>% 
    pull(taxon)
  taxa_praw <- daa_all[[i]]$resdf %>% pull(taxon)
  
  df2pca <- if(is.null(daa_all[[i]]$vst_counts_df)){ daa_all[[i]]$norm_counts_df}else{ daa_all[[i]]$vst_counts_df }
  
  all_pcas_adj <- makeAllPCAs(phobj, df2pca, taxa_padj, vars2pca, opt, "DiffTaxaPadj")
  all_pcas_praw <- makeAllPCAs(phobj, df2pca, taxa_praw, vars2pca, opt, "DiffTaxaPraw")
  
  taxa_padj01 <- daa_all[[i]]$all_contrasts$status_c2_Normal_vs_Excessive.gain$resdf %>% 
    dplyr::filter(padj <= 0.01 & abs(log2FoldChangeShrink) >= log2(opt$fc) ) %>% 
    pull(taxon)
  taxa_padj001 <-  daa_all[[i]]$all_contrasts$status_c2_Normal_vs_Excessive.gain$resdf %>% 
    dplyr::filter(padj <= 0.001 & abs(log2FoldChangeShrink) >= log2(opt$fc) ) %>% 
    pull(taxon)
  all_pcas_adj01 <- makeAllPCAs(phobj, df2pca, taxa_padj01, vars2pca, opt, "DiffTaxaPadj01")
  all_pcas_adj001 <- makeAllPCAs(phobj, df2pca, taxa_padj001, vars2pca, opt, "DiffTaxaPadj001")
  
  this_metadata <- sample_data(phobj) %>% data.frame %>% 
    dplyr::filter(sampleID %in% names(df2pca)) %>% 
    dplyr::filter(status_c2 != "Insufficient gain") %>% 
    dplyr::filter(status_c2 != "initially_overweight")
  
  
  food_PCA <- make_meta_PCA(this_metadata, food_variables, 
                            var2predict,
                            outdir,
                            make_log=TRUE,
                            name="foodVars")
  
  pcamat <- food_PCA$pca$x %>% 
    as.data.frame %>% 
    rownames_to_column("sample")
  names(pcamat)[2:ncol(pcamat)] <- paste("Food", names(pcamat)[2:ncol(pcamat)], sep="")
  byy <- join_by(sampleID == sample) 
  
  this_metadata <- this_metadata %>% inner_join(pcamat, by = byy)
  
  this_metadata <- this_metadata %>% 
    dplyr::mutate(edad_00_meses_logscale = scale(log(edad_00_meses+1)),
                  af_extraesc_m_00_logscale = scale(log(af_extraesc_m_00+1)))
  
  meta_predictvars1 <- c(names(pcamat)[2:ncol(pcamat)], extra_variables[1])
  meta_predictvars2 <- c(names(pcamat)[2:ncol(pcamat)], extra_variables[2])
  
  this_metadata <- this_metadata %>%
    mutate_at(food_variables, ~ scale(log(. +1)))
  
  all_model_results[[i]][["padj_taxa_taxa"]] <- taxa_padj
  all_model_results[[i]][["praw_taxa_taxa"]] <- taxa_praw
  all_model_results[[i]][["padj_taxa_pcas"]] <- all_pcas_adj
  all_model_results[[i]][["praw_taxa_pcas"]] <- all_pcas_praw
  all_model_results[[i]][["padj_taxa_res"]] <- callDoAllModelsFromALLPCAs(all_pcas_adj, name=paste0(i, "ConditionPadj"), 
                                                                          metadata=this_metadata, vars2pca=var2predict,
                                                                          variable_plim = 0.05, nfolds = NFOLDS)
  all_model_results[[i]][["praw_taxa_res"]] <- callDoAllModelsFromALLPCAs(all_pcas_praw, name=paste0(i, "ConditionPraw"), 
                                                                          metadata=this_metadata, vars2pca=var2predict,
                                                                          variable_plim = 0.05, nfolds = NFOLDS)
  all_model_results[[i]][["padj_taxa_res_foodAge"]] <- callDoAllModelsFromALLPCAs(all_pcas_adj, name=paste0(i, "ConditionPadjFoodAge"), 
                                                                          metadata=this_metadata, vars2pca=var2predict,
                                                                          variable_plim = 0.05,
                                                                          meta_vars = meta_predictvars1, nfolds = NFOLDS)
  all_model_results[[i]][["praw_taxa_res_foodAge"]] <- callDoAllModelsFromALLPCAs(all_pcas_praw, name=paste0(i, "ConditionPrawFoodAge"), 
                                                                          metadata=this_metadata, vars2pca=var2predict,
                                                                          variable_plim = 0.05,
                                                                          meta_vars = meta_predictvars1, nfolds = NFOLDS)
  all_model_results[[i]][["padj_taxa_res_foodExer"]] <- callDoAllModelsFromALLPCAs(all_pcas_adj, name=paste0(i, "ConditionPadjFoodAge"), 
                                                                               metadata=this_metadata, vars2pca=var2predict,
                                                                               variable_plim = 0.05,
                                                                               meta_vars = meta_predictvars2, nfolds = NFOLDS)
  all_model_results[[i]][["praw_taxa_res_foodExer"]] <- callDoAllModelsFromALLPCAs(all_pcas_praw, name=paste0(i, "ConditionPrawFoodAge"), 
                                                                               metadata=this_metadata, vars2pca=var2predict,
                                                                               variable_plim = 0.05,
                                                                               meta_vars = meta_predictvars2, nfolds = NFOLDS)
  all_model_results[[i]][["padj_taxa_res_food"]] <- callDoAllModelsFromALLPCAs(all_pcas_adj, name=paste0(i, "ConditionPadjFood"), 
                                                                               metadata=this_metadata, vars2pca=var2predict,
                                                                               variable_plim = 0.05,
                                                                               meta_vars = meta_predictvars[1:length(food_variables)], nfolds = NFOLDS)
  all_model_results[[i]][["praw_taxa_res_food"]] <- callDoAllModelsFromALLPCAs(all_pcas_praw, name=paste0(i, "ConditionPrawFood"), 
                                                                               metadata=this_metadata, vars2pca=var2predict,
                                                                               variable_plim = 0.05,
                                                                               meta_vars = meta_predictvars[1:length(food_variables)], nfolds = NFOLDS)
  
  all_model_results[[i]][["padj_taxa_res_foodraw1"]] <- callDoAllModelsFromALLPCAs(all_pcas_adj, name=paste0(i, "ConditionPadjFoodraw1"), 
                                                                               metadata=this_metadata, vars2pca=var2predict,
                                                                               variable_plim = 1,
                                                                               meta_vars = food_variables, nfolds = NFOLDS)
  all_model_results[[i]][["praw_taxa_res_foodraw1"]] <- callDoAllModelsFromALLPCAs(all_pcas_praw, name=paste0(i, "ConditionPrawFoodraw1"), 
                                                                               metadata=this_metadata, vars2pca=var2predict,
                                                                               variable_plim = 1,
                                                                               meta_vars = food_variables, nfolds = NFOLDS)
  all_model_results[[i]][["padj_taxa_res_foodraw05"]] <- callDoAllModelsFromALLPCAs(all_pcas_adj, name=paste0(i, "ConditionPadjFoodraw05"), 
                                                                                   metadata=this_metadata, vars2pca=var2predict,
                                                                                   variable_plim = 1,
                                                                                   meta_vars = food_variables, nfolds = NFOLDS)
  all_model_results[[i]][["praw_taxa_res_foodraw05"]] <- callDoAllModelsFromALLPCAs(all_pcas_praw, name=paste0(i, "ConditionPrawFoodraw05"), 
                                                                                   metadata=this_metadata, vars2pca=var2predict,
                                                                                   variable_plim = 1,
                                                                                   meta_vars = food_variables, nfolds = NFOLDS)
  
  #all_model_results[[i]][["padj_taxa_res01"]] <- callDoAllModelsFromALLPCAs(all_pcas_adj01, name=paste0(i, "ConditionPadj01"), 
  #                                                                          metadata=this_metadata, vars2pca=var2predict,
  #                                                                          variable_plim = 0.05)
  #all_model_results[[i]][["padj_taxa_res001"]] <- callDoAllModelsFromALLPCAs(all_pcas_adj001, name=paste0(i, "ConditionPadj001"),
  #                                                                           metadata=this_metadata, vars2pca=var2predict,
  #                                                                           variable_plim = 0.05)
  
  # PCs <- all_model_results[[i]][["padj_taxa_res"]]$varnames
  # modelo_svm <- all_model_results[[i]][["padj_taxa_res"]]$models$`SVM-linear`$mod_noscale
  # all_model_results[[i]][["padj_taxa_res_indiv"]] <- callDoAllModelsFromALLPCAsOriginalVars(all_pcas_adj, PCs, modelo_svm, 
  #                                                                                           df2pca, 
  #                                                                                           paste0(i, "_ConditionPadjIndiv"), 
  #                                                                                           vars2pca=c("Condition"), s_meta,
  #                                                                                           daa_all[[i]]$resdf, 
  #                                                                                           topns = c(5, 10, 20))
  
  all_model_results[[i]]$metadata <-this_metadata
  opt <- restaurar(opt)
}
opt <- restaurar(opt)

save(all_model_results, file=paste0(opt$out, "PredictDAA_onlyGain/all_model_results.RData"))
#load(file=paste0(opt$out, "PredictDAA/all_model_results.RData"))

#Integrate
opt$out <- paste0(opt$out, "PredictDAA_onlyGain/")

makeLinePlotComparingPhobjs(all_model_results, opt, models_name1 = "padj_taxa_res", models_name2 = "padj_taxa_res01")
## Compare with Bacteria in componets

n2plott <- names(all_model_results[[1]])
n2plott <- n2plott[grepl("_res", n2plott)]
  

walk(names(all_model_results), makeLinePlotComparingSamePhobjModels, 
     all_model_results, opt, 8, 12, n2plott)

allrestable <- data.frame()
for(nn in n2plott){ for(ph in names(all_model_results)){
  aux <- all_model_results[[ph]][[nn]]$modummary %>% 
    dplyr::mutate(Phobj = ph, Vars = nn, Cond = paste(ph, nn, sep="_") ) %>% 
    dplyr::select(Cond, Phobj, Vars, everything())
  allrestable <- rbind(allrestable, aux)
  
}
}
write_tsv(allrestable, file = paste0(opt$out, "All_predict_results.tsv"))


measures <- c("Accuracy_l1out", "Kappa_l1out", "Sensitivity_l1out", "Specificity_l1out", "Recall_l1out")

vv <- measures[1]




## Plot boxplot PCs
pcBoxplots <- map(names(all_model_results), makePCsBoxplot, all_model_results, opt, "padj_taxa_res", "padj_taxa_pcas", "status_c2", 6, 8)
names(pcBoxplots) <- names(all_model_results)

## Plot barplot PCs and LFC
pcBarplots <- map(names(all_model_results), makePCBarplot, all_model_results, pcBoxplots, daa_all, opt, "padj_taxa_res", "padj_taxa_pcas", "status_c2", w=10, h=10)
names(pcBarplots) <- names(all_model_results)

# Plot KNN (best model)

phname <- "remove_tanda2"
predplots <- map(names(all_model_results), plotAllModelPredictions, all_model_results, opt)

# Make sure that it works
opt <- restaurar(opt)



