# Predict  
source(opt$predictive_functions)
#opt$out <- "/home/carmoma/Desktop/202311_DEPRESION/results_rstudio_v2_1/"
var2predict <- "status_c2"
vars2pca <- c("status_c2", "hospital", "Sex", "age_months_t0")
phseq_to_use <- names(daa_all)[c(1,3,4,6,9,10)]
opt <- restaurar(opt)                  
opt$out <- paste0(opt$out, "PredictDAA_onlyGain")
if(!dir.exists(opt$out)) dir.create(opt$out)
opt <- restaurar(opt)

all_model_results <- list()
for(i in phseq_to_use){
  cat("Doing Predictive models for: ", i, "\n")
  all_model_results[[i]] <- list()
  phobj <- all_phyloseq[[i]]
  outdir <- paste0(opt$out, "PredictDAA_onlyGain/", i, "/")
  opt$reserva <- opt$out
  opt$out <- outdir
  if(!dir.exists(opt$out)) dir.create(opt$out)
  
  taxa_padj <- daa_all[[i]]$resdf %>% dplyr::filter(padj <= opt$pval & abs(log2FoldChangeShrink) >= log2(opt$fc) ) %>% 
    pull(taxon)
  taxa_praw <- daa_all[[i]]$resdf %>% dplyr::filter(pvalue <= opt$pval & abs(log2FoldChangeShrink) >= log2(opt$fc) ) %>% 
    pull(taxon)
  df2pca <- if(is.null(daa_all[[i]]$vst_counts_df)){ daa_all[[i]]$norm_counts_df}else{ daa_all[[i]]$vst_counts_df }
  
  all_pcas_adj <- makeAllPCAs(phobj, df2pca, taxa_padj, vars2pca, opt, "DiffTaxaPadj")
  all_pcas_praw <- makeAllPCAs(phobj, df2pca, taxa_praw, vars2pca, opt, "DiffTaxaPraw")
  
  taxa_padj01 <- daa_all[[i]]$resdf %>% dplyr::filter(padj <= 0.01 & abs(log2FoldChangeShrink) >= log2(opt$fc) ) %>% 
    pull(taxon)
  taxa_padj001 <- daa_all[[i]]$resdf %>% dplyr::filter(padj <= 0.001 & abs(log2FoldChangeShrink) >= log2(opt$fc) ) %>% 
    pull(taxon)
  all_pcas_adj01 <- makeAllPCAs(phobj, df2pca, taxa_padj01, vars2pca, opt, "DiffTaxaPadj01")
  #all_pcas_adj001 <- makeAllPCAs(phobj, df2pca, taxa_padj001, vars2pca, opt, "DiffTaxaPadj001")
  
  this_metadata <- sample_data(phobj) %>% data.frame %>% dplyr::filter(sampleID %in% names(df2pca))
  all_model_results[[i]][["padj_taxa_taxa"]] <- taxa_padj
  all_model_results[[i]][["praw_taxa_taxa"]] <- taxa_praw
  all_model_results[[i]][["padj_taxa_pcas"]] <- all_pcas_adj
  all_model_results[[i]][["praw_taxa_pcas"]] <- all_pcas_praw
  all_model_results[[i]][["padj_taxa_res"]] <- callDoAllModelsFromALLPCAs(all_pcas_adj, name=paste0(i, "ConditionPadj"), metadata=this_metadata, vars2pca=var2predict)
  #all_model_results[[i]][["praw_taxa_res"]] <- callDoAllModelsFromALLPCAs(all_pcas_praw, name=paste0(i, "ConditionPraw"), metadata=this_metadata, vars2pca=var2predict)
  all_model_results[[i]][["padj_taxa_res01"]] <- callDoAllModelsFromALLPCAs(all_pcas_adj01, name=paste0(i, "ConditionPadj01"), metadata=this_metadata, vars2pca=var2predict)
  #all_model_results[[i]][["padj_taxa_res001"]] <- callDoAllModelsFromALLPCAs(all_pcas_adj001, name=paste0(i, "ConditionPadj001"), metadata=this_metadata, vars2pca=var2predict)
  
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

walk(names(all_model_results), makeLinePlotComparingSamePhobjModels, 
     all_model_results, opt)

## Plot boxplot PCs
pcBoxplots <- map(names(all_model_results), makePCsBoxplot, all_model_results, opt)
names(pcBoxplots) <- names(all_model_results)

## Plot barplot PCs and LFC
pcBarplots <- map(names(all_model_results), makePCBarplot, all_model_results, pcBoxplots, daa_all, opt, w=8, h=14)
names(pcBarplots) <- names(all_model_results)

# Plot KNN (best model)

phname <- "remove_tanda2"
predplots <- map(names(all_model_results), plotAllModelPredictions, all_model_results, opt)

# Make sure that it works
opt <- restaurar(opt)



