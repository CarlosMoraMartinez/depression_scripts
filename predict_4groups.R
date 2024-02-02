#####################################
## Predict DEPR + OBESITY

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

opt$out <- paste0(opt$out, "PredictDAA_multiclass")
if(!dir.exists(opt$out)) dir.create(opt$out)
opt <- restaurar(opt)

all_model_multi_results <- list()
for(i in phseq_to_use){
  cat("Doing Predictive models for: ", i, "\n")
  all_model_results[[i]] <- list()
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
  taxa_cov_praw <- daa_all[[i]]$resdf %>% filter_taxa_praw
  taxa_padj01 <- daa_all[[i]]$resdf %>% filter_taxa_padj(plim=0.01)
  taxa_padj001 <- daa_all[[i]]$resdf %>% filter_taxa_padj(plim=0.001)
  taxa_cov_padj01 <- daa_all_corrected_only[[i]]$BMI_log %>% filter_taxa_padj(plim=0.01)
  taxa_cov_padj001 <- daa_all_corrected_only[[i]]$BMI_log %>% filter_taxa_padj(plim=0.001)
  
  taxa_padj_corr <- daa_all_corrected[[i]]$BMI_log %>% filter_taxa_padj 
  taxa_padj_corr_01 <- daa_all_corrected[[i]]$BMI_log %>% filter_taxa_padj(plim=0.01)
  taxa_padj_corr_001 <- daa_all_corrected[[i]]$BMI_log %>% filter_taxa_padj(plim=0.001)
  imcfname <- paste0(opt$reserva, "DESeq2_ControlVars/DeSEQ2/", i, "_BMI_log/",i, "_BMI_log_BMI_log_DAAshrinkNormal.tsv" )
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
  
  
  all_model_results[[i]] <- c(all_models_this, all_pcalists, taxa_list)
  all_model_results[[i]]$metadata <- sample_data(phobj_filt) %>% data.frame
  
  opt <- restaurar(opt)
}
opt <- restaurar(opt)
save(all_model_results, file=paste0(opt$out, "PredictDAA_multiclass/all_model_results.RData"))
#load(file=paste0(opt$out, "PredictDAA/all_model_results.RData"))

#for(k in names(all_model_results[[i]]$modresults)) {cat("\n\n", k);all_model_results[[i]]$modresults[[k]]$modummary %>% head(2) %>% print}

#Integrate
opt$out <- paste0(opt$out, "PredictDAA_multiclass/")
makeLinePlotComparingPhobjs( all_model_results, opt, models_name1 = "DiffTaxaPadjBase", models_name2 = "DiffTaxaPadjCov05")
## Compare with Bacteria in componets
condnames <- names(all_model_results$remove_tanda2)[c(1,6,7,9,10)]
makeLinePlotComparingSamePhobjModels_Cov(phname, condnames, all_model_results, "TaxaGroups_trim", opt, w=8, h=5)

walk(names(all_model_results), \(x)makeLinePlotComparingSamePhobjModels(x, all_model_results, 
                                                                        opt, w=6, h=8, get_pcnames_from = "DiffTaxaPadjCov05"))
## Plot boxplot PCs
pca2use <- "DiffTaxaPadjCov05"
pcBoxplots <- map(names(all_model_results), \(x){
  makePCsBoxplot(x, all_model_results, opt, 
                 get_pcnames_from = pca2use,
                 pca_name =  paste0("PCA_", pca2use), 
                 varname = "Depr_and_Ob",
                 w=12, h=8)
})
names(pcBoxplots) <- names(all_model_results)

## Plot barplot PCs and LFC
pcBarplots <- map(names(all_model_results), makePCBarplot, all_model_results, pcBoxplots, daa_all, opt,
                  get_pcnames_from = pca2use,
                  pca_name =  paste0("PCA_", pca2use), 
                  varname = "Depr_and_Ob",
                  w=12, h=20)
names(pcBarplots) <- names(all_model_results)

# Plot KNN (best model)

phname <- "remove_tanda2"
predplots <- map(names(all_model_results), plotAllModelPredictions, all_model_results, opt,
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
  all_model_results[[i]] <- list()
  phobj <- all_phyloseq[[i]]
  
  sample_data(phobj)[[var2add]] <- editObesity(sample_data(phobj)[[var2add]])
  phobj <- combineClasses(phobj, interestvar, var2add, newname)
  samples <- sample_data(phobj)$sampleID[! is.na(sample_data(phobj)[, var2add])]
  phobj_filt <- phyloseq::prune_samples(samples, phobj)
  
  outdir <- paste0(opt$out, "PredictDAA_multiclass/", i, "/")
  opt$reserva <- opt$out
  opt$out <- outdir
  if(!dir.exists(opt$out)) dir.create(opt$out)
  
  ## Get data for PCA 
  df2pca <- if(is.null(daa_all[[i]]$vst_counts_df)){ daa_all[[i]]$norm_counts_df}else{ daa_all[[i]]$vst_counts_df }
  
  ## Get taxa 
  taxa_padj <- daa_all[[i]]$resdf %>% filter_taxa_padj
  taxa_praw <- daa_all[[i]]$resdf %>% filter_taxa_praw
  taxa_cov_padj <- daa_all_corrected_only[[i]]$ob_o_sobrepeso %>% filter_taxa_padj
  taxa_cov_praw <- daa_all[[i]]$resdf %>% filter_taxa_praw
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
  
  
  all_model_results[[i]] <- c(all_models_this, all_pcalists, taxa_list)
  all_model_results[[i]]$metadata <- sample_data(phobj_filt) %>% data.frame
  
  opt <- restaurar(opt)
}
opt <- restaurar(opt)
save(all_model_results, file=paste0(opt$out, "PredictDAA_multiclass_fromOb/all_model_results.RData"))
#load(file=paste0(opt$out, "PredictDAA/all_model_results.RData"))

#for(k in names(all_model_results[[i]]$modresults)) {cat("\n\n", k);all_model_results[[i]]$modresults[[k]]$modummary %>% head(2) %>% print}

#Integrate
opt$out <- paste0(opt$out, "PredictDAA_multiclass/")
makeLinePlotComparingPhobjs( all_model_results, opt, models_name1 = "DiffTaxaPadjBase", models_name2 = "DiffTaxaPadjCov05")
## Compare with Bacteria in componets
condnames <- names(all_model_results$remove_tanda2)[c(1,6,7,9,10)]
makeLinePlotComparingSamePhobjModels_Cov(phname, condnames, all_model_results, "TaxaGroups_trim", opt, w=8, h=5)

walk(names(all_model_results), \(x)makeLinePlotComparingSamePhobjModels(x, all_model_results, 
                                                                        opt, w=6, h=8, get_pcnames_from = "DiffTaxaPadjCov05"))
## Plot boxplot PCs
pca2use <- "DiffTaxaPadjCov05"
pcBoxplots <- map(names(all_model_results), \(x){
  makePCsBoxplot(x, all_model_results, opt, 
                 get_pcnames_from = pca2use,
                 pca_name =  paste0("PCA_", pca2use), 
                 varname = "Depr_and_Ob",
                 w=12, h=8)
})
names(pcBoxplots) <- names(all_model_results)

## Plot barplot PCs and LFC
pcBarplots <- map(names(all_model_results), makePCBarplot, all_model_results, pcBoxplots, daa_all, opt,
                  get_pcnames_from = pca2use,
                  pca_name =  paste0("PCA_", pca2use), 
                  varname = "Depr_and_Ob",
                  w=12, h=20)
names(pcBarplots) <- names(all_model_results)

# Plot KNN (best model)

phname <- "remove_tanda2"
predplots <- map(names(all_model_results), plotAllModelPredictions, all_model_results, opt,
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
