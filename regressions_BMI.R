library(caret)

plotSignificantVars <- function(metadata_pcadiv, sig_vars, yvar, opt, w=8, h=6){
  metadata_sigvars <- metadata_pcadiv %>% dplyr::select(all_of(c("sampleID", varname, yvar, sig_vars))) %>% 
    gather("variable", "value", sig_vars)
  
  g1<-ggplot(metadata_sigvars, aes(x=value, y=!!sym(yvar), col=status_c2, fill=status_c2))+
    facet_wrap(.~variable, scales = "free")+
    geom_smooth(alpha=0.6, method="lm")+
    geom_point(alpha=0.4, size=0.3)+
    theme_bw()
  ggsave(filename = paste0(opt$out, "regression_byGroup.pdf"), g1, width = w, height = h)
  
  g2 <- ggplot(metadata_sigvars, aes(x=value, y=!!sym(yvar)))+
    facet_wrap(.~variable, scales = "free")+
    geom_smooth(alpha=0.6, method="lm")+
    geom_point(alpha=0.4, size=0.3)+
    theme_bw()
  ggsave(filename = paste0(opt$out, "regression_merged.pdf"), g2, width = w, height = h)
  g3 <- ggplot(metadata_sigvars, aes(x=status_c2, y=value, col=status_c2, fill=status_c2))+
    facet_wrap(.~variable, scales = "free")+
    geom_violin(alpha=0.5) + 
    geom_boxplot(width = 0.2, col="black")+
    geom_jitter(alpha=0.4, size=0.3)+
    theme_bw()
  ggsave(filename = paste0(opt$out, "boxplot_byGroup.pdf"), g3, width = w, height = h)
  
}

plotModelDiagnosticsAll <- function(models, opt, name="models"){
  pdf(paste0(opt$out, name, "_diag_plots.pdf"))
  map(names(models), \(x){
    x2 <- models[[x]] %>% summary() 
    titstr <- paste0(x, " - R^2=", as.character(round(x2$r.squared, 2)))
    par(mfrow=c(2,2))
    plot(models[[x]], pch=16, cex=0.5)
    mtext(titstr, side = 3, line = -2, outer = TRUE)
    
    
  })
  dev.off()
}

plotDataVsExpected <- function(metadata_pcadiv, opt){
  gpred <- ggplot(metadata_pcadiv, aes(x=imc_01, y=predicted, col=status_c2, shape=correct))+geom_point()+theme_bw() + ggtitle("With microbiome PCs")
  gpred2 <- ggplot(metadata_pcadiv, aes(x=imc_01, y=predicted_neutral, col=status_c2, shape=correct_neutral))+geom_point()+theme_bw() + ggtitle("Without microbiome PCs")
  cw <- cowplot::plot_grid(plotlist = list(gpred, gpred2))
  pdf( paste0(opt$out, "bmi_t1_real_vs_predicted.pdf"), width = 12, height = 5)
  print(cw)
  dev.off()
}

confmat2df<-function(results, opt, name="results_confussionmatrix_all"){
  confmat_df <- map(results, \(x){
    
    msum <- x$bestmodel %>% summary
    msum2 <- x$neutralmodel %>% summary
    confmat <- x$conf_mat_all
    confmat2 <- x$conf_mat_noloss
    confmat_sim <- x$conf_mat_all_neutral
    confmat2_sim <- x$conf_mat_noloss_neutral
    auxdf <- data.frame(  fullmodel=as.character(msum$call)[2], 
                          simplemodel =as.character(msum$call)[2],
                          r_squared_fullmodel=msum$r.squared,
                          r_squared_simplemodel=msum2$r.squared,
                          #summary(mod_glm)$coefficients[2, 4], 
                          Accuracy=confmat$overall["Accuracy"],
                          Sensitivity = confmat$byClass["Sensitivity"],
                          Specificity = confmat$byClass["Specificity"],
                          PPV = confmat$byClass["Pos Pred Value"],
                          NPV = confmat$byClass["Neg Pred Value"],
                          
                          Accuracy_simplemod=confmat$overall["Accuracy"],
                          Sensitivity_simplemod = confmat$byClass["Sensitivity"],
                          Specificity_simplemod = confmat$byClass["Specificity"],
                          PPV_simplemod = confmat$byClass["Pos Pred Value"],
                          NPV_simplemod = confmat$byClass["Neg Pred Value"],
                          
                          Accuracy_noWeightLoss=confmat$overall["Accuracy"],
                          Sensitivity_noWeightLoss = confmat$byClass["Sensitivity"],
                          Specificity_noWeightLoss = confmat$byClass["Specificity"],
                          PPV_noWeightLoss = confmat$byClass["Pos Pred Value"],
                          NPV_noWeightLoss = confmat$byClass["Neg Pred Value"],
                          
                          Accuracy_simplemod_noWeightLoss=confmat$overall["Accuracy"],
                          Sensitivity_simplemod_noWeightLoss = confmat$byClass["Sensitivity"],
                          Specificity_simplemod_noWeightLoss = confmat$byClass["Specificity"],
                          PPV_simplemod_noWeightLoss = confmat$byClass["Pos Pred Value"],
                          NPV_simplemod_noWeightLoss = confmat$byClass["Neg Pred Value"]
    ) 
    auxdf
  }) %>% bind_rows
  write_tsv(confmat_df, file = paste0(opt$out, name, "_confMat.tsv"))
  return(confmat_df)
}

opt <- restaurar(opt)
outdir_base <- opt$out

load(paste0(opt$out, "DeSEQ2/DESEQ2_all.RData"))
load(paste0(opt$out, "DeSEQ2/DESEQ2_all.RData"))
load(paste0(opt$out, "PredictDAA/all_model_results.RData"))
load(paste0(opt$out, "phyloseq/phyloseq_all_list.RData"))


opt$out <- paste0(opt$out, "regressions")
if(!dir.exists(opt$out)) dir.create(opt$out)

alpha_indices <- c("Observed", "Chao1", "Shannon", "InvSimpson")

phseq_to_use <- names(all_phyloseq)[c(1,9,10)]
varname <- "status_c2"
interestvar <- "imc_00_log"
imcdiff_var <- "imc_diff"
dependent_variables <- c("imc_01_log")
independent_variables <- c("age_months_t0_log", "Sex", "Observed", "Chao1", "Shannon", "InvSimpson") #hospital

results <- list()
for(i in phseq_to_use){
  opt <- restaurar(opt)
  opt$out <- paste0(opt$out, "regressions/", i, "_IMC_abs/")
  if(!dir.exists(opt$out)) dir.create(opt$out)
  
  phobj <- all_phyloseq[[i]]
  df2pca <- if(is.null(daa_all[[i]]$vst_counts_df)){ daa_all[[i]]$norm_counts_df}else{ daa_all[[i]]$vst_counts_df }
  taxa_padj <- daa_all[[i]]$resdf %>% 
                dplyr::filter(padj <= opt$pval & abs(log2FoldChangeShrink) >= log2(opt$fc) ) %>% 
                pull(taxon)
  pca_padj <- all_model_results[[i]][["padj_taxa_pcas"]][[varname]]$pca$x %>% 
    data.frame %>% 
    rownames_to_column("sampleID") #%>% 
    #dplyr::filter(sampleID %in% names(df2pca))
  assertthat::assert_that(all(pca_padj$sampleID %in% names(df2pca)))
  this_metadata <- sample_data(phobj) %>% data.frame %>% dplyr::filter(sampleID %in% pca_padj$sampleID)
  this_metadata$imc_diff <- this_metadata$imc_01 - this_metadata$imc_00
  
  assertthat::assert_that(nrow(this_metadata) == nrow(pca_padj))
  metadata_pca <- merge(this_metadata, pca_padj, by="sampleID", all=TRUE)
  assertthat::assert_that(nrow(this_metadata) == nrow(metadata_pca))
  
  divtab_f <- paste0(outdir_base, "AlphaDiversity/", i, "_AlphaDiv.tsv")
  divtab <- read_tsv(divtab_f) %>% dplyr::select(sampleID, all_of(alpha_indices))
  assertthat::assert_that(all(metadata_pca$sampleID %in% divtab$sampleID))
  
  metadata_pcadiv <- merge(metadata_pca, divtab, by="sampleID", all.x=T, all.y=F)
  
  metadata_pcadiv <- metadata_pcadiv %>% 
    dplyr::mutate(across(all_of(c("age_months_t0", "imc_00", "nreads_filtPhylum", "imc_01")), .fns=list( "log"=log)))
  write_tsv(metadata_pcadiv, file = paste0(opt$out, "metadata_merged_pca_diversity.tsv"))
  pcnames <- pca_padj %>% select(-sampleID) %>% names
  
  ## Model with IMC_t1
  models1 <- makeLinearModelsSingleVariable(metadata_pcadiv, interestvar = interestvar, 
                                            extravars = c(independent_variables, pcnames), 
                                            alpha_indices = dependent_variables, 
                                            combos=1,
                                            outdir = opt$out, name = paste0(i, "_Regress_linMod1var") )
  sig_vars <- models1$single_anovas %>% filter(`Pr(>F)` < 0.05) %>% rownames
  sig_vars <- sig_vars[sig_vars %in% names(metadata_pcadiv)]
  
  #write_tsv(models1$single_anovas, file = paste0(opt$out, i, "_single_anovas.tsv"))
  plotSignificantVars(metadata_pcadiv, sig_vars = sig_vars, yvar = dependent_variables[1], opt = opt, w=8, h=6)
  plotModelDiagnosticsAll(models1$models, opt, "single_var_model")
  
  models_all <- map(sig_vars, \(x){
    extravars <- sig_vars[sig_vars != x]
    models2 <- makeLinearModelsSingleVariable(metadata_pcadiv, interestvar = x, 
                                            extravars = extravars, 
                                            alpha_indices = dependent_variables, 
                                            combos=1:length(extravars),
                                            outdir = opt$out, name = paste0(i, "_", x, "_Regress_linModManyVars") )
    models2$anovas %>% mutate(tested_var = x) %>% dplyr::select(tested_var, everything())
  }
  ) %>% bind_rows()
  write_tsv(models_all, file = paste0(opt$out, i, "_all_anovas_signif_variables_withHospital.tsv"))
  #By hand I checked that IMC_t0, PC3, PC13 and PC17 are significant in the presence of each other
  # Observed diversity is also significant but not after adding PC3. Therefore, I adjust the following model:

  finalmod <- lm(imc_01_log ~ imc_00_log + PC3 + PC13 + PC17, data=metadata_pcadiv)
  neutralmod <- lm(imc_01_log ~ imc_00_log, data=metadata_pcadiv)
  plotModelDiagnosticsAll(list("imc_01_log ~ imc_00_log + PC3 + PC13 + PC17" = finalmod,
                               "imc_01_log ~ imc_00_log" = neutralmod
                               ), opt, "best_model")

  levs <- unique(metadata_pcadiv$status_c2)
  metadata_pcadiv <- metadata_pcadiv %>% dplyr::mutate(
    predicted_log = predict(finalmod),
    predicted = exp(predicted_log),
    predicted_status = ifelse(predicted < lower_bound_t1, "Insufficient gain",
                             ifelse(predicted > upper_bound_t1, "Excessive gain", "Normal") 
    ),
    predicted_neutral_log = predict(neutralmod),
    predicted_neutral = exp(predicted_neutral_log),
    predicted_neutral_status = ifelse(predicted_neutral < lower_bound_t1, "Insufficient gain",
                              ifelse(predicted_neutral > upper_bound_t1, "Excessive gain", "Normal") 
    ),
    correct =  ifelse(predicted_status == status_c2, "Right", "Wrong"),
    correct_neutral =  ifelse(predicted_neutral_status == status_c2, "Right", "Wrong")
  ) %>% 
    dplyr::mutate(across(all_of(c("predicted_status", "predicted_neutral_status", "status_c2")), \(x)factor(x, levels=levs)))
  write_tsv(metadata_pcadiv, file = paste0(opt$out, "metadata_merged_pca_diversity_predict.tsv"))
  plotDataVsExpected(metadata_pcadiv, opt)
  metadata_noloss <- metadata_pcadiv %>% filter(status_c2 != "Insufficient gain")
  results[[i]] <- list()
  results[[i]]$data <- metadata_pcadiv
  results[[i]]$bestmodel <-  finalmod
  results[[i]]$neutralmodel <-  neutralmod
  results[[i]]$conf_mat_all <- confusionMatrix(metadata_pcadiv$status_c2, metadata_pcadiv$predicted_status)
  results[[i]]$conf_mat_all_neutral <- confusionMatrix(metadata_pcadiv$status_c2, metadata_pcadiv$predicted_neutral_status)
  results[[i]]$conf_mat_noloss <- confusionMatrix(metadata_noloss$status_c2, metadata_noloss$predicted_status)
  results[[i]]$conf_mat_noloss_neutral <- confusionMatrix(metadata_noloss$status_c2, metadata_noloss$predicted_neutral_status)

}

opt <- restaurar(opt)
save(file = paste0(opt$out, "regressions/model_results_predRaBMI.RData"), results)

opt$out <- paste0(opt$out, "regressions/")
confmat_df <- confmat2df(results, opt, "BMI_raw_all")
opt <- restaurar(opt)


######################################################################
## Now predict with IMC Diff
results2 <- list()

independent_variables <- c(interestvar, "age_months_t0_log", "Sex", "Observed", "Chao1", "Shannon", "InvSimpson") 
for(i in phseq_to_use){
  opt <- restaurar(opt)
  opt$out <- paste0(opt$out, "regressions/", i, "_IMC_diff/")
  if(!dir.exists(opt$out)) dir.create(opt$out)
  results[[i]] <- list()
  phobj <- all_phyloseq[[i]]
  df2pca <- if(is.null(daa_all[[i]]$vst_counts_df)){ daa_all[[i]]$norm_counts_df}else{ daa_all[[i]]$vst_counts_df }
  taxa_padj <- daa_all[[i]]$resdf %>% 
    dplyr::filter(padj <= opt$pval & abs(log2FoldChangeShrink) >= log2(opt$fc) ) %>% 
    pull(taxon)
  pca_padj <- all_model_results[[i]][["padj_taxa_pcas"]][[varname]]$pca$x %>% 
    data.frame %>% 
    rownames_to_column("sampleID") #%>% 
  #dplyr::filter(sampleID %in% names(df2pca))
  assertthat::assert_that(all(pca_padj$sampleID %in% names(df2pca)))
  this_metadata <- sample_data(phobj) %>% data.frame %>% dplyr::filter(sampleID %in% pca_padj$sampleID)
  this_metadata$imc_diff <- this_metadata$imc_01 - this_metadata$imc_00
  
  assertthat::assert_that(nrow(this_metadata) == nrow(pca_padj))
  metadata_pca <- merge(this_metadata, pca_padj, by="sampleID", all=TRUE)
  assertthat::assert_that(nrow(this_metadata) == nrow(metadata_pca))
  
  divtab_f <- paste0(outdir_base, "AlphaDiversity/", i, "_AlphaDiv.tsv")
  divtab <- read_tsv(divtab_f) %>% dplyr::select(sampleID, all_of(alpha_indices))
  assertthat::assert_that(all(metadata_pca$sampleID %in% divtab$sampleID))
  
  metadata_pcadiv <- merge(metadata_pca, divtab, by="sampleID", all.x=T, all.y=F)
  
  metadata_pcadiv <- metadata_pcadiv %>% 
    dplyr::mutate(across(all_of(c("age_months_t0", "imc_00", "nreads_filtPhylum", "imc_01")), .fns=list( "log"=log)))
  write_tsv(metadata_pcadiv, file = paste0(opt$out, "metadata_merged_pca_diversity.tsv"))
  pcnames <- pca_padj %>% select(-sampleID) %>% names
  
  ## Model with IMC_t1
  models1 <- makeLinearModelsSingleVariable(metadata_pcadiv, interestvar = imcdiff_var, 
                                            extravars = c(independent_variables, pcnames), 
                                            alpha_indices = imcdiff_var, 
                                            combos=1,
                                            outdir = opt$out, name = paste0(i, "_Regress_linMod1var") )
  sig_vars <- models1$single_anovas %>% filter(`Pr(>F)` < 0.05) %>% rownames
  
  #write_tsv(models1$single_anovas, file = paste0(opt$out, i, "_single_anovas.tsv"))
  plotSignificantVars(metadata_pcadiv, sig_vars = sig_vars, yvar = imcdiff_var, opt = opt, w=8, h=6)
  plotModelDiagnosticsAll(models1$models, opt, "single_var_model")
  
  models_all <- map(sig_vars, \(x){
    extravars <- sig_vars[sig_vars != x]
    models2 <- makeLinearModelsSingleVariable(metadata_pcadiv, interestvar = x, 
                                              extravars = extravars, 
                                              alpha_indices = imcdiff_var, 
                                              combos=1:length(extravars),
                                              outdir = opt$out, name = paste0(i, "_", x, "_Regress_linModManyVars") )
    models2$anovas %>% mutate(tested_var = x) %>% dplyr::select(tested_var, everything())
  }
  ) %>% bind_rows()
  write_tsv(models_all, file = paste0(opt$out, i, "_all_anovas_signif_variables_withHospital.tsv"))
  #By hand I checked that IMC_t0, PC3, PC13 and PC17 are significant in the presence of each other
  # Observed diversity is also significant but not after adding PC3. Therefore, I adjust the following model:
  finalmod <- lm(imc_diff ~ age_months_t0_log + PC3 + PC13, data=metadata_pcadiv)
  neutralmod <- lm(imc_diff ~ age_months_t0_log, data=metadata_pcadiv)
  plotModelDiagnosticsAll(list("imc_diff ~ age_months_t0_log + PC3 + PC13" = finalmod,
                               "imc_diff ~ age_months_t0_log" = neutralmod
  ), opt, "best_model")
  
  levs <- unique(metadata_pcadiv$status_c2)
  metadata_pcadiv <- metadata_pcadiv %>% dplyr::mutate(
    predicted_log = predict(finalmod),
    predicted = exp(predicted_log),
    predicted_status = ifelse(predicted < lower_bound_t1, "Insufficient gain",
                              ifelse(predicted > upper_bound_t1, "Excessive gain", "Normal") 
    ),
    predicted_neutral_log = predict(neutralmod),
    predicted_neutral = exp(predicted_neutral_log),
    predicted_neutral_status = ifelse(predicted_neutral < lower_bound_t1, "Insufficient gain",
                                      ifelse(predicted_neutral > upper_bound_t1, "Excessive gain", "Normal") 
    ),
    correct =  ifelse(predicted_status == status_c2, "Right", "Wrong"),
    correct_neutral =  ifelse(predicted_neutral_status == status_c2, "Right", "Wrong")
  ) %>% 
    dplyr::mutate(across(all_of(c("predicted_status", "predicted_neutral_status", "status_c2")), \(x)factor(x, levels=levs)))
  write_tsv(metadata_pcadiv, file = paste0(opt$out, "metadata_merged_pca_diversity_predict.tsv"))
  plotDataVsExpected(metadata_pcadiv, opt)
  metadata_noloss <- metadata_pcadiv %>% filter(status_c2 != "Insufficient gain")
  results2[[i]]$data <- metadata_pcadiv
  results2[[i]]$bestmodel <-  finalmod
  results2[[i]]$neutralmodel <-  neutralmod
  results2[[i]]$conf_mat_all <- confusionMatrix(metadata_pcadiv$status_c2, metadata_pcadiv$predicted_status)
  results2[[i]]$conf_mat_all_neutral <- confusionMatrix(metadata_pcadiv$status_c2, metadata_pcadiv$predicted_neutral_status)
  results2[[i]]$conf_mat_noloss <- confusionMatrix(metadata_noloss$status_c2, metadata_noloss$predicted_status)
  results2[[i]]$conf_mat_noloss_neutral <- confusionMatrix(metadata_noloss$status_c2, metadata_noloss$predicted_neutral_status)
  
}

opt <- restaurar(opt)
save(file = paste0(opt$out, "regressions/model_results_predDiffBMI.RData"), results2)

opt$out <- paste0(opt$out, "regressions/")
confmat_df <- confmat2df(results2, opt, "BMI_diff_all")
opt <- restaurar(opt)
