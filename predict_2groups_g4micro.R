# Predict
detach("package:G4Micro", unload = TRUE)
library(G4Micro)
#source(opt$predictive_functions)
#opt$out <- "/home/carmoma/Desktop/202311_DEPRESION/results_rstudio_v2_1/"




load("/home/carlos/Documentos/CORALS/results_rstudio/results_Abril25_2/foodPCA/phyloseq_list_foodPCA_withNMF.RData")
load("/home/carlos/Documentos/CORALS/results_rstudio/results_Agosto25_1/DeSEQ2/DESEQ2_all.RData")

var2predict <- "status_c2"
levs2pred <- c("Normal", "Excessive gain")
vars2pca <- c("status_c2", "hospital", "Sex", "age_months_T0")
phseq_to_use <- "remove_tanda2" #names(daa_all) #[c(1,3,4,6,9,10)]
patnames <- c("Mediterranean", "Preprocessed", "Western")
foodnames <- sample_data(all_phyloseq$remove_tanda2) %>% data.frame() %>% names()
foodnames <- foodnames[grep("_clr", foodnames)]
foodnames_prop <- sample_data(all_phyloseq$remove_tanda2) %>% data.frame() %>% names()
foodnames_prop <- foodnames_prop[grep("_prop", foodnames_prop)]

opt <- restaurar(opt)
opt$out <- paste0(opt$out, "PredictDAA_onlyGain2")
if(!dir.exists(opt$out)) dir.create(opt$out)
opt <- restaurar(opt)

all_model_results <- list()
i <- phseq_to_use
## for(i in phseq_to_use){

  cat("Doing Predictive models for: ", i, "\n")
  all_model_results[[i]] <- list()
  phobj <- all_phyloseq[[i]]
  s_meta <- sample_data(phobj) %>% data.frame()
  s2use <- s_meta %>%
    filter(! is.na(status_c2)) %>%
    filter(status_c1=="normal") %>%
    filter(status_c2 != "NaN") %>%
    filter(status_c2 != "Insufficient gain") %>% pull(sampleID)
  phobj_filt <- phyloseq::prune_samples(s2use,phobj)
  outdir <- paste0(opt$out, "PredictDAA_onlyGain/", i, "/")
  opt$reserva <- opt$out
  opt$out <- outdir
  if(!dir.exists(opt$out)) dir.create(opt$out)

  taxa_padj <- daa_all[[i]]$resdf %>% dplyr::filter(padj <= opt$pval & abs(log2FoldChangeShrink) >= log2(opt$fc) ) %>%
    pull(taxon)
  taxa_praw <- daa_all[[i]]$resdf %>% dplyr::filter(pvalue <= opt$pval & abs(log2FoldChangeShrink) >= log2(opt$fc) ) %>%
    pull(taxon)
  df2pca <- if(is.null(daa_all[[i]]$vst_counts_df)){ daa_all[[i]]$norm_counts_df}else{ daa_all[[i]]$vst_counts_df }
  df2pca <- df2pca %>% dplyr::select(gene, all_of(s2use))
  all_pcas_adj <- G4Micro::makeAllPCAs(phobj_filt, df2pca, taxa_padj, vars2pca, opt, "DiffTaxaPadj")
  all_pcas_praw <- G4Micro::makeAllPCAs(phobj_filt, df2pca, taxa_praw, vars2pca, opt, "DiffTaxaPraw")

  taxa_padj01 <- daa_all[[i]]$resdf %>% dplyr::filter(padj <= 0.01 & abs(log2FoldChangeShrink) >= log2(opt$fc) ) %>%
    pull(taxon)
  taxa_padj001 <- daa_all[[i]]$resdf %>% dplyr::filter(padj <= 0.001 & abs(log2FoldChangeShrink) >= log2(opt$fc) ) %>%
    pull(taxon)
  all_pcas_adj01 <- G4Micro::makeAllPCAs(phobj_filt, df2pca, taxa_padj01, vars2pca, opt, "DiffTaxaPadj01")
  #all_pcas_adj001 <- makeAllPCAs(phobj, df2pca, taxa_padj001, vars2pca, opt, "DiffTaxaPadj001")

  this_metadata <- sample_data(phobj_filt) %>% data.frame %>% dplyr::filter(sampleID %in% names(df2pca))
  all_model_results[[i]][["padj_taxa_taxa"]] <- taxa_padj
  all_model_results[[i]][["praw_taxa_taxa"]] <- taxa_praw
  all_model_results[[i]][["padj_taxa_pcas"]] <- all_pcas_adj
  all_model_results[[i]][["praw_taxa_pcas"]] <- all_pcas_praw
  all_model_results[[i]][["padj_taxa_res"]] <- G4Micro::callDoAllModelsFromALLPCAs(all_pcas_adj, name=paste0(i, "ConditionPadj"),
                                                                                   variable_plim=0.05,
                                                                                   metadata=this_metadata, vars2pca=var2predict,
                                                                                   levs2predict = levs2pred,
                                                                                   nfolds = 10, opt=opt)
  all_model_results[[i]][["padj_taxa_res_p1"]] <- G4Micro::callDoAllModelsFromALLPCAs(all_pcas_adj, name=paste0(i, "ConditionPadjP1"),
                                                                                   variable_plim=1,
                                                                                   metadata=this_metadata, vars2pca=var2predict,
                                                                                   levs2predict = levs2pred,
                                                                                   nfolds = 10, opt=opt)
  all_model_results[[i]][["padj_taxa_res_p02"]] <- G4Micro::callDoAllModelsFromALLPCAs(all_pcas_adj, name=paste0(i, "ConditionPadjP02"),
                                                                                      variable_plim=0.2,
                                                                                      metadata=this_metadata, vars2pca=var2predict,
                                                                                      levs2predict = levs2pred,
                                                                                      nfolds = 10, opt=opt)

  this_metadata_imputed <- this_metadata %>%
    dplyr::mutate(across(all_of(c(patnames, foodnames, foodnames_prop)),
                         ~ ifelse(is.na(.x), median(.x, na.rm=TRUE), .x) ))
  all_model_results[[i]][["padj_taxa_res_patts"]] <- G4Micro::callDoAllModelsFromALLPCAs(all_pcas_adj, name=paste0(i, "ConditionPadjPatterns"),
                                                                                   variable_plim=0.05,
                                                                                   metadata=this_metadata_imputed, vars2pca=var2predict,
                                                                                   meta_vars = patnames,
                                                                                   levs2predict = levs2pred,
                                                                                   nfolds = 10, opt=opt)
  all_model_results[[i]][["padj_taxa_res_patts_p1"]] <- G4Micro::callDoAllModelsFromALLPCAs(all_pcas_adj, name=paste0(i, "ConditionPadjPatternsP1"),
                                                                                         variable_plim=1,
                                                                                         metadata=this_metadata_imputed, vars2pca=var2predict,
                                                                                         meta_vars = patnames,
                                                                                         levs2predict = levs2pred,
                                                                                         nfolds = 10, opt=opt)
  all_model_results[[i]][["padj_taxa_res_patts_p02"]] <- G4Micro::callDoAllModelsFromALLPCAs(all_pcas_adj, name=paste0(i, "ConditionPadjPatternsP02"),
                                                                                            variable_plim=0.2,
                                                                                            metadata=this_metadata_imputed, vars2pca=var2predict,
                                                                                            meta_vars = patnames,
                                                                                            levs2predict = levs2pred,
                                                                                            nfolds = 10, opt=opt)
  all_model_results[[i]][["padj_taxa_res_foodclr_p1"]] <- G4Micro::callDoAllModelsFromALLPCAs(all_pcas_adj, name=paste0(i, "ConditionPadjFoodCLRP1"),
                                                                                         variable_plim=1,
                                                                                         metadata=this_metadata_imputed, vars2pca=var2predict,
                                                                                         meta_vars = foodnames,
                                                                                         levs2predict = levs2pred,
                                                                                         nfolds = 10, opt=opt)
  all_model_results[[i]][["padj_taxa_res_foodclr_p02"]] <- G4Micro::callDoAllModelsFromALLPCAs(all_pcas_adj, name=paste0(i, "ConditionPadjFoodCLRP02"),
                                                                                           variable_plim=0.2,
                                                                                           metadata=this_metadata_imputed, vars2pca=var2predict,
                                                                                           meta_vars = foodnames,
                                                                                           levs2predict = levs2pred,
                                                                                           nfolds = 10, opt=opt)
  all_model_results[[i]][["padj_taxa_res_foodprop_p1"]] <- G4Micro::callDoAllModelsFromALLPCAs(all_pcas_adj, name=paste0(i, "ConditionPadjFoodPROPP1"),
                                                                                           variable_plim=1,
                                                                                           metadata=this_metadata_imputed, vars2pca=var2predict,
                                                                                           meta_vars = foodnames_prop,
                                                                                           levs2predict = levs2pred,
                                                                                           nfolds = 10, opt=opt)
  all_model_results[[i]][["padj_taxa_res_foodprop_p02"]] <- G4Micro::callDoAllModelsFromALLPCAs(all_pcas_adj, name=paste0(i, "ConditionPadjFoodPROPP02"),
                                                                                            variable_plim=0.2,
                                                                                            metadata=this_metadata_imputed, vars2pca=var2predict,
                                                                                            meta_vars = foodnames_prop,
                                                                                            levs2predict = levs2pred,
                                                                                            nfolds = 10, opt=opt)
  #all_model_results[[i]][["praw_taxa_res"]] <- callDoAllModelsFromALLPCAs(all_pcas_praw, name=paste0(i, "ConditionPraw"), metadata=this_metadata, vars2pca=var2predict)
  all_model_results[[i]][["padj_taxa_res01"]] <- G4Micro::callDoAllModelsFromALLPCAs(all_pcas_adj01, name=paste0(i, "ConditionPadj01"),
                                                                                     variable_plim=0.05,
                                                                                     metadata=this_metadata,
                                                                                     levs2predict = levs2pred,
                                                                                     vars2pca=var2predict, nfolds = 10, opt=opt)
  all_model_results[[i]][["padj_taxa_res01_p1"]] <- G4Micro::callDoAllModelsFromALLPCAs(all_pcas_adj01, name=paste0(i, "ConditionPadj01P1"),
                                                                                     variable_plim=1,
                                                                                     metadata=this_metadata,
                                                                                     levs2predict = levs2pred,
                                                                                     vars2pca=var2predict, nfolds = 10, opt=opt)

  all_model_results[[i]][["padj_taxa_res01_patts"]] <- G4Micro::callDoAllModelsFromALLPCAs(all_pcas_adj01, name=paste0(i, "ConditionPadj01Patterns"),
                                                                                         variable_plim=0.05,
                                                                                         metadata=this_metadata_imputed, vars2pca=var2predict,
                                                                                         meta_vars = patnames,
                                                                                         levs2predict = levs2pred,
                                                                                         nfolds = 10, opt=opt)
  all_model_results[[i]][["padj_taxa_res01_patts_p1"]] <- G4Micro::callDoAllModelsFromALLPCAs(all_pcas_adj01, name=paste0(i, "ConditionPadj01PatternsP1"),
                                                                                            variable_plim=1,
                                                                                            metadata=this_metadata_imputed, vars2pca=var2predict,
                                                                                            meta_vars = patnames,
                                                                                            levs2predict = levs2pred,
                                                                                            nfolds = 10, opt=opt)
  all_model_results[[i]][["padj_taxa_res01_patts_p02"]] <- G4Micro::callDoAllModelsFromALLPCAs(all_pcas_adj01, name=paste0(i, "ConditionPadj01PatternsP02"),
                                                                                             variable_plim=0.2,
                                                                                             metadata=this_metadata_imputed, vars2pca=var2predict,
                                                                                             meta_vars = patnames,
                                                                                             levs2predict = levs2pred,
                                                                                             nfolds = 10, opt=opt)
  all_model_results[[i]][["padj_taxa_res01_foodclr_p1"]] <- G4Micro::callDoAllModelsFromALLPCAs(all_pcas_adj01, name=paste0(i, "ConditionPadj01FoodCLRP1"),
                                                                                              variable_plim=1,
                                                                                              metadata=this_metadata_imputed, vars2pca=var2predict,
                                                                                              meta_vars = foodnames,
                                                                                              levs2predict = levs2pred,
                                                                                              nfolds = 10, opt=opt)
  all_model_results[[i]][["padj_taxa_res01_foodclr_p02"]] <- G4Micro::callDoAllModelsFromALLPCAs(all_pcas_adj01, name=paste0(i, "ConditionPadj01FoodCLRP02"),
                                                                                               variable_plim=0.2,
                                                                                               metadata=this_metadata_imputed, vars2pca=var2predict,
                                                                                               meta_vars = foodnames,
                                                                                               levs2predict = levs2pred,
                                                                                               nfolds = 10, opt=opt)
  all_model_results[[i]][["padj_taxa_res01_foodprop_p1"]] <- G4Micro::callDoAllModelsFromALLPCAs(all_pcas_adj01, name=paste0(i, "ConditionPadj01FoodPROPP1"),
                                                                                               variable_plim=1,
                                                                                               metadata=this_metadata_imputed, vars2pca=var2predict,
                                                                                               meta_vars = foodnames_prop,
                                                                                               levs2predict = levs2pred,
                                                                                               nfolds = 10, opt=opt)
  all_model_results[[i]][["padj_taxa_res01_foodprop_p02"]] <- G4Micro::callDoAllModelsFromALLPCAs(all_pcas_adj01, name=paste0(i, "ConditionPadj01FoodPROPP02"),
                                                                                                variable_plim=0.2,
                                                                                                metadata=this_metadata_imputed, vars2pca=var2predict,
                                                                                                meta_vars = foodnames_prop,
                                                                                                levs2predict = levs2pred,
                                                                                                nfolds = 10, opt=opt)
  #all_model_results[[i]][["padj_taxa_res001"]] <- callDoAllModelsFromALLPCAs(all_pcas_adj001, name=paste0(i, "ConditionPadj001"), metadata=this_metadata, vars2pca=var2predict)

   PCs <- all_model_results[[i]][["padj_taxa_res"]]$varnames
   modelo_svm <- all_model_results[[i]][["padj_taxa_res"]]$models$`SVM-linear`$mod_noscale
   all_model_results[[i]][["padj_taxa_res_indiv"]] <- callDoAllModelsFromALLPCAsOriginalVars(all_pcas_adj, PCs,
                                                                                             modelo_svm = modelo_svm,
                                                                                             vstdf = df2pca,
                                                                                             name=paste0(i, "_ConditionPadjIndiv"),
                                                                                             vars2pca=vars2pca[1], metadata = this_metadata,
                                                                                             daares = daa_all[[i]]$resdf,
                                                                                             topns = c(5, 10, 20, 50, 100))

  all_model_results[[i]]$metadata <-this_metadata
  opt <- restaurar(opt)
## }
opt <- restaurar(opt)

save(all_model_results, file=paste0(opt$out, "PredictDAA_onlyGain/all_model_results_withIndTaxa.RData"))
#load(file=paste0(opt$out, "PredictDAA/all_model_results.RData"))

#Integrate
opt$out <- paste0(opt$out, "PredictDAA_onlyGain/")
makeLinePlotComparingPhobjs(all_model_results, opt, models_name1 = "padj_taxa_res", models_name2 = "padj_taxa_res01")
## Compare with Bacteria in componets

walk(names(all_model_results), makeLinePlotComparingSamePhobjModels,
     all_model_results, opt)

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

#### plot again
library(colorspace)
library(scico)
library(ggnewscale)

n2plot <- names(all_model_results$remove_tanda2)
n2plot <- n2plot[grep("res", n2plot)]
n2plot <- n2plot[!grepl("res_indiv", n2plot)] # "padj_taxa_res_indiv"

merged_df <- data.frame()
for(nn in n2plot){
  aux <- all_model_results$remove_tanda2[[nn]]$modummary %>%
    dplyr::mutate(sel_method = nn,
                  varsused = paste(all_model_results$remove_tanda2[[nn]]$varnames, sep="|", collapse="|"))
  merged_df <- rbind(merged_df, aux)
}
merged_df <- rbind(merged_df, all_model_results$remove_tanda2$padj_taxa_res_indiv$allmodsum)
write_tsv(merged_df, file = paste0(opt$out, "250926_all_models_summary.tsv"))

METRIC <- "BalancedAccuracy_l1out"
METRIC <- "AUC_l1out"
merged_df_mod <- merged_df %>%
  dplyr::mutate(sel_method = gsub("_", " ", sel_method)) %>%
  dplyr::mutate(model = gsub("_", " ", model)) %>%
  dplyr::mutate(model = gsub("regression", "regr.", model)) %>%
  filter(model != "Ensemble2") %>%
  filter(model != "KMeans") %>%
  dplyr::mutate(sel_method_mod = sel_method) %>%
  dplyr::mutate(sel_method_mod = gsub("foodprop", "+ Food Groups proportion", sel_method_mod))  %>%
  dplyr::mutate(sel_method_mod = gsub("foodclr", "+ Food Groups CLR", sel_method_mod))  %>%
  dplyr::mutate(sel_method_mod = gsub("patts", "+ Food Factors", sel_method_mod))  %>%
  dplyr::mutate(sel_method_mod = gsub("DESeq pval top", "taxa selected by DESeq p-val top", sel_method_mod))  %>%
  dplyr::mutate(sel_method_mod = gsub("^PC", "taxa selected by loadings in PC", sel_method_mod, perl=T))  %>%
  dplyr::mutate(sel_method_mod = gsub("padj taxa res01", "PCA with taxa p<0.01", sel_method_mod))  %>%
  dplyr::mutate(sel_method_mod = gsub("padj taxa res", "PCA with taxa p<0.05", sel_method_mod)) %>%
  dplyr::mutate(sel_method_mod = gsub("top ", "- top ", sel_method_mod)) %>%
  dplyr::mutate(sel_method_mod = gsub(" p1$", " - all vars.", sel_method_mod)) %>%
  dplyr::mutate(sel_method_mod = gsub(" p02$", " - vars with p<0.2", sel_method_mod)) %>%
  dplyr::mutate(sel_method_mod = gsub("PCA with taxa p<0.01$", "PCA with taxa p<0.01 - vars with p<0.05", sel_method_mod, perl=T)) %>%
  dplyr::mutate(sel_method_mod = gsub("PCA with taxa p<0.05$", "PCA with taxa p<0.05 - vars with p<0.05", sel_method_mod, perl=T)) %>%
  dplyr::mutate(sel_method_mod = gsub("Factors$", "Factors - vars with p<0.05", sel_method_mod, perl=T)) %>%
  dplyr::mutate(type = "PCA",
                type = ifelse(grepl("taxa selected", sel_method_mod), "Taxa", type),
                type = ifelse(grepl("Food Groups", sel_method_mod), "PCA + Food Groups", type),
                type = ifelse(grepl("Food Factors", sel_method_mod), "PCA + Food Factors", type)
                )

order_sel <- merged_df_mod %>%
  group_by(sel_method_mod) %>%
  dplyr::summarise(best = max(!!sym(METRIC))) %>%
  arrange(desc(best)) %>% pull(sel_method_mod)

order_mod <-  merged_df_mod %>%
  group_by(model) %>%
  dplyr::summarise(best = max(!!sym(METRIC))) %>%
  arrange(desc(best)) %>% pull(model)

merged_df_mod <- merged_df_mod %>%
  dplyr::mutate(model = factor(model, levels=order_mod),
                sel_method_mod = factor(sel_method_mod, levels=order_sel))


(g0 <- ggplot(merged_df_mod, aes(x=model, y=sel_method_mod,
                                 fill=!!sym(METRIC),
                          col=!!sym(METRIC),
                          size = !!sym(METRIC)))+
  #geom_tile() +
  geom_point(shape = 21) +  #col="darkgray"
  scale_size(range = c(0.01, 4)) +
  theme_classic() +
  scale_fill_scico(palette = "vikO", direction = 1) +
  scale_color_scico(palette = "vikO", direction = 1) +
  theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1))+
  ylab("Variable Selection Method") +
  xlab("Model")

)

ggsave(filename = paste0(opt$out, "DOTPLOT_", METRIC, ".pdf"), g0, width = 12, height = 9)
write_tsv(merged_df_mod, file = paste0(opt$out, "DOTPLOT_", METRIC, ".tsv"))

(g1 <- ggplot(merged_df_mod, aes(y=sel_method_mod, fill = type)) +
  geom_tile(
    aes(x = -0.5),   # <- your new variable
    width = 0.1, height = 0.9) +
  #ggsci::scale_fill_uchicago() +
    scale_fill_manual(
      values = gray.colors(n = length(unique(merged_df_mod$type)),
                           start = 0.9, end = 0.2)  # light gray → near black
    ) +
    theme(axis.text.x = element_blank()) +
    theme(axis.title.x = element_blank()) +
    #theme(axis.title.y = element_blank()) +
    theme(axis.ticks.x = element_blank()) +
    theme(axis.ticks.y = element_blank()) +
    theme(legend.position = "left" ) +
    ylab("Variable Selection Method")

)

g1 <- g1 + theme(plot.margin = margin(5, 0, 5, 5))

g0b <- g0 +
  theme(axis.title.y = element_blank(), axis.text.y=element_blank()) +
  theme(axis.ticks.y = element_blank()) +
  theme(plot.margin = margin(5, 5, 5, 0.1)) +
  geom_hline(yintercept = as.numeric(merged_df_mod$sel_method_mod)-0.5, linetype=1, alpha=0.5, col="gray")

library(patchwork)
(g2 <- (g1 | g0b) +
  plot_layout(widths = c(0.5, 15))
)
ggsave(filename = paste0(opt$out, "DOTPLOT_", METRIC, "_annot2.pdf"), g2, width = 12, height = 9)
