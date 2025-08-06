


library(tidyverse)

outdir <- "/home/carlos/Escritorio/202311_DEPRESION/ReviewJune2025/Results_rstudio/results2/use/ML_tables/"
dir.create(outdir)

load("/home/carlos/Escritorio/202311_DEPRESION/ReviewJune2025/Results_rstudio/results2/PredictDAA_multiclass/all_model_multi_results.RData")

pp <- all_model_multi_results$remove_tanda2_3$DiffTaxaPadjCovOrSingle05$component_pvals %>% arrange(pval)
write_tsv(pp, paste0(outdir, "PCs_used_4groups.tsv"))

pp <- all_model_multi_results$remove_tanda2_3$DiffTaxaPadjCovOrSingle05$modummary
write_tsv(pp, paste0(outdir, "Summary_4groups.tsv"))

smotsum <- all_model_multi_results$remove_tanda2_3smote$DiffTaxaPadjCovOrSingle05$modummary
write_tsv(smotsum, paste0(outdir, "Summary_4groups_SMOTE.tsv"))

LINDAUSED <- "LinDADiffTaxaPadjCov"

pplin <- all_model_multi_results$remove_tanda2_3[[LINDAUSED]]$component_pvals %>% arrange(pval)
write_tsv(pplin, paste0(outdir, "PCs_used_4groups_LinDA.tsv"))

mm <- all_model_multi_results$remove_tanda2_3[[LINDAUSED]]$modummary
write_tsv(mm, paste0(outdir, "Summary_4groups_LinDA.tsv"))

mm <- all_model_multi_results$remove_tanda2_3smote[[LINDAUSED]]$modummary
write_tsv(mm, paste0(outdir, "Summary_4groups_LinDA_SMOTE.tsv"))


## 2 Classes
load("/home/carlos/Escritorio/202311_DEPRESION/ReviewJune2025/Results_rstudio/results2/PredictDAA_L1O_2/all_model_results_L1O.RData")

pp2 <- all_model_results$remove_tanda2$padj_taxa_res$modummary
write_tsv(pp2, paste0(outdir, "Summary_2groups.tsv"))

pp2 <- all_model_results$remove_tanda2$padj_taxa_res_SMOTE$modummary
write_tsv(pp2, paste0(outdir, "Summary_2groups_SMOTE.tsv"))

pp2 <- all_model_results$remove_tanda2$padj_taxa_res$component_pvals
write_tsv(pp2, paste0(outdir, "PCs_used_2groups.tsv"))

pp2 <- all_model_results$remove_tanda2$padj_taxa_res01linda$modummary
write_tsv(pp2, paste0(outdir, "Summary_2groups_LinDA.tsv"))

pp2 <- all_model_results$remove_tanda2$padj_taxa_res01linda_SMOTE$modummary
write_tsv(pp2, paste0(outdir, "Summary_2groups_LinDA_SMOTE.tsv"))

pp2 <- all_model_results$remove_tanda2$padj_taxa_res01linda$component_pvals %>% arrange(pval)
write_tsv(pp2, paste0(outdir, "PCs_used_2groups_LinDA.tsv"))

## All 
allsums <- read_tsv("/home/carlos/Escritorio/202311_DEPRESION/ReviewJune2025/Results_rstudio/results2/PredictDAA_L1O_2/remove_tanda2_orderBy_BalancedAccuracy_l1out_normal/remove_tanda2_modelSummariesWithIndividualSpecies.tsv")

allsums <- allsums %>% filter(model != "Ensemble2") %>% 
  filter(!grepl("combined 2 ", sel_method)) %>% 
filter(!grepl("PC2\\+PC11", sel_method))
write_tsv(allsums, paste0(outdir, "2groups_withIndividualSpeciesAlso.tsv"))

## Write to excel
library(openxlsx)

names2remove <- c("Accuracy", "Kappa", "Sensitivity", "Specificity", "PPV", "NPV", "Precision", "Recall", "BalancedAccuracy")
files <- list.files(outdir, pattern = "\\.tsv$", full.names = TRUE)
filedf <- data.frame(
  files = files, 
  basenames =  tools::file_path_sans_ext(basename(files))
) %>% 
  dplyr::mutate(
    program = ifelse(grepl("LinDA", basenames), 2, 1) %>% factor(),
    smote = ifelse(grepl("SMOTE", basenames), 2, 1) %>% factor(),
    type = ifelse(grepl("Summary", basenames), 2, 1),
    type = ifelse(grepl("IndividualSpeciesAlso", basenames), 3, type) %>% factor(),
    problem = ifelse(grepl("2groups", basenames), 1, 2) %>% factor(),
  ) %>% 
  arrange(program, problem, type, smote)

wb <- createWorkbook()

for (file in filedf$files) {
  sheet_name <- tools::file_path_sans_ext(basename(file))
  if(!grepl("LinDA", sheet_name)) sheet_name <- paste0(sheet_name, "_DESeq2")

  if(nchar(sheet_name) > 31) sheet_name <- paste(strsplit(sheet_name, "")[[1]][1:31], sep="", collapse="")
  
  data <- read_tsv(file) 
  if(grepl("Summary", sheet_name) | grepl("IndividualSpecies", sheet_name)) data <- data %>% select(-all_of(names2remove))
  
  addWorksheet(wb, sheet_name)
  writeData(wb, sheet_name, data)
}

saveWorkbook(wb, file = file.path(outdir, "TableS7_PredictionResults.xlsx"), overwrite = TRUE)

