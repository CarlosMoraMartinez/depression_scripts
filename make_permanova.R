
load(paste0(opt$out, "/phyloseq_original/phyloseq_all_list.RData"))

exclude_vars <- c("sampleID", "CODIGO", "CP")
dtypes <- c("bray", "jaccard")

outdir <- paste0(opt$out, "PERMANOVA/")
if(!dir.exists(outdir)) dir.create(outdir)
permaresults <- list()

phseq2use <- "remove_tanda2_rarefied_min"
for(i in phseq2use){
  cat("Doing PERMANOVA of: ", i)
  phobj <- all_phyloseq[[i]]
  permaresults[[i]] <- lapply(dtypes, FUN=function(dd, phobj, exclude_vars, SEED){
    oname <- paste0(outdir, "permanova_results_", i, "_", dd, ".tsv")
    makePermanova(phobj,dist_method = dd, 
                  seed = SEED, 
                  exclude_vars = exclude_vars, 
                  outname = oname) 
  }, phobj, exclude_vars, SEED)
  names(permaresults[[i]]) <- dtypes
}

permaresults$remove_tanda2_rarefied_min$bray %>% filter(variable %in% c("Condition", "Edad", "Sexo", "IMC"))
### ADONIS with multiple variables
permaformulas <- c(
  "braydist ~ Condition + Sexo",
  "braydist ~ Condition + BMI_log",
  "braydist ~ Condition + Edad_log",
  "braydist ~ Condition + IPAQ_act_fisica",
  "braydist ~ Condition + Mediterranean_diet_adherence2",
  "braydist ~ Condition + BMI_log + Sexo + Edad_log",
  "braydist ~ Condition + BMI_log + Sexo + Edad_log + IPAQ_act_fisica",
  "braydist ~ Condition + BMI_log + Sexo + Edad_log + Mediterranean_diet_adherence2",
  "braydist ~ Condition + BMI_log + Sexo + Edad_log + IPAQ_act_fisica + Mediterranean_diet_adherence2"
)

permaresults_mult <- list()
for(i in phseq2use){
  cat("Doing PERMANOVA of: ", i)
  phobj <- all_phyloseq[[i]]
  phobj <- updatePsWithLogs(phobj, c("Edad", "BMI"))
  permaresults_mult[[i]] <- lapply(dtypes, FUN=function(dd, phobj, exclude_vars, SEED){
    oname <- paste0(outdir, "permanova_resultsMult_", i, "_", dd, ".tsv")
    makePermanovaFormulas(phobj,
                  permaformulas,
                  dist_method = dd, 
                  seed = SEED, 
                  outname = oname) 
  }, phobj, exclude_vars, SEED)
  names(permaresults_mult[[i]]) <- dtypes
}
save(permaresults_mult, file = paste0(outdir, "PERMANOVA_MULT.RData"))

mm <- permaresults_mult$remove_tanda2_rarefied_min$bray$modelos

### get variables with significant dispersion tests

xx<- read_tsv(paste0(outdir, "/permanova_results_remove_tanda2_rarefied_min_bray.tsv"))

permanova_useful <- xx %>% arrange(perm_disp_P) %>% 
  filter(variable %in% c("Condition", "Sexo", "BMI", 
                         "Edad", "Edad_log", 
                         "IPAQ", "Mediterranean_diet_adherence2", 
                         "IPAQ_act_fisica", "Mediterranean_diet_adherence", 
                         "ob_o_sobrepeso",
                         "Smoking_status")) %>% 
  select(variable, DF_var, DF_Residual, DF_Total, R2_var, R2_Residual, F_statistic, P, perm_disp_P) %>% 
  dplyr::mutate(perm_disp_Padj = p.adjust(perm_disp_P, method = "BH"))
write_tsv(permanova_useful, file = paste0(outdir, "permanova_results_filtered_useful.tsv"))
