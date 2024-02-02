## DAA only with covariates, without depression condition
phseq_to_correct <- names(all_phyloseq)[4]
interestvar <- "Condition"
vars2test <- c("Tanda", "Edad_log", "Sexo", "ob_o_sobrepeso","obesidad", "BMI_log", "Estado.civil2", 
               "Educacion", "reads_log10", "Fumador", "DII",
               "TG", "TG_mayor_200",  "Colesterol", "Colesterol_mayor_200", 
               "PSS_estres", "Diabetes", "bristol_scale", "bristol_scale_cualitativo",
               "defecaciones_semana", "Euroqol", "IPAQ_act_fisica", "Mediterranean_diet_adherence",
               "tratamiento_ansioliticos", "tratamiento_anticonvulsivos",
               "tratamiento_ISRNs", "tratamiento_antidiabeticos", "tratamiento_coagul_betabloq_etc" 
)
vars2test <- c("BMI_log", "ob_o_sobrepeso", "Edad_log", "IPAQ_act_fisica", "Mediterranean_diet_adherence")
opt$reserva_0 <- opt$out
opt$out <- paste0(opt$out, "DESeq2_ControlVarsAlone/")
if(!dir.exists(opt$out)) dir.create(opt$out)
daa_all_corrected_only <- list()
for(phname in phseq_to_correct){
  cat("Doing DESeq2 Analysys with correction for: ", phname, "\n")
  phobj <- all_phyloseq[[phname]]
  phobj <- updatePsWithLogs(phobj, c("Edad", "BMI"))
  daa_all_corrected_only[[phname]] <- list()
  for(var in vars2test){
    cat("Doing DESeq2 Analysys with correction for: ", phname, '-', var, "\n")
    samples <- sample_data(phobj)$sampleID[! is.na(sample_data(phobj)[, var])]
    phobj_filt <- phyloseq::prune_samples(samples, phobj)
    cases <- sample_data(phobj_filt)[, var] %>% unlist %>%  table
    if(length(which(cases > 0)) < 2 ){cat("Only one level, skipping this variable");next}
    name <- paste0(phname, '_', var)
    res_tmp <-deseq_full_pipeline(phobj_filt, name, var, opt)
    daa_all_corrected_only[[phname]][[var]] <- res_tmp$resdf
  }}
opt$out <- opt$reserva_0
save(daa_all_corrected_only, file=paste0(opt$out, "DESeq2_ControlVarsAlone/DESEQ2_controlVarsAlone_all.RData"))
