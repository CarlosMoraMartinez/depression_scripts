# DAA correcting by covariates
phseq_to_correct <- names(all_phyloseq)[4]
interestvar <- "Condition"
vars2test <- c("Tanda", "Edad_log", "Sexo", "ob_o_sobrepeso","obesidad", "BMI_log",
               "ob_o_sobrepeso",
               "Estado.civil2", 
               "Educacion", "reads_log10", "Fumador", 
               "TG", "TG_mayor_200",  "Colesterol", "Colesterol_mayor_200", 
               "PSS_estres", "Diabetes", "bristol_scale", "bristol_scale_cualitativo",
               "defecaciones_semana", "Euroqol", "IPAQ_act_fisica", "Mediterranean_diet_adherence",
               "tratamiento_ansioliticos", "tratamiento_anticonvulsivos",
               "tratamiento_ISRNs", "tratamiento_antidiabeticos", "tratamiento_coagul_betabloq_etc", 
               "DII"
)
vars2test <- c("BMI_log", "ob_o_sobrepeso", "Edad_log",
               "Mediterranean_diet_adherence2", "Mediterranean_diet_adherence",
               "IPAQ_act_fisica")
opt$reserva_0 <- opt$out
opt$out <- paste0(opt$out, "DESeq2_ControlVars/")
if(!dir.exists(opt$out)) dir.create(opt$out)
daa_all_corrected <- list()
for(phname in phseq_to_correct){
  cat("Doing DESeq2 Analysys with correction for: ", phname, "\n")
  phobj <- all_phyloseq[[phname]]
  phobj <- updatePsWithLogs(phobj, c("Edad", "BMI"))
  daa_all_corrected[[phname]] <- list()
  for(var in vars2test){
    cat("Doing DESeq2 Analysys with correction for: ", phname, '-', var, "\n")
    samples <- sample_data(phobj)$sampleID[! is.na(sample_data(phobj)[, var])]
    phobj_filt <- phyloseq::prune_samples(samples, phobj)
    cases <- sample_data(phobj_filt)[, var] %>% unlist %>%  table
    if(length(which(cases > 0)) < 2 ){cat("Only one level, skipping this variable");next}
    
    vars2deseq <- c(interestvar, var)
    name <- paste0(phname, '_', var)
    res_tmp <-deseq_full_pipeline(phobj_filt, name, vars2deseq, opt)
    daa_all_corrected[[phname]][[var]] <- res_tmp$resdf
  }}
opt <- restaurar(opt)
save(daa_all_corrected, file=paste0(opt$out, "DESeq2_ControlVars/DESEQ2_controlVars_all.RData"))

## Plot heatmap with significant differences in all 
vargroups <- list(
  allvars=names(daa_all_corrected[[1]]),
  pobl1=c("Edad_log", "Sexo", "obesidad"),
  pobl2=c("Edad_log", "Sexo", "ob_o_sobrepeso"),
  pobl3=c("Edad_log", "Sexo", "IMC_log", "TG", "Colesterol"),
  pobl3=c("Edad_log", "Sexo", "ob_o_sobrepeso", "IMC_log", "TG", "Colesterol"),
  trats=vars2test[grep("trat", vars2test)],
  health=c("Diabetes", "TG_mayor_200", "Colesterol_mayor_200", "Mediterranean_diet_adherence", "IPAQ_act_fisica"),
  exer=c("Mediterranean_diet_adherence", "Euroqol","IPAQ_act_fisica", "bristol_scale_cualitativo", "PSS_estres")
)
vargroups <- list(pobl1=names(daa_all_corrected[[1]]),
                  pobl2 = c("BMI_log", "ob_o_sobrepeso"))
outdir <- paste0(opt$out, "DESeq2_ControlVars/SummaryHeatmaps/")
if(!dir.exists(outdir)) dir.create(outdir)
for(phname in phseq_to_correct){
  for(vgroup in names(vargroups)){
    oname <- paste0(phname, '_', vgroup, "_p05")
    makeHeatmapsFromMultipleDeseqResults(daa_all_corrected[[phname]], 
                                         daa_all[[phname]]$resdf, 
                                         main_name = interestvar,
                                         vars2plot = vargroups[[vgroup]],
                                         italics_rownames=T, 
                                         pfilt =0.01, pplot=0.05, name=oname, outdir=outdir)
    oname <- paste0(phname, '_', vgroup, "_p01")
    makeHeatmapsFromMultipleDeseqResults(daa_all_corrected[[phname]], 
                                         daa_all[[phname]]$resdf, 
                                         main_name = interestvar,
                                         vars2plot = vargroups[[vgroup]],
                                         italics_rownames=T, 
                                         pfilt =0.01, pplot=0.01,name=oname, outdir=outdir)
    
  }}
