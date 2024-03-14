## DAA only with covariates, without depression condition
phseq_to_correct <- names(all_phyloseq)
interestvar <- "status_c2"
vars2test <- c("status_c2", "Category_T0", "age_months_t0", "Sex"
)
opt$reserva_0 <- opt$out
opt$out <- paste0(opt$out, "DESeq2_ControlVarsAlone/")
if(!dir.exists(opt$out)) dir.create(opt$out)
daa_all_corrected_only <- list()
for(phname in phseq_to_correct){
  cat("Doing DESeq2 Analysys with correction for: ", phname, "\n")
  phobj <- all_phyloseq[[phname]]
  phobj <- updatePsWithLogs(phobj, c("age_months_t0"))
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
#load(paste0(opt$out, "DESeq2_ControlVarsAlone/DESEQ2_controlVarsAlone_all.RData"))

#daa_all = map(daa_all_corrected_only, \(x)x[["status_c2"]])
#save(daa_all, file=paste0(opt$out, "DESeq2_ControlVarsAlone/DESEQ2_all_mainCondition.RData"))
#load(paste0(opt$out, "DESeq2_ControlVarsAlone/DESEQ2_all_mainCondition.RData"))
