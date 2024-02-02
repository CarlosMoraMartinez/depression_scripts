
opt$reserva_0 <- opt$out
opt$out <- paste0(opt$out, "DESeq2_MainVariableScales/")
if(!dir.exists(opt$out)) dir.create(opt$out)

daa_all_scales <- list()
vars2deseq <- c(escalas_qual, escalas_quant)
opt$mincount <- 1
for(phname in names(all_phyloseq)){
  phobj <- all_phyloseq[[phname]]
  daa_all_scales[[phname]] <- list()
  for(var in vars2deseq){
    cat("Doing DESeq2 Analysys for: ", phname, ", with main variable ", var, "\n")
    samples <- sample_data(phobj)$sampleID[!is.na(sample_data(phobj)[, var]) ]
    phobj_filt <- phyloseq::prune_samples(samples, phobj)
    daa_all_scales[[phname]][[var]] <-deseq_full_pipeline(phobj_filt, paste0(phname, "_",var), var, opt)
  }}

save(daa_all_scales, file = paste0(opt$out, "DeSEQ2/DESEQ2_MainVariableScales.RData"))
opt$out <- opt$reserva_0