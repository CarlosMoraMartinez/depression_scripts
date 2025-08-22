opt <- restaurar(opt)
daa_all <- list()
vars2deseq <- c("edad_00")
opt$mincount <- 1
phseq_to_use <- names(all_phyloseq)#c("remove_tanda2", "rmbatch_tanda", "filt")
for(phname in phseq_to_use){
  cat("Doing DESeq2 Analysys for: ", phname, "\n")
  phobj <- all_phyloseq[[phname]]
  samples <- sample_data(phobj)$sampleID[! is.na(sample_data(phobj)[, vars2deseq[1] ])]
  phobj_filt <- phyloseq::prune_samples(samples, phobj)
  
  daa_all[[phname]] <-deseq_full_pipeline(phobj_filt, phname, vars2deseq, opt)
}
save(daa_all, file = paste0(opt$out, "DeSEQ2/DESEQ2_all_edad00.RData"))


### Age as class and covariates

opt <- restaurar(opt)
outdir <- paste0(opt$out, "DESeq2_ageclass_Covars1/")
if(!dir.exists(outdir)) dir.create(outdir)

daa_all <- list()

vars2deseq <- c("age_class1", "Sex", "edu_m_00", "z_imc_00")

phseq_to_use <- "remove_tanda2" #names(all_phyloseq)[-1]#c("remove_tanda2", "rmbatch_tanda", "filt")
for(phname in phseq_to_use){
  cat("Doing DESeq2 Analysys for: ", phname, "\n")
  phobj <- all_phyloseq[[phname]]
  
  samples_with_nas <- sample_data(phobj) %>% data.frame %>% dplyr::select(all_of(vars2deseq)) %>% mutate_all(is.na) %>% apply(MAR=1, any)
  samples <- sample_data(phobj)$sampleID[!samples_with_nas]
  phobj_filt <- phyloseq::prune_samples(samples, phobj)
  
  daa_all[[phname]] <-deseq_full_pipeline(phobj_filt, phname, vars2deseq, opt, deseqname = "DESeq2_ageclass_Covars1/")
}
save(daa_all, file = paste0(opt$out, "DESeq2_ageclass_Covars1/DESEQ2_all_edadCovs.RData"))

