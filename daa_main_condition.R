opt <- restaurar(opt)
daa_all <- list()
vars2deseq <- c("status_c2")
opt$mincount <- 1
phseq_to_use <- names(all_phyloseq)#c("remove_tanda2", "rmbatch_tanda", "filt")
for(phname in phseq_to_use){
  cat("Doing DESeq2 Analysys for: ", phname, "\n")
  phobj <- all_phyloseq[[phname]]
  samples <- sample_data(phobj)$sampleID[! is.na(sample_data(phobj)[, vars2deseq[1] ])]
  phobj_filt <- phyloseq::prune_samples(samples, phobj)
  
  daa_all[[phname]] <-deseq_full_pipeline(phobj_filt, phname, vars2deseq, opt)
}
save(daa_all, file = paste0(opt$out, "DeSEQ2/DESEQ2_all.RData"))

#load( paste0(opt$out, "DeSEQ2/DESEQ2_all.RData"))

daa_all <- list()
vars2deseq <- c("status_c2")
opt$mincount <- 1
levels2exclude <- c("Insufficient gain", "initially_overweight")

opt <- restaurar(opt)
opt$out <- paste0(opt$out, "DESeq2_MainConditionOnlyGain/")
if(!dir.exists(opt$out)) dir.create(opt$out)
#phseq_to_use <- c("remove_tanda2", "rmbatch_tanda", "filt")
for(phname in phseq_to_use){
  cat("Doing DESeq2 Analysys for: ", phname, "\n")
  phobj <- all_phyloseq[[phname]]
  samples <- sample_data(phobj) %>% data.frame %>% 
    dplyr::filter(! is.na( !!sym(vars2deseq[1]) )) %>% 
    dplyr::filter( !!sym(vars2deseq[1]) != levels2exclude[1] ) %>% 
    dplyr::filter( !!sym(vars2deseq[1]) != levels2exclude[2] ) %>% 
    pull(sampleID)
  phobj_filt <- phyloseq::prune_samples(samples, phobj)
  
  daa_all[[phname]] <-deseq_full_pipeline(phobj_filt, phname, vars2deseq, opt)
}
opt <- restaurar(opt)
save(daa_all, file = paste0(opt$out, "DESeq2_MainConditionOnlyGain/DESEQ2_all_removeWeightLoss.RData"))
