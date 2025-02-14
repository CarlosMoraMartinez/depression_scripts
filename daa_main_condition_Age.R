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

