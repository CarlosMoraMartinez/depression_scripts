daa_all <- list()
vars2deseq <- c("Condition")
opt$mincount <- 1
phseq_to_use <- c("remove_tanda2", "rmbatch_tanda", "filt")
for(phname in phseq_to_use){
  cat("Doing DESeq2 Analysys for: ", phname, "\n")
  daa_all[[phname]] <-deseq_full_pipeline(all_phyloseq[[phname]], phname, vars2deseq, opt)
}
save(daa_all, file = paste0(opt$out, "DeSEQ2/DESEQ2_all.RData"))