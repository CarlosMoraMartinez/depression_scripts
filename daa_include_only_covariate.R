## DAA only with covariates, without depression condition
phseq_to_correct <- names(all_phyloseq)
interestvar <- "status_c2"
vars2test <- c("status_c2", "age_months_t0", "Sex") #"Category_T0",
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

#Plot with age
for(phname in phseq_to_correct){
  phobj <- all_phyloseq[[phname]]
  age_orgs <- daa_all_corrected_only[[phname]]$age_months_t0 %>% filter(padj < 0.01) %>% arrange(log2FoldChangeShrink) %>% pull(taxon) %>% 
    gsub("-", ".", .) %>% 
    gsub("[\\[\\]]", "", ., perl=TRUE)
  counts_dir <- paste0(opt$out, "DESeq2_ControlVarsAlone/DeSEQ2/", phname, "_age_months_t0")
  counts_f <- list.files(counts_dir, pattern="norm_counts.tsv", full.names = T)
  counts <- read_tsv(counts_f)
  meta_otus <- sample_data(phobj) %>% data.frame
  otu_this <- counts %>% 
    dplyr::mutate(gene = gsub("-", ".", gene)) %>% 
    dplyr::mutate(gene = gsub("[\\[\\]]", "", gene, perl=TRUE)) %>% 
    column_to_rownames("gene") %>% as.matrix %>% t %>% data.frame %>% 
    dplyr::select(all_of(age_orgs)) %>% 
    rownames_to_column("sampleID")
  mergedtab <- merge(meta_otus, otu_this, by="sampleID", all=TRUE)
  merged_long <- mergedtab %>% gather("organism", "abundance", age_orgs)
  g1 <- ggplot(merged_long, aes(x=age_months_t0, y = abundance, group=organism, col=organism)) + 
    facet_wrap(organism ~ ., scales = "free", ncol=4)+
    #geom_point(show.legend = FALSE, size=0.2, alpha=0.5) + 
    geom_smooth(alpha=0.5, )+
    theme_bw()+
    theme(legend.position="none")

  ggsave(filename = paste0(counts_dir, "/age_vs_abundance_smooth.pdf"), g1, width = 12, height = 24, limitsize = F)
  
}

