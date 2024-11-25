opt <- restaurar(opt)
daa_all <- list()
vars2deseq <- c("Category_BinT1", "Category_BinT0")
opt$mincount <- 1
phseq_to_use <- c("remove_tanda2", "rmbatch_tanda", "filt", "phseq_rerefthenbatch_tanda")
for(phname in phseq_to_use){
  cat("Doing DESeq2 Analysys for: ", phname, "\n")
  phobj <- all_phyloseq[[phname]]
  samples <- sample_data(phobj)$sampleID[! is.na(sample_data(phobj)[, "Category_T1" ])]
  phobj_filt <- phyloseq::prune_samples(samples, phobj)
  
  meta_mod <- sample_data(phobj_filt) %>% data.frame %>% 
    dplyr::mutate(Category_BinT1 = ifelse(Category_T1 %in% c("Thinness", "Normal"), "Normal", "Overweight"),
                  Category_BinT0 = ifelse(Category_T0 %in% c("Thinness", "Normal"), "Normal", "Overweight"))
  sample_data(phobj_filt)$Category_BinT1 <- meta_mod$Category_BinT1
  sample_data(phobj_filt)$Category_BinT0 <- meta_mod$Category_BinT0
  
  daa_all[[phname]] <-deseq_full_pipeline(phobj_filt, phname, vars2deseq, opt)
}
save(daa_all, file = paste0(opt$out, "DeSEQ2/DESEQ2_all.RData"))

#load( paste0(opt$out, "DeSEQ2/DESEQ2_all.RData"))

############################# look at T0

daa_all <- list()
vars2deseq <- c("Category_T0")
opt$mincount <- 1

opt <- restaurar(opt)
opt$out <- paste0(opt$out, "DESeq2_statusAtT0_2/")
if(!dir.exists(opt$out)) dir.create(opt$out)
#phseq_to_use <- c("remove_tanda2", "rmbatch_tanda", "filt")
for(phname in phseq_to_use){
  cat("Doing DESeq2 Analysys for: ", phname, "\n")
  phobj <- all_phyloseq[[phname]]
  samples <- sample_data(phobj) %>% data.frame %>% 
    dplyr::filter(! is.na( !!sym(vars2deseq[1]) )) %>% 
  pull(sampleID)
  phobj_filt <- phyloseq::prune_samples(samples, phobj)
  daa_all[[phname]] <-deseq_full_pipeline(phobj_filt, phname, vars2deseq, opt)
}
opt <- restaurar(opt)
save(daa_all, file = paste0(opt$out, "DESeq2_statusAtT0_2/DESEQ2_all_statusAtT0.RData"))
load(paste0(opt$out, "DESeq2_statusAtT0_2/DESEQ2_all_statusAtT0.RData"))

### Several covariates

daa_all_covs <- list()
vars2deseq <- c("Category_T0", "edad_vs", "hospital")
opt$mincount <- 1

opt <- restaurar(opt)
opt$out <- paste0(opt$out, "DESeq2_statusAtT0_severalCovars/")
if(!dir.exists(opt$out)) dir.create(opt$out)
#phseq_to_use <- c("remove_tanda2", "rmbatch_tanda", "filt")
for(phname in phseq_to_use){
  cat("Doing DESeq2 Analysys for: ", phname, "\n")
  phobj <- all_phyloseq[[phname]]
  samples <- sample_data(phobj) %>% data.frame %>% 
    dplyr::filter(! is.na( !!sym(vars2deseq[1]) )) %>% 
    dplyr::filter(! is.na( !!sym(vars2deseq[2]) )) %>%
    dplyr::filter(! is.na( !!sym(vars2deseq[3]) )) %>% 
    pull(sampleID)
  phobj_filt <- phyloseq::prune_samples(samples, phobj)
  
  daa_all_covs[[phname]] <-deseq_full_pipeline(phobj_filt, phname, vars2deseq, opt)
}
opt <- restaurar(opt)
save(daa_all_covs, file = paste0(opt$out, "DESeq2_statusAtT0_severalCovars/DESEQ2_all_severalCovars.RData"))
load(paste0(opt$out, "DESeq2_statusAtT0_severalCovars/DESEQ2_all_severalCovars.RData"))

### Interaction with age 

phseq_to_use <- c("remove_tanda2", "rmbatch_tanda", "filt")
daa_all_covs_int <- list()
vars2deseq <- c("Category_T0", "edad_vs")
opt$mincount <- 1

opt <- restaurar(opt)
opt$out <- paste0(opt$out, "DESeq2_statusAtT0_severalCovarsInteraction/")
if(!dir.exists(opt$out)) dir.create(opt$out)
#phseq_to_use <- c("remove_tanda2", "rmbatch_tanda", "filt")
for(phname in phseq_to_use){
  cat("Doing DESeq2 Analysys for: ", phname, "\n")
  phobj <- all_phyloseq[[phname]]
  samples <- sample_data(phobj) %>% data.frame %>% 
    dplyr::filter(! is.na( !!sym(vars2deseq[1]) )) %>% 
    dplyr::filter(! is.na( !!sym(vars2deseq[2]) )) %>%
    pull(sampleID)
  phobj_filt <- phyloseq::prune_samples(samples, phobj)
  
  #daa_all_covs_int[[phname]] <- deseq_full_pipeline(phobj_filt, phname, vars2deseq, opt, interact = TRUE)
  daa_all_covs_int[[phname]] <- getDeseqResults(phobj_filt, opt, name=phname, variables = vars2deseq, interact=TRUE)
}
opt <- restaurar(opt)
save(daa_all_covs_int, file = paste0(opt$out, "DESeq2_statusAtT0_severalCovarsInteraction/DESEQ2_all_severalCovarsInteraction.RData"))
load(paste0(opt$out, "DESeq2_statusAtT0_severalCovarsInteraction/DESEQ2_all_severalCovarsInteraction.RData"))

#########################33
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
load(paste0(opt$out, "DESeq2_MainConditionOnlyGain/DESEQ2_all_removeWeightLoss.RData"))
