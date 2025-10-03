opt <- restaurar(opt)
daa_all <- list()
vars2deseq <- c("status_c2")
opt$mincount <- 1
phseq_to_use <- names(all_phyloseq)#c("remove_tanda2", "rmbatch_tanda", "filt")
for(phname in phseq_to_use){
  cat("Doing DESeq2 Analysys for: ", phname, "\n")
  phobj <- all_phyloseq[[phname]]
  aux <- sample_data(phobj) %>% data.frame()
  samples <- aux %>%
    filter(! is.na(status_c2)) %>%
    filter(status_c1=="normal") %>%
    filter(status_c2 != "NaN") %>%
    filter(status_c2 != "Insufficient gain") %>% pull(sampleID)
  phobj_filt <- phyloseq::prune_samples(samples, phobj)

  daa_all[[phname]] <-deseq_full_pipeline(phobj_filt, phname, vars2deseq, opt)
}
save(daa_all, file = paste0(opt$out, "DeSEQ2/DESEQ2_all.RData"))

#### To predict food
opt <- restaurar(opt)
opt$out <- paste0(opt$out, "/DAA_Diet_MedOrNot/")
if(!dir.exists(opt$out)) dir.create(opt$out)

daa_all <- list()
phseq_to_use <- "remove_tanda2"
vars2deseq <- c("Pattern_MedOrNot")
opt$mincount <- 1
#phseq_to_use <- names(all_phyloseq)#c("remove_tanda2", "rmbatch_tanda", "filt")
for(phname in phseq_to_use){
  cat("Doing DESeq2 Analysys for: ", phname, "\n")
  phobj <- all_phyloseq[[phname]]

  aux <- sample_data(phobj) %>% data.frame()
  sample_data(phobj)$Pattern_MedOrNot <- ifelse(aux$pattern == "Mediterranean", "Mediterranean", "Proc. or Western")
  samples <- aux %>%
    filter(! is.na(pattern)) %>% pull(sampleID)

  phobj_filt <- phyloseq::prune_samples(samples, phobj)

  daa_all[[phname]] <-deseq_full_pipeline(phobj_filt, phname, vars2deseq, opt)
}
all_phyloseq[[phname]] <- phobj
save(daa_all, file = paste0(opt$out, "DeSEQ2/DESEQ2_all_MedVsProcOrWes.RData"))
save(all_phyloseq, file = paste0(opt$out, "DeSEQ2/DESEQ2_all_MedVsProcOrWes_phyloseqList.RData"))
save(phobj_filt, file = paste0(opt$out, "DeSEQ2/DESEQ2_all_MedVsProcOrWes_phobj_filt.RData"))
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

## other covariates one by one

opt <- restaurar(opt)

all_vars2deseq <- c("z_t0", "z_t1", "status_c2", "edad_00", "bmi_t0", "bmi_t1")
opt$mincount <- 1
phseq_to_use <- names(all_phyloseq)#c("remove_tanda2", "rmbatch_tanda", "filt")
for(vars2deseq in all_vars2deseq){
  daa_all <- list()
  for(phname in phseq_to_use){

    cat("Doing DESeq2 Analysys for: ", phname, ", variable: ", var2use, "\n")
    phobj <- all_phyloseq[[phname]]
    samples <- sample_data(phobj)$sampleID[! is.na(sample_data(phobj)[, var2use])]
    phobj_filt <- phyloseq::prune_samples(samples, phobj)

    daa_all[[phname]] <-deseq_full_pipeline(phobj_filt, phname, vars2deseq, opt)
  }
  save(daa_all, file = paste0(opt$out, "DeSEQ2/", var2use, "DESEQ2_all.RData"))
}

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
