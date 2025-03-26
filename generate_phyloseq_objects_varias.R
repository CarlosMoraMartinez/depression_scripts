########################################
# Generate all Phyloseqs 
########################################

########################################
# Generate Phyloseq basic
########################################
filterPhyla <- c("Chloroplast", "Mitochondria", "Eukaryota", "Metazoa", "Viruses")

get_filtered_phyloseq <- function(pre_phyloseq1, phseqname="", filterPhyla, opt, outdir, path_phyloseq){

  #pre_phyloseq1 = subset_taxa(pre_phyloseq1, !Phylum %in% c(NA))
  
  pre_phyloseq1 <- subset_taxa(pre_phyloseq1, !Kingdom %in% filterPhyla)
  pre_phyloseq1 <- subset_taxa(pre_phyloseq1, !Phylum %in% filterPhyla)
  pre_phyloseq1 <- subset_taxa(pre_phyloseq1, !Class %in% filterPhyla)
  pre_phyloseq1 <- subset_taxa(pre_phyloseq1, !Order %in% filterPhyla) # 12 a nivel Order
  pre_phyloseq1 <- subset_taxa(pre_phyloseq1, !Family %in% filterPhyla) # 7 a nivel Family
  pre_phyloseq1 <- subset_taxa(pre_phyloseq1, !Genus %in% filterPhyla)
  
  save(file=paste0(path_phyloseq, "/phyloseq_object_raw_filt_by_Phylum_", 
                   janitor::make_clean_names(phseqname), ".RData"), pre_phyloseq1)

  ## Calculate prevalence
  ottmp <- phyloseq::otu_table(pre_phyloseq1)
  s_meta <- sample_data(pre_phyloseq1) %>% data.frame()
  
  pre_prevalence <- apply(X = ottmp,
                          MARGIN = ifelse(taxa_are_rows(pre_phyloseq1), yes = 1, no = 2),
                          FUN = function(x){sum(x > opt$mincountspersample)})
  pre_prevalence = data.frame(Prevalence = pre_prevalence,
                              TotalAbundance = phyloseq::taxa_sums(pre_phyloseq1),
                              tax_table(pre_phyloseq1), 
                              relative_prevalence = pre_prevalence/ nsamples(pre_phyloseq1)
  )
  write_tsv(pre_prevalence, paste0(path_phyloseq, "/raw_prevalence_", 
                                   janitor::make_clean_names(phseqname), ".tsv"))
  
  ## Filtered to frequency
  prevalenceThreshold = opt$minfreq * nsamples(pre_phyloseq1)
  keepTaxa = rownames(pre_prevalence)[(pre_prevalence$Prevalence >= prevalenceThreshold)]
  (pre_phyloseq_filt = prune_taxa(keepTaxa, pre_phyloseq1))
  filtered_phyloseq_filename <- paste0(path_phyloseq, 
                                       '/pre_phyloseq_filt_by_prevalence', 
                                       as.character(100*opt$minfreq), 
                                       janitor::make_clean_names(phseqname),'.RData')
  save(pre_phyloseq_filt, file = filtered_phyloseq_filename)
  
  #Reads before rarefeact
  nreads <- otu_table(pre_phyloseq_filt) %>% colSums()
  s_meta$nreads_filt <- nreads[s_meta$sampleID]
  write_tsv(s_meta, paste0(outdir, "/full_metadata2_",  janitor::make_clean_names(phseqname), ".tsv"))
  
  sample_data(pre_phyloseq_filt)$nreads_filt <- nreads[sample_data(pre_phyloseq_filt)$sampleID]
  
  return(pre_phyloseq_filt)
}

getRarefied <- function(phobj){
  pre_phyloseq_rarefied <-rarefy_even_depth(phobj, rngseed = SEED)
  return(pre_phyloseq_rarefied )
}

path_phyloseq <- paste0(opt$out, "/phyloseq")
if(! dir.exists(path_phyloseq)){dir.create(path_phyloseq)}


all_mpas <- all_mpas %>% mutate(
  Classification = map(Classification, \(x){rownames(x)<- x$Species; return(x)}),
  phyloseq_species = pmap(all_mpas, ~ phyloseq(sample_data(..9),
                                    otu_table(..7, taxa_are_rows = TRUE),
                                    tax_table(as.matrix(..6)))
                        ),
  phyloseq_filt = map2(phyloseq_species, Condition, .f=get_filtered_phyloseq, 
                       filterPhyla=filterPhyla, 
                       opt=opt, 
                       outdir=input_tabs_dir, 
                       path_phyloseq=path_phyloseq),
  phyloseq_rarefMin = map(phyloseq_filt, .f=getRarefied)
  )



walk2(all_mpas$Condition, all_mpas$phyloseq_species, \(condname, phobj) save(phobj, 
                                                                        file = paste0(path_phyloseq, 
                                                                               "/phyloseq_", 
                                                                               janitor::make_clean_names(condname), 
                                                                               ".RData")))
walk2(all_mpas$Condition, all_mpas$phyloseq_filt, \(condname, phobj) save(phobj, 
                                                                             file = paste0(path_phyloseq, 
                                                                                           "/phyloseq_filt_", 
                                                                                           janitor::make_clean_names(condname), 
                                                                                           ".RData")))


walk2(all_mpas$Condition, all_mpas$phyloseq_rarefMin, \(condname, phobj) save(phobj, 
                                                                          file = paste0(path_phyloseq, 
                                                                                        "/phyloseq_rarefMin_", 
                                                                                        janitor::make_clean_names(condname), 
                                                                                        ".RData")))



save(all_mpas, file=paste0(input_tabs_dir, "/all_otu_tables.RData"))


all_phyloseq <- list()
for(i in 1:nrow(all_mpas)){
  all_phyloseq[[paste0(all_mpas$Condition[i], "_raw")]]
  all_phyloseq[[paste0(all_mpas$Condition[i], "_filt")]]
  all_phyloseq[[paste0(all_mpas$Condition[i], "_rarefied_min")]]
}



