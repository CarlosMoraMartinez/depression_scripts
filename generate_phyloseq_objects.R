########################################
# Generate all Phyloseqs 
########################################

########################################
# Generate Phyloseq basic
########################################

path_phyloseq <- paste0(opt$out, "/phyloseq")
if(! dir.exists(path_phyloseq)){dir.create(path_phyloseq)}

ps_bracken_species <- phyloseq(sample_data(s_meta),
                               otu_table(s_otu_tab, taxa_are_rows = TRUE),
                               tax_table(as.matrix(classification)))
pre_phyloseq <- ps_bracken_species
save(file=paste0(path_phyloseq, "/phyloseq_object_analysis1.RData"), ps_bracken_species)

filterPhyla <- NA
(pre_phyloseq1 = subset_taxa(pre_phyloseq, !Phylum %in% filterPhyla))
filterPhyla <- c("Chloroplast", "Mitochondria", "Eukaryota", "Metazoa", "Viruses")
pre_phyloseq1 <- subset_taxa(pre_phyloseq1, !Kingdom %in% filterPhyla)
pre_phyloseq1 <- subset_taxa(pre_phyloseq1, !Phylum %in% filterPhyla)
pre_phyloseq1 <- subset_taxa(pre_phyloseq1, !Class %in% filterPhyla)
pre_phyloseq1 <- subset_taxa(pre_phyloseq1, !Order %in% filterPhyla) # 12 a nivel Order
pre_phyloseq1 <- subset_taxa(pre_phyloseq1, !Family %in% filterPhyla) # 7 a nivel Family
pre_phyloseq1 <- subset_taxa(pre_phyloseq1, !Genus %in% filterPhyla)

save(file=paste0(path_phyloseq, "/phyloseq_object_raw_filt_by_Phylum.RData"), pre_phyloseq1)

all_phyloseq <- list(raw = pre_phyloseq1)

## Calculate prevalence
ottmp <- phyloseq::otu_table(pre_phyloseq1)
pre_prevalence <- apply(X = ottmp,
                        MARGIN = ifelse(taxa_are_rows(pre_phyloseq1), yes = 1, no = 2),
                        FUN = function(x){sum(x > opt$mincountspersample)})
pre_prevalence = data.frame(Prevalence = pre_prevalence,
                            TotalAbundance = phyloseq::taxa_sums(pre_phyloseq1),
                            tax_table(pre_phyloseq1), 
                            relative_prevalence = pre_prevalence/ nsamples(pre_phyloseq1)
)
write_tsv(pre_prevalence, paste0(opt$out, "/raw_prevalence.tsv"))

## Filtered to frequency
prevalenceThreshold = opt$minfreq * nsamples(all_phyloseq$raw)
keepTaxa = rownames(pre_prevalence)[(pre_prevalence$Prevalence >= prevalenceThreshold)]
(pre_phyloseq_filt = prune_taxa(keepTaxa, pre_phyloseq1))
filtered_phyloseq_filename <- paste0(path_phyloseq,'/pre_phyloseq_filt_by_prevalence', as.character(100*opt$minfreq), '.RData')
save(pre_phyloseq_filt, file = filtered_phyloseq_filename)

#Reads before rarefeact
nreads <- otu_table(pre_phyloseq_filt) %>% colSums()
meta3$nreads_filt <- nreads[meta3$sample]
write_tsv(meta3, paste0(input_tabs_dir, "/full_metadata2.tsv"))

(greads <- ggplot(meta3, aes(x=Treatment,
                            y = log10(nreads), 
                            fill=Treatment,
                            col=Treatment))+
    facet_grid(Stress ~ Region_sequenced) +
  #geom_violin(alpha=0.6)+
  geom_boxplot(width=0.7, fill="white")+
  geom_point() +
  scale_color_npg() +
    theme_bw() +
    #scale_fill_npg() +
    theme(axis.text.x = element_text(vjust=1,hjust=1, angle = 45))
  
)
ggsave(filename = paste0(opt$out, "/reads_per_hospital_filtPhylum.pdf"), greads, width = 7, height = 4)

sample_data(pre_phyloseq_filt)$nreads_filt <- nreads[sample_data(pre_phyloseq_filt)$sampleID]
## Rarefaction min
raref_min_filename <- paste0(path_phyloseq,'/pre_phyloseq_filt', as.character(100*opt$minfreq), '_rarefMin.RData')
if(!file.exists(raref_min_filename) | opt$rewrite){
  pre_phyloseq_rarefied <-rarefy_even_depth(pre_phyloseq_filt, rngseed = SEED)
  save(pre_phyloseq_rarefied, file =raref_min_filename)
}else{
  load(raref_min_filename)
}

## Rarefaction 0.15
min_depth <- otu_table(pre_phyloseq_filt) %>% colSums() %>% quantile(opt$raref_quant)
raref_quant_filename <- paste0(path_phyloseq,'/pre_phyloseq_filt_raref_quant', as.character(100*opt$raref_quant), '.RData')
if(!file.exists(raref_quant_filename) | opt$rewrite){
  pre_phyloseq_rarefied2 <-rarefy_even_depth(pre_phyloseq_filt, sample.size = min_depth, rngseed = SEED)
  save(pre_phyloseq_rarefied2, file =raref_quant_filename)
}else{
  load(raref_quant_filename)
}
muestras_eliminadas <- sample_names(pre_phyloseq_filt)[!sample_names(pre_phyloseq_filt) %in% sample_names(pre_phyloseq_rarefied2)] 

eliminadas_df <- meta3 %>% 
  dplyr::filter(sampleID %in% muestras_eliminadas) %>% 
  dplyr::arrange(nreads_filt)
eliminadas_df %>% write_tsv(file=paste0(path_phyloseq, "/muestras_eliminadas_raref", as.character(opt$raref_quant), ".tsv"))


allphyloseqlist_fname <- paste0(path_phyloseq, "/phyloseq_all_list.RData")
if(!file.exists(allphyloseqlist_fname) | opt$rewrite){
  all_phyloseq <- list(#raw = pre_phyloseq1, 
    filt = pre_phyloseq_filt, 
    rarefied_min = pre_phyloseq_rarefied, 
    rarefied_quant = pre_phyloseq_rarefied2
  )
  save(all_phyloseq, file=allphyloseqlist_fname)
}else{
  load(allphyloseqlist_fname)
}
