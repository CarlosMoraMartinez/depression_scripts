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
s_meta$nreads_filt <- nreads[s_meta$sample]
write_tsv(meta3, paste0(input_tabs_dir, "/full_metadata2.tsv"))

(greads <- ggplot(s_meta, aes(x=Group,
                            y = log10(nreads), 
                            fill=Group,
                            col=Group))+
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

#### ONLY TRIBIOME
## Calculate prevalence
phseq_only_tribiome <- phyloseq::subset_samples(pre_phyloseq1, Group == "Control")

ottmp2 <- phyloseq::otu_table(phseq_only_tribiome)
pre_prevalence2 <- apply(X = ottmp2,
                        MARGIN = ifelse(taxa_are_rows(phseq_only_tribiome), yes = 1, no = 2),
                        FUN = function(x){sum(x > opt$mincountspersample)})
pre_prevalence2 = data.frame(Prevalence = pre_prevalence2,
                            TotalAbundance = phyloseq::taxa_sums(phseq_only_tribiome),
                            tax_table(phseq_only_tribiome), 
                            relative_prevalence = pre_prevalence2/ nsamples(phseq_only_tribiome)
)
write_tsv(pre_prevalence2, paste0(opt$out, "/raw_prevalence_onlyTribiome.tsv"))
prevalenceThreshold2 = opt$minfreq * nsamples(phseq_only_tribiome)
keepTaxa = rownames(pre_prevalence2)[(pre_prevalence2$Prevalence >= prevalenceThreshold2)]
(phseq_only_tribiome_filt = prune_taxa(keepTaxa, phseq_only_tribiome))
 
## Sibo + Tbm thin
phseq_only_sibonormw <- phyloseq::subset_samples(pre_phyloseq1, Weight != "Overweight")

ottmp2 <- phyloseq::otu_table(phseq_only_sibonormw)
pre_prevalence2 <- apply(X = ottmp2,
                         MARGIN = ifelse(taxa_are_rows(phseq_only_sibonormw), yes = 1, no = 2),
                         FUN = function(x){sum(x > opt$mincountspersample)})
pre_prevalence2 = data.frame(Prevalence = pre_prevalence2,
                             TotalAbundance = phyloseq::taxa_sums(phseq_only_sibonormw),
                             tax_table(phseq_only_sibonormw), 
                             relative_prevalence = pre_prevalence2/ nsamples(phseq_only_sibonormw)
)
write_tsv(pre_prevalence2, paste0(opt$out, "/raw_prevalence_onlySiboAndTribiomeNormalWeight.tsv"))
prevalenceThreshold2 = opt$minfreq * nsamples(phseq_only_sibonormw)
keepTaxa = rownames(pre_prevalence2)[(pre_prevalence2$Prevalence >= prevalenceThreshold2)]
(phseq_only_sibonormw_filt = prune_taxa(keepTaxa, phseq_only_sibonormw))


## Aggregated:


all_phyloseq <- list(all_raw = pre_phyloseq1, 
                     all_filt = pre_phyloseq_filt, 
                     tribiome_raw = phseq_only_tribiome, 
                     tribiome_filt = phseq_only_tribiome_filt,
                     normalw_raw = phseq_only_sibonormw,
                     normalw_raw_filt = phseq_only_sibonormw_filt
                     )


for(lev in c("Phylum", "Genus", "Species")){
  cat(lev, "\n")
  all_phyloseq[[paste0("all_raw_", lev)]] <- tax_glom_custom(pre_phyloseq1, taxrank = lev)
  all_phyloseq[[paste0("tribiome_raw_", lev)]] <- tax_glom_custom(phseq_only_tribiome, taxrank = lev)
  all_phyloseq[[paste0("normalw_raw_", lev)]] <- tax_glom_custom(phseq_only_sibonormw, taxrank = lev)
}

for(nn in names(all_phyloseq)){
  all_phyloseq[[paste0(nn, "_rarefied_min")]] <- rarefy_even_depth(all_phyloseq[[nn]], rngseed = SEED)
}


allphyloseqlist_fname <- paste0(path_phyloseq, "/phyloseq_all_list.RData")
save(all_phyloseq, file=allphyloseqlist_fname)
