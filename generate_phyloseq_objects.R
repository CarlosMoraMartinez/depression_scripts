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
metadata$nreads_filtPhylum <- nreads[rownames(metadata)]
greads <- ggplot(metadata, aes(x=hospital, y = log10(nreads), fill=hospital))+
  geom_violin(alpha=0.6)+
  geom_boxplot(width=0.2, fill="lightgray")+
  theme_bw()
ggsave(filename = paste0(outdir, "/reads_per_hospital_filtPhylum.pdf"), greads, width = 7, height = 4)

sample_data(pre_phyloseq_filt)$nreads_filtPhylum <- nreads[sample_data(pre_phyloseq_filt)$sampleID]
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

eliminadas_df <- metadata %>% 
  dplyr::filter(sampleID %in% muestras_eliminadas) %>% 
  dplyr::select(sampleID, hospital) %>% 
  dplyr::mutate(reads = nreads[ as.character(sampleID)]) %>% 
  dplyr::arrange(reads)
eliminadas_df %>% write_tsv(file=paste0(path_phyloseq, "/muestras_eliminadas_raref", as.character(opt$raref_quant), ".tsv"))

## Eliminar tanda 2
rmtanda2_fname <- paste0(path_phyloseq,'/pre_phyloseq_filt_noTandaZaragoza.RData')
standa1 <- metadata %>% dplyr::filter(hospital != "Zaragoza") %>% pull(sampleID) %>% as.character()
if(!file.exists(rmtanda2_fname) | opt$rewrite){
  pre_phyloseq_removet2 <- phyloseq::prune_samples(standa1, pre_phyloseq_filt) 
  save(pre_phyloseq_removet2, file = rmtanda2_fname)
}else{load(rmtanda2_fname)}

##  Eliminar tanda 2 - Rarefaction min
raref_min_filename_not2 <- paste0(path_phyloseq,'/pre_phyloseq_filt_noTandaZaragoza_rarefMin.RData')
if(!file.exists(raref_min_filename_not2) | opt$rewrite){
  pre_phyloseq_rarefied_not2 <-rarefy_even_depth(pre_phyloseq_removet2, rngseed = SEED)
  save(pre_phyloseq_rarefied_not2, file =raref_min_filename_not2)
}else{
  load(raref_min_filename_not2)
}


## Remove batch effect
# https://github.com/zhangyuqing/ComBat-seq
library(sva)
count_matrix <- otu_table(pre_phyloseq_filt)
data_matrix <- sample_data(pre_phyloseq_filt) %>% data.frame %>% dplyr::mutate(status_c2 = ifelse(is.na(status_c2), "initially_overweight", status_c2))
phseq_batch_tanda_fname <- paste0(path_phyloseq,'/pre_phyloseq_filt_Combat_Tanda2.RData')

if(!file.exists(phseq_batch_tanda_fname) | opt$rewrite){
  adjusted <- ComBat_seq(count_matrix, batch=data_matrix$hospital, group=data_matrix$status_c2)
  phseq_batch_tanda <- phyloseq(sample_data(data_matrix),
                                otu_table(adjusted, taxa_are_rows = TRUE),
                                tax_table(as.matrix(classification)))
  
  
  save(phseq_batch_tanda, file =phseq_batch_tanda_fname)
}else{
  load(phseq_batch_tanda_fname)
}

## Remove batch effect, with shrinkage
phseq_batch_tanda_shrink_fname <- paste0(path_phyloseq,'/pre_phyloseq_filt_Combat_Tanda2_shrink.RData')

if(!file.exists(phseq_batch_tanda_shrink_fname) | opt$rewrite){
  adjusted_shrink <- ComBat_seq(count_matrix, batch=data_matrix$hospital, group=data_matrix$status_c2, shrink = T)
  phseq_batch_tanda_shrink <- phyloseq(sample_data(data_matrix),
                                       otu_table(adjusted_shrink, taxa_are_rows = TRUE),
                                       tax_table(as.matrix(classification)))
  save(phseq_batch_tanda_shrink, file = phseq_batch_tanda_shrink_fname)
}else{
  load(phseq_batch_tanda_shrink_fname)
}

## Remove batch effect, rarefy
phseq_batch_tanda_raref_fname <- paste0(path_phyloseq,'/pre_phyloseq_filt_Combat_Tanda2_raref.RData')
if(!file.exists(phseq_batch_tanda_raref_fname) | opt$rewrite){
  phseq_batch_tanda_raref <-rarefy_even_depth(phseq_batch_tanda, rngseed = SEED)
  save(phseq_batch_tanda_raref, file=phseq_batch_tanda_raref_fname)
}else{
  load(phseq_batch_tanda_raref_fname)
}

## Remove batch effect with biological covariates --> DOES NOT WORK
# data_matrix2 <- data_matrix %>% dplyr::select(Sexo, Edad, IMC, Condition) %>% as.matrix
# s2remove <- is.na(data_matrix2) %>% rowSums 
# s2remove <- names(s2remove)[s2remove>0]
# data_matrix2 <- data_matrix2[! rownames(data_matrix2) %in% s2remove, ]
# data_matrix2 <- data_matrix2[! rownames(data_matrix2) %in% s2remove, ]
# count_matrix2 <- count_matrix[, !colnames(count_matrix) %in% s2remove]
# tanda2 <- data_matrix$Tanda[!data_matrix$sampleID %in% s2remove]
# 
# adjusted2 <- ComBat_seq(count_matrix2, batch=tanda2, group=data_matrix2)




## initially normal only
# https://github.com/zhangyuqing/ComBat-seq
library(sva)
count_matrix <- otu_table(pre_phyloseq_filt)
data_matrix <- sample_data(pre_phyloseq_filt) %>% data.frame %>% filter(!is.na(status_c2))
count_matrix <- count_matrix[, data_matrix$sampleID]
phseq_fname <- paste0(path_phyloseq,'/pre_phyloseq_filt_onlyNormalT0.RData')

if(!file.exists(phseq_fname) | opt$rewrite){
  phseq_onlyNorT0 <- phyloseq(sample_data(data_matrix),
                                          otu_table(count_matrix, taxa_are_rows = TRUE),
                                          tax_table(as.matrix(classification)))
  
  
  save(phseq_onlyNorT0, file =phseq_fname)
}else{
  load(phseq_fname)
}

## Remove batch effect AND initially normal only
phseq_batch_tanda_fname <- paste0(path_phyloseq,'/pre_phyloseq_filt_Combat_onlyNormalT0.RData')
if(!file.exists(phseq_batch_tanda_fname) | opt$rewrite){
  
  adjusted_onlyNorT0 <- ComBat_seq(count_matrix, batch=data_matrix$hospital, group=data_matrix$status_c2)
  phseq_batch_tanda_onlyNorT0 <- phyloseq(sample_data(data_matrix),
                                otu_table(adjusted_onlyNorT0, taxa_are_rows = TRUE),
                                tax_table(as.matrix(classification)))
  
  
  save(phseq_batch_tanda_onlyNorT0, file =phseq_batch_tanda_fname)
}else{
  load(phseq_batch_tanda_fname)
}

## Remove batch effect, with shrinkage, initially normal only
phseq_batch_tanda_shrink_fname <- paste0(path_phyloseq,'/pre_phyloseq_filt_Combat_onlyNormalT0_shrink.RData')

if(!file.exists(phseq_batch_tanda_shrink_fname) | opt$rewrite){
  adjusted_shrink_onlyNorT0 <- ComBat_seq(count_matrix, batch=data_matrix$hospital, group=data_matrix$status_c2, shrink = T)
  phseq_batch_tanda_shrink_onlyNorT0 <- phyloseq(sample_data(data_matrix),
                                       otu_table(adjusted_shrink_onlyNorT0, taxa_are_rows = TRUE),
                                       tax_table(as.matrix(classification)))
  save(phseq_batch_tanda_shrink_onlyNorT0, file = phseq_batch_tanda_shrink_fname)
}else{
  load(phseq_batch_tanda_shrink_fname)
}

## Remove batch effect, rarefy, initially normal only
phseq_batch_tanda_raref_fname <- paste0(path_phyloseq,'/pre_phyloseq_filt_Combat_onlyNormalT0_raref.RData')
if(!file.exists(phseq_batch_tanda_raref_fname) | opt$rewrite){
  phseq_batch_tanda_raref_onlyNorT0 <-rarefy_even_depth(phseq_batch_tanda_onlyNorT0, rngseed = SEED)
  save(phseq_batch_tanda_raref_onlyNorT0, file=phseq_batch_tanda_raref_fname)
}else{
  load(phseq_batch_tanda_raref_fname)
}

###############################################3
## Phyloseq list


allphyloseqlist_fname <- paste0(path_phyloseq, "/phyloseq_all_list.RData")
if(!file.exists(allphyloseqlist_fname) | opt$rewrite){
  all_phyloseq <- list(#raw = pre_phyloseq1, 
    filt = pre_phyloseq_filt, 
    rarefied_min = pre_phyloseq_rarefied, 
    rarefied_quant = pre_phyloseq_rarefied2, 
    remove_tanda2 = pre_phyloseq_removet2, 
    remove_tanda2_rarefied_min = pre_phyloseq_rarefied_not2,
    rmbatch_tanda =phseq_batch_tanda,
    rmbatch_tanda_shrink = phseq_batch_tanda_shrink,
    rmbatch_tanda_raref =phseq_batch_tanda_raref,
    
    rmbatch_onlyNorm0 =phseq_batch_tanda_onlyNorT0,
    rmbatch_onlyNorm0_shrink = phseq_batch_tanda_shrink_onlyNorT0,
    rmbatch_onlyNorm0_raref =phseq_batch_tanda_raref_onlyNorT0
    
  )
  save(all_phyloseq, file=allphyloseqlist_fname)
}else{
  load(allphyloseqlist_fname)
}
