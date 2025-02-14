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

ps_bracken_species_species <- phyloseq(sample_data(s_meta),
                                       otu_table(s_otu_tab_sp, taxa_are_rows = TRUE),
                                       tax_table(as.matrix(classification_sp)))
pre_phyloseq <- ps_bracken_species
save(file=paste0(path_phyloseq, "/phyloseq_object_analysis1_strain.RData"), ps_bracken_species)
save(file=paste0(path_phyloseq, "/phyloseq_object_analysis1_summedSpecies.RData"), ps_bracken_species_species)

filterPhyla <- c("Chloroplast", "Mitochondria", "Eukaryota", "Metazoa", "Viruses")

get_filtered_phyloseq <- function(pre_phyloseq1, phseqname="", filterPhyla){
  filterPhyla <- NA
  pre_phyloseq1 = subset_taxa(pre_phyloseq1, !Phylum %in% c(NA))
  
  pre_phyloseq1 <- subset_taxa(pre_phyloseq1, !Kingdom %in% filterPhyla)
  pre_phyloseq1 <- subset_taxa(pre_phyloseq1, !Phylum %in% filterPhyla)
  pre_phyloseq1 <- subset_taxa(pre_phyloseq1, !Class %in% filterPhyla)
  pre_phyloseq1 <- subset_taxa(pre_phyloseq1, !Order %in% filterPhyla) # 12 a nivel Order
  pre_phyloseq1 <- subset_taxa(pre_phyloseq1, !Family %in% filterPhyla) # 7 a nivel Family
  pre_phyloseq1 <- subset_taxa(pre_phyloseq1, !Genus %in% filterPhyla)

  save(file=paste0(path_phyloseq, "/phyloseq_object_raw_filt_by_Phylum_", phseqname,".RData"), pre_phyloseq1)
  
  all_phyloseq_tmp <- list()
  all_phyloseq_tmp[[paste0("raw_", phseqname)]] <- pre_phyloseq1
  
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
  write_tsv(pre_prevalence, paste0(opt$out, "/raw_prevalence_", phseqname, ".tsv"))

  ## Filtered to frequency
  prevalenceThreshold = opt$minfreq * nsamples(pre_phyloseq1)
  keepTaxa = rownames(pre_prevalence)[(pre_prevalence$Prevalence >= prevalenceThreshold)]
  (pre_phyloseq_filt = prune_taxa(keepTaxa, pre_phyloseq1))
  filtered_phyloseq_filename <- paste0(path_phyloseq,'/pre_phyloseq_filt_by_prevalence', as.character(100*opt$minfreq), phseqname,'.RData')
  save(pre_phyloseq_filt, file = filtered_phyloseq_filename)

  #Reads before rarefeact
  nreads <- otu_table(pre_phyloseq_filt) %>% colSums()
  s_meta$nreads_filt <- nreads[s_meta$sampleID]
  write_tsv(s_meta, paste0(outdir, "/", phseqname, "_full_metadata2.tsv"))


  sample_data(pre_phyloseq_filt)$nreads_filt <- nreads[sample_data(pre_phyloseq_filt)$sampleID]
  
  all_phyloseq_tmp[[paste0("filt_", phseqname)]] <- pre_phyloseq_filt
  
  return(all_phyloseq_tmp)

}

addRarefied <- function(phobj, phseqname, phlist){
  raref_min_filename <- paste0(path_phyloseq,'/pre_phyloseq_filt', as.character(100*opt$minfreq), '_rarefMin_', phseqname, '.RData')
  if(!file.exists(raref_min_filename) | opt$rewrite){
    pre_phyloseq_rarefied <-rarefy_even_depth(phobj, rngseed = SEED)
    save(pre_phyloseq_rarefied, file =raref_min_filename)
  }else{
    load(raref_min_filename)
  } 
  phlist[[paste0(phseqname, "_rarefied_min")]] <- pre_phyloseq_rarefied
  return(phlist)
}

addPrunedSamples <- function(phobj, phseqname, phlist, samples){
  rmtanda2_fname <- paste0(path_phyloseq,'/pre_phyloseq_filt_noTandaZaragoza', phseqname, '.RData')
  if(!file.exists(rmtanda2_fname) | opt$rewrite){
    pre_phyloseq_removet2 <- phyloseq::prune_samples(samples, phobj) 
    save(pre_phyloseq_removet2, file = rmtanda2_fname)
  }else{load(rmtanda2_fname)}
  
  phlist[[paste0(phseqname, "_rmTanda")]] <- pre_phyloseq_removet2
  return(phlist)
}

lista_strain <- get_filtered_phyloseq(ps_bracken_species, "strain", filterPhyla)
lista_spsum <- get_filtered_phyloseq(ps_bracken_species_species, "spsum", filterPhyla)

## Rarefaction min

lista_strain <- addRarefied(lista_strain$filt_strain, names(lista_strain)[2], lista_strain)
lista_spsum <- addRarefied(lista_spsum$filt_spsum, names(lista_spsum)[2], lista_spsum)

## Eliminar tanda 2

standa1 <- metadata %>% dplyr::filter(hospital != "Zaragoza") %>% pull(sampleID) %>% as.character()

lista_strain <- addPrunedSamples(lista_strain$filt_strain, names(lista_strain)[2], lista_strain, standa1)
lista_spsum <- addPrunedSamples(lista_spsum$filt_spsum, names(lista_spsum)[2], lista_spsum, standa1)


lista_strain <- addRarefied(lista_strain$filt_strain_rmTanda, names(lista_strain)[4], lista_strain)
lista_spsum <- addRarefied(lista_spsum$filt_spsum_rmTanda, names(lista_spsum)[4], lista_spsum)

table(names(lista_strain) %in% names(lista_spsum))
all_phyloseq <- append(lista_strain, lista_spsum)
allphyloseqlist_fname <- paste0(path_phyloseq, "/phyloseq_all_list.RData")
save(all_phyloseq, file=allphyloseqlist_fname)



######################################################################
#### de aqui hacia abajo: de momento no

## Remove batch effect
# https://github.com/zhangyuqing/ComBat-seq
library(sva)
count_matrix <- otu_table(pre_phyloseq_filt)
data_matrix <- sample_data(pre_phyloseq_filt) %>% data.frame %>% dplyr::mutate(status_c2 = ifelse(is.na(status_c2), "initially_overweight", status_c2))
phseq_batch_tanda_fname <- paste0(path_phyloseq,'/pre_phyloseq_filt_Combat_Tanda2.RData')
data_matrix$hospital[is.na(data_matrix$hospital)] <- "NA"
if(!file.exists(phseq_batch_tanda_fname) | opt$rewrite){
  adjusted <- ComBat_seq(count_matrix, batch=data_matrix$hospital, group=data_matrix$status_c2)
  phseq_batch_tanda <- phyloseq(sample_data(data_matrix),
                                otu_table(adjusted, taxa_are_rows = TRUE),
                                tax_table(as.matrix(classification)))
  
  
  save(phseq_batch_tanda, file =phseq_batch_tanda_fname)
}else{
  load(phseq_batch_tanda_fname)
}

## Remove batch effect with age as biological factor
count_matrix <- otu_table(pre_phyloseq_filt)
data_matrix <- sample_data(pre_phyloseq_filt) %>% data.frame %>% 
  dplyr::mutate(status_c2 = ifelse(is.na(status_c2), "initially_overweight", status_c2))
phseq_batch_tanda_fname2 <- paste0(path_phyloseq,'/pre_phyloseq_filt_Combat_Tanda2_AgeBF.RData')
data_matrix$hospital[is.na(data_matrix$hospital)] <- "NA"
if(!file.exists(phseq_batch_tanda_fname2) | opt$rewrite){
  adjusted2 <- ComBat_seq(count_matrix, batch=data_matrix$hospital, group=data_matrix$age_class1)
  phseq_batch_tanda_age <- phyloseq(sample_data(data_matrix),
                                otu_table(adjusted2, taxa_are_rows = TRUE),
                                tax_table(as.matrix(classification)))
  
  
  save(phseq_batch_tanda_age, file =phseq_batch_tanda_fname2)
}else{
  load(phseq_batch_tanda_fname2)
}

## Remove batch effect, with shrinkage
#phseq_batch_tanda_shrink_fname <- paste0(path_phyloseq,'/pre_phyloseq_filt_Combat_Tanda2_shrink.RData')
#
#if(!file.exists(phseq_batch_tanda_shrink_fname) | opt$rewrite){
#  adjusted_shrink <- ComBat_seq(count_matrix, batch=data_matrix$hospital, group=data_matrix$status_c2, shrink = T)
#  phseq_batch_tanda_shrink <- phyloseq(sample_data(data_matrix),
#                                       otu_table(adjusted_shrink, taxa_are_rows = TRUE),
#                                       tax_table(as.matrix(classification)))
#  save(phseq_batch_tanda_shrink, file = phseq_batch_tanda_shrink_fname)
#}else{
#  load(phseq_batch_tanda_shrink_fname)
#}

## Remove batch effect in rarefied samples
#phseq_rarthenbatch_tanda_fname <- paste0(path_phyloseq,'/pre_phyloseq_raref_then_Combat.RData')
#count_matrix <- otu_table(pre_phyloseq_rarefied)
#data_matrix <- sample_data(pre_phyloseq_rarefied) %>% data.frame %>% 
#  dplyr::mutate(status_c2 = ifelse(is.na(status_c2), "initially_overweight", status_c2),
#                hospital = ifelse(is.na(hospital), "NA", hospital))

#if(!file.exists(phseq_rarthenbatch_tanda_fname) | opt$rewrite){
#  adjusted <- ComBat_seq(count_matrix, batch=data_matrix$hospital, group=data_matrix$status_c2)
#  phseq_rerefthenbatch_tanda <- phyloseq(sample_data(data_matrix),
#                                otu_table(adjusted, taxa_are_rows = TRUE),
#                                tax_table(as.matrix(classification)))
  
  
#  save(phseq_batch_tanda, file =phseq_rarthenbatch_tanda_fname)
#}else{
#  load(phseq_batch_tanda_fname)
#}

## Remove batch effect, rarefy
#phseq_batch_tanda_raref_fname <- paste0(path_phyloseq,'/pre_phyloseq_filt_Combat_Tanda2_raref.RData')
#if(!file.exists(phseq_batch_tanda_raref_fname) | opt$rewrite){
#  phseq_batch_tanda_raref <-rarefy_even_depth(phseq_batch_tanda, rngseed = SEED)
#  save(phseq_batch_tanda_raref, file=phseq_batch_tanda_raref_fname)
#}else{
#  load(phseq_batch_tanda_raref_fname)
#}

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
#library(sva)
#count_matrix <- otu_table(pre_phyloseq_filt)
#data_matrix <- sample_data(pre_phyloseq_filt) %>% data.frame %>% filter(!is.na(status_c2))
#count_matrix <- count_matrix[, data_matrix$sampleID]
#phseq_fname <- paste0(path_phyloseq,'/pre_phyloseq_filt_onlyNormalT0.RData')
#
#if(!file.exists(phseq_fname) | opt$rewrite){
#  phseq_onlyNorT0 <- phyloseq(sample_data(data_matrix),
#                              otu_table(count_matrix, taxa_are_rows = TRUE),
#                              tax_table(as.matrix(classification)))
#  
#  
#  save(phseq_onlyNorT0, file =phseq_fname)
#}else{
#  load(phseq_fname)
#}

## Remove batch effect AND initially normal only
#phseq_batch_tanda_fname <- paste0(path_phyloseq,'/pre_phyloseq_filt_Combat_onlyNormalT0.RData')
#if(!file.exists(phseq_batch_tanda_fname) | opt$rewrite){
#  
#  adjusted_onlyNorT0 <- ComBat_seq(count_matrix, batch=data_matrix$hospital, group=data_matrix$status_c2)
#  phseq_batch_tanda_onlyNorT0 <- phyloseq(sample_data(data_matrix),
#                                          otu_table(adjusted_onlyNorT0, taxa_are_rows = TRUE),
#                                          tax_table(as.matrix(classification)))
#  
#  
#  save(phseq_batch_tanda_onlyNorT0, file =phseq_batch_tanda_fname)
#}else{
#  load(phseq_batch_tanda_fname)
#}
#
### Remove batch effect, with shrinkage, initially normal only
#phseq_batch_tanda_shrink_fname <- paste0(path_phyloseq,'/pre_phyloseq_filt_Combat_onlyNormalT0_shrink.RData')
#
#if(!file.exists(phseq_batch_tanda_shrink_fname) | opt$rewrite){
#  adjusted_shrink_onlyNorT0 <- ComBat_seq(count_matrix, batch=data_matrix$hospital, group=data_matrix$status_c2, shrink = T)
#  phseq_batch_tanda_shrink_onlyNorT0 <- phyloseq(sample_data(data_matrix),
#                                                 otu_table(adjusted_shrink_onlyNorT0, taxa_are_rows = TRUE),
#                                                 tax_table(as.matrix(classification)))
#  save(phseq_batch_tanda_shrink_onlyNorT0, file = phseq_batch_tanda_shrink_fname)
#}else{
#  load(phseq_batch_tanda_shrink_fname)
#}
#
### Remove batch effect, rarefy, initially normal only
#phseq_batch_tanda_raref_fname <- paste0(path_phyloseq,'/pre_phyloseq_filt_Combat_onlyNormalT0_raref.RData')
#if(!file.exists(phseq_batch_tanda_raref_fname) | opt$rewrite){
#  phseq_batch_tanda_raref_onlyNorT0 <-rarefy_even_depth(phseq_batch_tanda_onlyNorT0, rngseed = SEED)
#  save(phseq_batch_tanda_raref_onlyNorT0, file=phseq_batch_tanda_raref_fname)
#}else{
#  load(phseq_batch_tanda_raref_fname)
#}

###############################################3
## Phyloseq list

