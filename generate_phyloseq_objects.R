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
load(paste0(path_phyloseq, "/phyloseq_object_analysis1.RData"))

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
load(paste0(path_phyloseq, "/phyloseq_object_raw_filt_by_Phylum.RData"))

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
prevalenceThreshold = opt$minfreq * nsamples(pre_phyloseq1)
keepTaxa = rownames(pre_prevalence)[(pre_prevalence$Prevalence >= prevalenceThreshold)]
(pre_phyloseq_filt = prune_taxa(keepTaxa, pre_phyloseq1))
filtered_phyloseq_filename <- paste0(path_phyloseq,'/pre_phyloseq_filt_by_prevalence', as.character(100*opt$minfreq), '.RData')
save(pre_phyloseq_filt, file = filtered_phyloseq_filename)
load(filtered_phyloseq_filename)

#Reads before rarefeact
nreads <- otu_table(pre_phyloseq_filt) %>% colSums()
s_meta$nreads_filt <- nreads[s_meta$sampleID]
write_tsv(s_meta, paste0(outdir, "/full_metadata2.tsv"))


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
#min_depth <- otu_table(pre_phyloseq_filt) %>% colSums() %>% quantile(opt$raref_quant)
#raref_quant_filename <- paste0(path_phyloseq,'/pre_phyloseq_filt_raref_quant', as.character(100*opt$raref_quant), '.RData')
#if(!file.exists(raref_quant_filename) | opt$rewrite){
#  pre_phyloseq_rarefied2 <-rarefy_even_depth(pre_phyloseq_filt, sample.size = min_depth, rngseed = SEED)
#  save(pre_phyloseq_rarefied2, file =raref_quant_filename)
#}else{
#  load(raref_quant_filename)
#}
#muestras_eliminadas <- sample_names(pre_phyloseq_filt)[!sample_names(pre_phyloseq_filt) %in% sample_names(pre_phyloseq_rarefied2)] 

#eliminadas_df <- s_meta %>% 
#  dplyr::filter(sampleID %in% muestras_eliminadas) %>% 
#  dplyr::arrange(nreads_filt)
#eliminadas_df %>% write_tsv(file=paste0(path_phyloseq, "/muestras_eliminadas_raref", as.character(opt$raref_quant), ".tsv"))


## Eliminar tanda 2
rmtanda2_fname <- paste0(path_phyloseq,'/pre_phyloseq_filt_noTandaZaragoza.RData')
standa1 <- metadata %>% dplyr::filter(hospital != "Zaragoza") %>% pull(sampleID) %>% as.character()
if(!file.exists(rmtanda2_fname) | opt$rewrite){
  pre_phyloseq_removet2 <- phyloseq::prune_samples(standa1, pre_phyloseq_filt) 
  save(pre_phyloseq_removet2, file = rmtanda2_fname)
}else{
  load(rmtanda2_fname)
  }

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



allphyloseqlist_fname <- paste0(path_phyloseq, "/phyloseq_all_list.RData")
if(!file.exists(allphyloseqlist_fname) | opt$rewrite){
  all_phyloseq <- list(
    raw = pre_phyloseq1, 
    filt = pre_phyloseq_filt, 
    rarefied_min = pre_phyloseq_rarefied, 
    #rarefied_quant = pre_phyloseq_rarefied2,
    remove_tanda2 = pre_phyloseq_removet2,
    remove_tanda2_rarefied_min = pre_phyloseq_rarefied_not2
    #rmbatch_tanda =phseq_batch_tanda,
    #rmbatch_ageBF = phseq_batch_tanda_age
    #rmbatch_tanda_shrink = phseq_batch_tanda_shrink,
    #phseq_rerefthenbatch_tanda = phseq_rerefthenbatch_tanda
    #rmbatch_tanda_raref =phseq_batch_tanda_raref,
    
    #rmbatch_onlyNorm0 =phseq_batch_tanda_onlyNorT0,
    #rmbatch_onlyNorm0_shrink = phseq_batch_tanda_shrink_onlyNorT0,
    #rmbatch_onlyNorm0_raref =phseq_batch_tanda_raref_onlyNorT0
    
    
  )
  save(all_phyloseq, file=allphyloseqlist_fname)
}else{
  load(allphyloseqlist_fname)
}


## Modify exercise:



### Modify metadata 
#food_variable_names <- c(
#  "energy_kcal"              = "ffq_energia_00",
#  "carbohydrates_g"          = "ffq_h_carb_00",
#  "fiber_g"                  = "ffq_fibra_00",
#  "protein_g"                = "ffq_prot_00",
#  "total_fat_g"              = "ffq_grasa_00",
#  "dairy"              = "lacteos_00",
#  "dairy_derivatives"  = "derivalac_00",
#  "eggs"               = "huevos_00",
#  "meat"               = "carnes_00",
#  "fish"               = "pescados_00",
#  "vegetables"         = "vegetales_00",
#  "tubers"             = "tuberculos_00",
#  "fruits"             = "frutas_00",
#  "nuts"               = "frutosec_00",
#  "oleaginous_fruits"  = "frutoleo_00",
#  "refined_cereals"    = "cereref_00",
#  "whole_grain_cereals"= "cereint_00",
#  "legumes"            = "legum_00",
#  "fats_oils"          = "grasas_00",
#  "sweets_pastries"    = "dulces_bollpast_00",
#  "sugars_and_sweets"  = "azucdulc_00",
#  "snacks_savory"      = "snacks_00",
#  "prepared_foods"     = "alimprepa_00",
#  "sauces_condiments"  = "salscondi_00",
#  "water"              = "agua_00",
#  "juices_softdrinks"  = "refresc_00"
#)
#
#other_names <- c(
# "z_bmi_00" = "z_imc_00",
# "z_bmi_01" = "z_imc_01",
# "z_waist_00"="z_cintura_00",
# "z_waist_01" = "z_cintura_01",
# "mother_educ" = "educ_m_discrete",
# "age_months_T0" =  "edad_00_meses",
# "age_months_T1" =  "edad_01_meses"
#)
#
#
#metadata_t0cat <- read_csv("/home/carlos/Documentos/CORALS/METADATA/classified_kids_NEWDATA_PROVISIONAL_withZval_unfiltered.csv") %>% 
#  mutate(status_c1 = ifelse(status_c1 == "normal", status_c1, ifelse(Z_t0_ < 0, "low weight", "overweight"))) %>% 
#  dplyr::mutate(Status_c2 = ifelse(is.na(Z_t0_) | is.na(Z_t1_), NA, Status_c2))
#
#ggplot(metadata_t0cat, aes(x=status_c1, y=Z_t0_, col=status_c1)) +
#  facet_grid(. ~ edad_00_meses>=61)+
#  geom_point() + 
#  theme_bw()
#ss <- metadata_t0cat %>% mutate(agroup = ifelse(edad_00_meses>=61, "older", "younger")) %>% 
#  group_by(agroup, status_c1) %>% 
#  dplyr::summarise(minz = min(Z_t0_, na.rm = TRUE),
#            maxz = max(Z_t0_, na.rm = TRUE))
#
#metadata_t0cat
#for(phname in names(all_phyloseq)){
#  cat(phname, "\n")
#  metadata <- sample_data(all_phyloseq[[phname]]) %>% data.frame()
#  assertthat::assert_that(all(food_variable_names %in% names(metadata)))
#  sample_data(all_phyloseq[[phname]]) <- metadata %>%
#    dplyr::rename(!!!food_variable_names) %>% 
#    dplyr::rename(!!!other_names) %>% 
#    dplyr::mutate(age_T0 =  edad_00) %>% 
#    sample_data()
#}
#allphyloseqlist_fname <- paste0(path_phyloseq, "/phyloseq_all_list_modNames.RData")
#save(all_phyloseq, file=allphyloseqlist_fname)
