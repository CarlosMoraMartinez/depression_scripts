########################################
# Read MetaData
########################################
  
riga_tandas45 <- read_tsv(opt$metadata_riga_45)

metadata <- data.frame(read_tsv(opt$metadata))%>% 
  mutate(
         imc_00 = peso_00/(0.01*talla_00)^2,
         imc_01 = peso_01/(0.01*talla_01)^2,
         edad_00_round = round(edad_00) %>% as.factor()
  )


### Modify metadata 
food_variable_names <- c(
  "energy_kcal"              = "ffq_energia_00",
  "carbohydrates_g"          = "ffq_h_carb_00",
  "fiber_g"                  = "ffq_fibra_00",
  "protein_g"                = "ffq_prot_00",
  "total_fat_g"              = "ffq_grasa_00",
  "dairy"              = "lacteos_00",
  "dairy_derivatives"  = "derivalac_00",
  "eggs"               = "huevos_00",
  "meat"               = "carnes_00",
  "fish"               = "pescados_00",
  "vegetables"         = "vegetales_00",
  "tubers"             = "tuberculos_00",
  "fruits"             = "frutas_00",
  "nuts"               = "frutosec_00",
  "oleaginous_fruits"  = "frutoleo_00",
  "refined_cereals"    = "cereref_00",
  "whole_grain_cereals"= "cereint_00",
  "legumes"            = "legum_00",
  "fats_oils"          = "grasas_00",
  "sweets_pastries"    = "dulces_bollpast_00",
  "sugars_and_sweets"  = "azucdulc_00",
  "snacks_savory"      = "snacks_00",
  "prepared_foods"     = "alimprepa_00",
  "sauces_condiments"  = "salscondi_00",
  "water"              = "agua_00",
  "juices_softdrinks"  = "refresc_00"
)

other_names <- c(
  "z_bmi_00" = "z_imc_00",
  "z_bmi_01" = "z_imc_01",
  "z_waist_00"="z_cintura_00",
  "z_waist_01" = "z_cintura_01",
  "mother_educ" = "educ_m_discrete",
  "age_months_T0" =  "edad_00_meses",
  "age_months_T1" =  "edad_01_meses"
)


metadata_with_origin  <- read_xlsx(opt$metadata_with_origin) 

metadata_gain_labels <- read_csv(opt$metadata_class) %>%  
  clean_names() %>% data.frame() %>% 
  dplyr::mutate(status_c1 = ifelse(status_c1 == "normal", status_c1, ifelse(z_t0 < 0, "low weight", "overweight")))
rownames(metadata) <- paste0("C", metadata$id)

meta_mother <- read_xlsx(opt$metadata_mother)
rownames(meta_mother) <- paste0("C", meta_mother$id)



colnames(metadata)[2] <- "Sex"
metadata$Sex[metadata$Sex == 0] <- "Boy"
metadata$Sex[metadata$Sex == 1] <- "Girl"
metadata$Sex <- factor(metadata$Sex)

colnames(metadata)[33] <- "Category_T0"
metadata$Category_T0[metadata$Category_T0 == 0] <- "Thinness"
metadata$Category_T0[metadata$Category_T0 == 1] <- "Normal"
metadata$Category_T0[metadata$Category_T0 == 2] <- "Overweight"
metadata$Category_T0[metadata$Category_T0 == 3] <- "Obesity"
metadata$Category_T0 <- factor(metadata$Category_T0, levels = c("Thinness", "Normal", "Overweight", "Obesity"))

colnames(metadata)[38] <- "Category_T1"
metadata$Category_T1[metadata$Category_T1 == 0] <- "Thinness"
metadata$Category_T1[metadata$Category_T1 == 1] <- "Normal"
metadata$Category_T1[metadata$Category_T1 == 2] <- "Overweight"
metadata$Category_T1[metadata$Category_T1 == 3] <- "Obesity"
metadata$Category_T1 <- factor(metadata$Category_T1, levels = c("Thinness", "Normal", "Overweight", "Obesity"))

na_metadata_t0 <- rownames(metadata[is.na(metadata$Category_T0),])
na_metadata_t1 <- rownames(metadata[is.na(metadata$imc_01),])
metadata_faltan <- metadata[! rownames(metadata) %in% colnames(s_otu_tab), ]

#metadata$tanda <- 2
metadata$tanda <- riga_tandas45$N.Batch[match(as.character(metadata$id), as.character(riga_tandas45$codk2))]
metadata$tanda[is.na(metadata$tanda)] <- 2

colnames(metadata)[1] <- 'sampleID'
nreads <- s_otu_tab %>% colSums()

metadata$hospital <- metadata_with_origin$hospital[match(metadata$sampleID, metadata_with_origin$id)]

newnames <- names(metadata_gain_labels)[!names(metadata_gain_labels) %in% names(metadata)]
newnames <- newnames[!(newnames %in% c("x", "sexo_vs"))] ## , "cat_peso_00", "cat_peso_01"

metadata_gain_labels2 <- metadata_gain_labels %>% select(id, all_of(newnames))
table(metadata_gain_labels2$id %in% metadata$sampleID)

metadata <- merge(metadata, metadata_gain_labels2, by.x="sampleID", by.y="id", all.x = TRUE) %>% 
  dplyr::mutate(is_normal = !is.na(status_c2)) %>% 
  dplyr::mutate(status_c2 = dplyr::recode(status_c2, Normal = "Normal", `Ganancia excesiva`="Excessive gain", `Ganancia insuficiente`="Insufficient gain")) %>% 
  mutate(sampleID2 = sampleID,
         sampleID = paste("C", sampleID, sep="")) 
rownames(metadata) <- metadata$sampleID
  

assertthat::assert_that(all(meta_mother$id %in% metadata$sampleID2 ))
names(meta_mother) %in% names(metadata) %>% any
nrow(metadata)
nrow(meta_mother)
metadata <- metadata %>% merge(meta_mother, by.x="sampleID2", by.y = "id")
nrow(metadata)
rownames(metadata) <- metadata$sampleID


metadata <- metadata %>%
  dplyr::mutate(educ_m_discrete =ifelse(edu_m_00 <4, "1-3",
                                        ifelse(edu_m_00 <7, "4-6", "7-10"))
  ) %>% 
  dplyr::rename(!!!food_variable_names) %>% 
  dplyr::rename(!!!other_names) %>% 
  dplyr::mutate(age_T0 =  edad_00) 

write_tsv(metadata, file = paste0(outdir, "metadata_full.tsv"))
metadata_full <- metadata
s_otu_tab_unfilt <- s_otu_tab

#Filter samples in metadata
common_samples <- names(s_otu_tab)[names(s_otu_tab) %in% rownames(metadata)]
table(names(s_otu_tab) %in% metadata$sampleID)
names(s_otu_tab)[! names(s_otu_tab) %in% metadata$sampleID]
table(metadata$sampleID %in% names(s_otu_tab))

metadata <- metadata[common_samples, ]
s_otu_tab <- s_otu_tab_unfilt[, common_samples]

all(rownames(metadata) == colnames(s_otu_tab))

nreads <- s_otu_tab %>% colSums()
metadata$nreads <- nreads[rownames(metadata)]
metadata$hospital[is.na(metadata$hospital)] <- "Zaragoza"
greads <- ggplot(metadata, aes(x=hospital, y = log10(nreads), fill=hospital))+geom_violin(alpha=0.6)+geom_boxplot(width=0.2, fill="lightgray")+ theme_bw()
ggsave(filename = paste0(outdir, "/reads_per_hospital.pdf"), greads, width = 7, height = 4)

all(names(s_otu_tab) == rownames(metadata))
write_tsv(s_otu_tab %>% rownames_to_column("taxon") %>% select(taxon, everything()), 
          file = paste0(outdir, "otu_tab_names_presentInMetadata.tsv"))


metadata$age_class1 <- factor(ifelse(metadata$edad_00_round %in% c(7, 8), ">7", as.character(metadata$edad_00_round)),
                                 levels=c("3", "4", "5", "6", ">7"))
metadata$age_class2 <- factor(ifelse(metadata$edad_00_round %in% c(7, 8), "7-8",
                                     ifelse(metadata$edad_00_round %in% c(3, 4), "3-4",
                                     as.character(metadata$edad_00_round))),
                              levels=c("3-4", "5", "6", "7-8"))
write_tsv(metadata, file = paste0(outdir, "metadata_presentInOtus.tsv"))

#Filter only normal weight at T0
if(opt$only_normal_weight){
  metadata <- metadata %>% dplyr::filter(Category_T0 == "Normal")
  s_otu_tab <- s_otu_tab[, metadata$sampleID]
  write_tsv(metadata, file = paste0(outdir, "metadata_onlyNormalAtT0.tsv"))
}
s_meta <- metadata

