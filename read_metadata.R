########################################
# Read MetaData
########################################

metadata <- data.frame(read_excel(opt$metadata))
metadata_with_origin  <- read_xlsx(opt$metadata_with_origin)
metadata_gain_labels <- read_csv(opt$metadata_class) %>%  clean_names() %>% data.frame()
rownames(metadata) <- paste0("C", metadata$id)

colnames(metadata)[2] <- "Sex"
metadata$Sex[metadata$Sex == 0] <- "Boy"
metadata$Sex[metadata$Sex == 1] <- "Girl"
metadata$Sex <- factor(metadata$Sex)

colnames(metadata)[7] <- "Category_T0"
metadata$Category_T0[metadata$Category_T0 == 1] <- "Thinness"
metadata$Category_T0[metadata$Category_T0 == 2] <- "Normal"
metadata$Category_T0[metadata$Category_T0 == 3] <- "Overweight"
metadata$Category_T0[metadata$Category_T0 == 4] <- "Obesity"
metadata$Category_T0 <- factor(metadata$Category_T0, levels = c("Thinness", "Normal", "Overweight", "Obesity"))

metadata$cat_peso_01[metadata$cat_peso_01 == 1] <- "Thinness"
metadata$cat_peso_01[metadata$cat_peso_01 == 2] <- "Normal"
metadata$cat_peso_01[metadata$cat_peso_01 == 3] <- "Overweight"
metadata$cat_peso_01[metadata$cat_peso_01 == 4] <- "Obesity"
metadata$cat_peso_01 <- factor(metadata$cat_peso_01, levels = c("Thinness", "Normal", "Overweight", "Obesity"))

na_metadata_t0 <- rownames(metadata[is.na(metadata$Category_T0),])
na_metadata_t1 <- rownames(metadata[is.na(metadata$imc_01),])
metadata_faltan <- metadata[! rownames(metadata) %in% colnames(s_otu_tab), ]

metadata$tanda <- 2
metadata$tanda <- riga_tandas45$N.Batch[match(as.character(metadata$id), as.character(riga_tandas45$codk2))]


colnames(metadata)[1] <- 'sampleID'
nreads <- s_otu_tab %>% colSums()

metadata$hospital <- metadata_with_origin$hospital[match(metadata$sampleID, metadata_with_origin$id)]


metadata_gain_labels2 <- metadata_gain_labels %>% select(id_participante, age_months_t0, age_months_t1, status_c2, mean_t1, lower_bound_t1, upper_bound_t1, bmi_t1_2)
table(metadata_gain_labels$id_participante %in% metadata$sampleID)
metadata <- merge(metadata, metadata_gain_labels2, by.x="sampleID", by.y="id_participante", all.x = TRUE) %>% 
  dplyr::mutate(is_normal = !is.na(status_c2)) %>% 
  dplyr::mutate(status_c2 = recode(status_c2, Normal = "Normal", `Ganancia excesiva`="Excessive gain", `Ganancia insuficiente`="Insufficient gain")) %>% 
  mutate(sampleID2 = sampleID,
         sampleID = paste("C", sampleID, sep="")) 
rownames(metadata) <- metadata$sampleID
  


write_tsv(metadata, file = paste0(outdir, "metadata_full.tsv"))
metadata_full <- metadata
s_otu_tab_unfilt <- s_otu_tab

common_samples <- names(s_otu_tab)[names(s_otu_tab) %in% rownames(metadata)]
metadata <- metadata_full[common_samples, ]
s_otu_tab <- s_otu_tab_unfilt[, common_samples]

nreads <- s_otu_tab %>% colSums()
metadata$nreads <- nreads[rownames(metadata)]
greads <- ggplot(metadata, aes(x=hospital, y = log10(nreads), fill=hospital))+geom_violin(alpha=0.6)+geom_boxplot(width=0.2, fill="lightgray")+ theme_bw()
ggsave(filename = paste0(outdir, "/reads_per_hospital.pdf"), greads, width = 7, height = 4)



all(names(s_otu_tab) == rownames(metadata))
write_tsv(s_otu_tab %>% rownames_to_column("taxon") %>% select(taxon, everything()), file = paste0(outdir, "otu_tab_names_presentInMetadata.tsv"))
write_tsv(metadata, file = paste0(outdir, "metadata_presentInOtus.tsv"))
s_meta <- metadata
