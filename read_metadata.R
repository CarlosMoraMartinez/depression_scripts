########################################
# Read MetaData
########################################
escalas_qual <- c( "Beck_cualitativo", "Escala_Hamilton_cualitativo","Montgomery.Asberg_qual")
escalas_quant <- c( "Escala_depresión_Beck", "Escala_Hamilton", "Montgomery.Asberg", "DMSV_puntuacion_total")

metadata <- read_excel(opt$metadata,
                       sheet = 'Hoja1', na="NA") %>% 
  mutate(Tanda=as.factor(Tanda)) %>% 
  mutate( Condition = ifelse(CP =="DEPRESIVO", "Depression", "Control"), 
          tandacaso=paste0(Tanda,'-', Condition),
          procedenciacaso=paste(PROCEDENCIA, '-', Condition),
          obesidadcaso=paste(obesidad, '-', Condition),
          sexocaso=paste(Sexo, '-', Condition),
          BMI = IMC,
          Smoking_status = ifelse(is.na(Fumador), NA, ifelse(Fumador=="no", "No", "Yes")),
          Alcohol_abuse = ifelse(is.na(`Abuso Alcohol`), NA, ifelse(`Abuso Alcohol`=="no", "No", "Yes")),
          Beck_cualitativo = ifelse(Beck_cualitativo=="Ligero", "Leve", Beck_cualitativo),
          ob_o_sobrepeso = recode_factor(ob_o_sobrepeso, '1'="overweight", '2'="normal weight")) %>% 
  mutate_at(escalas_qual, \(x) dplyr::recode_factor(x, "No depresion"="Healthy", 
                                                    "Leve"="Mild", 
                                                    "Moderada"="Moderate",
                                                    "Severa"="Severe" ))

colnames(metadata)[3] <- 'sampleID'
nreads <- s_otu_tab %>% colSums()

tratamientos <- names(metadata)[grep("tratamiento", names(metadata))]
metadata <- metadata %>% 
  dplyr::mutate_at(tratamientos, as.factor) %>% 
  dplyr::mutate(reads = nreads[as.character(sampleID)], 
                reads_log10 = log10(nreads[as.character(sampleID)]))

names(metadata)[84] <- "Colesterol_ajuste"
metadata$Educacion[metadata$Educacion == "Estudios superiores:universidad y postgrado"] <- "universitarios" 
metadata$Educacion[metadata$Educacion == "fp"] <- "Bachiller y fp" 
metadata$Estado.civil2 <- ifelse(metadata$`Estado civil` == "pareja", "pareja", "soltero" )
metadata$DII.cualitativo[metadata$`DII cualitativo` == "Nuetra"] <- "Neutral"
metadata$Mediterranean_diet_adherence[metadata$Mediterranean_diet_adherence=="low"]<- "Low"
metadata$Mediterranean_diet_adherence2 <- ifelse(
  is.na(metadata$Mediterranean_diet_adherence), NA, 
  ifelse(metadata$Mediterranean_diet_adherence=="Low", "Low", "Good")
)
metadata$IPAQ_act_fisica <- with(metadata, recode(IPAQ_act_fisica, "Bajo"="Low", "Moderado"="Mid", "Alto"="High"))
metadata <- metadata %>% 
  mutate_if(is.character, str_to_title) %>% 
  mutate_if(~ (all(as.integer(.) == .) & length(unique(.)) < 5) | is.character(.),  as.factor)

omit_from_factors <- c("Condition", "CODIGO", "Presion.Sanguinea", "TRATAMIENTO" , "COMORBILIDADES", "tandacaso", "procedenciacaso", "obesidadcaso", "sexo", "sampleID")
#metad2 %>% select_if(is.numeric) %>% get_summary_stats

factors_all <- getValidFactors(metadata, max_nas = 1, min_pct = 0, max_classes = 5, exclude=omit_from_factors) %>% 
  filter(is_valid) %>% 
  pull(var)
numvars_all <- getValidNumeric(metadata)%>% filter(is_valid) %>% pull(var)
numvars_all <- numvars_all[numvars_all != "sampleID"]
## metadata table
s_meta <- data.frame(sampleID = names(s_otu_tab))
s_meta <- s_meta %>%
  dplyr::mutate(sampleNames_row = sampleID) %>%
  tibble::column_to_rownames(var = "sampleNames_row") 
# filter(! sampleID %in% c(82, 83)) #FALTAN ESTAS MUESTRAS EN EL EXCEL!
s_meta <- merge(s_meta, metadata, all.x = T)
row.names(s_meta) <- s_meta$sampleID