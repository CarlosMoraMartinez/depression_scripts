## DAA only in Depressive subjects
vars2test <- c("Edad_log", "Sexo", "BMI_log", "Estado.civil2", 
               "Educacion", "reads_log10", "Fumador", 
               "TG", "TG_mayor_200",  "Colesterol", "Colesterol_mayor_200", 
               "PSS_estres", "Diabetes", "bristol_scale", "bristol_scale_cualitativo",
               "defecaciones_semana", "Euroqol", "IPAQ_act_fisica", "Mediterranean_diet_adherence",
               "tratamiento_ansioliticos", "tratamiento_anticonvulsivos",
               "tratamiento_ISRNs", "tratamiento_antidiabeticos", "tratamiento_coagul_betabloq_etc",
               "DMSV_puntuacion_total",
               "Escala_depresión_Beck", "Beck_cualitativo", "Escala_Hamilton", "Escala_Hamilton_cualitativo",
               "Montgomery.Asberg", "Montgomery.Asberg_qual", "DII", "DII.cualitativo", "Indice_Alimentación_saludable"
)
vars2test <- c("Euroqol", "IPAQ_act_fisica", "Mediterranean_diet_adherence", "Mediterranean_diet_adherence2")
escalas_qual <- c( "Beck_cualitativo", "Escala_Hamilton_cualitativo","Montgomery.Asberg_qual")
escalas_quant <- c( "Escala_depresión_Beck", "Escala_Hamilton", "Montgomery.Asberg", "DMSV_puntuacion_total")

interestvar <- "Condition"
vlevs <- levels(metadata[, interestvar] %>% unlist)

opt$reserva_0 <- opt$out
opt$out <- paste0(opt$out, "DESeq2_ByGroup/")
if(!dir.exists(opt$out)) dir.create(opt$out)

load(allphyloseqlist_fname)
# for(phname in names(all_phyloseq)){
#   for(v in escalas_qual){
#      sample_data(all_phyloseq[[phname]])[, v] <- factor(unlist(sample_data(all_phyloseq[[phname]])[, v]),
#                                                         levels=c("No Depresion", "Leve", "Moderada",
#                                                                  "Severa" ))
#     sample_data(all_phyloseq[[phname]])[, v] <- dplyr::recode(unlist(sample_data(all_phyloseq[[phname]])[, v]),
#                                                      "No Depresion"="Healthy", 
#                                                      "Leve"="Mild", 
#                                                      "Moderada"="Moderate",
#                                                     "Severa"="Severe" )
#   }
# }

daa_all_bygroup <- list()
for(phname in phseq_to_correct){
  cat("Doing DESeq2 Analysys with correction for: ", phname, "\n")
  phobj <- all_phyloseq[[phname]]
  phobj <- updatePsWithLogs(phobj, c("Edad", "IMC"))
  daa_all_bygroup[[phname]] <- list()
  for(varlev in vlevs){
    daa_all_bygroup[[phname]][[varlev]] <- list()
    samples <- sample_data(phobj)$sampleID[unlist(sample_data(phobj)[, interestvar]) == varlev ]
    phobj_prefilt <- phyloseq::prune_samples(samples, phobj)
    for(var in vars2test){
      if(var %in% c(escalas_qual, escalas_quant) & varlev == "Control") next
      cat("Doing DESeq2 Analysys within group for: ", phname, '-', var,' inside ',interestvar, ':',varlev, "\n")
      samples <- sample_data(phobj_prefilt)$sampleID[!is.na(sample_data(phobj_prefilt)[, var]) ]
      phobj_filt <- phyloseq::prune_samples(samples, phobj_prefilt)
      if(var %in% escalas_qual){
        samples <- sample_data(phobj_filt)$sampleID[!unlist(sample_data(phobj_filt)[, var]) %in% c("No Depression", "Healthy", "No depresion") ]
        phobj_filt <- phyloseq::prune_samples(samples, phobj_filt)
      }
      cases <- sample_data(phobj_filt)[, var] %>% unlist %>%  table
      if(length(which(cases > 1)) < 2 & !is.numeric(unlist(sample_data(phobj_filt)[, var] ))){cat("Only one level, skipping this variable");next}
      
      name <- paste0(phname, '_only',varlev,'_', var)
      res_tmp <-deseq_full_pipeline(phobj_filt, name, var, opt)
      daa_all_bygroup[[phname]][[varlev]][[var]] <- res_tmp$resdf
    } # variables
  } # Levels in interestvar (case/control)
}# phyloseq objects

opt$out <- opt$reserva_0
save(daa_all_bygroup, file=paste0(opt$out, "DESeq2_ControlVars/DESEQ2_byGroup_all.RData"))

load(paste0(opt$out, "DESeq2_ControlVars/DESEQ2_byGroup_all.RData"))
