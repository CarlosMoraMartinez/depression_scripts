
opt$reserva_0 <- opt$out
opt$out <- paste0(opt$out, "DESeq2_MainVariableScales/")
if(!dir.exists(opt$out)) dir.create(opt$out)

daa_all_scales <- list()
vars2deseq <- c(escalas_qual, escalas_quant)
opt$mincount <- 1
for(phname in names(all_phyloseq)){
  phobj <- all_phyloseq[[phname]]
  daa_all_scales[[phname]] <- list()
  for(var in vars2deseq){
    cat("Doing DESeq2 Analysys for: ", phname, ", with main variable ", var, "\n")
    samples <- sample_data(phobj)$sampleID[!is.na(sample_data(phobj)[, var]) ]
    phobj_filt <- phyloseq::prune_samples(samples, phobj)
    daa_all_scales[[phname]][[var]] <-deseq_full_pipeline(phobj_filt, paste0(phname, "_",var), var, opt)
  }}

save(daa_all_scales, file = paste0(opt$out, "DeSEQ2/DESEQ2_MainVariableScales.RData"))
opt$out <- opt$reserva_0


## With scales, correcting for covariates
opt$reserva_0 <- opt$out
opt$out <- paste0(opt$out, "DESeq2_MainVariableScales_Corrected/")
if(!dir.exists(opt$out)) dir.create(opt$out)
interest_vars <- escalas_quant
daa_all_corrected_scales <- list()
vars2test <- c("Edad_log", "BMI_log", "Sexo")
for(phname in phseq_to_correct){
  cat("Doing DESeq2 Analysys with correction for: ", phname, "\n")
  phobj <- all_phyloseq[[phname]]
  phobj <- updatePsWithLogs(phobj, c("Edad", "IMC"))
  daa_all_corrected_scales[[phname]] <- list()
  for(interestvar_tmp in interest_vars){
    daa_all_corrected_scales[[phname]][[interestvar_tmp]] <- list()
    samples <- sample_data(phobj)$sampleID[! is.na(sample_data(phobj)[, interestvar_tmp])]
    phobj_prefilt <- phyloseq::prune_samples(samples, phobj)
    for(var in vars2test){
      cat("Doing DESeq2 Analysys with correction for: ", phname, '-',interestvar_tmp, ',correcting for', var, "\n")
      samples <- sample_data(phobj_prefilt)$sampleID[! is.na(sample_data(phobj_prefilt)[, var])]
      phobj_filt <- phyloseq::prune_samples(samples, phobj_prefilt)
      cases <- sample_data(phobj_filt)[, var] %>% unlist %>%  table
      if(length(which(cases > 0)) < 2 ){cat("Only one level, skipping this variable");next}
      
      vars2deseq <- c(interestvar_tmp, var)
      name <- paste0(phname, '_', var)
      res_tmp <-deseq_full_pipeline(phobj_filt, name, vars2deseq, opt)
      daa_all_corrected_scales[[phname]][[interestvar_tmp]][[var]] <- res_tmp
    }}}

save(daa_all_corrected_scales, file = paste0(opt$out, "DeSEQ2/DESEQ2_MainVariableScales.RData"))
opt$out <- opt$reserva_0

# Integrate
## Corr edad cualitativo vs corr edad cuantitativo
outdir <- paste0(opt$out, "DESeq2_MainVariableScales_Corrected/IntegrateContrasts/")
if(!dir.exists(outdir)) dir.create(outdir)
contrastlist2 <- list(daa_all_corrected_scales$remove_tanda2$Escala_depresión_Beck$Edad_log$all_contrasts[[1]],
                      daa_all_corrected_scales$remove_tanda2$Montgomery.Asberg$Edad_log$all_contrasts[[1]],
                      daa_all_corrected_scales$remove_tanda2$Escala_Hamilton$Edad_log$all_contrasts[[1]],
                      daa_all_corrected_scales$remove_tanda2$DMSV_puntuacion_total$Edad_log$all_contrasts[[1]]
)
name2remove2 <- "xxx"
names(contrastlist2) <- c("Beck", "Montgomery", "Hamiliton", "DMSV")
contrastNamesOrdered2 <- c("Depression vs Control", names(contrastlist2))
compareLFCContrats(contrastlist2, firstContrast_edad, 
                   contrastNamesOrdered2, mainContrastName, 
                   plim_select= 0.05, plim_plot=0.1,
                   name2remove = name2remove2,
                   resdfname="resdf", 
                   outdir = outdir, name="LFC_Comparison_CorrectEdad_AllNumCorrected_p05", w=12, h=8, scale_mode = "free")
compareLFCContrats(contrastlist2, firstContrast_edad, 
                   contrastNamesOrdered2, mainContrastName, 
                   plim_select= 0.001, plim_plot=0.1,
                   name2remove = name2remove2,
                   resdfname="resdf", 
                   outdir = outdir, name="LFC_Comparison_CorrectEdad_AllNumCorrected_pem3", w=12, h=8, 
                   scale_mode = "free")


## Corr IMC cualitativo vs corr IMC cuantitativo

contrastlist2 <- list(daa_all_corrected_scales$remove_tanda2$Escala_depresión_Beck$BMI_log$all_contrasts[[1]],
                      daa_all_corrected_scales$remove_tanda2$Montgomery.Asberg$BMI_log$all_contrasts[[1]],
                      daa_all_corrected_scales$remove_tanda2$Escala_Hamilton$BMI_log$all_contrasts[[1]],
                      daa_all_corrected_scales$remove_tanda2$DMSV_puntuacion_total$BMI_log$all_contrasts[[1]]
)
name2remove2 <- "xxx"
names(contrastlist2) <- c("Beck", "Montgomery", "Hamiliton", "DMSV")
contrastNamesOrdered2 <- c("Depression vs Control", names(contrastlist2))
compareLFCContrats(contrastlist2, firstContrast_IMC, 
                   contrastNamesOrdered2, mainContrastName, 
                   plim_select= 0.05, plim_plot=0.1,
                   name2remove = name2remove2,
                   resdfname="resdf", 
                   outdir = outdir, name="LFC_Comparison_CorrectIMC_AllNumCorrected_p05", w=12, h=8, scale_mode = "free")
compareLFCContrats(contrastlist2, firstContrast_IMC, 
                   contrastNamesOrdered2, mainContrastName, 
                   plim_select= 0.001, plim_plot=0.1,
                   name2remove = name2remove2,
                   resdfname="resdf", 
                   outdir = outdir, name="LFC_Comparison_CorrectIMC_AllNumCorrected_pem3", w=12, h=8, 
                   scale_mode = "free")