
outdir <- paste0(opt$out, "IntegrateWithAndWithoutCorrection/")
if(!dir.exists(outdir)) dir.create(outdir)
firstContrast <- daa_all$remove_tanda2
mainContrastName <- "Depression_vs_Control"
edad_correcting <- read_tsv(paste0(opt$out, "/DESeq2_ControlVars/DeSEQ2/remove_tanda2_Edad_log/remove_tanda2_Edad_log_Edad_log_DAAshrinkNormal.tsv"))
imc_correcting <- read_tsv(paste0(opt$out, "/DESeq2_ControlVars/DeSEQ2/remove_tanda2_BMI_log/remove_tanda2_BMI_log_BMI_log_DAAshrinkNormal.tsv"))
edad_corr2 <- read_tsv(paste0(opt$out,"/DESeq2_ControlVarsMany/DeSEQ2/remove_tanda2_EdaBmi/remove_tanda2_EdaBmi_Edad_log_DAAshrinkNormal.tsv"))
imc_corr2 <-  read_tsv(paste0(opt$out, "/DESeq2_ControlVarsMany/DeSEQ2/remove_tanda2_EdaBmi/remove_tanda2_EdaBmi_BMI_log_DAAshrinkNormal.tsv"))


#Sólo he guardado resdf
contrastlist2 <- list(
  daa_all_corrected$remove_tanda2$Edad_log,
  daa_all_corrected$remove_tanda2$BMI_log,
  daa_all_corrected_groups$remove_tanda2$EdaBmi,
  
  daa_all_corrected_only$remove_tanda2$Edad_log,
  edad_correcting,
  edad_corr2,
  
  daa_all_corrected_only$remove_tanda2$BMI_log,
  imc_correcting,
  imc_corr2
  
) %>% lapply(\(x)return(list(resdf=x)))
name2remove2 <- "xxx"
names(contrastlist2) <- c(
  "Condition_corrAge", 
  "Condition_corrBMI", 
  "Condition_corr2", 
  
  "Age_alone", 
  "Age_corrCond", 
  "Age_corr2",
  
  "BMI_alone", 
  "BMI_corrCond",
  "BMI_corr2")

contrastNamesOrdered2 <- c("Depression vs Control",  gsub("_", " ", names(contrastlist2)))
compareLFCContrats(contrastlist2, firstContrast, 
                   contrastNamesOrdered2, mainContrastName, 
                   plim_select= 0.05, plim_plot=0.1,
                   name2remove = name2remove2,
                   resdfname="resdf", 
                   outdir = outdir, 
                   name="LFC_Comparison_AgeAndBMI_allCombos_p05", w=12, h=12, scale_mode = "free")
compareLFCContrats(contrastlist2, firstContrast, 
                   contrastNamesOrdered2, mainContrastName, 
                   plim_select= 0.001, plim_plot=0.1,
                   name2remove = name2remove2,
                   resdfname="resdf", 
                   outdir = outdir, name="LFC_Comparison_AgeAndBMI_allCombos_pem3", w=12, h=12, 
                   scale_mode = "free")
dea2contrasts <- list(firstContrast = firstContrast, contrastlist2=contrastlist2)
save(dea2contrasts, file = paste0(outdir, "/LFC_Comparison_AgeAndBMI_allCombos.RData"))

### plots only with IMC

contrastlist2 <- list(
  #daa_all_corrected$remove_tanda2$Edad_log,
  daa_all_corrected$remove_tanda2$BMI_log,
  #daa_all_corrected_groups$remove_tanda2$EdaBmi,
  
  #daa_all_corrected_only$remove_tanda2$Edad_log,
  #edad_correcting,
  #edad_corr2,
  
  daa_all_corrected_only$remove_tanda2$BMI_log,
  imc_correcting
  #imc_corr2
  
) %>% lapply(\(x)return(list(resdf=x)))
name2remove2 <- "xxx"
names(contrastlist2) <- c(
  #"Condition_corrAge", 
  "Condition_corrIMC", 
  #"Condition_corr2", 
  
  #"Age_alone", 
  #"Age_corrCond", 
  #"Age_corr2",
  
  "BMI_alone", 
  "BMI_corrCond")
#"BMI_corr2")

contrastNamesOrdered2 <- c("Depression vs Control",  gsub("_", " ", names(contrastlist2)))
compareLFCContrats(contrastlist2, firstContrast, 
                   contrastNamesOrdered2, mainContrastName, 
                   plim_select= 0.05, plim_plot=0.1,
                   name2remove = name2remove2,
                   resdfname="resdf", 
                   outdir = outdir, name="LFC_Comparison_BMIcorrAndAlone_allCombos_p05", w=12, h=8, scale_mode = "free")
compareLFCContrats(contrastlist2, firstContrast, 
                   contrastNamesOrdered2, mainContrastName, 
                   plim_select= 0.001, plim_plot=0.1,
                   name2remove = name2remove2,
                   resdfname="resdf", 
                   outdir = outdir, name="LFC_Comparison_BMIcorrAndAlone_allCombos_pem3", w=12, h=8, 
                   scale_mode = "free")
dea2contrasts <- list(firstContrast = firstContrast, contrastlist2=contrastlist2)
save(dea2contrasts, file = paste0(outdir, "/LFC_Comparison_BMIcorrAndAlone_allCombos.RData"))

### Only Edad
contrastlist2 <- list(
  daa_all_corrected$remove_tanda2$Edad_log,
  #daa_all_corrected$remove_tanda2$BMI_log,
  #daa_all_corrected_groups$remove_tanda2$EdaBmi,
  
  daa_all_corrected_only$remove_tanda2$Edad_log,
  edad_correcting
  #edad_corr2,
  
  #daa_all_corrected_only$remove_tanda2$BMI_log,
  #imc_correcting
  #imc_corr2
  
) %>% lapply(\(x)return(list(resdf=x)))
name2remove2 <- "xxx"
names(contrastlist2) <- c(
  "Condition_corrAge", 
  #"Condition_corrIMC", 
  #"Condition_corr2", 
  
  "Age_alone", 
  "Age_corrCond") 
#"Age_corr2",

#  "BMI_alone", 
#  "BMI_corrCond")
#"BMI_corr2")

contrastNamesOrdered2 <- c("Depression vs Control",  gsub("_", " ", names(contrastlist2)))
compareLFCContrats(contrastlist2, firstContrast, 
                   contrastNamesOrdered2, mainContrastName, 
                   plim_select= 0.05, plim_plot=0.1,
                   name2remove = name2remove2,
                   resdfname="resdf", 
                   outdir = outdir, name="LFC_Comparison_AgecorrAndAlone_allCombos_p05", w=12, h=8, scale_mode = "free")
compareLFCContrats(contrastlist2, firstContrast, 
                   contrastNamesOrdered2, mainContrastName, 
                   plim_select= 0.001, plim_plot=0.1,
                   name2remove = name2remove2,
                   resdfname="resdf", 
                   outdir = outdir, name="LFC_Comparison_AgecorrAndAlone_allCombos_pem3", w=12, h=8, 
                   scale_mode = "free")
dea2contrasts <- list(firstContrast = firstContrast, contrastlist2=contrastlist2)
save(dea2contrasts, file = paste0(outdir, "/LFC_Comparison_AgecorrAndAlone_allCombos.RData"))



#### OBESITY/OVERWEIGHT AS CATEGORICAL

firstContrast <- daa_all$remove_tanda2
mainContrastName <- "Depression_vs_Control"

ob_correcting <- read_tsv(paste0(opt$out, "/DESeq2_ControlVars/DeSEQ2/remove_tanda2_ob_o_sobrepeso/remove_tanda2_ob_o_sobrepeso_ob_o_sobrepeso_normal.weight_vs_overweight_DAAshrinkNormal.tsv"))
edad_corr2_ob <- read_tsv(paste0(opt$out,"/DESeq2_ControlVarsMany/DeSEQ2/remove_tanda2_EdaObO/remove_tanda2_EdaObO_Edad_log_DAAshrinkNormal.tsv"))
ob_corr2 <-  read_tsv(paste0(opt$out, "/DESeq2_ControlVarsMany/DeSEQ2/remove_tanda2_EdaObO/remove_tanda2_EdaObO_ob_o_sobrepeso_normal.weight_vs_overweight_DAAshrinkNormal.tsv"))


contrastlist2_ob <- list(
  daa_all_corrected$remove_tanda2$Edad_log,
  daa_all_corrected$remove_tanda2$ob_o_sobrepeso,
  daa_all_corrected_groups$remove_tanda2$EdaObO,
  daa_all_corrected$remove_tanda2$BMI_log,
  daa_all_corrected_groups$remove_tanda2$EdaBmi,
  
  daa_all_corrected_only$remove_tanda2$Edad_log,
  edad_correcting,
  edad_corr2,
  edad_corr2_ob,
  
  daa_all_corrected_only$remove_tanda2$ob_o_sobrepeso,
  imc_correcting,
  imc_corr2,
  
  daa_all_corrected_only$remove_tanda2$ob_o_sobrepeso,
  ob_correcting,
  ob_corr2
  
) %>% lapply(\(x)return(list(resdf=x)))
name2remove2 <- "xxx"
names(contrastlist2) <- c(
  "Condition_corrAge", 
  "Condition_corrOverweight", 
  "Condition_corr2ob", 
  "Condition_corrBMI", 
  "Condition_corr2", 
  
  "Age_alone", 
  "Age_corrCond", 
  "Age_corr2",
  "Age_corr2ob",
  
  "BMI_alone", 
  "BMI_corrCond",
  "BMI_corr2",
  
  "Overweight_alone", 
  "Overweight_corrCond",
  "Overweight_corr2")

contrastNamesOrdered2 <- c("Depression vs Control",  gsub("_", " ", names(contrastlist2)))