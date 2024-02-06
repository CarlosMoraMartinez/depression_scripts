
outdir <- paste0(opt$out, "IntegrateAllContrasts/")
if(!dir.exists(outdir)) dir.create(outdir)

contrastlist <- daa_all_scales$remove_tanda2$Beck_cualitativo$all_contrasts
firstContrast <- daa_all$remove_tanda2
mainContrastName <- "Depression_vs_Control"
contrastNamesOrdered <- c("Depression vs Control", "Mild vs Healthy","Moderate vs Healthy", 
                          "Severe vs Healthy", "Moderate vs Mild","Severe vs Mild", "Severe vs Moderate")
plim_select= 0.000001
plim_plot=0.05
name2remove="Beck_cualitativo_"

## Against Beck levels
compareLFCContrats(contrastlist, firstContrast, 
                   contrastNamesOrdered, mainContrastName, 
                   plim_select= 0.05, plim_plot=0.1,
                   name2remove = name2remove,
                   resdfname="resdf", 
                   outdir = outdir, name="LFC_Comparison_BeckQuant_p05", w=12, h=8)
compareLFCContrats(contrastlist, firstContrast, 
                   contrastNamesOrdered, mainContrastName, 
                   plim_select= 0.000001, plim_plot=0.1,
                   name2remove = name2remove,
                   resdfname="resdf", 
                   outdir = outdir, name="LFC_Comparison_BeckQuant_pem6", w=12, h=8)
## Against Montgomery levels
contrastlist2 <- daa_all_scales$remove_tanda2$Montgomery.Asberg_qual$all_contrasts
name2remove2 <- "Montgomery.Asberg_qual_"
compareLFCContrats(contrastlist2, firstContrast, 
                   contrastNamesOrdered, mainContrastName, 
                   plim_select= 0.05, plim_plot=0.1,
                   name2remove = name2remove2,
                   resdfname="resdf", 
                   outdir = outdir, name="LFC_Comparison_MontgomQuant_p05", w=12, h=8)
compareLFCContrats(contrastlist2, firstContrast, 
                   contrastNamesOrdered, mainContrastName, 
                   plim_select= 0.000001, plim_plot=0.1,
                   name2remove = name2remove2,
                   resdfname="resdf", 
                   outdir = outdir, name="LFC_Comparison_MontgomQuant_pem6", w=12, h=8)
## Against Hamilton levels
contrastlist2 <- daa_all_scales$remove_tanda2$Escala_Hamilton_cualitativo$all_contrasts
name2remove2 <- "Escala_Hamilton_cualitativo_"
compareLFCContrats(contrastlist2, firstContrast, 
                   contrastNamesOrdered, mainContrastName, 
                   plim_select= 0.05, plim_plot=0.1,
                   name2remove = name2remove2,
                   resdfname="resdf", 
                   outdir = outdir, name="LFC_Comparison_HamiltonQuant_p05", w=12, h=8)
compareLFCContrats(contrastlist2, firstContrast, 
                   contrastNamesOrdered, mainContrastName, 
                   plim_select= 0.000001, plim_plot=0.1,
                   name2remove = name2remove2,
                   resdfname="resdf", 
                   outdir = outdir, name="LFC_Comparison_HamiltonQuant_pem6", w=12, h=8)

## Against Beck numerical
contrastlist2 <- daa_all_scales$remove_tanda2$Escala_depresión_Beck$all_contrasts
name2remove2 <- "xxx"
names(contrastlist2) <- c("Beck scale")
contrastNamesOrdered2 <- c("Depression vs Control", "Beck scale")
compareLFCContrats(contrastlist2, firstContrast, 
                   contrastNamesOrdered2, mainContrastName, 
                   plim_select= 0.05, plim_plot=0.1,
                   name2remove = name2remove2,
                   resdfname="resdf", 
                   outdir = outdir, name="LFC_Comparison_BeckNum_p05", w=12, h=8, scale_mode = "free")
compareLFCContrats(contrastlist2, firstContrast, 
                   contrastNamesOrdered2, mainContrastName, 
                   plim_select= 0.001, plim_plot=0.1,
                   name2remove = name2remove2,
                   resdfname="resdf", 
                   outdir = outdir, name="LFC_Comparison_BeckNum_pem3", w=12, h=8, 
                   scale_mode = "free")

### Quantitative scales against CtrlVsDepr
contrastlist2 <- list(daa_all_scales$remove_tanda2$Escala_depresión_Beck$all_contrasts[[1]],
                      daa_all_scales$remove_tanda2$Montgomery.Asberg$all_contrasts[[1]],
                      daa_all_scales$remove_tanda2$Escala_Hamilton$all_contrasts[[1]],
                      daa_all_scales$remove_tanda2$DMSV_puntuacion_total$all_contrasts[[1]]
)
name2remove2 <- "xxx"
names(contrastlist2) <- c("Beck", "Montgomery", "Hamiliton", "DMSV")
contrastNamesOrdered2 <- c("Depression vs Control", names(contrastlist2))
compareLFCContrats(contrastlist2, firstContrast, 
                   contrastNamesOrdered2, mainContrastName, 
                   plim_select= 0.05, plim_plot=0.1,
                   name2remove = name2remove2,
                   resdfname="resdf", 
                   outdir = outdir, name="LFC_Comparison_AllNum_p05", w=12, h=8, scale_mode = "free")
compareLFCContrats(contrastlist2, firstContrast, 
                   contrastNamesOrdered2, mainContrastName, 
                   plim_select= 0.001, plim_plot=0.1,
                   name2remove = name2remove2,
                   resdfname="resdf", 
                   outdir = outdir, name="LFC_Comparison_AllNum_pem3", w=12, h=8, 
                   scale_mode = "free")
compareLFCContrats2(contrastlist2, firstContrast, 
                    contrastNamesOrdered2, mainContrastName, 
                    plim_select= 0.05, plim_plot=0.1,
                    name2remove = name2remove2,
                    resdfname="resdf", 
                    outdir = outdir, name="LFC_Comparison_AllNum_p05", w=8, h=8, scale_mode = "free")
compareLFCContrats2(contrastlist2, firstContrast, 
                    contrastNamesOrdered2, mainContrastName, 
                    plim_select= 0.001, plim_plot=0.1,
                    name2remove = name2remove2,
                    resdfname="resdf", 
                    outdir = outdir, name="LFC_Comparison_AllNum_pem3", w=8, h=8, 
                    scale_mode = "free")

cont001 <- compareLFCContratsNumeric(contrastlist2, firstContrast, 
                                     contrastNamesOrdered2, mainContrastName, 
                                     plim_select= 0.001, plim_plot=0.1,
                                     name2remove = name2remove2,
                                     resdfname="resdf", 
                                     outdir = outdir, name="LFC_Comparison_AllNum_pem3", w=8, h=8, 
                                     scale_mode = "free")
cont01 <- compareLFCContratsNumeric(contrastlist2, firstContrast, 
                                    contrastNamesOrdered2, mainContrastName, 
                                    plim_select= 0.01, plim_plot=0.1,
                                    name2remove = name2remove2,
                                    resdfname="resdf", 
                                    outdir = outdir, name="LFC_Comparison_AllNum_pem2", w=8, h=8, 
                                    scale_mode = "free")
cont05 <- compareLFCContratsNumeric(contrastlist2, firstContrast, 
                                    contrastNamesOrdered2, mainContrastName, 
                                    plim_select= 0.05, plim_plot=0.1,
                                    name2remove = name2remove2,
                                    resdfname="resdf", 
                                    outdir = outdir, name="LFC_Comparison_AllNum_05", w=8, h=8, 
                                    scale_mode = "free")



## Quant depression scales Inside depressive subjects
#Sólo he guardado resdf
contrastlist2 <- list(daa_all_bygroup$remove_tanda2$Depression$Escala_depresión_Beck,
                      daa_all_bygroup$remove_tanda2$Depression$Montgomery.Asberg,
                      daa_all_bygroup$remove_tanda2$Depression$Escala_Hamilton,
                      daa_all_bygroup$remove_tanda2$Depression$DMSV_puntuacion_total,
                      daa_all_bygroup$remove_tanda2$Depression$PSS_estres
) %>% lapply(\(x)return(list(resdf=x)))
name2remove2 <- "xxx"
names(contrastlist2) <- c("Beck", "Montgomery", "Hamiliton", "DMSV", "PSS")
contrastNamesOrdered2 <- c("Depression vs Control", names(contrastlist2))
compareLFCContrats(contrastlist2, firstContrast, 
                   contrastNamesOrdered2, mainContrastName, 
                   plim_select= 0.05, plim_plot=0.1,
                   name2remove = name2remove2,
                   resdfname="resdf", 
                   outdir = outdir, name="LFC_Comparison_AllNum_OnlyDepr_p05", w=12, h=8, scale_mode = "free")
compareLFCContrats(contrastlist2, firstContrast, 
                   contrastNamesOrdered2, mainContrastName, 
                   plim_select= 0.001, plim_plot=0.1,
                   name2remove = name2remove2,
                   resdfname="resdf", 
                   outdir = outdir, name="LFC_Comparison_AllNum_OnlyDepr_pem3", w=12, h=8, 
                   scale_mode = "free")
compareLFCContrats2(contrastlist2, firstContrast, 
                    contrastNamesOrdered2, mainContrastName, 
                    plim_select= 0.05, plim_plot=0.05,
                    name2remove = name2remove2,
                    resdfname="resdf", 
                    outdir = outdir, name="LFC_Comparison_AllNum_OnlyDepr_p05", w=12, h=8, scale_mode = "free")
compareLFCContrats2(contrastlist2, firstContrast, 
                    contrastNamesOrdered2, mainContrastName, 
                    plim_select= 0.001, plim_plot=0.05,
                    name2remove = name2remove2,
                    resdfname="resdf", 
                    outdir = outdir, name="LFC_Comparison_AllNum_OnlyDepr_pem3", w=12, h=8, 
                    scale_mode = "free")

cont001 <- compareLFCContratsNumeric(contrastlist2, firstContrast, 
                                     contrastNamesOrdered2, mainContrastName, 
                                     plim_select= 0.001, plim_plot=0.1,
                                     name2remove = name2remove2,
                                     resdfname="resdf", 
                                     outdir = outdir, name="LFC_Comparison_AllNum_OnlyDepr_pem3", w=8, h=8, 
                                     scale_mode = "free")
cont01 <- compareLFCContratsNumeric(contrastlist2, firstContrast, 
                                    contrastNamesOrdered2, mainContrastName, 
                                    plim_select= 0.01, plim_plot=0.1,
                                    name2remove = name2remove2,
                                    resdfname="resdf", 
                                    outdir = outdir, name="LFC_Comparison_AllNum_OnlyDepr_pem2", w=8, h=8, 
                                    scale_mode = "free")
cont05 <- compareLFCContratsNumeric(contrastlist2, firstContrast, 
                                    contrastNamesOrdered2, mainContrastName, 
                                    plim_select= 0.05, plim_plot=0.1,
                                    name2remove = name2remove2,
                                    resdfname="resdf", 
                                    outdir = outdir, name="LFC_Comparison_AllNum_OnlyDepr_05", w=8, h=8, 
                                    scale_mode = "free")


## Beck, DSMV and stress Inside depressive subjects
#Sólo he guardado resdf
contrastlist2 <- list(daa_all_bygroup$remove_tanda2$Depression$Escala_depresión_Beck,
                      daa_all_bygroup$remove_tanda2$Depression$DMSV_puntuacion_total,
                      daa_all_bygroup$remove_tanda2$Depression$PSS_estres
                      
) %>% lapply(\(x)return(list(resdf=x)))
name2remove2 <- "xxx"
names(contrastlist2) <- c("Beck", "DMSV", "Stress")
contrastNamesOrdered2 <- c("Depression vs Control", names(contrastlist2))
compareLFCContrats(contrastlist2, firstContrast, 
                   contrastNamesOrdered2, mainContrastName, 
                   plim_select= 0.05, plim_plot=0.1,
                   name2remove = name2remove2,
                   resdfname="resdf", 
                   outdir = outdir, name="LFC_Comparison_BeckDMSStress_OnlyDepr_p05", w=12, h=8, scale_mode = "free")
compareLFCContrats(contrastlist2, firstContrast, 
                   contrastNamesOrdered2, mainContrastName, 
                   plim_select= 0.01, plim_plot=0.1,
                   name2remove = name2remove2,
                   resdfname="resdf", 
                   outdir = outdir, name="LFC_Comparison_BeckDMSStress_OnlyDepr_pem3", w=12, h=8, 
                   scale_mode = "free")

compareLFCContrats2(contrastlist2, firstContrast, 
                    contrastNamesOrdered2, mainContrastName, 
                    plim_select= 0.05, plim_plot=0.1,
                    name2remove = name2remove2,
                    resdfname="resdf", 
                    outdir = outdir, name="LFC_Comparison_BeckDMSStress_OnlyDepr_05", w=8, h=8, scale_mode = "free")
compareLFCContrats2(contrastlist2, firstContrast, 
                    contrastNamesOrdered2, mainContrastName, 
                    plim_select= 0.001, plim_plot=0.1,
                    name2remove = name2remove2,
                    resdfname="resdf", 
                    outdir = outdir, name="LFC_Comparison_BeckDMSStress_OnlyDepr_pem3", w=8, h=8, 
                    scale_mode = "free")

cont001 <- compareLFCContratsNumeric(contrastlist2, firstContrast, 
                                     contrastNamesOrdered2, mainContrastName, 
                                     plim_select= 0.001, plim_plot=0.1,
                                     name2remove = name2remove2,
                                     resdfname="resdf", 
                                     outdir = outdir, name="LFC_Comparison_BeckDMSStress_OnlyDepr_pem3", w=8, h=8, 
                                     scale_mode = "free")
cont01 <- compareLFCContratsNumeric(contrastlist2, firstContrast, 
                                    contrastNamesOrdered2, mainContrastName, 
                                    plim_select= 0.01, plim_plot=0.1,
                                    name2remove = name2remove2,
                                    resdfname="resdf", 
                                    outdir = outdir, name="LFC_Comparison_BeckDMSStress_OnlyDepr_pem2", w=8, h=8, 
                                    scale_mode = "free")
cont05 <- compareLFCContratsNumeric(contrastlist2, firstContrast, 
                                    contrastNamesOrdered2, mainContrastName, 
                                    plim_select= 0.05, plim_plot=0.1,
                                    name2remove = name2remove2,
                                    resdfname="resdf", 
                                    outdir = outdir, name="LFC_Comparison_BeckDMSStress_OnlyDepr_05", w=8, h=8, 
                                    scale_mode = "free")
#Beck levels inside Depression
fnames <- list.files(paste0(opt$out, "/DESeq2_ByGroup/DeSEQ2/remove_tanda2_onlyDepression_Beck_cualitativo"), 
                     pattern = "_vs_", full.names = T ) %>% 
  subset(grepl("Normal.tsv", .))
contrastlist2 <- map(fnames, read_tsv) %>% lapply(\(x)return(list(resdf=x)))
newnames <- basename(fnames) %>% sapply(\(x)strsplit(x, "Beck_cualitativo_")[[1]][3]) %>% 
  gsub("_DAAshrinkNormal.tsv", "", .) %>% gsub("_", " ", .)
names(contrastlist2) <- newnames
name2remove2 <- "xxx"

contrastNamesOrdered2 <- c("Depression vs Control", names(contrastlist2))
compareLFCContrats(contrastlist2, firstContrast, 
                   contrastNamesOrdered2, mainContrastName, 
                   plim_select= 0.05, plim_plot=0.1,
                   name2remove = name2remove2,
                   resdfname="resdf", 
                   outdir = outdir, name="LFC_Comparison_BeckQual_OnlyDepr_p05", w=12, h=8, scale_mode = "fixed")
compareLFCContrats(contrastlist2, firstContrast, 
                   contrastNamesOrdered2, mainContrastName, 
                   plim_select= 0.001, plim_plot=0.1,
                   name2remove = name2remove2,
                   resdfname="resdf", 
                   outdir = outdir, name="LFC_Comparison_BeckQual_OnlyDepr_pem3", w=12, h=8, 
                   scale_mode = "fixed")


#Exercise levels inside Depression
fnames <- list.files(paste0(opt$out, "/DESeq2_ByGroup/DeSEQ2/remove_tanda2_onlyDepression_IPAQ_act_fisica/"), 
                     pattern = "_vs_", full.names = T ) %>% 
  subset(grepl("Normal.tsv", .))
contrastlist2 <- map(fnames, read_tsv) %>% lapply(\(x)return(list(resdf=x)))
newnames <- basename(fnames) %>% sapply(\(x)strsplit(x, "fisica_")[[1]][3]) %>% 
  gsub("_DAAshrinkNormal.tsv", "", .) %>% gsub("_", " ", .)
names(contrastlist2) <- newnames
name2remove2 <- "xxx"

contrastNamesOrdered2 <- c("Depression vs Control", names(contrastlist2))
compareLFCContrats(contrastlist2, firstContrast, 
                   contrastNamesOrdered2, mainContrastName, 
                   plim_select= 0.05, plim_plot=0.1,
                   name2remove = name2remove2,
                   resdfname="resdf", 
                   outdir = outdir, name="LFC_Comparison_IPAQ_OnlyDepr_p05", w=12, h=8, scale_mode = "fixed")
compareLFCContrats(contrastlist2, firstContrast, 
                   contrastNamesOrdered2, mainContrastName, 
                   plim_select= 0.001, plim_plot=0.1,
                   name2remove = name2remove2,
                   resdfname="resdf", 
                   outdir = outdir, name="LFC_Comparison_IPAQ_OnlyDepr_pem3", w=12, h=8, 
                   scale_mode = "fixed")
compareLFCContrats2(contrastlist2, firstContrast, 
                    contrastNamesOrdered2, mainContrastName, 
                    plim_select= 0.05, plim_plot=0.1,
                    name2remove = name2remove2,
                    resdfname="resdf", 
                    outdir = outdir, name="LFC_Comparison_IPAQ_OnlyDepr_p05", w=10, h=12, scale_mode = "free")
compareLFCContrats2(contrastlist2, firstContrast, 
                    contrastNamesOrdered2, mainContrastName, 
                    plim_select= 0.001, plim_plot=0.1,
                    name2remove = name2remove2,
                    resdfname="resdf", 
                    outdir = outdir, name="LFC_Comparison_IPAQ_OnlyDepr_pem3", w=10, h=8, 
                    scale_mode = "free")
names(contrastlist2) <- gsub(" ", "_", names(contrastlist2))
contrastNamesOrdered2[2:length(contrastNamesOrdered2)] <- gsub(" ", "_",contrastNamesOrdered2[2:length(contrastNamesOrdered2)])
cont001 <- compareLFCContratsNumeric(contrastlist2, firstContrast, 
                                     contrastNamesOrdered2, mainContrastName, 
                                     plim_select= 0.001, plim_plot=0.1,
                                     name2remove = name2remove2,
                                     resdfname="resdf", 
                                     outdir = outdir, name="LFC_Comparison_IPAQ_OnlyDepr_pem3", w=12, h=8, 
                                     scale_mode = "free")
cont01 <- compareLFCContratsNumeric(contrastlist2, firstContrast, 
                                    contrastNamesOrdered2, mainContrastName, 
                                    plim_select= 0.01, plim_plot=0.1,
                                    name2remove = name2remove2,
                                    resdfname="resdf", 
                                    outdir = outdir, name="LFC_Comparison_IPAQ_OnlyDepr_pem2", w=8, h=8, 
                                    scale_mode = "free")
cont05 <- compareLFCContratsNumeric(contrastlist2, firstContrast, 
                                    contrastNamesOrdered2, mainContrastName, 
                                    plim_select= 0.05, plim_plot=0.1,
                                    name2remove = name2remove2,
                                    resdfname="resdf", 
                                    outdir = outdir, name="LFC_Comparison_IPAQ_OnlyDepr_p05", w=8, h=12, 
                                    scale_mode = "free")
cont05$correlations %>% select(1,2,3,5,9,10) %>% mutate_if(is.numeric, round,3) %>%  datatable()

#Exercise levels adjusted by Condition vs Condition (unadjusted)
fnames <- list.files(paste0(opt$out, "/DESeq2_ControlVars/DeSEQ2/remove_tanda2_IPAQ_act_fisica/"), 
                     pattern = "_vs_", full.names = T ) %>% 
  subset(grepl("Normal.tsv", .)) %>% subset(!grepl("Depression_vs_Control", .))
contrastlist2 <- map(fnames, read_tsv) %>% lapply(\(x)return(list(resdf=x)))
newnames <- basename(fnames) %>% sapply(\(x)strsplit(x, "fisica_")[[1]][3]) %>% 
  gsub("_DAAshrinkNormal.tsv", "", .) %>% gsub("_", " ", .)
names(contrastlist2) <- newnames
name2remove2 <- "xxx"

contrastNamesOrdered2 <- c("Depression vs Control", names(contrastlist2))
compareLFCContrats(contrastlist2, firstContrast, 
                   contrastNamesOrdered2, mainContrastName, 
                   plim_select= 0.05, plim_plot=0.1,
                   name2remove = name2remove2,
                   resdfname="resdf", 
                   outdir = outdir, name="LFC_Comparison_IPAQ_CorrectedByDepr_p05", w=12, h=8, scale_mode = "fixed")
compareLFCContrats(contrastlist2, firstContrast, 
                   contrastNamesOrdered2, mainContrastName, 
                   plim_select= 0.001, plim_plot=0.1,
                   name2remove = name2remove2,
                   resdfname="resdf", 
                   outdir = outdir, name="LFC_Comparison_IPAQ_CorrectedByDepr_pem3", w=12, h=8, 
                   scale_mode = "fixed")
compareLFCContrats2(contrastlist2, firstContrast, 
                    contrastNamesOrdered2, mainContrastName, 
                    plim_select= 0.05, plim_plot=0.1,
                    name2remove = name2remove2,
                    resdfname="resdf", 
                    outdir = outdir, name="LFC_Comparison2_IPAQ_CorrectedByDepr_p05", w=10, h=12, scale_mode = "free")
compareLFCContrats2(contrastlist2, firstContrast, 
                    contrastNamesOrdered2, mainContrastName, 
                    plim_select= 0.001, plim_plot=0.1,
                    name2remove = name2remove2,
                    resdfname="resdf", 
                    outdir = outdir, name="LFC_Comparison2_CorrectedByDepr_pem3", w=10, h=8, 
                    scale_mode = "free")

compareLFCContrats3(contrastlist2[1], firstContrast, 
                    contrastNamesOrdered2, mainContrastName, 
                    plim_select= 0.05, plim_plot=0.1,
                    name2remove = name2remove2,
                    resdfname="resdf", 
                    outdir = outdir, name="LFC_Comparison3_IPAQ_CorrectedByDepr_p05", w=10, h=4, scale_mode = "free")

names(contrastlist2) <- gsub(" ", "_", names(contrastlist2))
contrastNamesOrdered2[2:length(contrastNamesOrdered2)] <- gsub(" ", "_",contrastNamesOrdered2[2:length(contrastNamesOrdered2)])
cont001 <- compareLFCContratsNumeric(contrastlist2, firstContrast, 
                                     contrastNamesOrdered2, mainContrastName, 
                                     plim_select= 0.001, plim_plot=0.1,
                                     name2remove = name2remove2,
                                     resdfname="resdf", 
                                     outdir = outdir, name="LFC_Comparison_IPAQ_CorrectedByDepr_pem3", w=12, h=8, 
                                     scale_mode = "free")
cont01 <- compareLFCContratsNumeric(contrastlist2, firstContrast, 
                                    contrastNamesOrdered2, mainContrastName, 
                                    plim_select= 0.01, plim_plot=0.1,
                                    name2remove = name2remove2,
                                    resdfname="resdf", 
                                    outdir = outdir, name="LFC_Comparison_IPAQ_CorrectedByDepr_pem2", w=8, h=8, 
                                    scale_mode = "free")
cont05 <- compareLFCContratsNumeric(contrastlist2, firstContrast, 
                                    contrastNamesOrdered2, mainContrastName, 
                                    plim_select= 0.05, plim_plot=0.1,
                                    name2remove = name2remove2,
                                    resdfname="resdf", 
                                    outdir = outdir, name="LFC_Comparison_IPAQ_CorrectedByDepr_p05", w=8, h=12, 
                                    scale_mode = "free")
#cont05$correlations %>% select(1,2,3,5,9,10) %>% mutate_if(is.numeric, round,3) %>%  datatable()


## Beck inside Depr vs IPAQ inside depr
fnames <- list.files(paste0(opt$out, "/DESeq2_ByGroup/DeSEQ2/remove_tanda2_onlyDepression_IPAQ_act_fisica/"), 
                     pattern = "_vs_", full.names = T ) %>% 
  subset(grepl("Normal.tsv", .))
contrastlist2 <- map(fnames, read_tsv) %>% lapply(\(x)return(list(resdf=x)))
newnames <- basename(fnames) %>% sapply(\(x)strsplit(x, "fisica_")[[1]][3]) %>% 
  gsub("_DAAshrinkNormal.tsv", "", .) %>% gsub("_", " ", .)
names(contrastlist2) <- newnames
firstContrast2 <- list(resdf=daa_all_bygroup$remove_tanda2$Depression$Escala_depresión_Beck)
contrastNamesOrdered2 <- c("Beck", names(contrastlist2))
compareLFCContrats(contrastlist2, firstContrast2, 
                   contrastNamesOrdered2, mainContrastName, 
                   plim_select= 0.05, plim_plot=0.1,
                   name2remove = name2remove2,
                   resdfname="resdf", 
                   outdir = outdir, name="LFC_Comparison_BeckInside_IPAQ_OnlyDepr_p05", w=12, h=8, scale_mode = "free")
compareLFCContrats(contrastlist2, firstContrast, 
                   contrastNamesOrdered2, mainContrastName, 
                   plim_select= 0.001, plim_plot=0.1,
                   name2remove = name2remove2,
                   resdfname="resdf", 
                   outdir = outdir, name="LFC_Comparison_BeckInside_IPAQ_OnlyDepr_pem3", w=12, h=8, 
                   scale_mode = "free")


#Diet levels inside Depression
fnames <- list.files(paste0(opt$out, "/DESeq2_ByGroup/DeSEQ2/remove_tanda2_onlyControl_Indice_Alimentación_saludable//"), 
                     pattern = "_vs_", full.names = T ) %>% 
  subset(grepl("Normal.tsv", .))
contrastlist2 <- map(fnames, read_tsv) %>% lapply(\(x)return(list(resdf=x)))
newnames <- basename(fnames) %>% sapply(\(x)strsplit(x, "able_")[[1]][3]) %>% 
  gsub("_DAAshrinkNormal.tsv", "", .) %>% gsub("_", " ", .)
names(contrastlist2) <- newnames
name2remove2 <- "xxx"

contrastNamesOrdered2 <- c("Depression vs Control", names(contrastlist2))
compareLFCContrats(contrastlist2, firstContrast, 
                   contrastNamesOrdered2, mainContrastName, 
                   plim_select= 0.05, plim_plot=0.1,
                   name2remove = name2remove2,
                   resdfname="resdf", 
                   outdir = outdir, name="LFC_Comparison_Alim_OnlyDepr_p05", w=12, h=8, scale_mode = "fixed")
compareLFCContrats(contrastlist2, firstContrast, 
                   contrastNamesOrdered2, mainContrastName, 
                   plim_select= 0.001, plim_plot=0.1,
                   name2remove = name2remove2,
                   resdfname="resdf", 
                   outdir = outdir, name="LFC_Comparison_Alim_OnlyDepr_pem3", w=12, h=8, 
                   scale_mode = "fixed")

firstContrast2 <- list(resdf=daa_all_bygroup$remove_tanda2$Depression$Escala_depresión_Beck)
contrastNamesOrdered2 <- c("Beck", names(contrastlist2))
compareLFCContrats(contrastlist2, firstContrast2, 
                   contrastNamesOrdered2, mainContrastName, 
                   plim_select= 0.05, plim_plot=0.1,
                   name2remove = name2remove2,
                   resdfname="resdf", 
                   outdir = outdir, name="LFC_Comparison_BeckInside_Alim_OnlyDepr_p05", w=12, h=8, scale_mode = "free")
compareLFCContrats(contrastlist2, firstContrast, 
                   contrastNamesOrdered2, mainContrastName, 
                   plim_select= 0.001, plim_plot=0.1,
                   name2remove = name2remove2,
                   resdfname="resdf", 
                   outdir = outdir, name="LFC_Comparison_BeckInside_Alim_OnlyDepr_pem3", w=12, h=8, 
                   scale_mode = "free")


#Mediterranea
fnames <- list.files(paste0(opt$out, "/DESeq2_ByGroup/DeSEQ2/remove_tanda2_onlyControl_Mediterranean_diet_adherence//"), 
                     pattern = "_vs_", full.names = T ) %>% 
  subset(grepl("Normal.tsv", .))
contrastlist2 <- map(fnames, read_tsv) %>% lapply(\(x)return(list(resdf=x)))
newnames <- basename(fnames) %>% sapply(\(x)strsplit(x, "herence_")[[1]][3]) %>% 
  gsub("_DAAshrinkNormal.tsv", "", .) %>% gsub("_", " ", .)
names(contrastlist2) <- newnames
name2remove2 <- "xxx"

contrastNamesOrdered2 <- c("Depression vs Control", names(contrastlist2))
compareLFCContrats(contrastlist2, firstContrast, 
                   contrastNamesOrdered2, mainContrastName, 
                   plim_select= 0.05, plim_plot=0.1,
                   name2remove = name2remove2,
                   resdfname="resdf", 
                   outdir = outdir, name="LFC_Comparison_DMedit_OnlyDepr_p05", w=12, h=8, scale_mode = "fixed")
compareLFCContrats(contrastlist2, firstContrast, 
                   contrastNamesOrdered2, mainContrastName, 
                   plim_select= 0.001, plim_plot=0.1,
                   name2remove = name2remove2,
                   resdfname="resdf", 
                   outdir = outdir, name="LFC_Comparison_DMedit_OnlyDepr_pem3", w=12, h=8, 
                   scale_mode = "fixed")
compareLFCContrats2(contrastlist2, firstContrast, 
                    contrastNamesOrdered2, mainContrastName, 
                    plim_select= 0.05, plim_plot=0.1,
                    name2remove = name2remove2,
                    resdfname="resdf", 
                    outdir = outdir, name="LFC_Comparison_DMedit_OnlyDepr_p05", w=12, h=12, scale_mode = "fixed")
compareLFCContrats2(contrastlist2, firstContrast, 
                    contrastNamesOrdered2, mainContrastName, 
                    plim_select= 0.001, plim_plot=0.1,
                    name2remove = name2remove2,
                    resdfname="resdf", 
                    outdir = outdir, name="LFC_Comparison_DMedit_OnlyDepr_pem3", w=12, h=10, 
                    scale_mode = "fixed")

names(contrastlist2) <- gsub(" ", "_", names(contrastlist2))
contrastNamesOrdered2[2:length(contrastNamesOrdered2)] <- gsub(" ", "_",contrastNamesOrdered2[2:length(contrastNamesOrdered2)])
cont001 <- compareLFCContratsNumeric(contrastlist2, firstContrast, 
                                     contrastNamesOrdered2, mainContrastName, 
                                     plim_select= 0.001, plim_plot=0.1,
                                     name2remove = name2remove2,
                                     resdfname="resdf", 
                                     outdir = outdir, name="LFC_Comparison_DMedit_OnlyDepr_pem3", w=12, h=8, 
                                     scale_mode = "free")
cont01 <- compareLFCContratsNumeric(contrastlist2, firstContrast, 
                                    contrastNamesOrdered2, mainContrastName, 
                                    plim_select= 0.01, plim_plot=0.1,
                                    name2remove = name2remove2,
                                    resdfname="resdf", 
                                    outdir = outdir, name="LFC_Comparison_DMedit_OnlyDepr_pem2", w=8, h=8, 
                                    scale_mode = "free")
cont05 <- compareLFCContratsNumeric(contrastlist2, firstContrast, 
                                    contrastNamesOrdered2, mainContrastName, 
                                    plim_select= 0.05, plim_plot=0.1,
                                    name2remove = name2remove2,
                                    resdfname="resdf", 
                                    outdir = outdir, name="LFC_Comparison_DMedit_OnlyDepr_p05", w=8, h=12, 
                                    scale_mode = "free")

## Mediterranea ajustada por depresion
#Mediterranea
fnames <- list.files(paste0(opt$out, "/DESeq2_ControlVars/DeSEQ2/remove_tanda2_Mediterranean_diet_adherence"), 
                     pattern = "_vs_", full.names = T ) %>% 
  subset(grepl("Normal.tsv", .)) %>% subset(!grepl("Depression_vs_Control", .))
contrastlist2 <- map(fnames, read_tsv) %>% lapply(\(x)return(list(resdf=x)))
newnames <- basename(fnames) %>% sapply(\(x)strsplit(x, "herence_")[[1]][3]) %>% 
  gsub("_DAAshrinkNormal.tsv", "", .) %>% gsub("_", " ", .)
names(contrastlist2) <- newnames
name2remove2 <- "xxx"

contrastNamesOrdered2 <- c("Depression vs Control", names(contrastlist2))
compareLFCContrats(contrastlist2, firstContrast, 
                   contrastNamesOrdered2, mainContrastName, 
                   plim_select= 0.05, plim_plot=0.1,
                   name2remove = name2remove2,
                   resdfname="resdf", 
                   outdir = outdir, name="LFC_Comparison_DMedit_CorrectedByDepr_p05", w=12, h=8, scale_mode = "fixed")
compareLFCContrats(contrastlist2, firstContrast, 
                   contrastNamesOrdered2, mainContrastName, 
                   plim_select= 0.001, plim_plot=0.1,
                   name2remove = name2remove2,
                   resdfname="resdf", 
                   outdir = outdir, name="LFC_Comparison_DMedit_CorrectedByDepr_pem3", w=12, h=8, 
                   scale_mode = "fixed")
compareLFCContrats2(contrastlist2, firstContrast, 
                    contrastNamesOrdered2, mainContrastName, 
                    plim_select= 0.05, plim_plot=0.05,
                    name2remove = name2remove2,
                    resdfname="resdf", 
                    outdir = outdir, name="LFC_Comparison2_DMedit_CorrectedByDepr_p05", w=12, h=12, scale_mode = "fixed")

compareLFCContrats2(contrastlist2, firstContrast, 
                    contrastNamesOrdered2, mainContrastName, 
                    plim_select= 0.01, plim_plot=0.05,
                    name2remove = name2remove2,
                    resdfname="resdf", 
                    outdir = outdir, name="LFC_Comparison2_DMedit_CorrectedByDepr_p01", w=12, h=10, 
                    scale_mode = "fixed")

compareLFCContrats2(contrastlist2, firstContrast, 
                    contrastNamesOrdered2, mainContrastName, 
                    plim_select= 0.001, plim_plot=0.05,
                    name2remove = name2remove2,
                    resdfname="resdf", 
                    outdir = outdir, name="LFC_Comparison2_DMedit_CorrectedByDepr_pem3", w=12, h=10, 
                    scale_mode = "fixed")
contrastlist2_partial <- contrastlist2[2]
names(contrastlist2_partial) <- "Low vs Optimal"
contrastlist2_partial$`Low vs Optimal`$resdf$log2FoldChangeShrink <- -1*contrastlist2_partial$`Low vs Optimal`$resdf$log2FoldChangeShrink
contrastlist2_partial$`Low vs Optimal`$resdf$log2FoldChange <- -1*contrastlist2_partial$`Low vs Optimal`$resdf$log2FoldChange
contrastNamesOrdered2_partial <- c(contrastNamesOrdered2[1], "Low vs Optimal")
compareLFCContrats3(contrastlist2_partial, firstContrast, 
                    contrastNamesOrdered2_partial, mainContrastName, 
                    plim_select= 0.05, plim_plot=0.05,
                    name2remove = name2remove2,
                    resdfname="resdf", 
                    outdir = outdir, name="LFC_Comparison3_DMedit_CorrectedByDepr_p05", w=10, h=4, scale_mode = "fixed")

names(contrastlist2) <- gsub(" ", "_", names(contrastlist2))
contrastNamesOrdered2[2:length(contrastNamesOrdered2)] <- gsub(" ", "_",contrastNamesOrdered2[2:length(contrastNamesOrdered2)])
cont001 <- compareLFCContratsNumeric(contrastlist2, firstContrast, 
                                     contrastNamesOrdered2, mainContrastName, 
                                     plim_select= 0.001, plim_plot=0.1,
                                     name2remove = name2remove2,
                                     resdfname="resdf", 
                                     outdir = outdir, name="LFC_Comparison_CorrectedByDepr_OnlyDepr_pem3", w=12, h=8, 
                                     scale_mode = "free")
cont01 <- compareLFCContratsNumeric(contrastlist2, firstContrast, 
                                    contrastNamesOrdered2, mainContrastName, 
                                    plim_select= 0.01, plim_plot=0.1,
                                    name2remove = name2remove2,
                                    resdfname="resdf", 
                                    outdir = outdir, name="LFC_Comparison_CorrectedByDepr_OnlyDepr_pem2", w=8, h=8, 
                                    scale_mode = "free")
cont05 <- compareLFCContratsNumeric(contrastlist2, firstContrast, 
                                    contrastNamesOrdered2, mainContrastName, 
                                    plim_select= 0.05, plim_plot=0.1,
                                    name2remove = name2remove2,
                                    resdfname="resdf", 
                                    outdir = outdir, name="LFC_Comparison_CorrectedByDepr_OnlyDepr_p05", w=8, h=12, 
                                    scale_mode = "free")

## Mediterranea2 (2 clases) ajustada por depresion

fnames <- list.files(paste0(opt$out, "/DESeq2_ControlVars/DeSEQ2/remove_tanda2_Mediterranean_diet_adherence2"), 
                     pattern = "_vs_", full.names = T ) %>% 
  subset(grepl("Normal.tsv", .)) %>% subset(!grepl("Depression_vs_Control", .))
contrastlist2 <- map(fnames, read_tsv) %>% lapply(\(x)return(list(resdf=x)))
newnames <- basename(fnames) %>% sapply(\(x)strsplit(x, "herence2_")[[1]][3]) %>% 
  gsub("_DAAshrinkNormal.tsv", "", .) %>% gsub("_", " ", .)
names(contrastlist2) <- newnames
name2remove2 <- "xxx"

contrastNamesOrdered2 <- c("Depression vs Control", names(contrastlist2))
compareLFCContrats(contrastlist2, firstContrast, 
                   contrastNamesOrdered2, mainContrastName, 
                   plim_select= 0.05, plim_plot=0.1,
                   name2remove = name2remove2,
                   resdfname="resdf", 
                   outdir = outdir, name="LFC_Comparison_DMedit2_CorrectedByDepr_p05", w=12, h=8, scale_mode = "fixed")
compareLFCContrats(contrastlist2, firstContrast, 
                   contrastNamesOrdered2, mainContrastName, 
                   plim_select= 0.001, plim_plot=0.1,
                   name2remove = name2remove2,
                   resdfname="resdf", 
                   outdir = outdir, name="LFC_Comparison_DMedit2_CorrectedByDepr_pem3", w=12, h=8, 
                   scale_mode = "fixed")
compareLFCContrats2(contrastlist2, firstContrast, 
                    contrastNamesOrdered2, mainContrastName, 
                    plim_select= 0.05, plim_plot=0.1,
                    name2remove = name2remove2,
                    resdfname="resdf", 
                    outdir = outdir, name="LFC_Comparison_DMedit2_CorrectedByDepr_p05", w=12, h=12, scale_mode = "fixed")

compareLFCContrats2(contrastlist2, firstContrast, 
                    contrastNamesOrdered2, mainContrastName, 
                    plim_select= 0.01, plim_plot=0.01,
                    name2remove = name2remove2,
                    resdfname="resdf", 
                    outdir = outdir, name="LFC_Comparison_DMedit2_CorrectedByDepr_p01", w=12, h=10, 
                    scale_mode = "fixed")

compareLFCContrats2(contrastlist2, firstContrast, 
                    contrastNamesOrdered2, mainContrastName, 
                    plim_select= 0.001, plim_plot=0.1,
                    name2remove = name2remove2,
                    resdfname="resdf", 
                    outdir = outdir, name="LFC_Comparison_DMedit2_CorrectedByDepr_pem3", w=12, h=10, 
                    scale_mode = "fixed")

names(contrastlist2) <- gsub(" ", "_", names(contrastlist2))
contrastNamesOrdered2[2:length(contrastNamesOrdered2)] <- gsub(" ", "_",contrastNamesOrdered2[2:length(contrastNamesOrdered2)])
cont001 <- compareLFCContratsNumeric(contrastlist2, firstContrast, 
                                     contrastNamesOrdered2, mainContrastName, 
                                     plim_select= 0.001, plim_plot=0.1,
                                     name2remove = name2remove2,
                                     resdfname="resdf", 
                                     outdir = outdir, name="LFC_Comparison_CorrectedByDepr_OnlyDepr_pem3", w=12, h=8, 
                                     scale_mode = "free")
cont01 <- compareLFCContratsNumeric(contrastlist2, firstContrast, 
                                    contrastNamesOrdered2, mainContrastName, 
                                    plim_select= 0.01, plim_plot=0.1,
                                    name2remove = name2remove2,
                                    resdfname="resdf", 
                                    outdir = outdir, name="LFC_Comparison_CorrectedByDepr_OnlyDepr_pem2", w=8, h=8, 
                                    scale_mode = "free")
cont05 <- compareLFCContratsNumeric(contrastlist2, firstContrast, 
                                    contrastNamesOrdered2, mainContrastName, 
                                    plim_select= 0.05, plim_plot=0.1,
                                    name2remove = name2remove2,
                                    resdfname="resdf", 
                                    outdir = outdir, name="LFC_Comparison_CorrectedByDepr_OnlyDepr_p05", w=8, h=12, 
                                    scale_mode = "free")

## Escala de Beck dentro de Depr vs Mediterranea dentro de Depr
fnames <- list.files(paste0(opt$out, "/DESeq2_ByGroup/DeSEQ2/remove_tanda2_onlyControl_Mediterranean_diet_adherence//"), 
                     pattern = "_vs_", full.names = T ) %>% 
  subset(grepl("Normal.tsv", .))
firstContrast2 <- list(resdf=daa_all_bygroup$remove_tanda2$Depression$Escala_depresión_Beck)
contrastNamesOrdered2 <- c("Beck", names(contrastlist2))
compareLFCContrats(contrastlist2, firstContrast2, 
                   contrastNamesOrdered2, mainContrastName, 
                   plim_select= 0.05, plim_plot=0.1,
                   name2remove = name2remove2,
                   resdfname="resdf", 
                   outdir = outdir, name="LFC_Comparison_BeckInside_DMedit_OnlyDepr_p05", w=12, h=8, scale_mode = "free")
compareLFCContrats(contrastlist2, firstContrast, 
                   contrastNamesOrdered2, mainContrastName, 
                   plim_select= 0.001, plim_plot=0.1,
                   name2remove = name2remove2,
                   resdfname="resdf", 
                   outdir = outdir, name="LFC_Comparison_BeckInside_DMedit_OnlyDepr_pem3", w=12, h=8, 
                   scale_mode = "free")


# Mediterranea, Adjusted by Depression

#Euroqol 
fnames <- list.files(paste0(opt$out,"/DESeq2_ByGroup/DeSEQ2/remove_tanda2_onlyDepression_Euroqol/"), 
                     pattern = ".tsv", full.names = T ) %>% 
  subset(grepl("Normal.tsv", .))
contrastlist2 <- map(fnames, read_tsv) %>% lapply(\(x)return(list(resdf=x)))

names(contrastlist2) <- c("Euroqol")
name2remove2 <- "xxx"

contrastNamesOrdered2 <- c("Depression vs Control", names(contrastlist2))
compareLFCContrats(contrastlist2, firstContrast, 
                   contrastNamesOrdered2, mainContrastName, 
                   plim_select= 0.05, plim_plot=0.1,
                   name2remove = name2remove2,
                   resdfname="resdf", 
                   outdir = outdir, name="LFC_Comparison_Euroqol_OnlyDepr_p05", w=12, h=8, scale_mode = "fixed")
compareLFCContrats(contrastlist2, firstContrast, 
                   contrastNamesOrdered2, mainContrastName, 
                   plim_select= 0.001, plim_plot=0.1,
                   name2remove = name2remove2,
                   resdfname="resdf", 
                   outdir = outdir, name="LFC_Comparison_Euroqol_OnlyDepr_pem3", w=12, h=8, 
                   scale_mode = "fixed")

firstContrast2 <- list(resdf=daa_all_bygroup$remove_tanda2$Depression$Escala_depresión_Beck)
contrastNamesOrdered2 <- c("Beck", names(contrastlist2))
compareLFCContrats(contrastlist2, firstContrast2, 
                   contrastNamesOrdered2, mainContrastName, 
                   plim_select= 0.05, plim_plot=0.1,
                   name2remove = name2remove2,
                   resdfname="resdf", 
                   outdir = outdir, name="LFC_Comparison_BeckInside_Euroqol_OnlyDepr_p05", w=12, h=8, scale_mode = "free")
compareLFCContrats(contrastlist2, firstContrast, 
                   contrastNamesOrdered2, mainContrastName, 
                   plim_select= 0.001, plim_plot=0.1,
                   name2remove = name2remove2,
                   resdfname="resdf", 
                   outdir = outdir, name="LFC_Comparison_BeckInside_Euroqol_OnlyDepr_pem3", w=12, h=8, 
                   scale_mode = "free")


#DII 
fnames <- list.files(paste0(opt$out,"/DESeq2_ByGroup/DeSEQ2/remove_tanda2_onlyDepression_DII/"), 
                     pattern = ".tsv", full.names = T ) %>% 
  subset(grepl("Normal.tsv", .))
contrastlist2 <- map(fnames, read_tsv) %>% lapply(\(x)return(list(resdf=x)))

names(contrastlist2) <- c("DII")
name2remove2 <- "xxx"

contrastNamesOrdered2 <- c("Depression vs Control", names(contrastlist2))
compareLFCContrats(contrastlist2, firstContrast, 
                   contrastNamesOrdered2, mainContrastName, 
                   plim_select= 0.05, plim_plot=0.1,
                   name2remove = name2remove2,
                   resdfname="resdf", 
                   outdir = outdir, name="LFC_Comparison_DII_OnlyDepr_p05", w=12, h=8, scale_mode = "fixed")
compareLFCContrats(contrastlist2, firstContrast, 
                   contrastNamesOrdered2, mainContrastName, 
                   plim_select= 0.001, plim_plot=0.1,
                   name2remove = name2remove2,
                   resdfname="resdf", 
                   outdir = outdir, name="LFC_Comparison_DII_OnlyDepr_pem3", w=12, h=8, 
                   scale_mode = "fixed")

compareLFCContrats2(contrastlist2, firstContrast, 
                    contrastNamesOrdered2, mainContrastName, 
                    plim_select= 0.05, plim_plot=0.1,
                    name2remove = name2remove2,
                    resdfname="resdf", 
                    outdir = outdir, name="LFC_Comparison_DII_OnlyDepr_p05", w=12, h=12, scale_mode = "fixed")
compareLFCContrats2(contrastlist2, firstContrast, 
                    contrastNamesOrdered2, mainContrastName, 
                    plim_select= 0.001, plim_plot=0.1,
                    name2remove = name2remove2,
                    resdfname="resdf", 
                    outdir = outdir, name="LFC_Comparison_DII_OnlyDepr_pem3", w=12, h=10, 
                    scale_mode = "fixed")

names(contrastlist2) <- gsub(" ", "_", names(contrastlist2))
contrastNamesOrdered2[2:length(contrastNamesOrdered2)] <- gsub(" ", "_",contrastNamesOrdered2[2:length(contrastNamesOrdered2)])
cont001 <- compareLFCContratsNumeric(contrastlist2, firstContrast, 
                                     contrastNamesOrdered2, mainContrastName, 
                                     plim_select= 0.001, plim_plot=0.1,
                                     name2remove = name2remove2,
                                     resdfname="resdf", 
                                     outdir = outdir, name="LFC_Comparison_DII_OnlyDepr_pem3", w=12, h=8, 
                                     scale_mode = "free")
cont01 <- compareLFCContratsNumeric(contrastlist2, firstContrast, 
                                    contrastNamesOrdered2, mainContrastName, 
                                    plim_select= 0.01, plim_plot=0.1,
                                    name2remove = name2remove2,
                                    resdfname="resdf", 
                                    outdir = outdir, name="LFC_Comparison_DII_OnlyDepr_pem2", w=8, h=8, 
                                    scale_mode = "free")
cont05 <- compareLFCContratsNumeric(contrastlist2, firstContrast, 
                                    contrastNamesOrdered2, mainContrastName, 
                                    plim_select= 0.05, plim_plot=0.1,
                                    name2remove = name2remove2,
                                    resdfname="resdf", 
                                    outdir = outdir, name="LFC_Comparison_DII_OnlyDepr_p05", w=8, h=12, 
                                    scale_mode = "free")

firstContrast2 <- list(resdf=daa_all_bygroup$remove_tanda2$Depression$Escala_depresión_Beck)
contrastNamesOrdered2 <- c("Beck", names(contrastlist2))
compareLFCContrats(contrastlist2, firstContrast2, 
                   contrastNamesOrdered2, mainContrastName, 
                   plim_select= 0.05, plim_plot=0.1,
                   name2remove = name2remove2,
                   resdfname="resdf", 
                   outdir = outdir, name="LFC_Comparison_BeckInside_Euroqol_OnlyDepr_p05", w=12, h=8, scale_mode = "free")
compareLFCContrats(contrastlist2, firstContrast, 
                   contrastNamesOrdered2, mainContrastName, 
                   plim_select= 0.001, plim_plot=0.1,
                   name2remove = name2remove2,
                   resdfname="resdf", 
                   outdir = outdir, name="LFC_Comparison_BeckInside_Euroqol_OnlyDepr_pem3", w=12, h=8, 
                   scale_mode = "free")

## CorrectedByAge and BMI against Beck levels in depressed
firstContrast_edad <- list(resdf=daa_all_corrected$remove_tanda2$Edad_log)
firstContrast_IMC <- list(resdf=daa_all_corrected$remove_tanda2$BMI_log)

fnames <- list.files(paste0(opt$out,"/DESeq2_ByGroup/DeSEQ2/remove_tanda2_onlyDepression_Beck_cualitativo"), 
                     pattern = "_vs_", full.names = T ) %>% 
  subset(grepl("Normal.tsv", .))
contrastlist2 <- map(fnames, read_tsv) %>% lapply(\(x)return(list(resdf=x)))
newnames <- basename(fnames) %>% sapply(\(x)strsplit(x, "Beck_cualitativo_")[[1]][3]) %>% 
  gsub("_DAAshrinkNormal.tsv", "", .) %>% gsub("_", " ", .)
names(contrastlist2) <- newnames
name2remove2 <- "xxx"
contrastNamesOrdered2 <- c("Depression vs Control", names(contrastlist2))

compareLFCContrats(contrastlist2, firstContrast_edad, 
                   contrastNamesOrdered2, mainContrastName, 
                   plim_select= 0.05, plim_plot=0.1,
                   name2remove = name2remove2,
                   resdfname="resdf", 
                   outdir = outdir, name="LFC_Comparison_CorrectAge_BeckQual_OnlyDepr_p05", w=12, h=8, scale_mode = "fixed")
compareLFCContrats(contrastlist2, firstContrast_edad,
                   contrastNamesOrdered2, mainContrastName,
                   plim_select= 0.001, plim_plot=0.1,
                   name2remove = name2remove2,
                   resdfname="resdf", 
                   outdir = outdir, name="LFC_Comparison_CorrectAge_BeckQual_OnlyDepr_pem3", w=12, h=8, 
                   scale_mode = "fixed")
compareLFCContrats(contrastlist2, firstContrast_IMC, 
                   contrastNamesOrdered2, mainContrastName, 
                   plim_select= 0.05, plim_plot=0.1,
                   name2remove = name2remove2,
                   resdfname="resdf", 
                   outdir = outdir, name="LFC_Comparison_CorrectIMC_BeckQual_OnlyDepr_p05", w=12, h=8, scale_mode = "fixed")
compareLFCContrats(contrastlist2, firstContrast_IMC, 
                   contrastNamesOrdered2, mainContrastName, 
                   plim_select= 0.001, plim_plot=0.1,
                   name2remove = name2remove2,
                   resdfname="resdf", 
                   outdir = outdir, name="LFC_Comparison_CorrectIMC_BeckQual_OnlyDepr_pem3", w=12, h=8, 
                   scale_mode = "fixed")

#Correceted by Age and IMC, Again with quantitative scales inside depr
contrastlist2 <- list(daa_all_bygroup$remove_tanda2$Depression$Escala_depresión_Beck,
                      daa_all_bygroup$remove_tanda2$Depression$Montgomery.Asberg,
                      daa_all_bygroup$remove_tanda2$Depression$Escala_Hamilton,
                      daa_all_bygroup$remove_tanda2$Depression$DMSV_puntuacion_total
) %>% lapply(\(x)return(list(resdf=x)))
name2remove2 <- "xxx"
names(contrastlist2) <- c("Beck", "Montgomery", "Hamiliton", "DMSV")
contrastNamesOrdered2 <- c("Depression vs Control", names(contrastlist2))
compareLFCContrats(contrastlist2, firstContrast_edad, 
                   contrastNamesOrdered2, mainContrastName, 
                   plim_select= 0.05, plim_plot=0.1,
                   name2remove = name2remove2,
                   resdfname="resdf", 
                   outdir = outdir, name="LFC_Comparison_CorrAge_AllNum_OnlyDepr_p05", w=12, h=8, scale_mode = "free")
compareLFCContrats(contrastlist2, firstContrast_edad, 
                   contrastNamesOrdered2, mainContrastName, 
                   plim_select= 0.001, plim_plot=0.1,
                   name2remove = name2remove2,
                   resdfname="resdf", 
                   outdir = outdir, name="LFC_Comparison_CorrAge_AllNum_OnlyDepr_pem3", w=12, h=8, 
                   scale_mode = "free")

compareLFCContrats(contrastlist2, firstContrast_IMC, 
                   contrastNamesOrdered2, mainContrastName, 
                   plim_select= 0.05, plim_plot=0.1,
                   name2remove = name2remove2,
                   resdfname="resdf", 
                   outdir = outdir, name="LFC_Comparison_CorrIMC_AllNum_OnlyDepr_p05", w=12, h=8, scale_mode = "free")
compareLFCContrats(contrastlist2, firstContrast_IMC, 
                   contrastNamesOrdered2, mainContrastName, 
                   plim_select= 0.001, plim_plot=0.1,
                   name2remove = name2remove2,
                   resdfname="resdf", 
                   outdir = outdir, name="LFC_Comparison_CorrIMC_AllNum_OnlyDepr_pem3", w=12, h=8, 
                   scale_mode = "free")

# Compare Original with corrected, population variables
contrastlist2 <- list(daa_all_corrected$remove_tanda2$Edad_log,
                      daa_all_corrected$remove_tanda2$Sexo,
                      daa_all_corrected$remove_tanda2$ob_o_sobrepeso,
                      daa_all_corrected$remove_tanda2$BMI_log
                      #daa_all_corrected$remove_tanda2$Colesterol,
                      #daa_all_corrected$remove_tanda2$TG,
                      #daa_all_corrected$remove_tanda2$Euroqol
) %>% lapply(\(x)return(list(resdf=x)))
name2remove2 <- "xxx"
names(contrastlist2) <- c("Age", "Sex", "Overweight", "BMI")#, "Cholest", "TG", "Euroqol" )
contrastNamesOrdered2 <- c("Depression vs Control", names(contrastlist2))
compareLFCContrats(contrastlist2, firstContrast, 
                   contrastNamesOrdered2, mainContrastName, 
                   plim_select= 0.05, plim_plot=0.1,
                   name2remove = name2remove2,
                   resdfname="resdf", 
                   outdir = outdir, name="LFC_Comparison_withCorrected_p05", w=12, h=8, scale_mode = "free")
compareLFCContrats(contrastlist2, firstContrast, 
                   contrastNamesOrdered2, mainContrastName, 
                   plim_select= 0.001, plim_plot=0.1,
                   name2remove = name2remove2,
                   resdfname="resdf", 
                   outdir = outdir, name="LFC_Comparison_withCorrected_pem3", w=12, h=8, 
                   scale_mode = "free")

##Original vs corrected - treatments
contrastlist2 <- list(daa_all_corrected$remove_tanda2$tratamiento_ansioliticos,
                      daa_all_corrected$remove_tanda2$tratamiento_anticonvulsivos,
                      daa_all_corrected$remove_tanda2$tratamiento_ISRNs,
                      daa_all_corrected$remove_tanda2$tratamiento_antidiabeticos,
                      daa_all_corrected$remove_tanda2$tratamiento_coagul_betabloq_etc
                      #daa_all_corrected$remove_tanda2$TG,
                      #daa_all_corrected$remove_tanda2$Euroqol
) %>% lapply(\(x)return(list(resdf=x)))
name2remove2 <- "xxx"
names(contrastlist2) <- c("Ansiol", "Anticonv", "ISRN", "Antidiab", "Betab")#, "Cholest", "TG", "Euroqol" )
contrastNamesOrdered2 <- c("Depression vs Control", names(contrastlist2))
compareLFCContrats(contrastlist2, firstContrast, 
                   contrastNamesOrdered2, mainContrastName, 
                   plim_select= 0.05, plim_plot=0.1,
                   name2remove = name2remove2,
                   resdfname="resdf", 
                   outdir = outdir, name="LFC_Comparison_withCorrectedTreat_p05", w=12, h=8, scale_mode = "free")
compareLFCContrats(contrastlist2, firstContrast, 
                   contrastNamesOrdered2, mainContrastName, 
                   plim_select= 0.001, plim_plot=0.1,
                   name2remove = name2remove2,
                   resdfname="resdf", 
                   outdir = outdir, name="LFC_Comparison_withCorrectedTreat_pem3", w=12, h=8, 
                   scale_mode = "free")


