#library(plyr)
library(tidyverse)
library(gt)
library(phyloseq)
library(fantaxtic)
library(readxl)
library(ggpubr)
library(dendextend)
library(DESeq2)
library(EnhancedVolcano)
library(gridExtra)
library(cowplot)
library(pheatmap)
library(HMP)
library(knitr)
library(vegan)
library(janitor)
library(plyr)

SEED <- 123
MODE = "LOCAL"

if(MODE == "IATA"){
  opt <- list(out ="/home/ccarlos/Documentos/202309_DEPRESION/results_rstudio_v2_3/",
            indir = "/home/ccarlos/Documentos/202309_DEPRESION/results_cluster/mg09_combinempa/" ,
            r_functions="/home/ccarlos/repos/depression_analysis/metagenomics_core_functions.R",
            predictive_functions="/home/ccarlos/repos/depression_analysis/predictive_functions.R",
            metadata = "/home/ccarlos/Documentos/202309_DEPRESION/metadatos_MC_AL 12042023_CMcopy.xlsx",
            rewrite=TRUE,
            minfreq = 0.05,
            mincountspersample = 0,
            mincount= 1,
            minsampleswithcount = 0,
            raref_quant = 0.15,
            fc=1, 
            pval=0.05, 
            ptype="adjusted", 
            fctype="shrunk",
            num_genes_default=5
            )
}else{
  opt <- list(out ="/home/carmoma/Desktop/202311_DEPRESION/results_rstudio_9/",
              indir = "/home/carmoma/Desktop/202311_DEPRESION/mg09_combinempa/" ,
              r_functions="/home/carmoma/Desktop/202311_DEPRESION/depression_scripts/metagenomics_core_functions.R",
              predictive_functions="/home/carmoma/Desktop/202311_DEPRESION/depression_scripts/predictive_functions.R",
              read_metadata_script = "/home/carmoma/Desktop/202311_DEPRESION/depression_scripts/read_metadata.R",
              create_phyloseq_script = "/home/carmoma/Desktop/202311_DEPRESION/depression_scripts/generate_phyloseq_objects.R",
              read_otutable_script = "/home/carmoma/Desktop/202311_DEPRESION/depression_scripts/read_otu_table.R",
              alpha_beta_script = "/home/carmoma/Desktop/202311_DEPRESION/depression_scripts/alpha_beta_abund.R",
              daa_main_condition = "/home/carmoma/Desktop/202311_DEPRESION/depression_scripts/daa_main_condition.R",
              make_permanova_script = "/home/carmoma/Desktop/202311_DEPRESION/depression_scripts/make_permanova.R",
              daa_include_single_covariate = "/home/carmoma/Desktop/202311_DEPRESION/depression_scripts/daa_include_single_covariate.R",
              daa_only_covariates = "/home/carmoma/Desktop/202311_DEPRESION/depression_scripts/daa_include_only_covariate.R",
              daa_many_covariates = "/home/carmoma/Desktop/202311_DEPRESION/depression_scripts/daa_include_several_covariates.R",
              metadata = "/home/carmoma/Desktop/202311_DEPRESION/metadatos_MC_AL12042023_CM_corrected.xlsx",
              predict_2groups = "/home/carmoma/Desktop/202311_DEPRESION/depression_scripts/predict_2groups.R",
              daa_sep_by_group = "/home/carmoma/Desktop/202311_DEPRESION/depression_scripts/daa_sep_by_group.R",
              daa_with_scales = "/home/carmoma/Desktop/202311_DEPRESION/depression_scripts/daa_with_scales.R",
              daa_integrate_with_and_without_correction = "/home/carmoma/Desktop/202311_DEPRESION/depression_scripts/daa_integrate_with_and_without_correction.R",
              predict_4groups = "/home/carmoma/Desktop/202311_DEPRESION/depression_scripts/predict_4groups",
              rewrite=FALSE,
              minfreq = 0.05,
              mincountspersample = 0,
              mincount= 1,
              minsampleswithcount = 0,
              raref_quant = 0.15,
              fc=1, 
              pval=0.05, 
              ptype="adjusted", 
              fctype="shrunk",
              num_genes_default=5
  )
}
if(! dir.exists(opt$out)){dir.create(opt$out)}

source(opt$r_functions)

restaurar <- restauraropt_mk(opt)

# Read OTUs
source(opt$read_otutable_script)
# Read MetaData
source(opt$read_metadata_script)
# Create phyloseq objects
source(opt$create_phyloseq_script)
# Alpha and Beta diversity. Descriptive 
source(opt$alpha_beta_script)

# DESeq 4 each
source(opt$daa_main_condition)
#load(paste0(opt$out, "DeSEQ2/DESEQ2_all.RData"))

#  Predict
source(opt$predict_2groups)

# PERMANOVA
source(opt$make_permanova_script)

# DAA correcting by covariates
source(opt$daa_include_single_covariate)

## DAA correcting for several variables at the same time
source(opt$daa_many_covariates)

# Include only covariates
source(opt$daa_only_covariates)

##Integrate with and without correction
source(opt$daa_integrate_with_and_without_correction)

## DAA only in Depressive subjects
source(opt$daa_sep_by_group)

# Finally, all with scales
source(opt$daa_with_scales)

#Integrate all contrasts
source(opt$daa_integrate_with_and_without_correction)
#load("/home/carmoma/Desktop/202311_DEPRESION/results_rstudio_v2_4/DeSEQ2/remove_tanda2/DESEQ2_all_results_remove_tanda2.R")

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
                    outdir = outdir, name="LFC_Comparison_IPAQ_CorrectedByDepr_p05", w=10, h=12, scale_mode = "free")
compareLFCContrats2(contrastlist2, firstContrast, 
                    contrastNamesOrdered2, mainContrastName, 
                    plim_select= 0.001, plim_plot=0.1,
                    name2remove = name2remove2,
                    resdfname="resdf", 
                    outdir = outdir, name="LFC_Comparison_CorrectedByDepr_pem3", w=10, h=8, 
                    scale_mode = "free")
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
                    plim_select= 0.05, plim_plot=0.1,
                    name2remove = name2remove2,
                    resdfname="resdf", 
                    outdir = outdir, name="LFC_Comparison_DMedit_CorrectedByDepr_p05", w=12, h=12, scale_mode = "fixed")

compareLFCContrats2(contrastlist2, firstContrast, 
                    contrastNamesOrdered2, mainContrastName, 
                    plim_select= 0.01, plim_plot=0.01,
                    name2remove = name2remove2,
                    resdfname="resdf", 
                    outdir = outdir, name="LFC_Comparison_DMedit_CorrectedByDepr_p01", w=12, h=10, 
                    scale_mode = "fixed")

compareLFCContrats2(contrastlist2, firstContrast, 
                    contrastNamesOrdered2, mainContrastName, 
                    plim_select= 0.001, plim_plot=0.1,
                    name2remove = name2remove2,
                    resdfname="resdf", 
                    outdir = outdir, name="LFC_Comparison_DMedit_CorrectedByDepr_pem3", w=12, h=10, 
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


#####################################
## Predict DEPR + OBESITY

source(opt$predict_4groups)