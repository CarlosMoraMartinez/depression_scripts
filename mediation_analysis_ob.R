library(tidyverse)
library(DESeq2)
library(phyloseq)
library(janitor)
library(ggvenn)

## Mediation analysis based on tutorials:
## https://advstats.psychstat.org/book/mediation/index.php 
## https://rpubs.com/Momen/485122
## https://cran.r-project.org/web/packages/manymome/vignettes/med_lm.html
## https://stats.stackexchange.com/questions/19599/multiple-mediation-analysis-in-r

### Requirements:

### 1) Show that X is correlated with Y. Regress Y on X to estimate and test the path c
### 2) This step establishes that there is an effect that may be mediated.
### Show that X is correlated with M. Regress M on X to estimate and test path a
### 3) This step essentially involves treating the mediator as if it were an outcome variable.
### Show that M affects Y. Regress Y on both X and M to estimate and test path b
### 4) Note that it is not sufficient just to correlate the mediator with the outcome; the mediator and the outcome may be correlated because they are both caused by the input variable X. Thus, the input variable X must be controlled in establishing the effect of the mediator on the outcome.
### To establish that M completely mediates the X-Y relationship, the effect of X on Y controlling for M (path c′
### ) should be zero. The effects in both Steps 3 and 4 are estimated in the same equation.

### Only 2 and 3 are required. 

### How to meet them? 
### Our functional model is bacteria --> BMI --> depression AND bacteria --> depression. 
### 1) -> Differentially abundant bacteria (i.e, X bacteria are correlated with Y depression). We will select all bacteria in simple DAA
### 2) -> Mediator (BMI) is correlated to X (bacterial species). 
### 3) 

SEED <- 123
MODE = "LOCAL"

if(MODE == "IATA"){
  opt <- list()
}else{
  opt <- list(out ="/home/carmoma/Desktop/202311_DEPRESION/results_rstudio_9/mediation_analysis5_OB/",
              indir = "/home/carmoma/Desktop/202311_DEPRESION/results_rstudio_9/",
              phyloseq_list = "/home/carmoma/Desktop/202311_DEPRESION/results_rstudio_9/phyloseq/phyloseq_all_list.RData",
              phyloseq_name = "remove_tanda2",
              r_functions="/home/carmoma/Desktop/202311_DEPRESION/depression_scripts/metagenomics_core_functions.R",
              r_functions_mediation="/home/carmoma/Desktop/202311_DEPRESION/depression_scripts/mediation_functions.R",
              metadata = "/home/carmoma/Desktop/202311_DEPRESION/metadatos_MC_AL12042023_CMcopy.xlsx",
              rewrite=TRUE,
              fc=1, 
              pval=0.05,
              adjust_pvals = TRUE
  )
}
if(! dir.exists(opt$out)){dir.create(opt$out)}
AJUST_PVALS = opt$adjust_pvals
### LOAD DATA
source(opt$r_functions)
source(opt$r_functions_mediation)
restaurar <- restauraropt_mk(opt)

load(opt$phyloseq_list)
phobj <- all_phyloseq[[opt$phyloseq_name]]
phobj <- updatePsWithLogs(phobj, c("Edad", "IMC"))
metadata <- sample_data(phobj) %>% data.frame

alphadiv <- read_tsv(paste0(opt$indir, "AlphaDiversity/remove_tanda2_rarefied_min_AlphaDiv.tsv")) %>% 
  dplyr::select(sampleID, Observed, Chao1, Shannon, InvSimpson)

metadata2 <- merge(metadata, alphadiv, by="sampleID")

#load(paste0(opt$indir, "DESeq2_ControlVars/DeSEQ2/remove_tanda2_IMC_log/DESEQ2_all_results_remove_tanda2_IMC_log.R"))
load(paste0(opt$indir, "IntegrateAllContrasts/LFC_Comparison_AgeAndOb_allCombos.RData"))
vstdf <- read_tsv(paste0(opt$indir, "DeSEQ2/remove_tanda2/remove_tanda2_vst_counts.tsv"))
normdf <- read_tsv(paste0(opt$indir, "DeSEQ2/remove_tanda2/remove_tanda2_norm_counts.tsv"))

### Do tests


vars2venn <- list(
  "D vs C" = dea2contrasts$firstContrast$resdf %>% dplyr::filter(padj <= 0.05) %>% pull(taxon),
  "D vs C, adj OW" = dea2contrasts$contrastlist2$Condition_corr2ob$resdf %>% dplyr::filter(padj <= 0.05) %>% pull(taxon)
)
gv <- ggvenn(
  vars2venn, columns = names(vars2venn),
  stroke_size = 0.5,
  stroke_color = C_NS,
  fill_color = c(C_CASE, C_CTRL2, C_CTRL, C_CASE2),show_elements = F
)
ggsave(filename = paste0(opt$out, "VennDiagram_CvsD_adjOb.pdf"), gv, width = 6, height = 6)


#### Read data 

expr_df <- vstdf %>% column_to_rownames("gene") %>% 
  as.matrix() %>% t %>% data.frame %>% 
  rownames_to_column("sampleID")

df_all <- merge(metadata2, expr_df, by="sampleID")
df_all$Condition_bin <- ifelse(df_all$Condition=="Depression", TRUE, FALSE)
df_all$IMC_log <- ifelse(df_all$IMC > 25, TRUE, FALSE)

names(df_all) <- gsub("^X\\.", "", names(df_all)) %>% 
                 gsub("\\._", "_", .) %>% 
                 gsub("\\.$", "", .)

expr_df_norm <- normdf %>% column_to_rownames("gene") %>% 
  as.matrix() %>% t %>% data.frame %>% 
  rownames_to_column("sampleID")

df_all_norm <- merge(metadata2, expr_df_norm, by="sampleID")
df_all_norm$Condition_bin <- ifelse(df_all_norm$Condition=="Depression", TRUE, FALSE)
df_all_norm$IMC_log <- ifelse(df_all_norm$IMC > 25, TRUE, FALSE)

names(df_all_norm) <- gsub("^X\\.", "", names(df_all_norm)) %>% 
  gsub("\\._", "_", .) %>% 
  gsub("\\.$", "", .)


write_tsv(df_all, file=paste0(opt$out, "merged_data_vst.tsv"))
write_tsv(df_all_norm, file=paste0(opt$out, "merged_data_norm.tsv"))

### MEDIATION ANALYSIS WITH IMC

getvars_MAIN_ANALYSIS_CondAndIMC <- function(summary_df){
  vars2test <- summary_df %>% 
    dplyr::filter(depr_only_padj < plim &
                    imc_only_padj < plim  & 
                    imc_adjdepr_padj < plim) %>% 
    pull(variable)
  return(vars2test)
}

getvars_CondAdjAndIMC <- function(summary_df){
  vars2test <- summary_df %>% 
    dplyr::filter(depr_only_padj < plim &
                    depr_adjimc_padj < plim &
                    imc_only_padj < plim  & 
                    imc_adjdepr_padj < plim) %>% 
    pull(variable)
  return(vars2test)
}

getvars_onlyAfterAdjusting <- function(summary_df){
  vars2test <- summary_df %>% 
    dplyr::filter(depr_only_padj > plim &
                    depr_adjimc_padj < plim) %>% 
    pull(variable)
  return(vars2test)
}

getvars_onlyBeforeAdj <- function(summary_df){
  vars2test <- summary_df %>% 
    dplyr::filter(depr_only_padj < plim &
                    depr_adjimc_padj > plim) %>% 
    pull(variable) 
  return(vars2test)
}

getvars_onlyNotIMC <- function(summary_df){
  vars2test <- summary_df %>% 
    dplyr::filter(depr_only_padj < plim &
                    depr_adjimc_padj < plim &
                    imc_only_padj > plim  &
                    imc_adjdepr_padj > plim) %>% 
    pull(variable)
  return(vars2test)
}


#####################################################
mediator_name <- "IMC_log"
y_name <- "Condition_bin"
plim <- 0.05
plim_plot <- 0.05

custom_cols <- list(
  main_mediation_BMI_separate=list(total=TRUE, direct=NA, indirect=TRUE),
  main_mediation_BMI_separate2=list(total=TRUE, direct=TRUE, indirect=TRUE),
  mixed_mediation_BMI_CondSigBeforeAndAfterAdj=list(total=TRUE, direct=TRUE, indirect=TRUE),
  opos_mediation_BMI_OnlySigAfterAdjust_separate=list(total=FALSE, direct=TRUE, indirect=TRUE),
  indir_onlySigBeforeAdj=list(total=TRUE, direct=FALSE, indirect=TRUE),
  dir_onlyNotIMC=list(total=TRUE, direct=TRUE, indirect=FALSE)
)

getvars_funcs <- list(
  main_mediation_BMI_separate=getvars_MAIN_ANALYSIS_CondAndIMC,
  main_mediation_BMI_separate2=getvars_MAIN_ANALYSIS_CondAndIMC,
  mixed_mediation_BMI_CondSigBeforeAndAfterAdj=getvars_CondAdjAndIMC,
  opos_mediation_BMI_OnlySigAfterAdjust_separate=getvars_onlyAfterAdjusting,
  indir_onlySigBeforeAdj=getvars_onlyBeforeAdj,
  dir_onlyNotIMC=getvars_onlyNotIMC
)
heights <- list(
  main_mediation_BMI_separate=10,
  main_mediation_BMI_separate2=10,
  mixed_mediation_BMI_CondSigBeforeAndAfterAdj=10,
  opos_mediation_BMI_OnlySigAfterAdjust_separate=10,
  indir_onlySigBeforeAdj=10,
  dir_onlyNotIMC=15
)

allMedPlots <- map(c(0.05),\(plim_plot){
  map(names(getvars_funcs), \(x){
    opt <- restaurar(opt)
    opt$out <- paste0(opt$out, "/", x, "/")
    if(!dir.exists(opt$out)) dir.create(opt$out)
    makeFullMediationAnalysisOb(opt, 
                             getVarsFunction = getvars_funcs[[x]], 
                             mediator_name = mediator_name, 
                             y_name = y_name, 
                             plim = plim, 
                             plim_plot = plim_plot,
                             name = "analysis_IMC_separateModel_vjust",
                             wnet=14, hnet=heights[[x]], wbars=8, hbars=10, wbars2=10, hbars2=10, use_color_scale = FALSE,
                             fix_barplot_limits = TRUE, custom_colors=NULL, make_boxplots = FALSE) # make_boxplots will result in error because imc_log was transformed to binary
    if(plim_plot<1){
      makeFullMediationAnalysisOb(opt, 
                                 getVarsFunction = getvars_funcs[[x]], 
                                 mediator_name = mediator_name, 
                                 y_name = y_name, 
                                 plim = plim, 
                                 plim_plot = plim_plot,
                                 name = "analysis_IMC_separateModel_vjust",
                                 wnet=14, hnet=heights[[x]], wbars=8, hbars=10, wbars2=10, hbars2=10, use_color_scale = FALSE,
                                 fix_barplot_limits = TRUE, custom_colors=custom_cols[[x]], make_boxplots=FALSE)
    }
})})

