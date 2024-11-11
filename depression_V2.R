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
library(readxl)

SEED <- 123
MODE = "LOCAL"

if(MODE == "IATA"){
  opt <- list(out ="/home/ccarlos/Documentos/CLIMBOUT_CORALS/results_rstudio/240307_results5_kraken05_prev05_filtNormalBetter/",
              indir = "/home/ccarlos/Documentos/CLIMBOUT_CORALS/results_cluster3_Allsamples2_conf05rl100/mg09_combinempa/" , #results_cluster3_Allsamples2_conf05rl100 results_cluster3_Allsamples1_rl75
              
              r_functions="/home/ccarlos/Documentos/CLIMBOUT_CORALS/scripts_PAR/240307scripts/depression_scripts/metagenomics_core_functions.R",
              predictive_functions="/home/ccarlos/Documentos/CLIMBOUT_CORALS/scripts_PAR/240307scripts/depression_scripts/predictive_functions.R",
              read_metadata_script = "/home/ccarlos/Documentos/CLIMBOUT_CORALS/scripts_PAR/240307scripts/depression_scripts/read_metadata.R",
              create_phyloseq_script = "/home/ccarlos/Documentos/CLIMBOUT_CORALS/scripts_PAR/240307scripts/depression_scripts/generate_phyloseq_objects.R",
              read_otutable_script = "/home/ccarlos/Documentos/CLIMBOUT_CORALS/scripts_PAR/240307scripts/depression_scripts/read_otu_table.R",
              alpha_beta_script = "/home/ccarlos/Documentos/CLIMBOUT_CORALS/scripts_PAR/240307scripts/depression_scripts/alpha_beta_abund.R",
              daa_main_condition = "/home/ccarlos/Documentos/CLIMBOUT_CORALS/scripts_PAR/240307scripts/depression_scripts/daa_main_condition.R",
              make_permanova_script = "/home/ccarlos/Documentos/CLIMBOUT_CORALS/scripts_PAR/240307scripts/depression_scripts/make_permanova.R",
              daa_include_single_covariate = "/home/ccarlos/Documentos/CLIMBOUT_CORALS/scripts_PAR/240307scripts/depression_scripts/daa_include_single_covariate.R",
              daa_only_covariates = "/home/ccarlos/Documentos/CLIMBOUT_CORALS/scripts_PAR/240307scripts/depression_scripts/daa_include_only_covariate.R",
              daa_many_covariates = "/home/ccarlos/Documentos/CLIMBOUT_CORALS/scripts_PAR/240307scripts/depression_scripts/daa_include_several_covariates.R",
              
              metadata = "/home/ccarlos/Documentos/CLIMBOUT_CORALS/METADATA/CORALS_metagenomica_26.10.2022_metadataTeresa1.xlsx",
              metadata_class = "/home/ccarlos/Documentos/CLIMBOUT_CORALS/corals_scripts/classified_kids_CM.csv",
              metadata_riga_45 = "/home/ccarlos/Documentos/CLIMBOUT_CORALS/preprocess_data/muestras_Zaragoza_send.csv",
              metadata_with_origin = "/home/ccarlos/Documentos/CLIMBOUT_CORALS/climbout_sergio/ClimbOut/DataAnalysis/CORALS_metagenomica_26.10.2022.xlsx",
              oms_data = "/home/ccarlos/Documentos/CLIMBOUT_CORALS/corals_scripts/Referencia OMS/",
              predict_2groups = "/home/ccarlos/Documentos/CLIMBOUT_CORALS/scripts_PAR/240307scripts/depression_scripts/predict_2groups.R",
              daa_sep_by_group = "/home/ccarlos/Documentos/CLIMBOUT_CORALS/scripts_PAR/240307scripts/depression_scripts/daa_sep_by_group.R",
              daa_with_scales = "/home/ccarlos/Documentos/CLIMBOUT_CORALS/scripts_PAR/240307scripts/depression_scripts/daa_with_scales.R",
              daa_integrate_with_and_without_correction = "/home/ccarlos/Documentos/CLIMBOUT_CORALS/scripts_PAR/240307scripts/depression_scripts/daa_integrate_with_and_without_correction.R",
              daa_integrate_all_contrasts = "/home/ccarlos/Documentos/CLIMBOUT_CORALS/scripts_PAR/240307scripts/depression_scripts/daa_integrate_all_contrasts.R",
              predict_4groups = "/home/ccarlos/Documentos/CLIMBOUT_CORALS/scripts_PAR/240307scripts/depression_scripts/predict_4groups.R",
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
              num_genes_default=5,
              only_normal_weight=TRUE
  )
}else{
  CODEDIR = "/home/carmoma/Documents/CORALS/scripts_PAR/240806scripts/depression_scripts/"
  opt <- list(out ="/home/carmoma/Documents/CORALS/results_rstudio/results_predictOctober_1/",
              indir = "/home/carmoma/Documents/CORALS/results_cluster_240924/names_changed_all/" , #results_cluster3_Allsamples2_conf05rl100 results_cluster3_Allsamples1_rl75
              
              r_functions=paste0(CODEDIR, "metagenomics_core_functions.R"),
              predictive_functions=paste0(CODEDIR, "predictive_functions.R"),
              read_metadata_script = paste0(CODEDIR, "read_metadata.R"),
              create_phyloseq_script = paste0(CODEDIR, "generate_phyloseq_objects.R"),
              read_otutable_script = paste0(CODEDIR, "read_otu_table.R"),
              alpha_beta_script = paste0(CODEDIR, "alpha_beta_abund.R"),
              daa_main_condition = paste0(CODEDIR, "daa_main_condition.R"),
              make_permanova_script = paste0(CODEDIR, "make_permanova.R"),
              daa_include_single_covariate = paste0(CODEDIR, "daa_include_single_covariate.R"),
              daa_only_covariates = paste0(CODEDIR, "daa_include_only_covariate.R"),
              daa_many_covariates = paste0(CODEDIR, "daa_include_several_covariates.R"),
              predict_2groups = paste0(CODEDIR, "predict_2groups.R"),
              daa_sep_by_group = paste0(CODEDIR, "daa_sep_by_group.R"),
              daa_with_scales = paste0(CODEDIR, "daa_with_scales.R"),
              daa_integrate_with_and_without_correction = paste0(CODEDIR, "daa_integrate_with_and_without_correction.R"),
              daa_integrate_all_contrasts = paste0(CODEDIR, "daa_integrate_all_contrasts.R"),
              predict_4groups = paste0(CODEDIR, "predict_4groups.R"),
              
              metadata = "/home/carmoma/Documents/CORALS/METADATA/metadata_NEW_only_CORALS.csv",
              metadata_class = "/home/carmoma/Documents/CORALS/METADATA/classified_kids_NEWDATA_PROVISIONAL.csv",
              metadata_riga_45 = "/home/carmoma/Documents/CORALS/preprocess_data/muestras_Zaragoza_send.csv",
              metadata_with_origin = "/home/carmoma/Documents/CORALS//climbout_sergio/ClimbOut/DataAnalysis/CORALS_metagenomica_26.10.2022.xlsx",
              oms_data = "/home/carmoma/Documents/CORALS/METADATA/Referencia OMS/",
              functional_dir = "/home/carmoma/Documents/CORALS/results_cluster_240924/Humann3_analisis_funcional/MERGED_changed_names/",
              
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
              num_genes_default=5,
              only_normal_weight=FALSE
  )
}
if(! dir.exists(opt$out)){dir.create(opt$out)}
outdir <- paste0(opt$out, "inputdata/")
if(! dir.exists(outdir)){dir.create(outdir)}

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

# Include only covariates
source(opt$daa_only_covariates)

#  Predict
source(opt$predict_2groups)

# PERMANOVA
source(opt$make_permanova_script)

# DAA correcting by covariates
source(opt$daa_include_single_covariate)

## DAA correcting for several variables at the same time
source(opt$daa_many_covariates)



##Integrate with and without correction
source(opt$daa_integrate_with_and_without_correction)

## DAA only in Depressive subjects
source(opt$daa_sep_by_group)

# Finally, all with scales
source(opt$daa_with_scales)

#Integrate all contrasts
source(opt$daa_integrate_with_and_without_correction)
#load("/home/carmoma/Desktop/202311_DEPRESION/results_rstudio_v2_4/DeSEQ2/remove_tanda2/DESEQ2_all_results_remove_tanda2.R")
source(opt$daa_integrate_all_contrasts)

#####################################
## Predict DEPR + OBESITY

source(opt$predict_4groups)