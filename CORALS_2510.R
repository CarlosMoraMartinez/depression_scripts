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
  opt <- list()
}else{
  CODEDIR = "/home/carlos/Documentos/CORALS/scripts_PAR/240806scripts/depression_scripts/"
  opt <- list(out ="/home/carlos/Documentos/CORALS/results_rstudio/results_Agosto25_1/",
              indir = "/home/carlos/Documentos/CORALS/results_cluster_240924/names_changed_all/" , #results_cluster3_Allsamples2_conf05rl100 results_cluster3_Allsamples1_rl75
              input_funcional = "/home/carlos/Documentos/CORALS/results_cluster_240924/Humann3_analisis_funcional/MERGED_changed_names/",

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
              functional_functions=paste0(CODEDIR, "functional_auxiliary_functions.R"),
              functional_script=paste0(CODEDIR, "functional_daa.R"),


              metadata = "/home/carlos/Documentos/CORALS/METADATA/metadata_NEW_only_CORALS.csv",
              metadata_class = "/home/carlos/Documentos/CORALS/METADATA/classified_kids_NEWDATA_PROVISIONAL_withZval_unfiltered.csv",
              metadata_riga_45 = "/home/carlos/Documentos/CORALS/preprocess_data/muestras_Zaragoza_send.csv",
              metadata_with_origin = "/home/carlos/Documentos/CORALS//climbout_sergio/ClimbOut/DataAnalysis/CORALS_metagenomica_26.10.2022.xlsx",
              oms_data = "/home/carlos/Documentos/CORALS/METADATA/Referencia OMS/",
              metadata_mother = "/home/carlos/Documentos/CORALS/METADATA/CORALS_pabdompgrasanivelsocio_2025.xlsx",
              functional_dir = "/home/carlos/Documentos/CORALS/results_cluster_240924/Humann3_analisis_funcional/MERGED_changed_names/",

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
#load("/home/carlos/Desktop/202311_DEPRESION/results_rstudio_v2_4/DeSEQ2/remove_tanda2/DESEQ2_all_results_remove_tanda2.R")
source(opt$daa_integrate_all_contrasts)

#####################################
## Predict DEPR + OBESITY

source(opt$predict_4groups)
