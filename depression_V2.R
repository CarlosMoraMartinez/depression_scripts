#library(plyr)
library(tidyverse)
library(phyloseq)
library(readxl)
library(janitor)
library(devtools)

devtools::load_all("/home/carlos/Escritorio/202311_DEPRESION/ReviewJune2025/scripts/G4Micro")


SEED <- 123
MODE = "LOCAL"

if(MODE == "IATA"){
  opt <- list()
}else{
  opt <- list(out ="/home/carlos/Escritorio/202311_DEPRESION/ReviewJune2025/Results_rstudio/results3_pckg/",
              indir = "/home/carlos/Escritorio/202311_DEPRESION/202311_DEPRESION/mg09_combinempa/" ,
              read_metadata_script = "/home/carlos/Escritorio/202311_DEPRESION/ReviewJune2025/scripts/depression_scripts/read_metadata.R",
              create_phyloseq_script = "/home/carlos/Escritorio/202311_DEPRESION/ReviewJune2025/scripts/depression_scripts/generate_phyloseq_objects.R",
              read_otutable_script = "/home/carlos/Escritorio/202311_DEPRESION/ReviewJune2025/scripts/depression_scripts/read_otu_table.R",
              alpha_beta_script = "/home/carlos/Escritorio/202311_DEPRESION/ReviewJune2025/scripts/depression_scripts/alpha_beta_abund.R",
              daa_main_condition = "/home/carlos/Escritorio/202311_DEPRESION/ReviewJune2025/scripts/depression_scripts/daa_main_condition.R",
              make_permanova_script = "/home/carlos/Escritorio/202311_DEPRESION/ReviewJune2025/scripts/depression_scripts/make_permanova.R",
              daa_include_single_covariate = "/home/carlos/Escritorio/202311_DEPRESION/ReviewJune2025/scripts/depression_scripts/daa_include_single_covariate.R",
              daa_only_covariates = "/home/carlos/Escritorio/202311_DEPRESION/ReviewJune2025/scripts/depression_scripts/daa_include_only_covariate.R",
              daa_many_covariates = "/home/carlos/Escritorio/202311_DEPRESION/ReviewJune2025/scripts/depression_scripts/daa_include_several_covariates.R",
              metadata = "/home/carlos/Escritorio/202311_DEPRESION/202311_DEPRESION/metadatos_MC_AL12042023_CM_corrected.xlsx",
              predict_2groups = "/home/carlos/Escritorio/202311_DEPRESION/ReviewJune2025/scripts/depression_scripts/predict_2groups.R",
              daa_sep_by_group = "/home/carlos/Escritorio/202311_DEPRESION/ReviewJune2025/scripts/depression_scripts/daa_sep_by_group.R",
              daa_with_scales = "/home/carlos/Escritorio/202311_DEPRESION/ReviewJune2025/scripts/depression_scripts/daa_with_scales.R",
              daa_integrate_with_and_without_correction = "/home/carlos/Escritorio/202311_DEPRESION/ReviewJune2025/scripts/depression_scripts/daa_integrate_with_and_without_correction.R",
              daa_integrate_all_contrasts = "/home/carlos/Escritorio/202311_DEPRESION/ReviewJune2025/scripts/depression_scripts/daa_integrate_all_contrasts.R",
              predict_4groups = "/home/carlos/Escritorio/202311_DEPRESION/ReviewJune2025/scripts/depression_scripts/predict_4groups.R",
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
              num_genes_default=5 # meaning genes or taxa, depending on the context
  )
}
if(! dir.exists(opt$out)){dir.create(opt$out)}

restaurar <- restauraropt_mk(opt)

# Read OTUs
source(opt$read_otutable_script)
# Read MetaData
source(opt$read_metadata_script)
# Create phyloseq objects
source(opt$create_phyloseq_script)

#Alternatively, simply load the data:
data("all_phyloseq")

# Alpha and Beta diversity. Descriptive
source(opt$alpha_beta_script)

# PERMANOVA
source(opt$make_permanova_script)

# DESeq 4 each
source(opt$daa_main_condition)
#load(paste0(opt$out, "DeSEQ2/DESEQ2_all.RData"))

#  Predict 2 groups
source(opt$predict_2groups)

## Predict DEPR + OBESITY
source(opt$predict_4groups)

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

#load("/home/carlos/Escritorio/202311_DEPRESION/results_rstudio_v2_4/DeSEQ2/remove_tanda2/DESEQ2_all_results_remove_tanda2.R")
source(opt$daa_integrate_all_contrasts)

#####################################

