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
#library(plyr)

SEED <- 123
MODE = "LOCAL"

if(MODE == "IATA"){
  opt <- list()
}else{
  BASEDIR = "/home/carmoma/projects/EvaMouse_Otcubre2024/"
  opt <- list(out =paste0(BASEDIR, "Results_rstudio/results1/"),
              indir = paste0(BASEDIR, "/Results_cluster/Results1/mg09_combinempa/") ,
              r_functions= paste0(BASEDIR, "/code/depression_scripts/metagenomics_core_functions.R"),
              predictive_functions=paste0(BASEDIR, "/code//depression_scripts/predictive_functions.R"),
              read_metadata_script = paste0(BASEDIR, "/code//depression_scripts/read_metadata_EvaMouse.R"),
              create_phyloseq_script = paste0(BASEDIR, "/code/depression_scripts/generate_phyloseq_objects.R"),
              read_otutable_script = paste0(BASEDIR, "/code/depression_scripts/read_otu_table.R"),
              alpha_beta_script = paste0(BASEDIR, "/code/depression_scripts/alpha_beta_abund.R"),
              daa_main_condition = paste0(BASEDIR, "/code/depression_scripts/daa_main_condition.R"),
              make_permanova_script = paste0(BASEDIR, "/code/depression_scripts/make_permanova.R"),
              daa_include_single_covariate = paste0(BASEDIR, "/code/depression_scripts/daa_include_single_covariate.R"),
              daa_only_covariates = paste0(BASEDIR, "/code/depression_scripts/daa_include_only_covariate.R"),
              daa_many_covariates =paste0(BASEDIR, "/code/depression_scripts/daa_include_several_covariates.R"),
              predict_2groups = paste0(BASEDIR, "/code/depression_scripts/predict_2groups.R"),
              daa_sep_by_group = paste0(BASEDIR, "/code/depression_scripts/daa_sep_by_group.R"),
              daa_with_scales = paste0(BASEDIR, "/code/depression_scripts/daa_with_scales.R"),
              daa_integrate_with_and_without_correction = paste0(BASEDIR, "/code//depression_scripts/daa_integrate_with_and_without_correction.R"),
              daa_integrate_all_contrasts =paste0(BASEDIR, "/code/depression_scripts/daa_integrate_all_contrasts.R"),
              predict_4groups = paste0(BASEDIR, "/code/depression_scripts/predict_4groups"),
                                                                 
              metadata = paste0(BASEDIR, "/Metadata/metadata_CM.xlsx"),
              flowcelldata =  paste0(BASEDIR, "/Metadata/flowcell.txt"),

              rewrite=FALSE,
              minfreq = 0.1,
              mincountspersample = 0,
              mincount= 20,
              minsampleswithcount = 0,
              raref_quant = 0.15,
              fc=1, 
              pval=0.01, 
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

# #  Predict
# source(opt$predict_2groups)
# 
# # PERMANOVA
# source(opt$make_permanova_script)
# 
# # DAA correcting by covariates
# source(opt$daa_include_single_covariate)
# 
# ## DAA correcting for several variables at the same time
# source(opt$daa_many_covariates)
# 
# # Include only covariates
# source(opt$daa_only_covariates)
# 
# ##Integrate with and without correction
# source(opt$daa_integrate_with_and_without_correction)
# 
# ## DAA only in Depressive subjects
# source(opt$daa_sep_by_group)
# 
# # Finally, all with scales
# source(opt$daa_with_scales)
# 
# #Integrate all contrasts
# source(opt$daa_integrate_with_and_without_correction)
# #load("/home/carmoma/Desktop/202311_DEPRESION/results_rstudio_v2_4/DeSEQ2/remove_tanda2/DESEQ2_all_results_remove_tanda2.R")
# source(opt$daa_integrate_all_contrasts)
# 
# #####################################
# ## Predict DEPR + OBESITY
# 
# source(opt$predict_4groups)