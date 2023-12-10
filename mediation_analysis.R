library(tidyverse)
library(bmem)
library(DESeq2)
library(phyloseq)

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
  opt <- list(out ="/home/carmoma/Desktop/202311_DEPRESION/results_rstudio_v2_4/mediation_analysis",
              indir = "/home/carmoma/Desktop/202311_DEPRESION/results_rstudio_v2_4/",
              phyloseq_list = "/home/carmoma/Desktop/202311_DEPRESION/results_rstudio_v2_4/phyloseq/phyloseq_all_list.RData",
              phyloseq_name = "remove_tanda2",
              r_functions="/home/carmoma/Desktop/202311_DEPRESION/depression_scripts/metagenomics_core_functions.R",
              metadata = "/home/carmoma/Desktop/202311_DEPRESION/metadatos_MC_AL12042023_CMcopy.xlsx",
              rewrite=TRUE,
              fc=1, 
              pval=0.05
  )
}
if(! dir.exists(opt$out)){dir.create(opt$out)}

### LOAD DATA
load(opt$phyloseq_list)
phobj <- all_phyloseq[[opt$phyloseq_name]]
phobj <- updatePsWithLogs(phobj, c("Edad", "IMC"))
metadata <- sample_data(phobj)


