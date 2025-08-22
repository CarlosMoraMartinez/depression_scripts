library(tidyverse)
library(wesanderson)
library(LinDA)
library(phyloseq)
library(ggvenn)
library(ggVennDiagram)


make_linDA<- function(phobj, var2test,
                                confounders = c(),
                     outdir,
                     alpha = 0.05,
                     prev.cut = 0.05,
                     lib.cut = 1000,
                     winsor.quan = 0.97, base_level=""){
    otu.tab <- otu_table(phobj)
    meta <- phyloseq::sample_data(phobj) %>% data.frame()
    if(base_level != ""){

      levs <- meta %>% filter(!is.na(!!sym(var2test))) %>%
        pull(!!sym(var2test)) %>% unique
      levs <- c(base_level, levs[levs != base_level])

      meta <- meta %>%
        dplyr::mutate(!!sym(var2test) := factor(!!sym(var2test), levels=levs))
    }

    form <- paste0("~", var2test)
    if(length(confounders) > 0){
      form <- paste0(form, "+", paste(confounders, sep="+", collapse=""))
    }
    linda.obj_cond <- LinDA::linda(otu.tab, meta, formula = form,
                        alpha = 0.05,
                        prev.cut = 0.05,
                        lib.cut = 1000,
                        winsor.quan = 0.97)
    return(linda.obj_cond)

}

pal <- wes_palette("AsteroidCity1", 2)

outdir <- "/home/carlos/Documentos/CORALS/results_rstudio/results_Agosto25_1/DAA_linda/"
if(!dir.exists(outdir)) dir.create(outdir)

opt <- restaurar(opt)

load("/home/carlos/Documentos/CORALS/results_rstudio/results_Agosto25_1//foodPCA/phyloseq_list_foodPCA_withNMF_withInc.RData")


s_meta <- sample_data(all_phyloseq$remove_tanda2) %>% data.frame

food_variables<-c(names(s_meta)[5], names(s_meta)[grep("_clr$", names(s_meta))]) # names(sample_data(all_phyloseq$remove_tanda2))[5:30]
patnames <- c("Preprocessed", "Mediterranean", "Western")
quant_vars <- c("age_T0", "af_extraesc_m_00",
                "z_bmi_00", "z_bmi_01", "inc_z_bmiCOR", # calc by collabs
                "z_t0", "z_t1","inc_z_bmi", # calc by us
                "mg_p_00", "z_waist_00", "z_waist_01","inc_z_waist",
                "nreads"
)
cat_vars <-  c("pattern", "Sex", "hospital", "age_class1", "mother_educ", "exercise_cat",
               "Category_T0", "Category_T1",  # calc by collabs
               "status_c2", "status_c1" # calc by us
)
vars2deseq <- c(food_variables, patnames, cat_vars, quant_vars)

vars_basename <- rep("", length(vars2deseq))
names(vars_basename) <- vars2deseq
vars_basename["status_c2"] <- "Normal"
vars_basename["status_c1"] <- "normal"
vars_basename["Category_T0"] <- "Normal"
vars_basename["Category_T1"] <- "Normal"
vars_basename["mother_educ"] <- "1-3"
vars_basename["exercise_cat"] <- "2h"

daa_all <- list()
daa_all_onlyNormalT0 <- list()

opt$mincount <- 1
phseq_to_use <-  "remove_tanda2" #names(all_phyloseq)#c("remove_tanda2", "rmbatch_tanda", "filt")
opt <- restaurar(opt)
for(phname in phseq_to_use){
  daa_all[[phname]] <- list()
  daa_all_onlyNormalT0[[phname]] <- list()
  for(var in vars2deseq){
    cat("Doing LinDA Analysys for: ", phname, ", var=",var, "(", which(var==vars2deseq), " of ", length(vars2deseq), "), all data\n")
    phobj <- all_phyloseq[[phname]]
    samples <- sample_data(phobj)$sampleID[! is.na(sample_data(phobj)[, var ])]
    phobj_filt <- phyloseq::prune_samples(samples, phobj)

    daa_all[[phname]][[var]] <- make_linDA(phobj_filt, var, outdir, base_level="")

    cat("Doing LinDA Analysys for: ", phname, ", var=",var, "(", which(var==vars2deseq), " of ", length(vars2deseq), "), only normal at T0 \n")
    samples_normal <- sample_data(phobj_filt) %>% data.frame() %>%
      filter(status_c1 == "normal") %>% pull(sampleID)
    phobj_filt_normal <- phyloseq::prune_samples(samples_normal, phobj_filt)

    daa_all_onlyNormalT0[[phname]][[var]] <- make_linDA(phobj_filt, var, outdir, base_level="")


  }
}
save(daa_all, file = paste0(outdir, "/LinDA_all_food_variables_Patterns_250814.RData"))
save(daa_all_onlyNormalT0, file = paste0(outdir, "/LinDA_all_food_variables_Patterns_onlyNormalT0_250814.RData"))

names(daa_all_onlyNormalT0)

### with confounders
confvars <- c( "Sex", "hospital", "age_T0", "mother_educ")
vars2deseq_conf <- vars2deseq[!vars2deseq %in% confvars]
vars2deseq_conf <- vars2deseq_conf[vars2deseq_conf != "age_class1"]

daa_all_adj <- list()
daa_all_onlyNormalT0_adj <- list()

for(phname in phseq_to_use){
  daa_all_adj[[phname]] <- list()
  daa_all_onlyNormalT0_adj[[phname]] <- list()
  for(var in vars2deseq){
    cat("Doing LinDA Analysys for: ", phname, ", var=",var, "(", which(var==vars2deseq), " of ", length(vars2deseq), "), all data\n")
    phobj <- all_phyloseq[[phname]]
    samples <- sample_data(phobj)$sampleID[! is.na(sample_data(phobj)[, var ])]
    phobj_filt <- phyloseq::prune_samples(samples, phobj)

    daa_all_adj[[phname]][[var]] <- make_linDA(phobj_filt, var, confvars, outdir, base_level="")

    cat("Doing LinDA Analysys for: ", phname, ", var=",var, "(", which(var==vars2deseq), " of ", length(vars2deseq), "), only normal at T0 \n")
    samples_normal <- sample_data(phobj_filt) %>% data.frame() %>%
      filter(status_c1 == "normal") %>% pull(sampleID)
    phobj_filt_normal <- phyloseq::prune_samples(samples_normal, phobj_filt)

    daa_all_onlyNormalT0_adj[[phname]][[var]] <- make_linDA(phobj_filt, var, confvars, outdir, base_level="")


  }
}
save(daa_all, file = paste0(outdir, "/LinDA_all_food_variables_Patterns_AdjSexHosAgeEdu_250814.RData"))
save(daa_all_onlyNormalT0, file = paste0(outdir, "/LinDA_all_food_variables_Patterns_onlyNormalT0_AdjSexHosAgeEdu_250814.RData"))
