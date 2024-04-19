# Alpha 4 each


## Cualitativas

alpha_indices <- c("Observed", "Chao1", "Shannon", "InvSimpson")
vars2test <- c("Condition", "Sexo", "PROCEDENCIA", 
               "Estado.civil2", "Educacion", 
               "Fumador", "Colesterol_mayor_200",
               "Mediterranean_diet_adherence",
               "obesidad", "ob_o_sobrepeso", 
               "defecaciones_semana", 
               "bristol_scale_cualitativo",
               "sexocaso", 
               "Tanda",
               "tandacaso", 
               "procedenciacaso", 
               "IPAQ_act_fisica")
vars2test_ampl <- c(vars2test, escalas_qual)

quant_vars <- c("Edad", "Edad_log", "BMI","BMI_log",  "TG", "Colesterol", "Glucosa.ayunas", "DII", "PSS_estres")
quant_vars_ext <- c(quant_vars, escalas_quant)
interestvar <- "Condition"
extravars <- c(quant_vars, vars2test)
extravars <- extravars[extravars != interestvar]

outdir <- paste0(opt$out, "/AlphaDiversity/")
if(!dir.exists(outdir)) dir.create(outdir)


interestvar <- "Condition"
extravars <- c(quant_vars, vars2test)
extravars <- extravars[extravars != interestvar]
extravars2 <- c("Sexo", "Edad_log", "BMI_log", "PSS_estres", "IPAQ_act_fisica", 
                "Mediterranean_diet_adherence", "bristol_scale_cualitativo", 
                "defecaciones_semana")

phseq_to_use <- c("remove_tanda2_rarefied_min")
#load(allphyloseqlist_fname)

for(phname in phseq_to_use){
  cat("Alpha diversity in ", phname, "\n")
  phobj <- all_phyloseq[[phname]]
  phobj <- updatePsWithLogs(phobj, c("Edad", "BMI"))
  sample_data(phobj)$IPAQ_act_fisica <- factor(sample_data(phobj)$IPAQ_act_fisica, levels=c("Low", "Mid", "High"))
  divtab <- calculateAlphaDiversityTable(phobj, outdir, alpha_indices, paste0(phname, "_AlphaDiv") )
  divtab$IPAQ_act_fisica <- factor(divtab$IPAQ_act_fisica, levels=c("Low", "Mid", "High"))
  models1 <- makeLinearModelsSingleVariable(divtab, interestvar, 
                                            extravars, 
                                            alpha_indices, 
                                            combos=1,
                                            outdir = outdir, name = paste0(phname, "_AlphaDiv_linMod1var") )
  
  models2 <- makeLinearModelsSingleVariable(divtab, interestvar, 
                                            extravars2, 
                                            alpha_indices, 
                                            combos=1:3,
                                            outdir = outdir, name = paste0(phname, "_AlphaDiv_linModManyVars") )
   models3 <- makeLinearModelsSingleVariable(divtab %>% dplyr::filter(!is.na(PSS_estres)),"PSS_estres", 
                                             c(interestvar, "BMI_log"), 
                                             alpha_indices, 
                                             combos=1,
                                             outdir = outdir, name = paste0(phname, "_AlphaDiv_linModPSS") )
   
   models4 <- makeLinearModelsSingleVariable(divtab %>% dplyr::filter(!is.na(IPAQ_act_fisica)),"IPAQ_act_fisica", 
                                             c(interestvar, "BMI_log"), 
                                             alpha_indices, 
                                             combos=1,
                                             outdir = outdir, name = paste0(phname, "_AlphaDiv_linModIPAQ") )
   
   models5 <- makeLinearModelsSingleVariable(divtab %>% dplyr::filter(!is.na(BMI_log)),"BMI_log", 
                                             c(interestvar, "IPAQ_act_fisica", "Mediterranean_diet_adherence"), 
                                             alpha_indices, 
                                             combos=1,
                                             outdir = outdir, name = paste0(phname, "_AlphaDiv_linModBMI") )
   
   models6 <- makeLinearModelsSingleVariable(divtab %>% dplyr::filter(!is.na(Mediterranean_diet_adherence)),"Mediterranean_diet_adherence", 
                                             c(interestvar, "BMI_log", "IPAQ_act_fisica"), 
                                             alpha_indices, 
                                             combos=1,
                                             outdir = outdir, name = paste0(phname, "_AlphaDiv_linMEDDIET") )
   
   
  alphadif <- testDiversityDifferences(divtab, alpha_indices, vars2test, outdir, "AlphaDiv_rawdata")
  # Ya se hace dentro de la siguiente funcion
  gipaq <- make_IPAQ_Boxplot(divtab, "IPAQ_act_fisica", test2show = "wilcox.test", 
                             alpha_indices = alpha_indices, outdir = outdir, 
                             name=phname,correct_pvalues = TRUE)
  gipaq <- make_IPAQ_Boxplot(divtab, "IPAQ_act_fisica", test2show = "wilcox.test", 
                             alpha_indices = alpha_indices, outdir = outdir, 
                             name=paste0(phname, "_unadj"),correct_pvalues = FALSE)
  
  divplots <- getAlphaDiversity(phobj, vars2test_ampl, quant_vars_ext,
                                opt,
                                indices= alpha_indices,
                                correct_pvalues = T, correct_pvalues_indices = F,
                                name = paste0(phname, "_AlphaDiv"), w = 10, h = 4)
  divplots <- getAlphaDiversity(phobj, vars2test_ampl, quant_vars_ext,
                                opt,
                                indices= alpha_indices,
                                correct_pvalues = T, correct_pvalues_indices = T,
                                name = paste0(phname, "_AlphaDivAdjInd"), w = 10, h = 4)
}

## REMOVING OUTLIERS
outliers <- c(108, 113) %>% as.character

for(phname in phseq_to_use){
  cat("Alpha diversity in ", phname, "\n")
  phobj <- all_phyloseq[[phname]]
  phobj <- updatePsWithLogs(phobj, c("Edad", "BMI"))
  samples <- sample_data(phobj) %>% data.frame %>% pull(sampleID)
  samples <- samples[! samples %in% outliers ]
  phobj_filt <- phyloseq::prune_samples(samples, phobj)
  sample_data(phobj_filt)$IPAQ_act_fisica <- factor(sample_data(phobj_filt)$IPAQ_act_fisica, levels=c("Low", "Mid", "High"))
  divtab <- calculateAlphaDiversityTable(phobj_filt, outdir, alpha_indices, paste0(phname, "_AlphaDiv_RMOL") )
  divtab$IPAQ_act_fisica <- factor(divtab$IPAQ_act_fisica, levels=c("Low", "Mid", "High"))
  models1 <- makeLinearModelsSingleVariable(divtab, interestvar, 
                                            extravars, 
                                            alpha_indices, 
                                            combos=1,
                                            outdir = outdir, name = paste0(phname, "_AlphaDiv_linMod1var_RMOL") )
  
  models2 <- makeLinearModelsSingleVariable(divtab, interestvar, 
                                            extravars2, 
                                            alpha_indices, 
                                            combos=1:3,
                                            outdir = outdir, name = paste0(phname, "_AlphaDiv_linModManyVars_RMOL") )
  models3 <- makeLinearModelsSingleVariable(divtab %>% dplyr::filter(!is.na(PSS_estres)),"PSS_estres", 
                                            c(interestvar, "BMI_log"), 
                                            alpha_indices, 
                                            combos=1,
                                            outdir = outdir, name = paste0(phname, "_AlphaDiv_linModPSS_RMOL") )
  models4 <- makeLinearModelsSingleVariable(divtab %>% dplyr::filter(!is.na(IPAQ_act_fisica)),"IPAQ_act_fisica", 
                                            c(interestvar, "BMI_log"), 
                                            alpha_indices, 
                                            combos=1,
                                            outdir = outdir, name = paste0(phname, "_AlphaDiv_linModIPAQ_RMOL") )
  
  alphadif <- testDiversityDifferences(divtab, alpha_indices, vars2test, outdir, "AlphaDiv_rawdata_RMOL")
  # Ya se hace dentro de la siguiente funcion
  gipaq <- make_IPAQ_Boxplot(divtab, "IPAQ_act_fisica", test2show = "wilcox.test", 
                             alpha_indices = alpha_indices, outdir = outdir, 
                             name=paste0(phname, "_RMOL"),correct_pvalues = TRUE)
  gipaq <- make_IPAQ_Boxplot(divtab, "IPAQ_act_fisica", test2show = "wilcox.test", 
                             alpha_indices = alpha_indices, outdir = outdir, 
                             name=paste0(phname, "_unadj_RMOL"),correct_pvalues = FALSE)
  
  divplots <- getAlphaDiversity(phobj_filt, vars2test_ampl, quant_vars_ext,
                                opt,
                                indices= alpha_indices,
                                correct_pvalues = T, correct_pvalues_indices = F,
                                name = paste0(phname, "_AlphaDiv_RMOL"), w = 10, h = 4)
  divplots <- getAlphaDiversity(phobj_filt, vars2test_ampl, quant_vars_ext,
                                opt,
                                indices= alpha_indices,
                                correct_pvalues = T, correct_pvalues_indices = T,
                                name = paste0(phname, "_AlphaDivAdjInd_RMOL"), w = 10, h = 4)
}

# Alpha inside depression
outdir <- paste0(opt$out, "/AlphaDiversityInsideDepr/")
if(!dir.exists(outdir)) dir.create(outdir)


interestvar <- "Condition"
quant_vars_onlydepr <-  c(escalas_quant, "PSS_estres", "DII", "BMI_log")
qual_vars_onlydepr <- escalas_qual

phseq_to_use <- c("rarefied_min", "remove_tanda2_rarefied_min", "rmbatch_tanda_raref")
#load(allphyloseqlist_fname)

for(phname in phseq_to_use){
  cat("Alpha diversity inside depression ", phname, "\n")
  phobj <- all_phyloseq[[phname]]
  phobj <- updatePsWithLogs(phobj, c("Edad", "BMI"))
  samples <- sample_data(phobj)$sampleID[unlist(sample_data(phobj)[, interestvar]) == "Depression" ]
  phobj_filt <- phyloseq::prune_samples(samples, phobj)
  #alphadif <- testDiversityDifferences(divtab, alpha_indices, vars2test, outdir, "AlphaDiv_rawdata")
  # Ya se hace dentro de la siguiente funcion
  divplots <- getAlphaDiversity(phobj_filt, c(interestvar, "IPAQ_act_fisica"), quant_vars_onlydepr,
                                opt,
                                indices= alpha_indices,
                                correct_pvalues = T,
                                name = paste0(phname, "_AlphaDivOnlyDepr"), w = 10, h = 4)
}

outliers <- c(108, 113) %>% as.character
for(phname in phseq_to_use){
  cat("Alpha diversity inside depression ", phname, "\n")
  phobj <- all_phyloseq[[phname]]
  phobj <- updatePsWithLogs(phobj, c("Edad", "BMI"))
  samples <- sample_data(phobj)$sampleID[unlist(sample_data(phobj)[, interestvar]) == "Depression" ]
  samples <- samples[! samples %in% outliers ]
  phobj_filt <- phyloseq::prune_samples(samples, phobj)
  #alphadif <- testDiversityDifferences(divtab, alpha_indices, vars2test, outdir, "AlphaDiv_rawdata")
  # Ya se hace dentro de la siguiente funcion
  divplots <- getAlphaDiversity(phobj_filt, c(interestvar, "IPAQ_act_fisica"), quant_vars_onlydepr,
                                opt,
                                indices= alpha_indices,
                                correct_pvalues = T,
                                name = paste0(phname, "_AlphaDivOnlyDepr_RMOUTLIERS"), w = 10, h = 4)
}

quant_vars_onlyctr <-  c("Escala_depresión_Beck", "PSS_estres", "DII", "BMI_log")
qual_vars_onlyctr <- escalas_qual
for(phname in phseq_to_use){
  cat("Alpha diversity inside control ", phname, "\n")
  phobj <- all_phyloseq[[phname]]
  phobj <- updatePsWithLogs(phobj, c("Edad", "BMI"))
  samples <- sample_data(phobj)$sampleID[unlist(sample_data(phobj)[, interestvar]) == "Control" ]
  phobj_filt <- phyloseq::prune_samples(samples, phobj)
  #alphadif <- testDiversityDifferences(divtab, alpha_indices, vars2test, outdir, "AlphaDiv_rawdata")
  # Ya se hace dentro de la siguiente funcion
  divplots <- getAlphaDiversity(phobj_filt, c(interestvar, "IPAQ_act_fisica"), quant_vars_onlyctr,
                                opt,
                                indices= alpha_indices,
                                correct_pvalues = T,
                                name = paste0(phname, "_AlphaDivOnlyCtrl"), w = 10, h = 4)
}
# Beta 4 each
outdir <- paste0(opt$out, "/BetaDiversity/")
if(!dir.exists(outdir)) dir.create(outdir)

dists <- c("bray", "jaccard")
vars2pcoa <- c(quant_vars_ext, vars2test_ampl)
#Repeat some plots changing size
vars2pcoa_long_plots <- c("ob_o_sobrepeso", 
                          "Educacion", 
                          "defecaciones_semana",
                          "bristol_scale_cualitativo",
                          "IPAQ_act_fisica", 
                          "Mediterranean_diet_adherence", 
                          "ob_o_sobrepeso", 
                          "DMSV_puntuacion_total", 
                          escalas_qual)

ccaplots <- list()
for(phname in phseq_to_use){
  for(method in c("PCoA", "NMDS")){
    for(dist in dists){
      name <- paste0(phname, "_", dist, "_", method)
      cat("Beta diversity for ", name, "\n")
      
      if(phname %in% c("remove_tanda2", "remove_tanda2_rarefied_min")){
        reserva1 <- vars2pcoa
        vars2pcoa <- vars2pcoa[vars2pcoa != "Tanda"]
      }
      ccaplots[[name]] <- makeAllPCoAs(all_phyloseq[[phname]], outdir,
                                       method = method,
                                       name = name, 
                                       dist_type = dist, 
                                       dist_name = dist,
                                       vars2plot = vars2pcoa, 
                                       extradims = 2:3, 
                                       create_pdfs = T)
      #Just repeat a few plots with adapted sizes
      xx <- makeAllPCoAs(all_phyloseq[[phname]], outdir,
                                       method = method,
                                       name = paste0(name, "_size2"), 
                                       dist_type = dist, 
                                       dist_name = dist,
                                       vars2plot = vars2pcoa_long_plots, 
                                       extradims = 2:3, 
                                       create_pdfs = T, w=16)
      if(phname %in% c("remove_tanda2", "remove_tanda2_rarefied_min")){
        vars2pcoa <- reserva1
      }
    }}}

# Composition 4 each

outdir <- paste0(opt$out, "/DescriptiveAbundances/")
if(!dir.exists(outdir)) dir.create(outdir)
tops <- c(5, 10)
n_species <- 20
interestvar <- "Condition"

for(phname in phseq_to_use){
  cat("Doing Abundance Plots for: ", phname, "\n")
  abund_plots <- plotAbundanceFullPipeline(all_phyloseq[[phname]], interestvar, outdir, phname, c("Control", "Depression"), tops)
}
