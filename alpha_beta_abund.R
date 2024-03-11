# Alpha 4 each


## Cualitativas

alpha_indices <- c("Observed", "Chao1", "Shannon", "InvSimpson")
vars2test <- c("status_c2", "Sex", "hospital", "tanda")

quant_vars <- c("age_months_t0", "imc_00", "imc_01", "nreads_filtPhylum", "bmi_t1_2")
vars2log <- c( "age_months_t0", "imc_00", "imc_01")

quant_vars_ext <- c(quant_vars, paste(vars2log, "_log", sep=""))
interestvar <- "status_c2"
extravars <- c(quant_vars, vars2test)
extravars <- extravars[extravars != interestvar]

outdir <- paste0(opt$out, "/AlphaDiversity/")
if(!dir.exists(outdir)) dir.create(outdir)

extravars2 <- c("Sex", "age_months_t0_log", "imc_00_log", "imc_01_log")


phseq_to_use <- names(all_phyloseq)
#load(allphyloseqlist_fname)

for(phname in phseq_to_use){
  cat("Alpha diversity in ", phname, "\n")
  phobj <- all_phyloseq[[phname]]
  phobj <- updatePsWithLogs(phobj, vars2log)
  sample_data(phobj)$tanda[is.na(sample_data(phobj)$tanda)] <- 2
  sample_data(phobj)$tanda <- as.character(sample_data(phobj)$tanda)

  divtab <- calculateAlphaDiversityTable(phobj, outdir, alpha_indices, paste0(phname, "_AlphaDiv") )

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

   
  alphadif <- testDiversityDifferences(divtab, alpha_indices, vars2test, outdir, "AlphaDiv_rawdata")

  
  divplots <- getAlphaDiversity(phobj, vars2test, quant_vars_ext,
                                opt,
                                indices= alpha_indices,
                                correct_pvalues = T, correct_pvalues_indices = F,
                                name = paste0(phname, "_AlphaDiv"), w = 12, h = 6)
  divplots <- getAlphaDiversity(phobj, vars2test, quant_vars_ext,
                                opt,
                                indices= alpha_indices,
                                correct_pvalues = T, correct_pvalues_indices = T,
                                name = paste0(phname, "_AlphaDivAdjInd"), w = 12, h = 6)
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
