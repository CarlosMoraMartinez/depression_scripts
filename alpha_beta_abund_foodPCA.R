# Alpha 4 each

load("/home/carlos/Documentos/CORALS/results_rstudio/results_Agosto25_1//foodPCA/phyloseq_list_foodPCA_withNMF_withInc.RData")
## Cualitativas

alpha_indices <- c("Observed", "Chao1", "Shannon", "InvSimpson")
vars2test <- c("educ_m_discrete")


food_variables<- names(sample_data(all_phyloseq$remove_tanda2_rarefied_min))[5:30]
pcnames <- names(sample_data(all_phyloseq$remove_tanda2_rarefied_min))
pcnames <- pcnames[grepl("^PC", pcnames)]

quant_vars <- c(pcnames, food_variables)
vars2log <- food_variables

quant_vars_ext <- c(quant_vars, paste(vars2log, "_log", sep=""))
interestvar <- "edad_00" # "status_c2"
extravars <- c(quant_vars, vars2test)
extravars <- extravars[extravars != interestvar]
extravars <- extravars[extravars != "age_class1"]
extravars <- extravars[extravars != "edad_00_meses"]

opt <- restaurar(opt)
opt$out <- paste0(opt$out, "AlphaBetaFood/")
if(!dir.exists(opt$out)) dir.create(opt$out)
outdir <- paste0(opt$out, "/AlphaDiversity/")
if(!dir.exists(outdir)) dir.create(outdir)

extravars2 <- c() #"edad_00_meses_log"

phseq_to_use <- c("remove_tanda2_rarefied_min") # names(all_phyloseq)
#load(allphyloseqlist_fname)

for(phname in phseq_to_use){
  cat("Alpha diversity in ", phname, "\n")
  phobj <- all_phyloseq[[phname]]
  phobj <- updatePsWithLogs(phobj, vars2log)

  divtab <- calculateAlphaDiversityTable(phobj, outdir, alpha_indices, paste0(phname, "_AlphaDiv") )

  models1 <- makeLinearModelsSingleVariable(divtab, interestvar,
                                            extravars,
                                            alpha_indices,
                                            combos=1,
                                            outdir = outdir, name = paste0(phname, "_AlphaDiv_linMod1var") )


  alphadif <- testDiversityDifferences(divtab, alpha_indices, vars2test, outdir, "AlphaDiv_rawdata")


  divplots <- getAlphaDiversityCustomPlot(phobj, vars2test, quant_vars_ext,
                                opt,
                                indices= alpha_indices,
                                correct_pvalues = T, correct_pvalues_indices = F,
                                name = paste0(phname, "_AlphaDivCusP"), w = 10, h = 4)
  divplots <- getAlphaDiversityCustomPlot(phobj, vars2test, quant_vars_ext,
                                opt,
                                indices= alpha_indices,
                                correct_pvalues = T, correct_pvalues_indices = T,
                                name = paste0(phname, "_AlphaDivCusPAdjInd"), w = 10, h = 4)
  divplots <- getAlphaDiversityCustomPlot(phobj, vars2test, quant_vars_ext,
                                          opt,
                                          indices= alpha_indices,
                                          correct_pvalues = T, correct_pvalues_indices = F,
                                          name = paste0(phname, "_AlphaDivCusPraw"), w = 10, h = 3) #3 better for regressions
}


# Beta 4 each
outdir <- paste0(opt$out, "/BetaDiversity/")
if(!dir.exists(outdir)) dir.create(outdir)

phseq_to_use <- names(all_phyloseq)[3:5]  #[c(9,10,2,7)] # [c(2,3,7,9)]

dists <- c("bray") # "jaccard"
METHODS <- c("PCoA", "NMDS") #, "NMDS"
vars2pcoa <- c(vars2test, quant_vars_ext)
ccaplots <- list()
for(phname in phseq_to_use){
  for(method in METHODS){
    for(dist in dists){
      name <- paste0(phname, "_", dist, "_", method)
      cat("Beta diversity for ", name, "\n")
      if(method != "NMDS"){
        extradims_use <- 2:3
        w <- 12
        h <- 4
      }else{
        extradims_use <- c(2)
        w <- 6
        h <- 4
      }
      logvars <- vars2pcoa[grepl("_log$", vars2pcoa, perl=T)]
      origvars <- gsub("_log", "", logvars)
      phobj <- updatePsWithLogs(all_phyloseq[[phname]], origvars)

      ccaplots[[name]] <- makeAllPCoAs(phobj, outdir,
                                       method = method,
                                       name = name,
                                       dist_type = dist,
                                       dist_name = dist,
                                       vars2plot = vars2pcoa,
                                       extradims = extradims_use,
                                       create_pdfs = T, w=w, h=h) #w=16, h=12
    }}}


