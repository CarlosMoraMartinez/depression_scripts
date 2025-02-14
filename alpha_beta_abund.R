# Alpha 4 each


## Cualitativas

alpha_indices <- c("Observed", "Chao1", "Shannon", "InvSimpson")
vars2test <- c("status_c2", "Sex", "Category_T0", "Category_T1", "hospital", "tanda", "age_class1")

quant_vars <- c("edad_00_meses", "edad_00", "imc_00", "imc_01", "nreads", "edu_m_00", "mg_p_00", "cintura_00", "z_cintura_00")
vars2log <- c( "edad_00_meses", "imc_00", "imc_01")

quant_vars_ext <- c(quant_vars, paste(vars2log, "_log", sep=""))
interestvar <- "edad_00" # "status_c2"
extravars <- c(quant_vars, vars2test)
extravars <- extravars[extravars != interestvar]
extravars <- extravars[extravars != "age_class1"]
extravars <- extravars[extravars != "edad_00_meses"]

outdir <- paste0(opt$out, "/AlphaDiversity/")
if(!dir.exists(outdir)) dir.create(outdir)

extravars2 <- c("Sex", "imc_00_log", "imc_01_log", "edu_m_00", "mg_p_00", "cintura_00") #"edad_00_meses_log"


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

# Composition 4 each

outdir <- paste0(opt$out, "/DescriptiveAbundances/")
if(!dir.exists(outdir)) dir.create(outdir)
tops <- c(5, 10)
interestvar <- "status_c2"

for(phname in phseq_to_use){
  cat("Doing Abundance Plots for: ", phname, "\n")
  abund_plots <- plotAbundanceFullPipeline(all_phyloseq[[phname]], interestvar, outdir, phname, c("Control", "Depression"), tops)
}
