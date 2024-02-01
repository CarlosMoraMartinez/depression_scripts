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

SEED <- 123
MODE = "LOCAL"

if(MODE == "IATA"){
  opt <- list(out ="/home/ccarlos/Documentos/202309_DEPRESION/results_rstudio_v2_3/",
            indir = "/home/ccarlos/Documentos/202309_DEPRESION/results_cluster/mg09_combinempa/" ,
            r_functions="/home/ccarlos/repos/depression_analysis/metagenomics_core_functions.R",
            predictive_functions="/home/ccarlos/repos/depression_analysis/predictive_functions.R",
            metadata = "/home/ccarlos/Documentos/202309_DEPRESION/metadatos_MC_AL 12042023_CMcopy.xlsx",
            rewrite=TRUE,
            minfreq = 0.05,
            mincountspersample = 0,
            mincount= 1,
            minsampleswithcount = 0,
            raref_quant = 0.15,
            fc=1, 
            pval=0.05, 
            ptype="adjusted", 
            fctype="shrunk",
            num_genes_default=5
            )
}else{
  opt <- list(out ="/home/carmoma/Desktop/202311_DEPRESION/results_rstudio_9/",
              indir = "/home/carmoma/Desktop/202311_DEPRESION/mg09_combinempa/" ,
              r_functions="/home/carmoma/Desktop/202311_DEPRESION/depression_scripts/metagenomics_core_functions.R",
              predictive_functions="/home/carmoma/Desktop/202311_DEPRESION/depression_scripts/predictive_functions.R",
              read_metadata_script = "/home/carmoma/Desktop/202311_DEPRESION/depression_scripts/read_metadata.R",
              create_phyloseq_script = "/home/carmoma/Desktop/202311_DEPRESION/depression_scripts/generate_phyloseq_objects.R",
              read_otutable_script = "/home/carmoma/Desktop/202311_DEPRESION/depression_scripts/read_otu_table.R",
              metadata = "/home/carmoma/Desktop/202311_DEPRESION/metadatos_MC_AL12042023_CM_corrected.xlsx",
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

# Alpha 4 each

## Cualitativas

alpha_indices <- c("Observed", "Chao1", "Shannon", "InvSimpson")
vars2test <- c("Condition", "Sexo", "PROCEDENCIA", "Estado.civil2", "Educacion", 
               "Fumador", "Colesterol_mayor_200",
               "Mediterranean_diet_adherence",
               "obesidad", "ob_o_sobrepeso", "defecaciones_semana",  "bristol_scale_cualitativo",
               "sexocaso", "Tanda", "tandacaso", "procedenciacaso", "IPAQ_act_fisica")
vars2test_ampl <- c(vars2test, escalas_qual)

quant_vars <- c("Edad", "IMC",  "TG", "Colesterol", "Glucosa.ayunas", "DII")
quant_vars_ext <- c(quant_vars, escalas_quant)
interestvar <- "Condition"
extravars <- c(quant_vars, vars2test)
extravars <- extravars[extravars != interestvar]

outdir <- paste0(opt$out, "/AlphaDiversity/")
if(!dir.exists(outdir)) dir.create(outdir)


interestvar <- "Condition"
extravars <- c(quant_vars, vars2test)
extravars <- extravars[extravars != interestvar]
extravars2 <- c("Sexo", "Edad", "IMC")

phseq_to_use <- c("remove_tanda2_rarefied_min")
#load(allphyloseqlist_fname)

for(phname in phseq_to_use){
  cat("Alpha diversity in ", phname, "\n")
  phobj <- all_phyloseq[[phname]]
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
#alphadif <- testDiversityDifferences(divtab, alpha_indices, vars2test, outdir, "AlphaDiv_rawdata")
# Ya se hace dentro de la siguiente funcion
  divplots <- getAlphaDiversity(phobj, vars2test_ampl, quant_vars_ext,
                              opt,
                              indices= alpha_indices,
                              correct_pvalues = T,
                              name = paste0(phname, "_AlphaDiv"))
}

# Alpha inside depression
outdir <- paste0(opt$out, "/AlphaDiversityInsideDepr/")
if(!dir.exists(outdir)) dir.create(outdir)


interestvar <- "Condition"
quant_vars_onlydepr <-  c(escalas_quant, "PSS_estres", "DII", "IMC")
qual_vars_onlydepr <- escalas_qual

phseq_to_use <- c("rarefied_min", "remove_tanda2_rarefied_min", "rmbatch_tanda_raref")
#load(allphyloseqlist_fname)

for(phname in phseq_to_use){
  cat("Alpha diversity inside depression ", phname, "\n")
  phobj <- all_phyloseq[[phname]]
  samples <- sample_data(phobj)$sampleID[unlist(sample_data(phobj)[, interestvar]) == "Depression" ]
  phobj_filt <- phyloseq::prune_samples(samples, phobj)
  #alphadif <- testDiversityDifferences(divtab, alpha_indices, vars2test, outdir, "AlphaDiv_rawdata")
  # Ya se hace dentro de la siguiente funcion
  divplots <- getAlphaDiversity(phobj_filt, interestvar, quant_vars_onlydepr,
                                opt,
                                indices= alpha_indices,
                                correct_pvalues = T,
                                name = paste0(phname, "_AlphaDivOnlyDepr"))
}

# Beta 4 each
outdir <- paste0(opt$out, "/BetaDiversity/")
if(!dir.exists(outdir)) dir.create(outdir)

dists <- c("bray", "jaccard")
vars2pcoa <- c(quant_vars_ext, vars2test_ampl)

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


# DESeq 4 each
daa_all <- list()
vars2deseq <- c("Condition")
opt$mincount <- 1
phseq_to_use <- c("remove_tanda2", "rmbatch_tanda", "filt")
for(phname in phseq_to_use){
  cat("Doing DESeq2 Analysys for: ", phname, "\n")
  daa_all[[phname]] <-deseq_full_pipeline(all_phyloseq[[phname]], phname, vars2deseq, opt)
}
save(daa_all, file = paste0(opt$out, "DeSEQ2/DESEQ2_all.RData"))
#load(paste0(opt$out, "DeSEQ2/DESEQ2_all.RData"))

# Predict  
source(opt$predictive_functions)
#opt$out <- "/home/carmoma/Desktop/202311_DEPRESION/results_rstudio_v2_1/"
vars2pca <- c("Condition", "Sexo", "Edad")
opt$out <- paste0(opt$out, "PredictDAA")
if(!dir.exists(opt$out)) dir.create(opt$out)
opt <- restaurar(opt)

all_model_results <- list()
for(i in phseq_to_use){
  cat("Doing Predictive models for: ", i, "\n")
  all_model_results[[i]] <- list()
  phobj <- all_phyloseq[[i]]
  outdir <- paste0(opt$out, "PredictDAA/", i, "/")
  opt$reserva <- opt$out
  opt$out <- outdir
  if(!dir.exists(opt$out)) dir.create(opt$out)
  
  taxa_padj <- daa_all[[i]]$resdf %>% dplyr::filter(padj <= opt$pval & abs(log2FoldChangeShrink) >= log2(opt$fc) ) %>% 
    pull(taxon)
  taxa_praw <- daa_all[[i]]$resdf %>% dplyr::filter(pvalue <= opt$pval & abs(log2FoldChangeShrink) >= log2(opt$fc) ) %>% 
    pull(taxon)
  df2pca <- if(is.null(daa_all[[i]]$vst_counts_df)){ daa_all[[i]]$norm_counts_df}else{ daa_all[[i]]$vst_counts_df }
  all_pcas_adj <- makeAllPCAs(phobj, df2pca, taxa_padj, vars2pca, opt, "DiffTaxaPadj")
  all_pcas_praw <- makeAllPCAs(phobj, df2pca, taxa_praw, vars2pca, opt, "DiffTaxaPraw")
  
  taxa_padj01 <- daa_all[[i]]$resdf %>% dplyr::filter(padj <= 0.01 & abs(log2FoldChangeShrink) >= log2(opt$fc) ) %>% 
    pull(taxon)
  taxa_padj001 <- daa_all[[i]]$resdf %>% dplyr::filter(padj <= 0.001 & abs(log2FoldChangeShrink) >= log2(opt$fc) ) %>% 
    pull(taxon)
  all_pcas_adj01 <- makeAllPCAs(phobj, df2pca, taxa_padj01, vars2pca, opt, "DiffTaxaPadj01")
  all_pcas_adj001 <- makeAllPCAs(phobj, df2pca, taxa_padj001, vars2pca, opt, "DiffTaxaPadj001")
  
  this_metadata <- sample_data(phobj) %>% data.frame
  all_model_results[[i]][["padj_taxa_taxa"]] <- taxa_padj
  all_model_results[[i]][["praw_taxa_taxa"]] <- taxa_praw
  all_model_results[[i]][["padj_taxa_pcas"]] <- all_pcas_adj
  all_model_results[[i]][["praw_taxa_pcas"]] <- all_pcas_praw
  all_model_results[[i]][["padj_taxa_res"]] <- callDoAllModelsFromALLPCAs(all_pcas_adj, name=paste0(i, "ConditionPadj"), metadata=this_metadata, vars2pca=c("Condition"))
  all_model_results[[i]][["praw_taxa_res"]] <- callDoAllModelsFromALLPCAs(all_pcas_praw, name=paste0(i, "ConditionPraw"), metadata=this_metadata, vars2pca=c("Condition"))
  all_model_results[[i]][["padj_taxa_res01"]] <- callDoAllModelsFromALLPCAs(all_pcas_adj01, name=paste0(i, "ConditionPadj01"), metadata=this_metadata, vars2pca=c("Condition"))
  all_model_results[[i]][["padj_taxa_res001"]] <- callDoAllModelsFromALLPCAs(all_pcas_adj001, name=paste0(i, "ConditionPadj001"), metadata=this_metadata, vars2pca=c("Condition"))
  
  PCs <- all_model_results[[i]][["padj_taxa_res"]]$varnames
  modelo_svm <- all_model_results[[i]][["padj_taxa_res"]]$models$`SVM-linear`$mod_noscale
  all_model_results[[i]][["padj_taxa_res_indiv"]] <- callDoAllModelsFromALLPCAsOriginalVars(all_pcas_adj, PCs, modelo_svm, 
                                                                                            df2pca, 
                                                                                            paste0(i, "_ConditionPadjIndiv"), 
                                                                                            vars2pca=c("Condition"), s_meta,
                                                                                            daa_all[[i]]$resdf, 
                                                                                            topns = c(5, 10, 20))
  
  all_model_results[[i]]$metadata <- sample_data(phobj) %>% data.frame
  opt <- restaurar(opt)
}
opt <- restaurar(opt)
save(all_model_results, file=paste0(opt$out, "PredictDAA/all_model_results.RData"))
#load(file=paste0(opt$out, "PredictDAA/all_model_results.RData"))

#Integrate
opt$out <- paste0(opt$out, "PredictDAA/")
makeLinePlotComparingPhobjs(all_model_results, opt)
## Compare with Bacteria in componets

walk(names(all_model_results), makeLinePlotComparingSamePhobjModels, 
     all_model_results, opt)

## Plot boxplot PCs
pcBoxplots <- map(names(all_model_results), makePCsBoxplot, all_model_results, opt)
names(pcBoxplots) <- names(all_model_results)

## Plot barplot PCs and LFC
pcBarplots <- map(names(all_model_results), makePCBarplot, all_model_results, pcBoxplots, daa_all, opt, w=8, h=14)
names(pcBarplots) <- names(all_model_results)

# Plot KNN (best model)

phname <- "remove_tanda2"
predplots <- map(names(all_model_results), plotAllModelPredictions, all_model_results, opt)

 # Make sure that it works
opt <- restaurar(opt)

library(VennDiagram)
outdir <- paste0(opt$out, "PredictDAA/VennDiagrams/")
if(!dir.exists(outdir)) dir.create(outdir)
fname <- paste0(outdir, "/venndiagram_3noraref_padj.png")
genes2compare <- list(filtered=all_model_results$filt$padj_taxa_taxa,
                      remove=all_model_results$remove_tanda2$padj_taxa_taxa,
                      correct_batch=all_model_results$phseq_batch_tanda$padj_taxa_taxa)
xx <-venn.diagram(genes2compare, output=TRUE, filename=fname, na="remove")

fname <- paste0(outdir, "/venndiagram_3noraref_praw.png")
genes2compare <- list(filtered=all_model_results$filt$praw_taxa_taxa,
                      remove=all_model_results$remove_tanda2$praw_taxa_taxa,
                      correct_batch=all_model_results$phseq_batch_tanda$praw_taxa_taxa)
xx <-venn.diagram(genes2compare, output=TRUE, filename=fname, na="remove")

### 

fname <- paste0(outdir, "/venndiagram_3raref_padj.png")
genes2compare <- list(filtered_min=all_model_results$rarefied_min$padj_taxa_taxa,
                      #filtered_quant=all_model_results$rarefied_quant$padj_taxa_taxa,
                      remove=all_model_results$remove_tanda2_rarefied_min$padj_taxa_taxa,
                      correct_batch=all_model_results$phseq_batch_tanda_raref$padj_taxa_taxa)
xx <-venn.diagram(genes2compare, output=TRUE, filename=fname, na="remove")

fname <- paste0(outdir, "/venndiagram_3raref_praw.png")
genes2compare <- list(filtered_min=all_model_results$rarefied_min$praw_taxa_taxa,
                      #filtered_quant=all_model_results$rarefied_quant$praw_taxa_taxa,
                      remove=all_model_results$remove_tanda2_rarefied_min$praw_taxa_taxa,
                      correct_batch=all_model_results$phseq_batch_tanda_raref$praw_taxa_taxa)
xx <-venn.diagram(genes2compare, output=TRUE, filename=fname, na="remove")

### 

fname <- paste0(outdir, "/venndiagram_filt_rawVSraref_padj.png")
genes2compare <- list(filtered=all_model_results$filt$padj_taxa_taxa,
                      raref_min=all_model_results$rarefied_min$padj_taxa_taxa,
                      raref_quant=all_model_results$rarefied_quant$padj_taxa_taxa)
xx <-venn.diagram(genes2compare, output=TRUE, filename=fname, na="remove")

fname <- paste0(outdir, "/venndiagram_filt_rawVSraref_praw.png")
genes2compare <- list(filtered=all_model_results$filt$praw_taxa_taxa,
                      raref_min=all_model_results$rarefied_min$praw_taxa_taxa,
                      raref_quant=all_model_results$rarefied_quant$praw_taxa_taxa)
xx <-venn.diagram(genes2compare, output=TRUE, filename=fname, na="remove")

#### 

fname <- paste0(outdir, "/venndiagram_remTanda2_rawVSraref_padj.png")
genes2compare <- list(rm_tanda2=all_model_results$remove_tanda2$padj_taxa_taxa,
                      rm_tanda2_raref=all_model_results$remove_tanda2_rarefied_min$padj_taxa_taxa)
xx <-venn.diagram(genes2compare, output=TRUE, filename=fname, na="remove")

fname <- paste0(outdir, "/venndiagram_remTanda2_rawVSraref_praw.png")
genes2compare <- list(rm_tanda2=all_model_results$remove_tanda2$praw_taxa_taxa,
                      rm_tanda2_raref=all_model_results$remove_tanda2_rarefied_min$praw_taxa_taxa)
xx <-venn.diagram(genes2compare, output=TRUE, filename=fname, na="remove")

####
fname <- paste0(outdir, "/venndiagram_batchrm_rawVSraref_padj.png")
genes2compare <- list(batch=all_model_results$rmbatch_tanda$padj_taxa_taxa,
                      shrink=all_model_results$rmbatch_tanda_shrink$padj_taxa_taxa,
                      batch_raref=all_model_results$rmbatch_tanda_raref$padj_taxa_taxa)
xx <-venn.diagram(genes2compare, output=TRUE, filename=fname, na="remove")

fname <- paste0(outdir, "/venndiagram_batchrm_rawVSraref_praw.png")
genes2compare <- list(batch=all_model_results$rmbatch_tanda$praw_taxa_taxa,
                      shrink=all_model_results$rmbatch_tanda_shrink$praw_taxa_taxa,
                      batch_raref=all_model_results$rmbatch_tanda_raref$praw_taxa_taxa)
xx <-venn.diagram(genes2compare, output=TRUE, filename=fname, na="remove")

####

fname <- paste0(outdir, "/venndiagram_remTanda2_rmComorb_padj.png")
genes2compare <- list(rm_tanda2=all_model_results$remove_tanda2$padj_taxa_taxa,
                      rm_tanda2_raref=all_model_results$remove_t2_and_comorb$padj_taxa_taxa)
xx <-venn.diagram(genes2compare, output=TRUE, filename=fname, na="remove")

fname <- paste0(outdir, "/venndiagram_remTanda2_rmComorb_praw.png")
genes2compare <- list(rm_tanda2=all_model_results$remove_tanda2$praw_taxa_taxa,
                      rm_tanda2_raref=all_model_results$remove_t2_and_comorb$praw_taxa_taxa)
xx <-venn.diagram(genes2compare, output=TRUE, filename=fname, na="remove")


# PERMANOVA

exclude_vars <- c("sampleID", "CODIGO", "CP")
dtypes <- c("bray", "jaccard")

outdir <- paste0(opt$out, "PERMANOVA/")
if(!dir.exists(outdir)) dir.create(outdir)
permaresults <- list()

for(i in names(all_phyloseq)){
  cat("Doing PERMANOVA of: ", i)
  phobj <- all_phyloseq[[i]]
  permaresults[[i]] <- lapply(dtypes, FUN=function(dd, phobj, exclude_vars, SEED){
    oname <- paste0(outdir, "permanova_results_", i, "_", dd, ".tsv")
    makePermanova(phobj,dist_method = dd, 
                seed = SEED, 
                exclude_vars = exclude_vars, 
                outname = oname) 
  }, phobj, exclude_vars, SEED)
  names(permaresults[[i]]) <- dtypes
}


# DAA correcting by covariates
phseq_to_correct <- names(all_phyloseq)[4]
interestvar <- "Condition"
vars2test <- c("Tanda", "Edad_log", "Sexo", "ob_o_sobrepeso","obesidad", "IMC_log", "Estado.civil2", 
               "Educacion", "reads_log10", "Fumador", 
               "TG", "TG_mayor_200",  "Colesterol", "Colesterol_mayor_200", 
                "PSS_estres", "Diabetes", "bristol_scale", "bristol_scale_cualitativo",
               "defecaciones_semana", "Euroqol", "IPAQ_act_fisica", "Mediterranean_diet_adherence",
               "tratamiento_ansioliticos", "tratamiento_anticonvulsivos",
               "tratamiento_ISRNs", "tratamiento_antidiabeticos", "tratamiento_coagul_betabloq_etc", 
               "DII"
)
vars2test <- c("IMC_log", "Mediterranean_diet_adherence2", "Mediterranean_diet_adherence", "IPAQ_act_fisica")
opt$reserva_0 <- opt$out
opt$out <- paste0(opt$out, "DESeq2_ControlVars/")
if(!dir.exists(opt$out)) dir.create(opt$out)
daa_all_corrected <- list()
for(phname in phseq_to_correct){
  cat("Doing DESeq2 Analysys with correction for: ", phname, "\n")
  phobj <- all_phyloseq[[phname]]
  phobj <- updatePsWithLogs(phobj, c("Edad", "IMC"))
  daa_all_corrected[[phname]] <- list()
  for(var in vars2test){
    cat("Doing DESeq2 Analysys with correction for: ", phname, '-', var, "\n")
    samples <- sample_data(phobj)$sampleID[! is.na(sample_data(phobj)[, var])]
    phobj_filt <- phyloseq::prune_samples(samples, phobj)
    cases <- sample_data(phobj_filt)[, var] %>% unlist %>%  table
    if(length(which(cases > 0)) < 2 ){cat("Only one level, skipping this variable");next}
    
    vars2deseq <- c(interestvar, var)
    name <- paste0(phname, '_', var)
    res_tmp <-deseq_full_pipeline(phobj_filt, name, vars2deseq, opt)
    daa_all_corrected[[phname]][[var]] <- res_tmp$resdf
  }}
opt <- restaurar(opt)
save(daa_all_corrected, file=paste0(opt$out, "DESeq2_ControlVars/DESEQ2_controlVars_all.RData"))

## Plot heatmap with significant differences in all 
vargroups <- list(
  allvars=names(daa_all_corrected[[1]]),
  pobl1=c("Edad_log", "Sexo", "obesidad"),
  pobl2=c("Edad_log", "Sexo", "ob_o_sobrepeso"),
  pobl3=c("Edad_log", "Sexo", "IMC_log", "TG", "Colesterol"),
  pobl3=c("Edad_log", "Sexo", "ob_o_sobrepeso", "IMC_log", "TG", "Colesterol"),
  trats=vars2test[grep("trat", vars2test)],
  health=c("Diabetes", "TG_mayor_200", "Colesterol_mayor_200", "Mediterranean_diet_adherence", "IPAQ_act_fisica"),
  exer=c("Mediterranean_diet_adherence", "Euroqol","IPAQ_act_fisica", "bristol_scale_cualitativo", "PSS_estres")
)
vargroups <- list(pobl1=names(daa_all_corrected[[1]]))
outdir <- paste0(opt$out, "DESeq2_ControlVars/SummaryHeatmaps/")
if(!dir.exists(outdir)) dir.create(outdir)
for(phname in phseq_to_correct){
  for(vgroup in names(vargroups)){
    oname <- paste0(phname, '_', vgroup, "_p05")
    makeHeatmapsFromMultipleDeseqResults(daa_all_corrected[[phname]], 
                                       daa_all[[phname]]$resdf, 
                                        main_name = interestvar,
                                        vars2plot = vargroups[[vgroup]],
                                        italics_rownames=T, 
                                        pfilt =0.01, pplot=0.05, name=oname, outdir=outdir)
    oname <- paste0(phname, '_', vgroup, "_p01")
    makeHeatmapsFromMultipleDeseqResults(daa_all_corrected[[phname]], 
                                         daa_all[[phname]]$resdf, 
                                         main_name = interestvar,
                                         vars2plot = vargroups[[vgroup]],
                                         italics_rownames=T, 
                                         pfilt =0.01, pplot=0.01,name=oname, outdir=outdir)
  
}}


## DAA correcting for several variables at the same time
phseq_to_correct <- names(all_phyloseq)[4]
interestvar <- "Condition"

groups2test <- list(
  c("Edad_log", "IMC_log"),
  c("Sexo", "Diabetes", "IMC_log"),
  c("Edad_log", "Sexo", "IMC_log"),
  c("Edad_log", "Sexo", "IMC_log", "Diabetes"),
  c("Edad_log", "Sexo", "IMC_log", "Diabetes", "Fumador"),
  c("tratamiento_ansioliticos", "tratamiento_anticonvulsivos", "tratamiento_ISRNs"),
  c("tratamiento_ansioliticos", "tratamiento_anticonvulsivos"),
  c("tratamiento_ansioliticos", "tratamiento_ISRNs"),
  c("tratamiento_ansioliticos", "tratamiento_anticonvulsivos",
    "tratamiento_ISRNs", "tratamiento_antidiabeticos", "tratamiento_coagul_betabloq_etc"),
  c("Edad_log", "Sexo", "IMC_log", "Diabetes",
    "tratamiento_ansioliticos", "tratamiento_anticonvulsivos")
)
getGroupName <- function(groupvars){
  gsub("_log|tratamiento_", "", groupvars, perl=T) %>% 
    gsub(" ", "", .) %>% 
    gsub("anti", "anti_", .) %>% 
    make_clean_names(case="big_camel") %>% 
    gsub("Anti", "A", .) %>% 
    substr(1, 3) %>% 
    paste(collapse="")
}
names(groups2test) <- lapply(groups2test, getGroupName)

opt$reserva_0 <- opt$out
opt$out <- paste0(opt$out, "DESeq2_ControlVarsMany/")
if(!dir.exists(opt$out)) dir.create(opt$out)
daa_all_corrected_groups <- list()
for(phname in phseq_to_correct){
  cat("Doing DESeq2 Analysys with correction for: ", phname, "\n")
  phobj <- all_phyloseq[[phname]]
  phobj <- updatePsWithLogs(phobj, c("Edad", "IMC"))
  daa_all_corrected_groups[[phname]] <- list()
  for(vargroup in groups2test){
    cat("Doing DESeq2 Analysys with correction for: ", phname, '-', paste(vargroup, collapse=", "), "\n")
    samples_with_nas <- sample_data(phobj) %>% data.frame %>% dplyr::select(all_of(c(vargroup))) %>% mutate_all(is.na) %>% apply(MAR=1, any)
    samples <- sample_data(phobj)$sampleID[!samples_with_nas]
    phobj_filt <- phyloseq::prune_samples(samples, phobj)
    
    vars2deseq <- c(interestvar, vargroup)
    vgroupname <- getGroupName(vargroup)
    name <- paste0(phname, '_', vgroupname)
    res_tmp <-deseq_full_pipeline(phobj_filt, name, vars2deseq, opt)
    daa_all_corrected_groups[[phname]][[vgroupname]] <- res_tmp$resdf
  }}
opt$out <- opt$reserva_0
save(daa_all_corrected_groups, file=paste0(opt$out, "DESeq2_ControlVarsMany/DESEQ2_controlVarsMany_all.RData"))

## Make heatmap correcting by group (many variables at a time)
outdir <- paste0(opt$out, "DESeq2_ControlVarsMany/SummaryHeatmaps/")
if(!dir.exists(outdir)) dir.create(outdir)

for(phname in phseq_to_correct){
  oname <- paste0(phname, '_', vgroup, "_p01")
  makeHeatmapsFromMultipleDeseqResults(daa_all_corrected_groups[[phname]], 
                                         daa_all[[phname]]$resdf, 
                                         main_name = interestvar,
                                         vars2plot = names(groups2test),
                                         italics_rownames=T, 
                                         pfilt =0.01, pplot=0.05, name=oname, outdir=outdir)
  oname <- paste0(phname, '_', vgroup, "_p05")
  makeHeatmapsFromMultipleDeseqResults(daa_all_corrected_groups[[phname]], 
                                         daa_all[[phname]]$resdf, 
                                         main_name = interestvar,
                                         vars2plot = names(groups2test),
                                         italics_rownames=T, 
                                         pfilt =0.01, pplot=0.01,name=oname, outdir=outdir)
    
}

## DAA only with covariates, without depression condition
phseq_to_correct <- names(all_phyloseq)[4]
interestvar <- "Condition"
vars2test <- c("Tanda", "Edad_log", "Sexo", "ob_o_sobrepeso","obesidad", "IMC_log", "Estado.civil2", 
               "Educacion", "reads_log10", "Fumador", "DII",
               "TG", "TG_mayor_200",  "Colesterol", "Colesterol_mayor_200", 
               "PSS_estres", "Diabetes", "bristol_scale", "bristol_scale_cualitativo",
               "defecaciones_semana", "Euroqol", "IPAQ_act_fisica", "Mediterranean_diet_adherence",
               "tratamiento_ansioliticos", "tratamiento_anticonvulsivos",
               "tratamiento_ISRNs", "tratamiento_antidiabeticos", "tratamiento_coagul_betabloq_etc" 
)
vars2test <- "IMC_log"
opt$reserva_0 <- opt$out
opt$out <- paste0(opt$out, "DESeq2_ControlVarsAlone/")
if(!dir.exists(opt$out)) dir.create(opt$out)
daa_all_corrected_only <- list()
for(phname in phseq_to_correct){
  cat("Doing DESeq2 Analysys with correction for: ", phname, "\n")
  phobj <- all_phyloseq[[phname]]
  phobj <- updatePsWithLogs(phobj, c("Edad", "IMC"))
  daa_all_corrected_only[[phname]] <- list()
  for(var in vars2test){
    cat("Doing DESeq2 Analysys with correction for: ", phname, '-', var, "\n")
    samples <- sample_data(phobj)$sampleID[! is.na(sample_data(phobj)[, var])]
    phobj_filt <- phyloseq::prune_samples(samples, phobj)
    cases <- sample_data(phobj_filt)[, var] %>% unlist %>%  table
    if(length(which(cases > 0)) < 2 ){cat("Only one level, skipping this variable");next}
    name <- paste0(phname, '_', var)
    res_tmp <-deseq_full_pipeline(phobj_filt, name, var, opt)
    daa_all_corrected_only[[phname]][[var]] <- res_tmp$resdf
  }}
opt$out <- opt$reserva_0
save(daa_all_corrected_only, file=paste0(opt$out, "DESeq2_ControlVarsAlone/DESEQ2_controlVarsAlone_all.RData"))

##Integrate with and without correction
outdir <- paste0(opt$out, "IntegrateWithAndWithoutCorrection/")
if(!dir.exists(outdir)) dir.create(outdir)
firstContrast <- daa_all$remove_tanda2
mainContrastName <- "Depression_vs_Control"
edad_correcting <- read_tsv(paste0(opt$out, "/DESeq2_ControlVars/DeSEQ2/remove_tanda2_Edad_log/remove_tanda2_Edad_log_Edad_log_DAAshrinkNormal.tsv"))
imc_correcting <- read_tsv(paste0(opt$out, "/DESeq2_ControlVars/DeSEQ2/remove_tanda2_IMC_log/remove_tanda2_IMC_log_IMC_log_DAAshrinkNormal.tsv"))
edad_corr2 <- read_tsv(paste0(opt$out,"/DESeq2_ControlVarsMany/DeSEQ2/remove_tanda2_EdaImc/remove_tanda2_EdaImc_Edad_log_DAAshrinkNormal.tsv"))
imc_corr2 <-  read_tsv(paste0(opt$out, "/DESeq2_ControlVarsMany/DeSEQ2/remove_tanda2_EdaImc/remove_tanda2_EdaImc_IMC_log_DAAshrinkNormal.tsv"))

#Sólo he guardado resdf
contrastlist2 <- list(
                      daa_all_corrected$remove_tanda2$Edad_log,
                      daa_all_corrected$remove_tanda2$IMC_log,
                      daa_all_corrected_groups$remove_tanda2$EdaImc,
                      
                      daa_all_corrected_only$remove_tanda2$Edad_log,
                      edad_correcting,
                      edad_corr2,
                      
                      daa_all_corrected_only$remove_tanda2$IMC_log,
                      imc_correcting,
                      imc_corr2
                      
) %>% lapply(\(x)return(list(resdf=x)))
name2remove2 <- "xxx"
names(contrastlist2) <- c(
                          "Condition_corrAge", 
                          "Condition_corrIMC", 
                          "Condition_corr2", 
                          
                          "Age_alone", 
                          "Age_corrCond", 
                          "Age_corr2",
                          
                          "BMI_alone", 
                          "BMI_corrCond",
                          "BMI_corr2")

contrastNamesOrdered2 <- c("Depression vs Control",  gsub("_", " ", names(contrastlist2)))
compareLFCContrats(contrastlist2, firstContrast, 
                   contrastNamesOrdered2, mainContrastName, 
                   plim_select= 0.05, plim_plot=0.1,
                   name2remove = name2remove2,
                   resdfname="resdf", 
                   outdir = outdir, 
                   name="LFC_Comparison_AgeAndBMI_allCombos_p05", w=12, h=12, scale_mode = "free")
compareLFCContrats(contrastlist2, firstContrast, 
                   contrastNamesOrdered2, mainContrastName, 
                   plim_select= 0.001, plim_plot=0.1,
                   name2remove = name2remove2,
                   resdfname="resdf", 
                   outdir = outdir, name="LFC_Comparison_AgeAndBMI_allCombos_pem3", w=12, h=12, 
                   scale_mode = "free")
dea2contrasts <- list(firstContrast = firstContrast, contrastlist2=contrastlist2)
save(dea2contrasts, file = paste0(outdir, "/LFC_Comparison_AgeAndBMI_allCombos.RData"))

### plots only with IMC

contrastlist2 <- list(
  #daa_all_corrected$remove_tanda2$Edad_log,
  daa_all_corrected$remove_tanda2$IMC_log,
  #daa_all_corrected_groups$remove_tanda2$EdaImc,
  
  #daa_all_corrected_only$remove_tanda2$Edad_log,
  #edad_correcting,
  #edad_corr2,
  
  daa_all_corrected_only$remove_tanda2$IMC_log,
  imc_correcting
  #imc_corr2
  
) %>% lapply(\(x)return(list(resdf=x)))
name2remove2 <- "xxx"
names(contrastlist2) <- c(
  #"Condition_corrAge", 
  "Condition_corrIMC", 
  #"Condition_corr2", 
  
  #"Age_alone", 
  #"Age_corrCond", 
  #"Age_corr2",
  
  "BMI_alone", 
  "BMI_corrCond")
  #"BMI_corr2")

contrastNamesOrdered2 <- c("Depression vs Control",  gsub("_", " ", names(contrastlist2)))
compareLFCContrats(contrastlist2, firstContrast, 
                   contrastNamesOrdered2, mainContrastName, 
                   plim_select= 0.05, plim_plot=0.1,
                   name2remove = name2remove2,
                   resdfname="resdf", 
                   outdir = outdir, name="LFC_Comparison_BMIcorrAndAlone_allCombos_p05", w=12, h=8, scale_mode = "free")
compareLFCContrats(contrastlist2, firstContrast, 
                   contrastNamesOrdered2, mainContrastName, 
                   plim_select= 0.001, plim_plot=0.1,
                   name2remove = name2remove2,
                   resdfname="resdf", 
                   outdir = outdir, name="LFC_Comparison_BMIcorrAndAlone_allCombos_pem3", w=12, h=8, 
                   scale_mode = "free")
dea2contrasts <- list(firstContrast = firstContrast, contrastlist2=contrastlist2)
save(dea2contrasts, file = paste0(outdir, "/LFC_Comparison_BMIcorrAndAlone_allCombos.RData"))

### Only Edad
contrastlist2 <- list(
  daa_all_corrected$remove_tanda2$Edad_log,
  #daa_all_corrected$remove_tanda2$IMC_log,
  #daa_all_corrected_groups$remove_tanda2$EdaImc,
  
  daa_all_corrected_only$remove_tanda2$Edad_log,
  edad_correcting
  #edad_corr2,
  
  #daa_all_corrected_only$remove_tanda2$IMC_log,
  #imc_correcting
  #imc_corr2
  
) %>% lapply(\(x)return(list(resdf=x)))
name2remove2 <- "xxx"
names(contrastlist2) <- c(
  "Condition_corrAge", 
  #"Condition_corrIMC", 
  #"Condition_corr2", 
  
  "Age_alone", 
  "Age_corrCond") 
  #"Age_corr2",
  
#  "BMI_alone", 
#  "BMI_corrCond")
#"BMI_corr2")

contrastNamesOrdered2 <- c("Depression vs Control",  gsub("_", " ", names(contrastlist2)))
compareLFCContrats(contrastlist2, firstContrast, 
                   contrastNamesOrdered2, mainContrastName, 
                   plim_select= 0.05, plim_plot=0.1,
                   name2remove = name2remove2,
                   resdfname="resdf", 
                   outdir = outdir, name="LFC_Comparison_AgecorrAndAlone_allCombos_p05", w=12, h=8, scale_mode = "free")
compareLFCContrats(contrastlist2, firstContrast, 
                   contrastNamesOrdered2, mainContrastName, 
                   plim_select= 0.001, plim_plot=0.1,
                   name2remove = name2remove2,
                   resdfname="resdf", 
                   outdir = outdir, name="LFC_Comparison_AgecorrAndAlone_allCombos_pem3", w=12, h=8, 
                   scale_mode = "free")
dea2contrasts <- list(firstContrast = firstContrast, contrastlist2=contrastlist2)
save(dea2contrasts, file = paste0(outdir, "/LFC_Comparison_AgecorrAndAlone_allCombos.RData"))


## DAA only in Depressive subjects
vars2test <- c("Edad_log", "Sexo", "IMC_log", "Estado.civil2", 
               "Educacion", "reads_log10", "Fumador", 
               "TG", "TG_mayor_200",  "Colesterol", "Colesterol_mayor_200", 
               "PSS_estres", "Diabetes", "bristol_scale", "bristol_scale_cualitativo",
               "defecaciones_semana", "Euroqol", "IPAQ_act_fisica", "Mediterranean_diet_adherence",
               "tratamiento_ansioliticos", "tratamiento_anticonvulsivos",
               "tratamiento_ISRNs", "tratamiento_antidiabeticos", "tratamiento_coagul_betabloq_etc",
               "DMSV_puntuacion_total",
               "Escala_depresión_Beck", "Beck_cualitativo", "Escala_Hamilton", "Escala_Hamilton_cualitativo",
               "Montgomery.Asberg", "Montgomery.Asberg_qual", "DII", "DII.cualitativo", "Indice_Alimentación_saludable"
)
vars2test <- c("Euroqol", "IPAQ_act_fisica", "Mediterranean_diet_adherence", "Mediterranean_diet_adherence2")
escalas_qual <- c( "Beck_cualitativo", "Escala_Hamilton_cualitativo","Montgomery.Asberg_qual")
escalas_quant <- c( "Escala_depresión_Beck", "Escala_Hamilton", "Montgomery.Asberg", "DMSV_puntuacion_total")

interestvar <- "Condition"
vlevs <- levels(metadata[, interestvar] %>% unlist)

opt$reserva_0 <- opt$out
opt$out <- paste0(opt$out, "DESeq2_ByGroup/")
if(!dir.exists(opt$out)) dir.create(opt$out)

load(allphyloseqlist_fname)
# for(phname in names(all_phyloseq)){
#   for(v in escalas_qual){
#      sample_data(all_phyloseq[[phname]])[, v] <- factor(unlist(sample_data(all_phyloseq[[phname]])[, v]),
#                                                         levels=c("No Depresion", "Leve", "Moderada",
#                                                                  "Severa" ))
#     sample_data(all_phyloseq[[phname]])[, v] <- dplyr::recode(unlist(sample_data(all_phyloseq[[phname]])[, v]),
#                                                      "No Depresion"="Healthy", 
#                                                      "Leve"="Mild", 
#                                                      "Moderada"="Moderate",
#                                                     "Severa"="Severe" )
#   }
# }

daa_all_bygroup <- list()
for(phname in phseq_to_correct){
  cat("Doing DESeq2 Analysys with correction for: ", phname, "\n")
  phobj <- all_phyloseq[[phname]]
  phobj <- updatePsWithLogs(phobj, c("Edad", "IMC"))
  daa_all_bygroup[[phname]] <- list()
  for(varlev in vlevs){
    daa_all_bygroup[[phname]][[varlev]] <- list()
    samples <- sample_data(phobj)$sampleID[unlist(sample_data(phobj)[, interestvar]) == varlev ]
    phobj_prefilt <- phyloseq::prune_samples(samples, phobj)
    for(var in vars2test){
      if(var %in% c(escalas_qual, escalas_quant) & varlev == "Control") next
      cat("Doing DESeq2 Analysys within group for: ", phname, '-', var,' inside ',interestvar, ':',varlev, "\n")
      samples <- sample_data(phobj_prefilt)$sampleID[!is.na(sample_data(phobj_prefilt)[, var]) ]
      phobj_filt <- phyloseq::prune_samples(samples, phobj_prefilt)
      if(var %in% escalas_qual){
        samples <- sample_data(phobj_filt)$sampleID[!unlist(sample_data(phobj_filt)[, var]) %in% c("No Depression", "Healthy", "No depresion") ]
        phobj_filt <- phyloseq::prune_samples(samples, phobj_filt)
      }
      cases <- sample_data(phobj_filt)[, var] %>% unlist %>%  table
      if(length(which(cases > 1)) < 2 & !is.numeric(unlist(sample_data(phobj_filt)[, var] ))){cat("Only one level, skipping this variable");next}
    
    name <- paste0(phname, '_only',varlev,'_', var)
    res_tmp <-deseq_full_pipeline(phobj_filt, name, var, opt)
    daa_all_bygroup[[phname]][[varlev]][[var]] <- res_tmp$resdf
    } # variables
  } # Levels in interestvar (case/control)
}# phyloseq objects

opt$out <- opt$reserva_0
save(daa_all_bygroup, file=paste0(opt$out, "DESeq2_ControlVars/DESEQ2_byGroup_all.RData"))

# Finally, all with scales

opt$reserva_0 <- opt$out
opt$out <- paste0(opt$out, "DESeq2_MainVariableScales/")
if(!dir.exists(opt$out)) dir.create(opt$out)

daa_all_scales <- list()
vars2deseq <- c(escalas_qual, escalas_quant)
opt$mincount <- 1
for(phname in names(all_phyloseq)){
  phobj <- all_phyloseq[[phname]]
  daa_all_scales[[phname]] <- list()
  for(var in vars2deseq){
    cat("Doing DESeq2 Analysys for: ", phname, ", with main variable ", var, "\n")
    samples <- sample_data(phobj)$sampleID[!is.na(sample_data(phobj)[, var]) ]
    phobj_filt <- phyloseq::prune_samples(samples, phobj)
    daa_all_scales[[phname]][[var]] <-deseq_full_pipeline(phobj_filt, paste0(phname, "_",var), var, opt)
  }}

save(daa_all_scales, file = paste0(opt$out, "DeSEQ2/DESEQ2_MainVariableScales.RData"))
opt$out <- opt$reserva_0

#Integrate all contrasts

#load("/home/carmoma/Desktop/202311_DEPRESION/results_rstudio_v2_4/DeSEQ2/remove_tanda2/DESEQ2_all_results_remove_tanda2.R")

outdir <- paste0(opt$out, "IntegrateAllContrasts/")
if(!dir.exists(outdir)) dir.create(outdir)

contrastlist <- daa_all_scales$remove_tanda2$Beck_cualitativo$all_contrasts
firstContrast <- daa_all$remove_tanda2
mainContrastName <- "Depression_vs_Control"
contrastNamesOrdered <- c("Depression vs Control", "Mild vs Healthy","Moderate vs Healthy", 
                          "Severe vs Healthy", "Moderate vs Mild","Severe vs Mild", "Severe vs Moderate")
plim_select= 0.000001
plim_plot=0.05
name2remove="Beck_cualitativo_"

## Against Beck levels
compareLFCContrats(contrastlist, firstContrast, 
                   contrastNamesOrdered, mainContrastName, 
                   plim_select= 0.05, plim_plot=0.1,
                   name2remove = name2remove,
                   resdfname="resdf", 
                   outdir = outdir, name="LFC_Comparison_BeckQuant_p05", w=12, h=8)
compareLFCContrats(contrastlist, firstContrast, 
                   contrastNamesOrdered, mainContrastName, 
                   plim_select= 0.000001, plim_plot=0.1,
                   name2remove = name2remove,
                   resdfname="resdf", 
                   outdir = outdir, name="LFC_Comparison_BeckQuant_pem6", w=12, h=8)
## Against Montgomery levels
contrastlist2 <- daa_all_scales$remove_tanda2$Montgomery.Asberg_qual$all_contrasts
name2remove2 <- "Montgomery.Asberg_qual_"
compareLFCContrats(contrastlist2, firstContrast, 
                   contrastNamesOrdered, mainContrastName, 
                   plim_select= 0.05, plim_plot=0.1,
                   name2remove = name2remove2,
                   resdfname="resdf", 
                   outdir = outdir, name="LFC_Comparison_MontgomQuant_p05", w=12, h=8)
compareLFCContrats(contrastlist2, firstContrast, 
                   contrastNamesOrdered, mainContrastName, 
                   plim_select= 0.000001, plim_plot=0.1,
                   name2remove = name2remove2,
                   resdfname="resdf", 
                   outdir = outdir, name="LFC_Comparison_MontgomQuant_pem6", w=12, h=8)
## Against Hamilton levels
contrastlist2 <- daa_all_scales$remove_tanda2$Escala_Hamilton_cualitativo$all_contrasts
name2remove2 <- "Escala_Hamilton_cualitativo_"
compareLFCContrats(contrastlist2, firstContrast, 
                   contrastNamesOrdered, mainContrastName, 
                   plim_select= 0.05, plim_plot=0.1,
                   name2remove = name2remove2,
                   resdfname="resdf", 
                   outdir = outdir, name="LFC_Comparison_HamiltonQuant_p05", w=12, h=8)
compareLFCContrats(contrastlist2, firstContrast, 
                   contrastNamesOrdered, mainContrastName, 
                   plim_select= 0.000001, plim_plot=0.1,
                   name2remove = name2remove2,
                   resdfname="resdf", 
                   outdir = outdir, name="LFC_Comparison_HamiltonQuant_pem6", w=12, h=8)

## Against Beck numerical
contrastlist2 <- daa_all_scales$remove_tanda2$Escala_depresión_Beck$all_contrasts
name2remove2 <- "xxx"
names(contrastlist2) <- c("Beck scale")
contrastNamesOrdered2 <- c("Depression vs Control", "Beck scale")
compareLFCContrats(contrastlist2, firstContrast, 
                   contrastNamesOrdered2, mainContrastName, 
                   plim_select= 0.05, plim_plot=0.1,
                   name2remove = name2remove2,
                   resdfname="resdf", 
                   outdir = outdir, name="LFC_Comparison_BeckNum_p05", w=12, h=8, scale_mode = "free")
compareLFCContrats(contrastlist2, firstContrast, 
                   contrastNamesOrdered2, mainContrastName, 
                   plim_select= 0.001, plim_plot=0.1,
                   name2remove = name2remove2,
                   resdfname="resdf", 
                   outdir = outdir, name="LFC_Comparison_BeckNum_pem3", w=12, h=8, 
                   scale_mode = "free")

### Quantitative scales against CtrlVsDepr
contrastlist2 <- list(daa_all_scales$remove_tanda2$Escala_depresión_Beck$all_contrasts[[1]],
                      daa_all_scales$remove_tanda2$Montgomery.Asberg$all_contrasts[[1]],
                      daa_all_scales$remove_tanda2$Escala_Hamilton$all_contrasts[[1]],
                      daa_all_scales$remove_tanda2$DMSV_puntuacion_total$all_contrasts[[1]]
)
name2remove2 <- "xxx"
names(contrastlist2) <- c("Beck", "Montgomery", "Hamiliton", "DMSV")
contrastNamesOrdered2 <- c("Depression vs Control", names(contrastlist2))
compareLFCContrats(contrastlist2, firstContrast, 
                   contrastNamesOrdered2, mainContrastName, 
                   plim_select= 0.05, plim_plot=0.1,
                   name2remove = name2remove2,
                   resdfname="resdf", 
                   outdir = outdir, name="LFC_Comparison_AllNum_p05", w=12, h=8, scale_mode = "free")
compareLFCContrats(contrastlist2, firstContrast, 
                   contrastNamesOrdered2, mainContrastName, 
                   plim_select= 0.001, plim_plot=0.1,
                   name2remove = name2remove2,
                   resdfname="resdf", 
                   outdir = outdir, name="LFC_Comparison_AllNum_pem3", w=12, h=8, 
                   scale_mode = "free")
compareLFCContrats2(contrastlist2, firstContrast, 
                   contrastNamesOrdered2, mainContrastName, 
                   plim_select= 0.05, plim_plot=0.1,
                   name2remove = name2remove2,
                   resdfname="resdf", 
                   outdir = outdir, name="LFC_Comparison_AllNum_p05", w=8, h=8, scale_mode = "free")
compareLFCContrats2(contrastlist2, firstContrast, 
                   contrastNamesOrdered2, mainContrastName, 
                   plim_select= 0.001, plim_plot=0.1,
                   name2remove = name2remove2,
                   resdfname="resdf", 
                   outdir = outdir, name="LFC_Comparison_AllNum_pem3", w=8, h=8, 
                   scale_mode = "free")

cont001 <- compareLFCContratsNumeric(contrastlist2, firstContrast, 
                    contrastNamesOrdered2, mainContrastName, 
                    plim_select= 0.001, plim_plot=0.1,
                    name2remove = name2remove2,
                    resdfname="resdf", 
                    outdir = outdir, name="LFC_Comparison_AllNum_pem3", w=8, h=8, 
                    scale_mode = "free")
cont01 <- compareLFCContratsNumeric(contrastlist2, firstContrast, 
                                     contrastNamesOrdered2, mainContrastName, 
                                     plim_select= 0.01, plim_plot=0.1,
                                     name2remove = name2remove2,
                                     resdfname="resdf", 
                                     outdir = outdir, name="LFC_Comparison_AllNum_pem2", w=8, h=8, 
                                     scale_mode = "free")
cont05 <- compareLFCContratsNumeric(contrastlist2, firstContrast, 
                                    contrastNamesOrdered2, mainContrastName, 
                                    plim_select= 0.05, plim_plot=0.1,
                                    name2remove = name2remove2,
                                    resdfname="resdf", 
                                    outdir = outdir, name="LFC_Comparison_AllNum_05", w=8, h=8, 
                                    scale_mode = "free")



## Quant depression scales Inside depressive subjects
#Sólo he guardado resdf
contrastlist2 <- list(daa_all_bygroup$remove_tanda2$Depression$Escala_depresión_Beck,
                      daa_all_bygroup$remove_tanda2$Depression$Montgomery.Asberg,
                      daa_all_bygroup$remove_tanda2$Depression$Escala_Hamilton,
                      daa_all_bygroup$remove_tanda2$Depression$DMSV_puntuacion_total,
                      daa_all_bygroup$remove_tanda2$Depression$PSS_estres
) %>% lapply(\(x)return(list(resdf=x)))
name2remove2 <- "xxx"
names(contrastlist2) <- c("Beck", "Montgomery", "Hamiliton", "DMSV", "PSS")
contrastNamesOrdered2 <- c("Depression vs Control", names(contrastlist2))
compareLFCContrats(contrastlist2, firstContrast, 
                   contrastNamesOrdered2, mainContrastName, 
                   plim_select= 0.05, plim_plot=0.1,
                   name2remove = name2remove2,
                   resdfname="resdf", 
                   outdir = outdir, name="LFC_Comparison_AllNum_OnlyDepr_p05", w=12, h=8, scale_mode = "free")
compareLFCContrats(contrastlist2, firstContrast, 
                   contrastNamesOrdered2, mainContrastName, 
                   plim_select= 0.001, plim_plot=0.1,
                   name2remove = name2remove2,
                   resdfname="resdf", 
                   outdir = outdir, name="LFC_Comparison_AllNum_OnlyDepr_pem3", w=12, h=8, 
                   scale_mode = "free")
compareLFCContrats2(contrastlist2, firstContrast, 
                    contrastNamesOrdered2, mainContrastName, 
                    plim_select= 0.05, plim_plot=0.05,
                    name2remove = name2remove2,
                    resdfname="resdf", 
                    outdir = outdir, name="LFC_Comparison_AllNum_OnlyDepr_p05", w=12, h=8, scale_mode = "free")
compareLFCContrats2(contrastlist2, firstContrast, 
                    contrastNamesOrdered2, mainContrastName, 
                    plim_select= 0.001, plim_plot=0.05,
                    name2remove = name2remove2,
                    resdfname="resdf", 
                    outdir = outdir, name="LFC_Comparison_AllNum_OnlyDepr_pem3", w=12, h=8, 
                    scale_mode = "free")

cont001 <- compareLFCContratsNumeric(contrastlist2, firstContrast, 
                                     contrastNamesOrdered2, mainContrastName, 
                                     plim_select= 0.001, plim_plot=0.1,
                                     name2remove = name2remove2,
                                     resdfname="resdf", 
                                     outdir = outdir, name="LFC_Comparison_AllNum_OnlyDepr_pem3", w=8, h=8, 
                                     scale_mode = "free")
cont01 <- compareLFCContratsNumeric(contrastlist2, firstContrast, 
                                    contrastNamesOrdered2, mainContrastName, 
                                    plim_select= 0.01, plim_plot=0.1,
                                    name2remove = name2remove2,
                                    resdfname="resdf", 
                                    outdir = outdir, name="LFC_Comparison_AllNum_OnlyDepr_pem2", w=8, h=8, 
                                    scale_mode = "free")
cont05 <- compareLFCContratsNumeric(contrastlist2, firstContrast, 
                                    contrastNamesOrdered2, mainContrastName, 
                                    plim_select= 0.05, plim_plot=0.1,
                                    name2remove = name2remove2,
                                    resdfname="resdf", 
                                    outdir = outdir, name="LFC_Comparison_AllNum_OnlyDepr_05", w=8, h=8, 
                                    scale_mode = "free")


## Beck, DSMV and stress Inside depressive subjects
#Sólo he guardado resdf
contrastlist2 <- list(daa_all_bygroup$remove_tanda2$Depression$Escala_depresión_Beck,
                      daa_all_bygroup$remove_tanda2$Depression$DMSV_puntuacion_total,
                      daa_all_bygroup$remove_tanda2$Depression$PSS_estres
                      
) %>% lapply(\(x)return(list(resdf=x)))
name2remove2 <- "xxx"
names(contrastlist2) <- c("Beck", "DMSV", "Stress")
contrastNamesOrdered2 <- c("Depression vs Control", names(contrastlist2))
compareLFCContrats(contrastlist2, firstContrast, 
                   contrastNamesOrdered2, mainContrastName, 
                   plim_select= 0.05, plim_plot=0.1,
                   name2remove = name2remove2,
                   resdfname="resdf", 
                   outdir = outdir, name="LFC_Comparison_BeckDMSStress_OnlyDepr_p05", w=12, h=8, scale_mode = "free")
compareLFCContrats(contrastlist2, firstContrast, 
                   contrastNamesOrdered2, mainContrastName, 
                   plim_select= 0.01, plim_plot=0.1,
                   name2remove = name2remove2,
                   resdfname="resdf", 
                   outdir = outdir, name="LFC_Comparison_BeckDMSStress_OnlyDepr_pem3", w=12, h=8, 
                   scale_mode = "free")

compareLFCContrats2(contrastlist2, firstContrast, 
                    contrastNamesOrdered2, mainContrastName, 
                    plim_select= 0.05, plim_plot=0.1,
                    name2remove = name2remove2,
                    resdfname="resdf", 
                    outdir = outdir, name="LFC_Comparison_BeckDMSStress_OnlyDepr_05", w=8, h=8, scale_mode = "free")
compareLFCContrats2(contrastlist2, firstContrast, 
                    contrastNamesOrdered2, mainContrastName, 
                    plim_select= 0.001, plim_plot=0.1,
                    name2remove = name2remove2,
                    resdfname="resdf", 
                    outdir = outdir, name="LFC_Comparison_BeckDMSStress_OnlyDepr_pem3", w=8, h=8, 
                    scale_mode = "free")

cont001 <- compareLFCContratsNumeric(contrastlist2, firstContrast, 
                                     contrastNamesOrdered2, mainContrastName, 
                                     plim_select= 0.001, plim_plot=0.1,
                                     name2remove = name2remove2,
                                     resdfname="resdf", 
                                     outdir = outdir, name="LFC_Comparison_BeckDMSStress_OnlyDepr_pem3", w=8, h=8, 
                                     scale_mode = "free")
cont01 <- compareLFCContratsNumeric(contrastlist2, firstContrast, 
                                    contrastNamesOrdered2, mainContrastName, 
                                    plim_select= 0.01, plim_plot=0.1,
                                    name2remove = name2remove2,
                                    resdfname="resdf", 
                                    outdir = outdir, name="LFC_Comparison_BeckDMSStress_OnlyDepr_pem2", w=8, h=8, 
                                    scale_mode = "free")
cont05 <- compareLFCContratsNumeric(contrastlist2, firstContrast, 
                                    contrastNamesOrdered2, mainContrastName, 
                                    plim_select= 0.05, plim_plot=0.1,
                                    name2remove = name2remove2,
                                    resdfname="resdf", 
                                    outdir = outdir, name="LFC_Comparison_BeckDMSStress_OnlyDepr_05", w=8, h=8, 
                                    scale_mode = "free")
#Beck levels inside Depression
fnames <- list.files(paste0(opt$out, "/DESeq2_ByGroup/DeSEQ2/remove_tanda2_onlyDepression_Beck_cualitativo"), 
                     pattern = "_vs_", full.names = T ) %>% 
  subset(grepl("Normal.tsv", .))
contrastlist2 <- map(fnames, read_tsv) %>% lapply(\(x)return(list(resdf=x)))
newnames <- basename(fnames) %>% sapply(\(x)strsplit(x, "Beck_cualitativo_")[[1]][3]) %>% 
  gsub("_DAAshrinkNormal.tsv", "", .) %>% gsub("_", " ", .)
names(contrastlist2) <- newnames
name2remove2 <- "xxx"

contrastNamesOrdered2 <- c("Depression vs Control", names(contrastlist2))
compareLFCContrats(contrastlist2, firstContrast, 
                   contrastNamesOrdered2, mainContrastName, 
                   plim_select= 0.05, plim_plot=0.1,
                   name2remove = name2remove2,
                   resdfname="resdf", 
                   outdir = outdir, name="LFC_Comparison_BeckQual_OnlyDepr_p05", w=12, h=8, scale_mode = "fixed")
compareLFCContrats(contrastlist2, firstContrast, 
                   contrastNamesOrdered2, mainContrastName, 
                   plim_select= 0.001, plim_plot=0.1,
                   name2remove = name2remove2,
                   resdfname="resdf", 
                   outdir = outdir, name="LFC_Comparison_BeckQual_OnlyDepr_pem3", w=12, h=8, 
                   scale_mode = "fixed")


#Exercise levels inside Depression
fnames <- list.files(paste0(opt$out, "/DESeq2_ByGroup/DeSEQ2/remove_tanda2_onlyDepression_IPAQ_act_fisica/"), 
                     pattern = "_vs_", full.names = T ) %>% 
  subset(grepl("Normal.tsv", .))
contrastlist2 <- map(fnames, read_tsv) %>% lapply(\(x)return(list(resdf=x)))
newnames <- basename(fnames) %>% sapply(\(x)strsplit(x, "fisica_")[[1]][3]) %>% 
  gsub("_DAAshrinkNormal.tsv", "", .) %>% gsub("_", " ", .)
names(contrastlist2) <- newnames
name2remove2 <- "xxx"

contrastNamesOrdered2 <- c("Depression vs Control", names(contrastlist2))
compareLFCContrats(contrastlist2, firstContrast, 
                   contrastNamesOrdered2, mainContrastName, 
                   plim_select= 0.05, plim_plot=0.1,
                   name2remove = name2remove2,
                   resdfname="resdf", 
                   outdir = outdir, name="LFC_Comparison_IPAQ_OnlyDepr_p05", w=12, h=8, scale_mode = "fixed")
compareLFCContrats(contrastlist2, firstContrast, 
                   contrastNamesOrdered2, mainContrastName, 
                   plim_select= 0.001, plim_plot=0.1,
                   name2remove = name2remove2,
                   resdfname="resdf", 
                   outdir = outdir, name="LFC_Comparison_IPAQ_OnlyDepr_pem3", w=12, h=8, 
                   scale_mode = "fixed")
compareLFCContrats2(contrastlist2, firstContrast, 
                    contrastNamesOrdered2, mainContrastName, 
                    plim_select= 0.05, plim_plot=0.1,
                    name2remove = name2remove2,
                    resdfname="resdf", 
                    outdir = outdir, name="LFC_Comparison_IPAQ_OnlyDepr_p05", w=10, h=12, scale_mode = "free")
compareLFCContrats2(contrastlist2, firstContrast, 
                    contrastNamesOrdered2, mainContrastName, 
                    plim_select= 0.001, plim_plot=0.1,
                    name2remove = name2remove2,
                    resdfname="resdf", 
                    outdir = outdir, name="LFC_Comparison_IPAQ_OnlyDepr_pem3", w=10, h=8, 
                    scale_mode = "free")
names(contrastlist2) <- gsub(" ", "_", names(contrastlist2))
contrastNamesOrdered2[2:length(contrastNamesOrdered2)] <- gsub(" ", "_",contrastNamesOrdered2[2:length(contrastNamesOrdered2)])
cont001 <- compareLFCContratsNumeric(contrastlist2, firstContrast, 
                                     contrastNamesOrdered2, mainContrastName, 
                                     plim_select= 0.001, plim_plot=0.1,
                                     name2remove = name2remove2,
                                     resdfname="resdf", 
                                     outdir = outdir, name="LFC_Comparison_IPAQ_OnlyDepr_pem3", w=12, h=8, 
                                     scale_mode = "free")
cont01 <- compareLFCContratsNumeric(contrastlist2, firstContrast, 
                                    contrastNamesOrdered2, mainContrastName, 
                                    plim_select= 0.01, plim_plot=0.1,
                                    name2remove = name2remove2,
                                    resdfname="resdf", 
                                    outdir = outdir, name="LFC_Comparison_IPAQ_OnlyDepr_pem2", w=8, h=8, 
                                    scale_mode = "free")
cont05 <- compareLFCContratsNumeric(contrastlist2, firstContrast, 
                                    contrastNamesOrdered2, mainContrastName, 
                                    plim_select= 0.05, plim_plot=0.1,
                                    name2remove = name2remove2,
                                    resdfname="resdf", 
                                    outdir = outdir, name="LFC_Comparison_IPAQ_OnlyDepr_p05", w=8, h=12, 
                                    scale_mode = "free")
cont05$correlations %>% select(1,2,3,5,9,10) %>% mutate_if(is.numeric, round,3) %>%  datatable()

#Exercise levels adjusted by Condition vs Condition (unadjusted)
fnames <- list.files(paste0(opt$out, "/DESeq2_ControlVars/DeSEQ2/remove_tanda2_IPAQ_act_fisica/"), 
                     pattern = "_vs_", full.names = T ) %>% 
  subset(grepl("Normal.tsv", .)) %>% subset(!grepl("Depression_vs_Control", .))
contrastlist2 <- map(fnames, read_tsv) %>% lapply(\(x)return(list(resdf=x)))
newnames <- basename(fnames) %>% sapply(\(x)strsplit(x, "fisica_")[[1]][3]) %>% 
  gsub("_DAAshrinkNormal.tsv", "", .) %>% gsub("_", " ", .)
names(contrastlist2) <- newnames
name2remove2 <- "xxx"

contrastNamesOrdered2 <- c("Depression vs Control", names(contrastlist2))
compareLFCContrats(contrastlist2, firstContrast, 
                   contrastNamesOrdered2, mainContrastName, 
                   plim_select= 0.05, plim_plot=0.1,
                   name2remove = name2remove2,
                   resdfname="resdf", 
                   outdir = outdir, name="LFC_Comparison_IPAQ_CorrectedByDepr_p05", w=12, h=8, scale_mode = "fixed")
compareLFCContrats(contrastlist2, firstContrast, 
                   contrastNamesOrdered2, mainContrastName, 
                   plim_select= 0.001, plim_plot=0.1,
                   name2remove = name2remove2,
                   resdfname="resdf", 
                   outdir = outdir, name="LFC_Comparison_IPAQ_CorrectedByDepr_pem3", w=12, h=8, 
                   scale_mode = "fixed")
compareLFCContrats2(contrastlist2, firstContrast, 
                    contrastNamesOrdered2, mainContrastName, 
                    plim_select= 0.05, plim_plot=0.1,
                    name2remove = name2remove2,
                    resdfname="resdf", 
                    outdir = outdir, name="LFC_Comparison_IPAQ_CorrectedByDepr_p05", w=10, h=12, scale_mode = "free")
compareLFCContrats2(contrastlist2, firstContrast, 
                    contrastNamesOrdered2, mainContrastName, 
                    plim_select= 0.001, plim_plot=0.1,
                    name2remove = name2remove2,
                    resdfname="resdf", 
                    outdir = outdir, name="LFC_Comparison_CorrectedByDepr_pem3", w=10, h=8, 
                    scale_mode = "free")
names(contrastlist2) <- gsub(" ", "_", names(contrastlist2))
contrastNamesOrdered2[2:length(contrastNamesOrdered2)] <- gsub(" ", "_",contrastNamesOrdered2[2:length(contrastNamesOrdered2)])
cont001 <- compareLFCContratsNumeric(contrastlist2, firstContrast, 
                                     contrastNamesOrdered2, mainContrastName, 
                                     plim_select= 0.001, plim_plot=0.1,
                                     name2remove = name2remove2,
                                     resdfname="resdf", 
                                     outdir = outdir, name="LFC_Comparison_IPAQ_CorrectedByDepr_pem3", w=12, h=8, 
                                     scale_mode = "free")
cont01 <- compareLFCContratsNumeric(contrastlist2, firstContrast, 
                                    contrastNamesOrdered2, mainContrastName, 
                                    plim_select= 0.01, plim_plot=0.1,
                                    name2remove = name2remove2,
                                    resdfname="resdf", 
                                    outdir = outdir, name="LFC_Comparison_IPAQ_CorrectedByDepr_pem2", w=8, h=8, 
                                    scale_mode = "free")
cont05 <- compareLFCContratsNumeric(contrastlist2, firstContrast, 
                                    contrastNamesOrdered2, mainContrastName, 
                                    plim_select= 0.05, plim_plot=0.1,
                                    name2remove = name2remove2,
                                    resdfname="resdf", 
                                    outdir = outdir, name="LFC_Comparison_IPAQ_CorrectedByDepr_p05", w=8, h=12, 
                                    scale_mode = "free")
#cont05$correlations %>% select(1,2,3,5,9,10) %>% mutate_if(is.numeric, round,3) %>%  datatable()


## Beck inside Depr vs IPAQ inside depr
fnames <- list.files(paste0(opt$out, "/DESeq2_ByGroup/DeSEQ2/remove_tanda2_onlyDepression_IPAQ_act_fisica/"), 
                     pattern = "_vs_", full.names = T ) %>% 
  subset(grepl("Normal.tsv", .))
contrastlist2 <- map(fnames, read_tsv) %>% lapply(\(x)return(list(resdf=x)))
newnames <- basename(fnames) %>% sapply(\(x)strsplit(x, "fisica_")[[1]][3]) %>% 
  gsub("_DAAshrinkNormal.tsv", "", .) %>% gsub("_", " ", .)
names(contrastlist2) <- newnames
firstContrast2 <- list(resdf=daa_all_bygroup$remove_tanda2$Depression$Escala_depresión_Beck)
contrastNamesOrdered2 <- c("Beck", names(contrastlist2))
compareLFCContrats(contrastlist2, firstContrast2, 
                   contrastNamesOrdered2, mainContrastName, 
                   plim_select= 0.05, plim_plot=0.1,
                   name2remove = name2remove2,
                   resdfname="resdf", 
                   outdir = outdir, name="LFC_Comparison_BeckInside_IPAQ_OnlyDepr_p05", w=12, h=8, scale_mode = "free")
compareLFCContrats(contrastlist2, firstContrast, 
                   contrastNamesOrdered2, mainContrastName, 
                   plim_select= 0.001, plim_plot=0.1,
                   name2remove = name2remove2,
                   resdfname="resdf", 
                   outdir = outdir, name="LFC_Comparison_BeckInside_IPAQ_OnlyDepr_pem3", w=12, h=8, 
                   scale_mode = "free")


#Diet levels inside Depression
fnames <- list.files(paste0(opt$out, "/DESeq2_ByGroup/DeSEQ2/remove_tanda2_onlyControl_Indice_Alimentación_saludable//"), 
                     pattern = "_vs_", full.names = T ) %>% 
  subset(grepl("Normal.tsv", .))
contrastlist2 <- map(fnames, read_tsv) %>% lapply(\(x)return(list(resdf=x)))
newnames <- basename(fnames) %>% sapply(\(x)strsplit(x, "able_")[[1]][3]) %>% 
  gsub("_DAAshrinkNormal.tsv", "", .) %>% gsub("_", " ", .)
names(contrastlist2) <- newnames
name2remove2 <- "xxx"

contrastNamesOrdered2 <- c("Depression vs Control", names(contrastlist2))
compareLFCContrats(contrastlist2, firstContrast, 
                   contrastNamesOrdered2, mainContrastName, 
                   plim_select= 0.05, plim_plot=0.1,
                   name2remove = name2remove2,
                   resdfname="resdf", 
                   outdir = outdir, name="LFC_Comparison_Alim_OnlyDepr_p05", w=12, h=8, scale_mode = "fixed")
compareLFCContrats(contrastlist2, firstContrast, 
                   contrastNamesOrdered2, mainContrastName, 
                   plim_select= 0.001, plim_plot=0.1,
                   name2remove = name2remove2,
                   resdfname="resdf", 
                   outdir = outdir, name="LFC_Comparison_Alim_OnlyDepr_pem3", w=12, h=8, 
                   scale_mode = "fixed")

firstContrast2 <- list(resdf=daa_all_bygroup$remove_tanda2$Depression$Escala_depresión_Beck)
contrastNamesOrdered2 <- c("Beck", names(contrastlist2))
compareLFCContrats(contrastlist2, firstContrast2, 
                   contrastNamesOrdered2, mainContrastName, 
                   plim_select= 0.05, plim_plot=0.1,
                   name2remove = name2remove2,
                   resdfname="resdf", 
                   outdir = outdir, name="LFC_Comparison_BeckInside_Alim_OnlyDepr_p05", w=12, h=8, scale_mode = "free")
compareLFCContrats(contrastlist2, firstContrast, 
                   contrastNamesOrdered2, mainContrastName, 
                   plim_select= 0.001, plim_plot=0.1,
                   name2remove = name2remove2,
                   resdfname="resdf", 
                   outdir = outdir, name="LFC_Comparison_BeckInside_Alim_OnlyDepr_pem3", w=12, h=8, 
                   scale_mode = "free")


#Mediterranea
fnames <- list.files(paste0(opt$out, "/DESeq2_ByGroup/DeSEQ2/remove_tanda2_onlyControl_Mediterranean_diet_adherence//"), 
                     pattern = "_vs_", full.names = T ) %>% 
  subset(grepl("Normal.tsv", .))
contrastlist2 <- map(fnames, read_tsv) %>% lapply(\(x)return(list(resdf=x)))
newnames <- basename(fnames) %>% sapply(\(x)strsplit(x, "herence_")[[1]][3]) %>% 
  gsub("_DAAshrinkNormal.tsv", "", .) %>% gsub("_", " ", .)
names(contrastlist2) <- newnames
name2remove2 <- "xxx"

contrastNamesOrdered2 <- c("Depression vs Control", names(contrastlist2))
compareLFCContrats(contrastlist2, firstContrast, 
                   contrastNamesOrdered2, mainContrastName, 
                   plim_select= 0.05, plim_plot=0.1,
                   name2remove = name2remove2,
                   resdfname="resdf", 
                   outdir = outdir, name="LFC_Comparison_DMedit_OnlyDepr_p05", w=12, h=8, scale_mode = "fixed")
compareLFCContrats(contrastlist2, firstContrast, 
                   contrastNamesOrdered2, mainContrastName, 
                   plim_select= 0.001, plim_plot=0.1,
                   name2remove = name2remove2,
                   resdfname="resdf", 
                   outdir = outdir, name="LFC_Comparison_DMedit_OnlyDepr_pem3", w=12, h=8, 
                   scale_mode = "fixed")
compareLFCContrats2(contrastlist2, firstContrast, 
                   contrastNamesOrdered2, mainContrastName, 
                   plim_select= 0.05, plim_plot=0.1,
                   name2remove = name2remove2,
                   resdfname="resdf", 
                   outdir = outdir, name="LFC_Comparison_DMedit_OnlyDepr_p05", w=12, h=12, scale_mode = "fixed")
compareLFCContrats2(contrastlist2, firstContrast, 
                   contrastNamesOrdered2, mainContrastName, 
                   plim_select= 0.001, plim_plot=0.1,
                   name2remove = name2remove2,
                   resdfname="resdf", 
                   outdir = outdir, name="LFC_Comparison_DMedit_OnlyDepr_pem3", w=12, h=10, 
                   scale_mode = "fixed")

names(contrastlist2) <- gsub(" ", "_", names(contrastlist2))
contrastNamesOrdered2[2:length(contrastNamesOrdered2)] <- gsub(" ", "_",contrastNamesOrdered2[2:length(contrastNamesOrdered2)])
cont001 <- compareLFCContratsNumeric(contrastlist2, firstContrast, 
                                     contrastNamesOrdered2, mainContrastName, 
                                     plim_select= 0.001, plim_plot=0.1,
                                     name2remove = name2remove2,
                                     resdfname="resdf", 
                                     outdir = outdir, name="LFC_Comparison_DMedit_OnlyDepr_pem3", w=12, h=8, 
                                     scale_mode = "free")
cont01 <- compareLFCContratsNumeric(contrastlist2, firstContrast, 
                                    contrastNamesOrdered2, mainContrastName, 
                                    plim_select= 0.01, plim_plot=0.1,
                                    name2remove = name2remove2,
                                    resdfname="resdf", 
                                    outdir = outdir, name="LFC_Comparison_DMedit_OnlyDepr_pem2", w=8, h=8, 
                                    scale_mode = "free")
cont05 <- compareLFCContratsNumeric(contrastlist2, firstContrast, 
                                    contrastNamesOrdered2, mainContrastName, 
                                    plim_select= 0.05, plim_plot=0.1,
                                    name2remove = name2remove2,
                                    resdfname="resdf", 
                                    outdir = outdir, name="LFC_Comparison_DMedit_OnlyDepr_p05", w=8, h=12, 
                                    scale_mode = "free")

## Mediterranea ajustada por depresion
#Mediterranea
fnames <- list.files(paste0(opt$out, "/DESeq2_ControlVars/DeSEQ2/remove_tanda2_Mediterranean_diet_adherence"), 
                     pattern = "_vs_", full.names = T ) %>% 
  subset(grepl("Normal.tsv", .)) %>% subset(!grepl("Depression_vs_Control", .))
contrastlist2 <- map(fnames, read_tsv) %>% lapply(\(x)return(list(resdf=x)))
newnames <- basename(fnames) %>% sapply(\(x)strsplit(x, "herence_")[[1]][3]) %>% 
  gsub("_DAAshrinkNormal.tsv", "", .) %>% gsub("_", " ", .)
names(contrastlist2) <- newnames
name2remove2 <- "xxx"

contrastNamesOrdered2 <- c("Depression vs Control", names(contrastlist2))
compareLFCContrats(contrastlist2, firstContrast, 
                   contrastNamesOrdered2, mainContrastName, 
                   plim_select= 0.05, plim_plot=0.1,
                   name2remove = name2remove2,
                   resdfname="resdf", 
                   outdir = outdir, name="LFC_Comparison_DMedit_CorrectedByDepr_p05", w=12, h=8, scale_mode = "fixed")
compareLFCContrats(contrastlist2, firstContrast, 
                   contrastNamesOrdered2, mainContrastName, 
                   plim_select= 0.001, plim_plot=0.1,
                   name2remove = name2remove2,
                   resdfname="resdf", 
                   outdir = outdir, name="LFC_Comparison_DMedit_CorrectedByDepr_pem3", w=12, h=8, 
                   scale_mode = "fixed")
compareLFCContrats2(contrastlist2, firstContrast, 
                    contrastNamesOrdered2, mainContrastName, 
                    plim_select= 0.05, plim_plot=0.1,
                    name2remove = name2remove2,
                    resdfname="resdf", 
                    outdir = outdir, name="LFC_Comparison_DMedit_CorrectedByDepr_p05", w=12, h=12, scale_mode = "fixed")

compareLFCContrats2(contrastlist2, firstContrast, 
                    contrastNamesOrdered2, mainContrastName, 
                    plim_select= 0.01, plim_plot=0.01,
                    name2remove = name2remove2,
                    resdfname="resdf", 
                    outdir = outdir, name="LFC_Comparison_DMedit_CorrectedByDepr_p01", w=12, h=10, 
                    scale_mode = "fixed")

compareLFCContrats2(contrastlist2, firstContrast, 
                    contrastNamesOrdered2, mainContrastName, 
                    plim_select= 0.001, plim_plot=0.1,
                    name2remove = name2remove2,
                    resdfname="resdf", 
                    outdir = outdir, name="LFC_Comparison_DMedit_CorrectedByDepr_pem3", w=12, h=10, 
                    scale_mode = "fixed")

names(contrastlist2) <- gsub(" ", "_", names(contrastlist2))
contrastNamesOrdered2[2:length(contrastNamesOrdered2)] <- gsub(" ", "_",contrastNamesOrdered2[2:length(contrastNamesOrdered2)])
cont001 <- compareLFCContratsNumeric(contrastlist2, firstContrast, 
                                     contrastNamesOrdered2, mainContrastName, 
                                     plim_select= 0.001, plim_plot=0.1,
                                     name2remove = name2remove2,
                                     resdfname="resdf", 
                                     outdir = outdir, name="LFC_Comparison_CorrectedByDepr_OnlyDepr_pem3", w=12, h=8, 
                                     scale_mode = "free")
cont01 <- compareLFCContratsNumeric(contrastlist2, firstContrast, 
                                    contrastNamesOrdered2, mainContrastName, 
                                    plim_select= 0.01, plim_plot=0.1,
                                    name2remove = name2remove2,
                                    resdfname="resdf", 
                                    outdir = outdir, name="LFC_Comparison_CorrectedByDepr_OnlyDepr_pem2", w=8, h=8, 
                                    scale_mode = "free")
cont05 <- compareLFCContratsNumeric(contrastlist2, firstContrast, 
                                    contrastNamesOrdered2, mainContrastName, 
                                    plim_select= 0.05, plim_plot=0.1,
                                    name2remove = name2remove2,
                                    resdfname="resdf", 
                                    outdir = outdir, name="LFC_Comparison_CorrectedByDepr_OnlyDepr_p05", w=8, h=12, 
                                    scale_mode = "free")

## Mediterranea2 (2 clases) ajustada por depresion

fnames <- list.files(paste0(opt$out, "/DESeq2_ControlVars/DeSEQ2/remove_tanda2_Mediterranean_diet_adherence2"), 
                     pattern = "_vs_", full.names = T ) %>% 
  subset(grepl("Normal.tsv", .)) %>% subset(!grepl("Depression_vs_Control", .))
contrastlist2 <- map(fnames, read_tsv) %>% lapply(\(x)return(list(resdf=x)))
newnames <- basename(fnames) %>% sapply(\(x)strsplit(x, "herence2_")[[1]][3]) %>% 
  gsub("_DAAshrinkNormal.tsv", "", .) %>% gsub("_", " ", .)
names(contrastlist2) <- newnames
name2remove2 <- "xxx"

contrastNamesOrdered2 <- c("Depression vs Control", names(contrastlist2))
compareLFCContrats(contrastlist2, firstContrast, 
                   contrastNamesOrdered2, mainContrastName, 
                   plim_select= 0.05, plim_plot=0.1,
                   name2remove = name2remove2,
                   resdfname="resdf", 
                   outdir = outdir, name="LFC_Comparison_DMedit2_CorrectedByDepr_p05", w=12, h=8, scale_mode = "fixed")
compareLFCContrats(contrastlist2, firstContrast, 
                   contrastNamesOrdered2, mainContrastName, 
                   plim_select= 0.001, plim_plot=0.1,
                   name2remove = name2remove2,
                   resdfname="resdf", 
                   outdir = outdir, name="LFC_Comparison_DMedit2_CorrectedByDepr_pem3", w=12, h=8, 
                   scale_mode = "fixed")
compareLFCContrats2(contrastlist2, firstContrast, 
                    contrastNamesOrdered2, mainContrastName, 
                    plim_select= 0.05, plim_plot=0.1,
                    name2remove = name2remove2,
                    resdfname="resdf", 
                    outdir = outdir, name="LFC_Comparison_DMedit2_CorrectedByDepr_p05", w=12, h=12, scale_mode = "fixed")

compareLFCContrats2(contrastlist2, firstContrast, 
                    contrastNamesOrdered2, mainContrastName, 
                    plim_select= 0.01, plim_plot=0.01,
                    name2remove = name2remove2,
                    resdfname="resdf", 
                    outdir = outdir, name="LFC_Comparison_DMedit2_CorrectedByDepr_p01", w=12, h=10, 
                    scale_mode = "fixed")

compareLFCContrats2(contrastlist2, firstContrast, 
                    contrastNamesOrdered2, mainContrastName, 
                    plim_select= 0.001, plim_plot=0.1,
                    name2remove = name2remove2,
                    resdfname="resdf", 
                    outdir = outdir, name="LFC_Comparison_DMedit2_CorrectedByDepr_pem3", w=12, h=10, 
                    scale_mode = "fixed")

names(contrastlist2) <- gsub(" ", "_", names(contrastlist2))
contrastNamesOrdered2[2:length(contrastNamesOrdered2)] <- gsub(" ", "_",contrastNamesOrdered2[2:length(contrastNamesOrdered2)])
cont001 <- compareLFCContratsNumeric(contrastlist2, firstContrast, 
                                     contrastNamesOrdered2, mainContrastName, 
                                     plim_select= 0.001, plim_plot=0.1,
                                     name2remove = name2remove2,
                                     resdfname="resdf", 
                                     outdir = outdir, name="LFC_Comparison_CorrectedByDepr_OnlyDepr_pem3", w=12, h=8, 
                                     scale_mode = "free")
cont01 <- compareLFCContratsNumeric(contrastlist2, firstContrast, 
                                    contrastNamesOrdered2, mainContrastName, 
                                    plim_select= 0.01, plim_plot=0.1,
                                    name2remove = name2remove2,
                                    resdfname="resdf", 
                                    outdir = outdir, name="LFC_Comparison_CorrectedByDepr_OnlyDepr_pem2", w=8, h=8, 
                                    scale_mode = "free")
cont05 <- compareLFCContratsNumeric(contrastlist2, firstContrast, 
                                    contrastNamesOrdered2, mainContrastName, 
                                    plim_select= 0.05, plim_plot=0.1,
                                    name2remove = name2remove2,
                                    resdfname="resdf", 
                                    outdir = outdir, name="LFC_Comparison_CorrectedByDepr_OnlyDepr_p05", w=8, h=12, 
                                    scale_mode = "free")

## Escala de Beck dentro de Depr vs Mediterranea dentro de Depr
fnames <- list.files(paste0(opt$out, "/DESeq2_ByGroup/DeSEQ2/remove_tanda2_onlyControl_Mediterranean_diet_adherence//"), 
                     pattern = "_vs_", full.names = T ) %>% 
  subset(grepl("Normal.tsv", .))
firstContrast2 <- list(resdf=daa_all_bygroup$remove_tanda2$Depression$Escala_depresión_Beck)
contrastNamesOrdered2 <- c("Beck", names(contrastlist2))
compareLFCContrats(contrastlist2, firstContrast2, 
                   contrastNamesOrdered2, mainContrastName, 
                   plim_select= 0.05, plim_plot=0.1,
                   name2remove = name2remove2,
                   resdfname="resdf", 
                   outdir = outdir, name="LFC_Comparison_BeckInside_DMedit_OnlyDepr_p05", w=12, h=8, scale_mode = "free")
compareLFCContrats(contrastlist2, firstContrast, 
                   contrastNamesOrdered2, mainContrastName, 
                   plim_select= 0.001, plim_plot=0.1,
                   name2remove = name2remove2,
                   resdfname="resdf", 
                   outdir = outdir, name="LFC_Comparison_BeckInside_DMedit_OnlyDepr_pem3", w=12, h=8, 
                   scale_mode = "free")


# Mediterranea, Adjusted by Depression

#Euroqol 
fnames <- list.files(paste0(opt$out,"/DESeq2_ByGroup/DeSEQ2/remove_tanda2_onlyDepression_Euroqol/"), 
                     pattern = ".tsv", full.names = T ) %>% 
  subset(grepl("Normal.tsv", .))
contrastlist2 <- map(fnames, read_tsv) %>% lapply(\(x)return(list(resdf=x)))

names(contrastlist2) <- c("Euroqol")
name2remove2 <- "xxx"

contrastNamesOrdered2 <- c("Depression vs Control", names(contrastlist2))
compareLFCContrats(contrastlist2, firstContrast, 
                   contrastNamesOrdered2, mainContrastName, 
                   plim_select= 0.05, plim_plot=0.1,
                   name2remove = name2remove2,
                   resdfname="resdf", 
                   outdir = outdir, name="LFC_Comparison_Euroqol_OnlyDepr_p05", w=12, h=8, scale_mode = "fixed")
compareLFCContrats(contrastlist2, firstContrast, 
                   contrastNamesOrdered2, mainContrastName, 
                   plim_select= 0.001, plim_plot=0.1,
                   name2remove = name2remove2,
                   resdfname="resdf", 
                   outdir = outdir, name="LFC_Comparison_Euroqol_OnlyDepr_pem3", w=12, h=8, 
                   scale_mode = "fixed")

firstContrast2 <- list(resdf=daa_all_bygroup$remove_tanda2$Depression$Escala_depresión_Beck)
contrastNamesOrdered2 <- c("Beck", names(contrastlist2))
compareLFCContrats(contrastlist2, firstContrast2, 
                   contrastNamesOrdered2, mainContrastName, 
                   plim_select= 0.05, plim_plot=0.1,
                   name2remove = name2remove2,
                   resdfname="resdf", 
                   outdir = outdir, name="LFC_Comparison_BeckInside_Euroqol_OnlyDepr_p05", w=12, h=8, scale_mode = "free")
compareLFCContrats(contrastlist2, firstContrast, 
                   contrastNamesOrdered2, mainContrastName, 
                   plim_select= 0.001, plim_plot=0.1,
                   name2remove = name2remove2,
                   resdfname="resdf", 
                   outdir = outdir, name="LFC_Comparison_BeckInside_Euroqol_OnlyDepr_pem3", w=12, h=8, 
                   scale_mode = "free")


#DII 
fnames <- list.files(paste0(opt$out,"/DESeq2_ByGroup/DeSEQ2/remove_tanda2_onlyDepression_DII/"), 
                     pattern = ".tsv", full.names = T ) %>% 
  subset(grepl("Normal.tsv", .))
contrastlist2 <- map(fnames, read_tsv) %>% lapply(\(x)return(list(resdf=x)))

names(contrastlist2) <- c("DII")
name2remove2 <- "xxx"

contrastNamesOrdered2 <- c("Depression vs Control", names(contrastlist2))
compareLFCContrats(contrastlist2, firstContrast, 
                   contrastNamesOrdered2, mainContrastName, 
                   plim_select= 0.05, plim_plot=0.1,
                   name2remove = name2remove2,
                   resdfname="resdf", 
                   outdir = outdir, name="LFC_Comparison_DII_OnlyDepr_p05", w=12, h=8, scale_mode = "fixed")
compareLFCContrats(contrastlist2, firstContrast, 
                   contrastNamesOrdered2, mainContrastName, 
                   plim_select= 0.001, plim_plot=0.1,
                   name2remove = name2remove2,
                   resdfname="resdf", 
                   outdir = outdir, name="LFC_Comparison_DII_OnlyDepr_pem3", w=12, h=8, 
                   scale_mode = "fixed")

compareLFCContrats2(contrastlist2, firstContrast, 
                    contrastNamesOrdered2, mainContrastName, 
                    plim_select= 0.05, plim_plot=0.1,
                    name2remove = name2remove2,
                    resdfname="resdf", 
                    outdir = outdir, name="LFC_Comparison_DII_OnlyDepr_p05", w=12, h=12, scale_mode = "fixed")
compareLFCContrats2(contrastlist2, firstContrast, 
                    contrastNamesOrdered2, mainContrastName, 
                    plim_select= 0.001, plim_plot=0.1,
                    name2remove = name2remove2,
                    resdfname="resdf", 
                    outdir = outdir, name="LFC_Comparison_DII_OnlyDepr_pem3", w=12, h=10, 
                    scale_mode = "fixed")

names(contrastlist2) <- gsub(" ", "_", names(contrastlist2))
contrastNamesOrdered2[2:length(contrastNamesOrdered2)] <- gsub(" ", "_",contrastNamesOrdered2[2:length(contrastNamesOrdered2)])
cont001 <- compareLFCContratsNumeric(contrastlist2, firstContrast, 
                                     contrastNamesOrdered2, mainContrastName, 
                                     plim_select= 0.001, plim_plot=0.1,
                                     name2remove = name2remove2,
                                     resdfname="resdf", 
                                     outdir = outdir, name="LFC_Comparison_DII_OnlyDepr_pem3", w=12, h=8, 
                                     scale_mode = "free")
cont01 <- compareLFCContratsNumeric(contrastlist2, firstContrast, 
                                    contrastNamesOrdered2, mainContrastName, 
                                    plim_select= 0.01, plim_plot=0.1,
                                    name2remove = name2remove2,
                                    resdfname="resdf", 
                                    outdir = outdir, name="LFC_Comparison_DII_OnlyDepr_pem2", w=8, h=8, 
                                    scale_mode = "free")
cont05 <- compareLFCContratsNumeric(contrastlist2, firstContrast, 
                                    contrastNamesOrdered2, mainContrastName, 
                                    plim_select= 0.05, plim_plot=0.1,
                                    name2remove = name2remove2,
                                    resdfname="resdf", 
                                    outdir = outdir, name="LFC_Comparison_DII_OnlyDepr_p05", w=8, h=12, 
                                    scale_mode = "free")

firstContrast2 <- list(resdf=daa_all_bygroup$remove_tanda2$Depression$Escala_depresión_Beck)
contrastNamesOrdered2 <- c("Beck", names(contrastlist2))
compareLFCContrats(contrastlist2, firstContrast2, 
                   contrastNamesOrdered2, mainContrastName, 
                   plim_select= 0.05, plim_plot=0.1,
                   name2remove = name2remove2,
                   resdfname="resdf", 
                   outdir = outdir, name="LFC_Comparison_BeckInside_Euroqol_OnlyDepr_p05", w=12, h=8, scale_mode = "free")
compareLFCContrats(contrastlist2, firstContrast, 
                   contrastNamesOrdered2, mainContrastName, 
                   plim_select= 0.001, plim_plot=0.1,
                   name2remove = name2remove2,
                   resdfname="resdf", 
                   outdir = outdir, name="LFC_Comparison_BeckInside_Euroqol_OnlyDepr_pem3", w=12, h=8, 
                   scale_mode = "free")

## CorrectedByAge and BMI against Beck levels in depressed
firstContrast_edad <- list(resdf=daa_all_corrected$remove_tanda2$Edad_log)
firstContrast_IMC <- list(resdf=daa_all_corrected$remove_tanda2$IMC_log)

fnames <- list.files(paste0(opt$out,"/DESeq2_ByGroup/DeSEQ2/remove_tanda2_onlyDepression_Beck_cualitativo"), 
                     pattern = "_vs_", full.names = T ) %>% 
  subset(grepl("Normal.tsv", .))
contrastlist2 <- map(fnames, read_tsv) %>% lapply(\(x)return(list(resdf=x)))
newnames <- basename(fnames) %>% sapply(\(x)strsplit(x, "Beck_cualitativo_")[[1]][3]) %>% 
  gsub("_DAAshrinkNormal.tsv", "", .) %>% gsub("_", " ", .)
names(contrastlist2) <- newnames
name2remove2 <- "xxx"
contrastNamesOrdered2 <- c("Depression vs Control", names(contrastlist2))

compareLFCContrats(contrastlist2, firstContrast_edad, 
                   contrastNamesOrdered2, mainContrastName, 
                   plim_select= 0.05, plim_plot=0.1,
                   name2remove = name2remove2,
                   resdfname="resdf", 
                   outdir = outdir, name="LFC_Comparison_CorrectAge_BeckQual_OnlyDepr_p05", w=12, h=8, scale_mode = "fixed")
compareLFCContrats(contrastlist2, firstContrast_edad,
                   contrastNamesOrdered2, mainContrastName,
                   plim_select= 0.001, plim_plot=0.1,
                   name2remove = name2remove2,
                   resdfname="resdf", 
                   outdir = outdir, name="LFC_Comparison_CorrectAge_BeckQual_OnlyDepr_pem3", w=12, h=8, 
                   scale_mode = "fixed")
compareLFCContrats(contrastlist2, firstContrast_IMC, 
                   contrastNamesOrdered2, mainContrastName, 
                   plim_select= 0.05, plim_plot=0.1,
                   name2remove = name2remove2,
                   resdfname="resdf", 
                   outdir = outdir, name="LFC_Comparison_CorrectIMC_BeckQual_OnlyDepr_p05", w=12, h=8, scale_mode = "fixed")
compareLFCContrats(contrastlist2, firstContrast_IMC, 
                   contrastNamesOrdered2, mainContrastName, 
                   plim_select= 0.001, plim_plot=0.1,
                   name2remove = name2remove2,
                   resdfname="resdf", 
                   outdir = outdir, name="LFC_Comparison_CorrectIMC_BeckQual_OnlyDepr_pem3", w=12, h=8, 
                   scale_mode = "fixed")

#Correceted by Age and IMC, Again with quantitative scales inside depr
contrastlist2 <- list(daa_all_bygroup$remove_tanda2$Depression$Escala_depresión_Beck,
                      daa_all_bygroup$remove_tanda2$Depression$Montgomery.Asberg,
                      daa_all_bygroup$remove_tanda2$Depression$Escala_Hamilton,
                      daa_all_bygroup$remove_tanda2$Depression$DMSV_puntuacion_total
) %>% lapply(\(x)return(list(resdf=x)))
name2remove2 <- "xxx"
names(contrastlist2) <- c("Beck", "Montgomery", "Hamiliton", "DMSV")
contrastNamesOrdered2 <- c("Depression vs Control", names(contrastlist2))
compareLFCContrats(contrastlist2, firstContrast_edad, 
                   contrastNamesOrdered2, mainContrastName, 
                   plim_select= 0.05, plim_plot=0.1,
                   name2remove = name2remove2,
                   resdfname="resdf", 
                   outdir = outdir, name="LFC_Comparison_CorrAge_AllNum_OnlyDepr_p05", w=12, h=8, scale_mode = "free")
compareLFCContrats(contrastlist2, firstContrast_edad, 
                   contrastNamesOrdered2, mainContrastName, 
                   plim_select= 0.001, plim_plot=0.1,
                   name2remove = name2remove2,
                   resdfname="resdf", 
                   outdir = outdir, name="LFC_Comparison_CorrAge_AllNum_OnlyDepr_pem3", w=12, h=8, 
                   scale_mode = "free")

compareLFCContrats(contrastlist2, firstContrast_IMC, 
                   contrastNamesOrdered2, mainContrastName, 
                   plim_select= 0.05, plim_plot=0.1,
                   name2remove = name2remove2,
                   resdfname="resdf", 
                   outdir = outdir, name="LFC_Comparison_CorrIMC_AllNum_OnlyDepr_p05", w=12, h=8, scale_mode = "free")
compareLFCContrats(contrastlist2, firstContrast_IMC, 
                   contrastNamesOrdered2, mainContrastName, 
                   plim_select= 0.001, plim_plot=0.1,
                   name2remove = name2remove2,
                   resdfname="resdf", 
                   outdir = outdir, name="LFC_Comparison_CorrIMC_AllNum_OnlyDepr_pem3", w=12, h=8, 
                   scale_mode = "free")

# Compare Original with corrected, population variables
contrastlist2 <- list(daa_all_corrected$remove_tanda2$Edad_log,
                      daa_all_corrected$remove_tanda2$Sexo,
                      daa_all_corrected$remove_tanda2$ob_o_sobrepeso,
                      daa_all_corrected$remove_tanda2$IMC_log
                      #daa_all_corrected$remove_tanda2$Colesterol,
                      #daa_all_corrected$remove_tanda2$TG,
                      #daa_all_corrected$remove_tanda2$Euroqol
) %>% lapply(\(x)return(list(resdf=x)))
name2remove2 <- "xxx"
names(contrastlist2) <- c("Age", "Sex", "Overweight", "BMI")#, "Cholest", "TG", "Euroqol" )
contrastNamesOrdered2 <- c("Depression vs Control", names(contrastlist2))
compareLFCContrats(contrastlist2, firstContrast, 
                   contrastNamesOrdered2, mainContrastName, 
                   plim_select= 0.05, plim_plot=0.1,
                   name2remove = name2remove2,
                   resdfname="resdf", 
                   outdir = outdir, name="LFC_Comparison_withCorrected_p05", w=12, h=8, scale_mode = "free")
compareLFCContrats(contrastlist2, firstContrast, 
                   contrastNamesOrdered2, mainContrastName, 
                   plim_select= 0.001, plim_plot=0.1,
                   name2remove = name2remove2,
                   resdfname="resdf", 
                   outdir = outdir, name="LFC_Comparison_withCorrected_pem3", w=12, h=8, 
                   scale_mode = "free")

##Original vs corrected - treatments
contrastlist2 <- list(daa_all_corrected$remove_tanda2$tratamiento_ansioliticos,
                      daa_all_corrected$remove_tanda2$tratamiento_anticonvulsivos,
                      daa_all_corrected$remove_tanda2$tratamiento_ISRNs,
                      daa_all_corrected$remove_tanda2$tratamiento_antidiabeticos,
                      daa_all_corrected$remove_tanda2$tratamiento_coagul_betabloq_etc
                      #daa_all_corrected$remove_tanda2$TG,
                      #daa_all_corrected$remove_tanda2$Euroqol
) %>% lapply(\(x)return(list(resdf=x)))
name2remove2 <- "xxx"
names(contrastlist2) <- c("Ansiol", "Anticonv", "ISRN", "Antidiab", "Betab")#, "Cholest", "TG", "Euroqol" )
contrastNamesOrdered2 <- c("Depression vs Control", names(contrastlist2))
compareLFCContrats(contrastlist2, firstContrast, 
                   contrastNamesOrdered2, mainContrastName, 
                   plim_select= 0.05, plim_plot=0.1,
                   name2remove = name2remove2,
                   resdfname="resdf", 
                   outdir = outdir, name="LFC_Comparison_withCorrectedTreat_p05", w=12, h=8, scale_mode = "free")
compareLFCContrats(contrastlist2, firstContrast, 
                   contrastNamesOrdered2, mainContrastName, 
                   plim_select= 0.001, plim_plot=0.1,
                   name2remove = name2remove2,
                   resdfname="resdf", 
                   outdir = outdir, name="LFC_Comparison_withCorrectedTreat_pem3", w=12, h=8, 
                   scale_mode = "free")

## With scales, correcting for covariates
opt$reserva_0 <- opt$out
opt$out <- paste0(opt$out, "DESeq2_MainVariableScales_Corrected/")
if(!dir.exists(opt$out)) dir.create(opt$out)
interest_vars <- escalas_quant
daa_all_corrected_scales <- list()
vars2test <- c("Edad_log", "IMC_log", "Sexo")
for(phname in phseq_to_correct){
  cat("Doing DESeq2 Analysys with correction for: ", phname, "\n")
  phobj <- all_phyloseq[[phname]]
  phobj <- updatePsWithLogs(phobj, c("Edad", "IMC"))
  daa_all_corrected_scales[[phname]] <- list()
  for(interestvar_tmp in interest_vars){
    daa_all_corrected_scales[[phname]][[interestvar_tmp]] <- list()
    samples <- sample_data(phobj)$sampleID[! is.na(sample_data(phobj)[, interestvar_tmp])]
    phobj_prefilt <- phyloseq::prune_samples(samples, phobj)
    for(var in vars2test){
      cat("Doing DESeq2 Analysys with correction for: ", phname, '-',interestvar_tmp, ',correcting for', var, "\n")
      samples <- sample_data(phobj_prefilt)$sampleID[! is.na(sample_data(phobj_prefilt)[, var])]
      phobj_filt <- phyloseq::prune_samples(samples, phobj_prefilt)
      cases <- sample_data(phobj_filt)[, var] %>% unlist %>%  table
      if(length(which(cases > 0)) < 2 ){cat("Only one level, skipping this variable");next}
    
      vars2deseq <- c(interestvar_tmp, var)
      name <- paste0(phname, '_', var)
      res_tmp <-deseq_full_pipeline(phobj_filt, name, vars2deseq, opt)
      daa_all_corrected_scales[[phname]][[interestvar_tmp]][[var]] <- res_tmp
  }}}

save(daa_all_corrected_scales, file = paste0(opt$out, "DeSEQ2/DESEQ2_MainVariableScales.RData"))
opt$out <- opt$reserva_0

# Integrate
## Corr edad cualitativo vs corr edad cuantitativo
outdir <- paste0(opt$out, "DESeq2_MainVariableScales_Corrected/IntegrateContrasts/")
if(!dir.exists(outdir)) dir.create(outdir)
contrastlist2 <- list(daa_all_corrected_scales$remove_tanda2$Escala_depresión_Beck$Edad_log$all_contrasts[[1]],
                      daa_all_corrected_scales$remove_tanda2$Montgomery.Asberg$Edad_log$all_contrasts[[1]],
                      daa_all_corrected_scales$remove_tanda2$Escala_Hamilton$Edad_log$all_contrasts[[1]],
                      daa_all_corrected_scales$remove_tanda2$DMSV_puntuacion_total$Edad_log$all_contrasts[[1]]
)
name2remove2 <- "xxx"
names(contrastlist2) <- c("Beck", "Montgomery", "Hamiliton", "DMSV")
contrastNamesOrdered2 <- c("Depression vs Control", names(contrastlist2))
compareLFCContrats(contrastlist2, firstContrast_edad, 
                   contrastNamesOrdered2, mainContrastName, 
                   plim_select= 0.05, plim_plot=0.1,
                   name2remove = name2remove2,
                   resdfname="resdf", 
                   outdir = outdir, name="LFC_Comparison_CorrectEdad_AllNumCorrected_p05", w=12, h=8, scale_mode = "free")
compareLFCContrats(contrastlist2, firstContrast_edad, 
                   contrastNamesOrdered2, mainContrastName, 
                   plim_select= 0.001, plim_plot=0.1,
                   name2remove = name2remove2,
                   resdfname="resdf", 
                   outdir = outdir, name="LFC_Comparison_CorrectEdad_AllNumCorrected_pem3", w=12, h=8, 
                   scale_mode = "free")


## Corr IMC cualitativo vs corr IMC cuantitativo

contrastlist2 <- list(daa_all_corrected_scales$remove_tanda2$Escala_depresión_Beck$IMC_log$all_contrasts[[1]],
                      daa_all_corrected_scales$remove_tanda2$Montgomery.Asberg$IMC_log$all_contrasts[[1]],
                      daa_all_corrected_scales$remove_tanda2$Escala_Hamilton$IMC_log$all_contrasts[[1]],
                      daa_all_corrected_scales$remove_tanda2$DMSV_puntuacion_total$IMC_log$all_contrasts[[1]]
)
name2remove2 <- "xxx"
names(contrastlist2) <- c("Beck", "Montgomery", "Hamiliton", "DMSV")
contrastNamesOrdered2 <- c("Depression vs Control", names(contrastlist2))
compareLFCContrats(contrastlist2, firstContrast_IMC, 
                   contrastNamesOrdered2, mainContrastName, 
                   plim_select= 0.05, plim_plot=0.1,
                   name2remove = name2remove2,
                   resdfname="resdf", 
                   outdir = outdir, name="LFC_Comparison_CorrectIMC_AllNumCorrected_p05", w=12, h=8, scale_mode = "free")
compareLFCContrats(contrastlist2, firstContrast_IMC, 
                   contrastNamesOrdered2, mainContrastName, 
                   plim_select= 0.001, plim_plot=0.1,
                   name2remove = name2remove2,
                   resdfname="resdf", 
                   outdir = outdir, name="LFC_Comparison_CorrectIMC_AllNumCorrected_pem3", w=12, h=8, 
                   scale_mode = "free")


#####################################
## Predict DEPR + OBESITY

editObesity <- function(vec){
  return(ifelse(is.na(vec), NA, ifelse(vec == "Normopeso", "normal weight", "overweight")))
}
combineClasses <- function(phobj, v1, v2, newname, sep=":"){
  v1vec <- sample_data(phobj)[[v1]] %>% gsub("Depression", "D", .) %>% gsub("Control", "C", .)
  sample_data(phobj)[[newname]] <- paste(v1vec, sample_data(phobj)[[v2]], sep=sep)
  sample_data(phobj)[[newname]][is.na(sample_data(phobj)[[v1]]) | is.na(sample_data(phobj)[[v2]])] <- NA
  
  return(phobj)
}


interestvar <- "Condition"
var2add <- "obesidad"
newname <- "Depr_and_Ob"
vars2pca <- c(interestvar, var2add, newname)
phseq_to_use <- c("remove_tanda2")

opt$out <- paste0(opt$out, "PredictDAA_multiclass")
if(!dir.exists(opt$out)) dir.create(opt$out)
opt <- restaurar(opt)

all_model_multi_results <- list()
for(i in phseq_to_use){
  cat("Doing Predictive models for: ", i, "\n")
  all_model_results[[i]] <- list()
  phobj <- all_phyloseq[[i]]
  
  sample_data(phobj)[[var2add]] <- editObesity(sample_data(phobj)[[var2add]])
  phobj <- combineClasses(phobj, interestvar, var2add, newname)
  samples <- sample_data(phobj)$sampleID[! is.na(sample_data(phobj)[, var2add])]
  phobj_filt <- phyloseq::prune_samples(samples, phobj)
  
  outdir <- paste0(opt$out, "PredictDAA_multiclass/", i, "/")
  opt$reserva <- opt$out
  opt$out <- outdir
  if(!dir.exists(opt$out)) dir.create(opt$out)
  
  ## Get data for PCA 
  df2pca <- if(is.null(daa_all[[i]]$vst_counts_df)){ daa_all[[i]]$norm_counts_df}else{ daa_all[[i]]$vst_counts_df }
  
  ## Get taxa 
  taxa_padj <- daa_all[[i]]$resdf %>% filter_taxa_padj
  taxa_praw <- daa_all[[i]]$resdf %>% filter_taxa_praw
  taxa_cov_padj <- daa_all_corrected_only[[i]]$IMC_log %>% filter_taxa_padj
  taxa_cov_praw <- daa_all[[i]]$resdf %>% filter_taxa_praw
  taxa_padj01 <- daa_all[[i]]$resdf %>% filter_taxa_padj(plim=0.01)
  taxa_padj001 <- daa_all[[i]]$resdf %>% filter_taxa_padj(plim=0.001)
  taxa_cov_padj01 <- daa_all_corrected_only[[i]]$IMC_log %>% filter_taxa_padj(plim=0.01)
  taxa_cov_padj001 <- daa_all_corrected_only[[i]]$IMC_log %>% filter_taxa_padj(plim=0.001)
  
  taxa_padj_corr <- daa_all_corrected[[i]]$IMC_log %>% filter_taxa_padj 
  taxa_padj_corr_01 <- daa_all_corrected[[i]]$IMC_log %>% filter_taxa_padj(plim=0.01)
  taxa_padj_corr_001 <- daa_all_corrected[[i]]$IMC_log %>% filter_taxa_padj(plim=0.001)
  imcfname <- paste0(opt$reserva, "DESeq2_ControlVars/DeSEQ2/", i, "_IMC_log/",i, "_IMC_log_IMC_log_DAAshrinkNormal.tsv" )
  imctab <- read_tsv(imcfname)
  taxa_padj_cov_corr <- imctab %>% filter_taxa_padj
  taxa_padj_cov_corr_01 <- imctab %>% filter_taxa_padj(plim = 0.01)
  taxa_padj_cov_corr_001 <- imctab %>% filter_taxa_padj(plim = 0.001)
  
  taxa_padj_05 <- unique(c(taxa_padj, taxa_cov_padj))
  taxa_praw_05 <- unique(c(taxa_praw, taxa_cov_praw))
  taxa_padj_01 <- unique(c(taxa_padj01, taxa_cov_padj01))
  taxa_padj_001 <- unique(c(taxa_cov_padj001, taxa_padj001))
  taxa_padj_cov_05 <- unique(c(taxa_padj_corr, taxa_padj_cov_corr))
  taxa_padj_cov_01 <- unique(c(taxa_padj_corr_01, taxa_padj_cov_corr_01))
  taxa_padj_cov_001 <- unique(c(taxa_padj_corr_001, taxa_padj_cov_corr_001))
  
  taxa_padj_covorsingle_05 <- unique(c(taxa_padj_cov_05, taxa_padj_05))
  taxa_padj_covorsingle_01 <- unique(c(taxa_padj_01, taxa_padj_cov_01))
  taxa_padj_covorsingle_001 <- unique(c(taxa_padj_001, taxa_padj_cov_001))
  ## Make PCAs
  
  taxa_list <- list(
    "DiffTaxaPadjBase" = taxa_padj,
    "DiffTaxaPraw" = taxa_praw_05,
    "DiffTaxaPadj" = taxa_padj_05,
    "DiffTaxaPadj01" = taxa_padj_01,
    "DiffTaxaPadj001" = taxa_padj_001,
    "DiffTaxaPadjCov05" = taxa_padj_cov_05,
    "DiffTaxaPadjCov01" = taxa_padj_cov_01,
    "DiffTaxaPadjCov001" = taxa_padj_cov_001,
    "DiffTaxaPadjCovOrSingle05" = taxa_padj_covorsingle_05,
    "DiffTaxaPadjCovOrSingle01" = taxa_padj_covorsingle_01,
    "DiffTaxaPadjCovOrSingle001" = taxa_padj_covorsingle_001
  )
  all_pcalists <- map(names(taxa_list), 
                      \(x)makeAllPCAs(phobj_filt, df2pca, taxa_list[[x]], vars2pca, opt, x))
  names(all_pcalists) <- paste0("PCA_", names(taxa_list))
  
  all_models_this <- map2(all_pcalists, names(taxa_list),
                                            \(x, y)callDoAllModelsFromALLPCAs(x, name=paste0(i, y), 
                                                                              metadata=this_metadata, 
                                                                              vars2pca=c("Depr_and_Ob")))
  names(all_models_this) <- names(taxa_list)
  names(taxa_list) <- paste0("taxa_", names(taxa_list))
  
  
  all_model_results[[i]] <- c(all_models_this, all_pcalists, taxa_list)
  all_model_results[[i]]$metadata <- sample_data(phobj_filt) %>% data.frame
  
  opt <- restaurar(opt)
}
opt <- restaurar(opt)
save(all_model_results, file=paste0(opt$out, "PredictDAA_multiclass/all_model_results.RData"))
#load(file=paste0(opt$out, "PredictDAA/all_model_results.RData"))

#for(k in names(all_model_results[[i]]$modresults)) {cat("\n\n", k);all_model_results[[i]]$modresults[[k]]$modummary %>% head(2) %>% print}

#Integrate
opt$out <- paste0(opt$out, "PredictDAA_multiclass/")
makeLinePlotComparingPhobjs( all_model_results, opt, models_name1 = "DiffTaxaPadjBase", models_name2 = "DiffTaxaPadjCov05")
## Compare with Bacteria in componets
condnames <- names(all_model_results$remove_tanda2)[c(1,6,7,9,10)]
makeLinePlotComparingSamePhobjModels_Cov(phname, condnames, all_model_results, "TaxaGroups_trim", opt, w=8, h=5)

walk(names(all_model_results), \(x)makeLinePlotComparingSamePhobjModels(x, all_model_results, 
     opt, w=6, h=8, get_pcnames_from = "DiffTaxaPadjCov05"))
## Plot boxplot PCs
pca2use <- "DiffTaxaPadjCov05"
pcBoxplots <- map(names(all_model_results), \(x){
  makePCsBoxplot(x, all_model_results, opt, 
                 get_pcnames_from = pca2use,
                 pca_name =  paste0("PCA_", pca2use), 
                 varname = "Depr_and_Ob",
                 w=12, h=8)
  })
names(pcBoxplots) <- names(all_model_results)

## Plot barplot PCs and LFC
pcBarplots <- map(names(all_model_results), makePCBarplot, all_model_results, pcBoxplots, daa_all, opt,
                  get_pcnames_from = pca2use,
                  pca_name =  paste0("PCA_", pca2use), 
                  varname = "Depr_and_Ob",
                  w=12, h=20)
names(pcBarplots) <- names(all_model_results)

# Plot KNN (best model)

phname <- "remove_tanda2"
predplots <- map(names(all_model_results), plotAllModelPredictions, all_model_results, opt,
                 get_pcnames_from = pca2use,
                 pca_name =  paste0("PCA_", pca2use), 
                 varname = "Depr_and_Ob",
                 pred_mode="l1o")


opt <- restaurar(opt)

#Contingency table
library(ggmosaic)
df <- sample_data(phobj) %>% data.frame %>% 
  dplyr::filter(!is.na(IMC)) %>% 
  dplyr::mutate(overweight = ifelse(IMC>25, "overweight", "normal weight")) %>% 
  dplyr::mutate(Depr_ow = paste(Condition, overweight, sep=":"))
var1 <- "Condition"
var2 <- "overweight"

makeContingencyPlot(df, var1, var2, opt$out, "ObesityVsDepression", w=6, h=5)

