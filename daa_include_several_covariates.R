
opt <- restaurar(opt) 

phseq_to_correct <- names(all_phyloseq)[4]
interestvar <- "Condition"

groups2test <- list(
  c("Edad_log", "BMI_log"),
  c("Edad_log", "ob_o_sobrepeso"),
  c("Edad_log", "Sexo", "BMI_log"),
  c("Edad_log", "Sexo", "ob_o_sobrepeso"),
  c("Sexo", "Diabetes", "BMI_log"),
  c("Edad_log", "Sexo", "BMI_log", "Diabetes"),
  c("Edad_log", "Sexo", "BMI_log", "Diabetes", "Fumador"),
  c("tratamiento_ansioliticos", "tratamiento_anticonvulsivos", "tratamiento_ISRNs"),
  c("tratamiento_ansioliticos", "tratamiento_anticonvulsivos"),
  c("tratamiento_ansioliticos", "tratamiento_ISRNs"),
  c("tratamiento_ansioliticos", "tratamiento_anticonvulsivos",
    "tratamiento_ISRNs", "tratamiento_antidiabeticos", "tratamiento_coagul_betabloq_etc"),
  c("Edad_log", "Sexo", "BMI_log", "Diabetes",
    "tratamiento_ansioliticos", "tratamiento_anticonvulsivos")
)

groups2test <- groups2test[1:4]
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
  phobj <- updatePsWithLogs(phobj, c("Edad", "BMI"))
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
#load(paste0(opt$out, "DESeq2_ControlVarsMany/DESEQ2_controlVarsMany_all.RData"))
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