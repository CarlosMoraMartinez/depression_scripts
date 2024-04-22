
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


### ADONIS with multiple variables
permaformulas <- c(
  "braydist ~ Condition + Sexo",
  "braydist ~ Condition + BMI_log",
  "braydist ~ Condition + Edad_log",
  "braydist ~ Condition + IPAQ_act_fisica",
  "braydist ~ Condition + Mediterranean_diet_adherence2",
  "braydist ~ Condition + BMI_log + Sexo + Edad_log",
  "braydist ~ Condition + BMI_log + Sexo + Edad_log + IPAQ_act_fisica",
  "braydist ~ Condition + BMI_log + Sexo + Edad_log + Mediterranean_diet_adherence2",
  "braydist ~ Condition + BMI_log + Sexo + Edad_log + IPAQ_act_fisica + Mediterranean_diet_adherence2"
)

permaresults_mult <- list()
for(i in names(all_phyloseq)){
  cat("Doing PERMANOVA of: ", i)
  phobj <- all_phyloseq[[i]]
  phobj <- updatePsWithLogs(phobj, c("Edad", "BMI"))
  permaresults_mult[[i]] <- lapply(dtypes, FUN=function(dd, phobj, exclude_vars, SEED){
    oname <- paste0(outdir, "permanova_resultsMult_", i, "_", dd, ".tsv")
    makePermanovaFormulas(phobj,
                  permaformulas,
                  dist_method = dd, 
                  seed = SEED, 
                  outname = oname) 
  }, phobj, exclude_vars, SEED)
  names(permaresults_mult[[i]]) <- dtypes
}
save(permaresults_mult, file = paste0(outdir, "PERMANOVA_MULT.RData"))

mm <- permaresults_mult$remove_tanda2_rarefied_min$bray$modelos
