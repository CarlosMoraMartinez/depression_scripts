
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