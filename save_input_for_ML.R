library(tidyverse)
library(janitor)

fooddir <- paste0(opt$out, "foodPCA/")
fname <- paste0(fooddir, "phyloseq_list_foodPCA.RData")
load(fname)

deadir <- paste0(opt$out, "DeSEQ2_age00/DeSEQ2/")
fname <- paste0(deadir, "DESEQ2_all_edad00.RData")
load(fname)

functdir <-  paste0(opt$out, "Functional_age/")
functdf <- read_tsv(paste0(functdir, "functTabInput_MetaCyc_filteredByProcess.tsv"))


outdir <- paste0(opt$out, "input_ML_prepared/")
if(!dir.exists(outdir)) dir.create(outdir)

abnmat <- function(df, genecolumn="gene", extraname=""){
  df <- df %>% column_to_rownames(genecolumn) %>% 
    as.matrix %>% t %>% 
    as.data.frame() %>% 
    rownames_to_column("sampleID") %>% 
    clean_names() %>% 
    dplyr::rename(sampleID = sample_id)
  if(extraname != ""){
    names(df)[2:ncol(df)] <- paste(extraname,  names(df)[2:ncol(df)], sep="__")
  }
  return(df)
}

ph2use <- names(all_phyloseq)[c(4, 5)]
for(nn in ph2use){
  meta <- sample_data(all_phyloseq[[nn]]) %>% data.frame()
  foodonly <- meta %>% select(sampleID, 5:30) %>% clean_names() %>% dplyr::rename(sampleID = sample_id)
  names(meta)[grep("PC", names(meta))] <- paste("food", names(meta)[grep("PC", names(meta))], sep="__")
  vstdf <- daa_all[[nn]]$vst_counts_df %>% abnmat(extraname = "vstabn")
  normdf <- daa_all[[nn]]$norm_counts_df %>% abnmat(extraname = "normabn")
  rawdf <- daa_all[[nn]]$raw_df %>% abnmat(extraname = "rawabn")
  funcfilt <- functdf %>% dplyr::select(Pathway, all_of(meta$sampleID)) %>% abnmat(genecolumn = "Pathway", extraname = "cpmfun")
  
  write_tsv(meta, file = paste0(outdir, nn, "_metadata.tsv"))
  write_tsv(foodonly, file = paste0(outdir, nn, "_foodonly.tsv"))
  write_tsv(vstdf, file = paste0(outdir, nn, "_vstabund.tsv"))
  write_tsv(normdf, file = paste0(outdir, nn, "_normabund.tsv"))
  write_tsv(rawdf, file = paste0(outdir, nn, "_rawabund.tsv"))
  write_tsv(funcfilt, file = paste0(outdir, nn, "_functabund.tsv"))
  
  dim(meta); dim(normdf); dim(vstdf); dim(rawdf); dim(funcfilt) 
  
  fulldf <- meta %>% 
    merge(vstdf, by.x="sampleID", by.y="sampleID", all.x=TRUE) %>% 
    merge(normdf, by.x="sampleID", by.y="sampleID", all.x=TRUE) %>% 
    merge(rawdf, by.x="sampleID", by.y="sampleID", all.x=TRUE) %>% 
    merge(funcfilt, by.x="sampleID", by.y="sampleID", all.x=TRUE)
    
  dim(fulldf)
  
  write_tsv(fulldf, file = paste0(outdir, nn, "_fulldf.tsv"))
  
  food_variables<- names(meta)[5:30]
  fndf <- data.frame(ind=5:30, variable=food_variables)
  write_tsv(fndf, paste0(outdir, nn, "_food_names.tsv"))
  
  bac_variables <- daa_all[[nn]]$raw_df %>% abnmat(extraname = "") %>% 
    select(-sampleID) %>% names
  bac_original_names <-  daa_all[[nn]]$raw_df %>% pull(gene)
  fndf <- data.frame(ind=1:length(bac_variables), variable=bac_variables,
                     original_name = bac_original_names)
  write_tsv(fndf, paste0(outdir, nn, "_tax_names.tsv"))
  
  fun_variables <- funcfilt %>% select(-sampleID) %>% names()
  fun_variables_original <- functdf %>% pull(Pathway)
  fndf <- data.frame(ind=1:length(fun_variables), variable=fun_variables,
                     original_name = fun_variables_original)
  write_tsv(fndf, paste0(outdir, nn, "_functional_names.tsv"))

}
