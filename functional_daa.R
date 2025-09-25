library(tidyverse)
library(phyloseq)
library(G4Micro)
library(Maaslin2)
library(foreach)
library(doParallel)

opt <- list(fc=1,
            pval=0.05,
            minfreq = 0.05,
            num_genes_default=5,
            ptype="adjusted",
            fctype="shrunk",
            mincount= 1,
            minsampleswithcount= 6,
            out= "/home/carlos/Documentos/CORALS/results_rstudio/results_Agosto25_1/functional_1/",
            input=  "/home/carlos/Documentos/CORALS/results_cluster_240924/Humann3_analisis_funcional/MERGED_changed_names/",
            modules_path=  "/home/carlos/Escritorio//202311_DEPRESION/202311_DEPRESION/funcional1/kegg_modules/",
            cazy="/home/carlos/Escritorio//202311_DEPRESION/202311_DEPRESION/cazy/",
            input_funcional = "/home/carlos/Documentos/CORALS/results_cluster_240924/Humann3_analisis_funcional/MERGED_changed_names/",
            #dearesult = "/home/carmoma/Desktop/202311_DEPRESION/results_rstudio_v2_4/DeSEQ2/remove_tanda2/DESEQ2_all_results_remove_tanda2.R",
            #metadata =  "/home/carmoma/Desktop/202311_DEPRESION/metadatos_MC_AL12042023_CMcopy.xlsx",
            phyloseq_list = "/home/carlos/Documentos/CORALS/results_rstudio/results_Agosto25_1//foodPCA/phyloseq_list_foodPCA_withNMF_withInc.RData"
)

restaurar <- restauraropt_mk(opt)
opt <- restaurar(opt)
load(opt$phyloseq_list)

if(!dir.exists(opt$out)) dir.create(opt$out)


samples2keep <- sample_data(all_phyloseq$remove_tanda2) %>% data.frame %>% pull(sampleID)

metacyc_ab <-  readFunctionalMatrix(opt, "humann3_merged.renamed.tsv") %>% dplyr::select(all_of(c("Pathway",samples2keep)))
#metacyc_rxn <- readFunctionalMatrix(opt, "humann3_merged_genetables_RXN_CPM.renamed.tsv") %>% dplyr::select(all_of(c("Pathway",samples2keep)))
#ko <- readFunctionalMatrix(opt, "humann3_merged_genetables_KO_CPM.renamed.tsv") %>% dplyr::select(all_of(c("Pathway",samples2keep)))
#go <- readFunctionalMatrix(opt, "humann3_merged_genetables_GO_CPM.renamed.tsv") %>% dplyr::select(all_of(c("Pathway",samples2keep)))
#cazy <- readFunctionalMatrix(opt, "humann3_merged_genetables_CAZY_CPM.renamed.tsv") %>% dplyr::select(all_of(c("Pathway",samples2keep)))

tabs2test <- list("MetaCyc"=metacyc_ab) #, "MetaCyc Reactions"=metacyc_rxn, "KEGG"=ko, "Gene Ontology"=go, "CAZy" = cazy ) #"GOmixerModules"=modules_nosp2
walk2(tabs2test, paste0(opt$out, "/functTab2Use_",names(tabs2test), ".tsv"), \(df, name)if(!is.null(df))write_tsv(df, file=name))
save(tabs2test, file = paste0(opt$out, "functional_tabs.RData"))


tabs2annot <- list("MetaCyc"=NULL) #, "MetaCyc Reactions"=metacyc_mapnames, "KEGG"=ko_mapnames, "Gene Ontology"=go_mapnames, "CAZy" = cazy_mapnames, "GOmixerModules"=modules_mapnames_nosp )

walk2(tabs2test, paste0(opt$out, "/functTabInput_",names(tabs2test), ".tsv"), \(df, name)if(!is.null(df))write_tsv(df, file=name))


#### Limma
#met2use <- sample_data(all_phyloseq$filt) %>% data.frame %>%
#  dplyr::filter(status_c2 != "Insufficient gain") %>%
#  dplyr::filter(status_c2 != "initially_overweight")

s_meta <- sample_data(all_phyloseq$remove_tanda2) %>% data.frame

opt$minsampleswithcount <- opt$minfreq*nrow(s_meta)
filtered_samples <- lapply(tabs2test, \(x) x %>% dplyr::select(Pathway, all_of(s_meta$sampleID)))
filtered_byproc <- lapply(filtered_samples, FUN=filterWholeProcessesAndFreq, opt)
walk2(filtered_byproc, paste0(opt$out, "/functTabInput_",names(tabs2test), "_filteredByProcess.tsv"), \(df, name)if(!is.null(df))write_tsv(df, file=name))

filtered_byspecies <- lapply(tabs2test, FUN=filterBySpeciesAndFreq, opt)
walk2(filtered_byproc, paste0(opt$out, "/functTabInput_",names(tabs2test), "_filteredBySpecies.tsv"), \(df, name)if(!is.null(df))write_tsv(df, file=name))


vars2test1 <-  c("status_c1",
                 "z_t0", "z_t1", "inc_z_bmi",
                 "z_waist_00", "z_waist_01", "inc_z_waist",
                 "exercise_cat",  "mother_educ",
                 "status_c2", "nreads")

covars <- c("age_T0", "Sex")
covars_random <- c("hospital")

patnames <- c("Preprocessed",  "Mediterranean", "Western")
quant_vars_diet <- c(names(s_meta)[5], names(s_meta)[grepl("_clr", names(s_meta))])
diet_aggr_vars <- quant_vars_diet[1:5]

#check absent variants
vars2test1[!vars2test1 %in% names(s_meta)]

allvars <- c(vars2test1, patnames, quant_vars_diet)
proc <- filtered_byproc$MetaCyc %>% column_to_rownames("Pathway")

# for use with only normal at T0
s_meta_normT0 <- s_meta %>% filter(status_c1 == "normal")
proc_normT0 <- proc[, s_meta_normT0$sampleID]
allvars_normT0 <- allvars[allvars != "status_c1"]

######################################
#### Run with all data
ncores <- 4 #parallel::detectCores() - 1
cl <- makeCluster(ncores)
registerDoParallel(cl)
#for(var in allvars){
res <- foreach(var = allvars, .packages = c("tidyverse", "phyloseq", "G4Micro", "Maaslin2"),
               .export   = c("s_meta", "proc", "covars", "covars_random", "opt")) %dopar% {

    this_vars <- c(var, covars)
    var_class <- s_meta %>% dplyr::select(all_of(this_vars)) %>% map_vec(class)
    ref <- map_vec(this_vars[var_class %in% c("character", "factor")], \(x) paste(x, unique(s_meta[!is.na(s_meta[, x]), x])[1], sep=",", collapse=","))
    cat(ref, "\n")

    this_meta <- this_vars %>% reduce(~ .x %>% filter(!is.na(!!sym(.y))) , .init=s_meta)
    this_proc <- proc[, this_meta$sampleID]
    #tryCatch({
    fit_data = Maaslin2(input_data     = this_proc,
                    input_metadata = this_meta,
                    min_prevalence = 0,
                    normalization  = "NONE",
                    transform = "LOG",
                    output         = paste0(opt$out,  "Maaslin2_AllSexHosAge_", var),
                    fixed_effects  = c(var, covars),
                    random_effects = covars_random,
                    plot_scatter=FALSE,
                    max_pngs = 5,
                    reference      = ref)
    rm(fit_data)
    #}
   #, error = function(e) {
  #    message("Error in variable ", var, ": ", conditionMessage(e))
   #   return(NA)  # return NA or any placeholder value
    #})
}

stopCluster(cl)


####################################
##Now only normal at T0


cl <- makeCluster(ncores)
registerDoParallel(cl)
#for(var in allvars){
res <- foreach(var = allvars_normT0, .packages = c("tidyverse", "phyloseq", "G4Micro", "Maaslin2"),
               .export   = c("s_meta_normT0", "proc_normT0", "covars", "covars_random", "opt")) %dopar% {

                 this_vars <- c(var, covars)
                 var_class <- s_meta %>% dplyr::select(all_of(this_vars)) %>% map_vec(class)
                 ref <- map_vec(this_vars[var_class %in% c("character", "factor")], \(x) paste(x, unique(s_meta_normT0[!is.na(s_meta_normT0[, x]), x])[1], sep=",", collapse=","))
                 cat(ref, "\n")

                 this_meta <- this_vars %>% reduce(~ .x %>% filter(!is.na(!!sym(.y))) , .init=s_meta_normT0)
                 this_proc <- proc_normT0[, this_meta$sampleID]
                 #tryCatch({
                 fit_data = Maaslin2(input_data     = this_proc,
                                     input_metadata = this_meta,
                                     min_prevalence = 0,
                                     normalization  = "NONE",
                                     transform = "LOG",
                                     output         = paste0(opt$out,  "Maaslin2_NormT0SexHosAge_", var),
                                     fixed_effects  = c(var, covars),
                                     random_effects = covars_random,
                                     plot_scatter=FALSE,
                                     max_pngs = 5,
                                     reference      = ref)
                 rm(fit_data)
                 #}
                 #, error = function(e) {
                 #    message("Error in variable ", var, ": ", conditionMessage(e))
                 #   return(NA)  # return NA or any placeholder value
                 #})
               }

stopCluster(cl)

#### Run with all data by SPECIES
proc_sp <- filtered_byspecies$MetaCyc %>% column_to_rownames("Pathway")
proc_sp_normT0 <- proc_sp[, s_meta_normT0$sampleID]

ncores <- 6 #parallel::detectCores() - 1
cl <- makeCluster(ncores)
registerDoParallel(cl)
#for(var in allvars){
res <- foreach(var = allvars, .packages = c("tidyverse", "phyloseq", "G4Micro", "Maaslin2"),
               .export   = c("s_meta", "proc_sp", "covars", "covars_random", "opt")) %dopar% {

                 this_vars <- c(var, covars)
                 var_class <- s_meta %>% dplyr::select(all_of(this_vars)) %>% map_vec(class)
                 ref <- map_vec(this_vars[var_class %in% c("character", "factor")], \(x) paste(x, unique(s_meta[!is.na(s_meta[, x]), x])[1], sep=",", collapse=","))
                 cat(ref, "\n")

                 this_meta <- this_vars %>% reduce(~ .x %>% filter(!is.na(!!sym(.y))) , .init=s_meta)
                 this_proc <- proc_sp[, this_meta$sampleID]
                 #tryCatch({
                 fit_data = Maaslin2(input_data     = this_proc,
                                     input_metadata = this_meta,
                                     min_prevalence = 0,
                                     normalization  = "NONE",
                                     transform = "LOG",
                                     output         = paste0(opt$out,  "Maaslin2_BySpeciesAllSexHosAge_", var),
                                     fixed_effects  = c(var, covars),
                                     random_effects = covars_random,
                                     plot_scatter=FALSE,
                                     max_pngs = 5,
                                     reference      = ref)
                 rm(fit_data)
                 #}
                 #, error = function(e) {
                 #    message("Error in variable ", var, ": ", conditionMessage(e))
                 #   return(NA)  # return NA or any placeholder value
                 #})
               }

stopCluster(cl)

##Now only normal at T0


cl <- makeCluster(ncores)
registerDoParallel(cl)
#for(var in allvars){
res <- foreach(var = allvars_normT0, .packages = c("tidyverse", "phyloseq", "G4Micro", "Maaslin2"),
               .export   = c("s_meta_normT0", "proc_sp_normT0", "covars", "covars_random", "opt")) %dopar% {

                 this_vars <- c(var, covars)
                 var_class <- s_meta %>% dplyr::select(all_of(this_vars)) %>% map_vec(class)
                 ref <- map_vec(this_vars[var_class %in% c("character", "factor")], \(x) paste(x, unique(s_meta_normT0[!is.na(s_meta_normT0[, x]), x])[1], sep=",", collapse=","))
                 cat(ref, "\n")

                 this_meta <- this_vars %>% reduce(~ .x %>% filter(!is.na(!!sym(.y))) , .init=s_meta_normT0)
                 this_proc <- proc_normT0[, this_meta$sampleID]
                 #tryCatch({
                 fit_data = Maaslin2(input_data     = this_proc,
                                     input_metadata = this_meta,
                                     min_prevalence = 0,
                                     normalization  = "NONE",
                                     transform = "LOG",
                                     output         = paste0(opt$out,  "Maaslin2_NormT0SexHosAge_", var),
                                     fixed_effects  = c(var, covars),
                                     random_effects = covars_random,
                                     plot_scatter=FALSE,
                                     max_pngs = 5,
                                     reference      = ref)
                 rm(fit_data)
                 #}
                 #, error = function(e) {
                 #    message("Error in variable ", var, ": ", conditionMessage(e))
                 #   return(NA)  # return NA or any placeholder value
                 #})
               }

stopCluster(cl)





####################################################################################
dirlist <- list.dirs(opt$out, recursive = FALSE)
dirlist_all <- dirlist[grep("_AllSexHosAge", dirlist)]
dirlist_t0norm <- dirlist[grep("NormT0SexHosAge", dirlist)]

funcres_all <- map(allvars, \(x){
  dd <- dirlist_all[grep(x, dirlist_all)]
  if(length(dd) > 1) {stringr::str_glue("Error! more than 1 dir matching {x}: {dd}"): return(data.frame())}
  df <- read_tsv(paste0(dd, "/all_results.tsv")) %>%
    dplyr::mutate(variable = x, dataset="all")
  return(df)

}) %>% bind_rows()

funcres_all <- map(allvars_normT0, \(x){
  dd <- dirlist_t0norm[grep(x, dirlist_t0norm)]
  if(length(dd) > 1) {stringr::str_glue("Error! more than 1 dir matching {x}: {dd}"): return(data.frame())}
  df <- read_tsv(paste0(dd, "/all_results.tsv")) %>%
    dplyr::mutate(variable = x, dataset="Normal T0")
  return(df)

}) %>% bind_rows() %>% rbind(funcres_all)

funcres_all <- funcres_all %>%
  dplyr::mutate(path_code = map_vec(feature, \(x) strsplit(x, "\\.\\.")[[1]][1])) %>%
  dplyr::mutate(path_code = gsub("\\.", "-", path_code))

funcres_all %>% dim
funcres_all$dataset %>% table

write_tsv(funcres_all, file = paste0(opt$out, "Maaslin2_allMerged.tsv"))

funcres_all <- read_tsv(paste0(opt$out, "Maaslin2_allMerged.tsv"))
proc2use <- funcres_all %>% filter(dataset=="all") %>% filter(qval < 0.01) %>% pull(path_code) %>% unique


mat1 <- funcres_all %>% filter(path_code %in% proc2use) %>%
  dplyr::mutate(var_full = paste(variable, value, sep=":")) %>%
  filter(! value %in% c("age_T0", "Girl") ) %>%
  filter(dataset=="all") %>%
  #dplyr::mutate(coef_sig = ifelse(qval <=0.05 & !is.na(coef), coef, 0)) %>%
  dplyr::mutate(coef_sig = coef) %>%
  dplyr::select(value, path_code, coef_sig) %>%
  pivot_wider(names_from = path_code, values_from = coef_sig)

mat2<- mat1 %>% column_to_rownames("value") %>% as.matrix %>%
  t %>% scale
mat2[is.infinite(mat2)] <- 0
mat2[is.na(mat2)] <- 0

library(pheatmap)

pheatmap(mat2)


######## plot lolipop
clean_names <- function(tax){
  gsub("_", " ", tax) %>%
    gsub("[\\[\\]]", "", .)
}
makeLoliplotFromMaaslin2 <- function(daa_df,
                         vars2loliplot = c(), # Variables which LFC include in plot (all in column named 'Contrast')
                         vars2sort = c(), # Variables to sort taxa, ordered
                         vars2filter = c(), # Variables to select significant taxa
                         outdir = "./",
                         name="test",
                         plim = 0.01, plim_col = 0.05, lfclim = 0,
                         strip_fontsize=9,
                         w=12, h=12,
                         contrast_name = "Contrast",
                         pval_name =  "padj", # "qval"
                         lfc_name = "log2FoldChangeShrink", #"coef"
                         feature_name = "taxon", # feature
                         axis_face = "plain",
                         add_lines = TRUE
){

  if(length(vars2sort) == 0) vars2sort <- vars2loliplot
  if(length(vars2filter) == 0) vars2filter <- vars2loliplot

  usedf <- daa_df %>%
    dplyr::filter(!!sym(contrast_name) %in% vars2loliplot) %>%
    dplyr::mutate(
      Sig = ifelse(!is.na(!!sym(pval_name)) & !!sym(pval_name) <= plim_col, ifelse(!!sym(lfc_name) < 0 , "Down", "Up"), "NS"),
      !!sym(contrast_name) := factor(!!sym(contrast_name), levels=vars2loliplot)
    )

  tax2use <- daa_df %>%
    dplyr::filter(!!sym(contrast_name) %in% vars2filter) %>%
    dplyr::filter(!is.na(!!sym(pval_name)) & !!sym(pval_name) <= plim & abs(!!sym(lfc_name)) > lfclim) %>%
    dplyr::pull(!!sym(feature_name)) %>% unique
  tax2use %>% length

  taxorder <- usedf %>%
    dplyr::filter(!!sym(feature_name) %in% tax2use) %>%
    dplyr::select(all_of(c(feature_name, contrast_name, lfc_name, pval_name))) %>%
    tidyr::gather("vart", "valt", !!sym((lfc_name)), !!sym(pval_name)) %>%
    unite("vart2", !!sym(contrast_name), vart, sep="__") %>%
    tidyr::spread(vart2, valt) %>%
    dplyr::mutate(bmisig = "NS")

  taxorder_b <- purrr::reduce(vars2sort, ~ .y %>% dplyr::mutate(bmisig := ifelse(!!sym(.x) < plim_col,
                                                                                 .x, bmisig )),
                              .init= taxorder, .dir="backward")

  newlevs <- if("NS" %in% taxorder_b$bmisig){c("NS", rev(vars2sort))}else{ c("NS", rev(vars2sort))}

  taxorder_b <- taxorder_b %>%
    dplyr::mutate(bmisig = factor(bmisig, levels = newlevs)) %>%
    dplyr::group_by(bmisig) %>%
    dplyr::arrange(bmisig) %>%
    group_split()

  sortlevs <- purrr::map_vec(taxorder_b, ~ unique(.x[["bmisig"]])) %>% as.character
  if("NS" %in% sortlevs){
    sortlevs[sortlevs=="NS"] <- sortlevs[length(sortlevs)]
  }
  taxorder_b <- purrr::map2(taxorder_b, sortlevs, ~ .x %>%
                              dplyr::arrange( desc(.data[[gsub(glue::glue("__{pval_name}"), glue::glue("__{lfc_name}"), .y)]]) )) %>%
    bind_rows() %>%
    dplyr::mutate(!!sym(feature_name) :=clean_names(!!sym(feature_name))) %>%
    dplyr::mutate(!!sym(feature_name) :=factor(!!sym(feature_name), levels=!!sym(feature_name)))

  TAXORDER <- taxorder_b %>% pull(!!sym(feature_name))

  usedf2 <- usedf %>%
    filter(!!sym(feature_name) %in% tax2use) %>%
    dplyr::mutate(!!sym(feature_name) := clean_names(!!sym(feature_name))) %>%
    dplyr::mutate(!!sym(feature_name) := factor(!!sym(feature_name), levels=TAXORDER))

  linedf <- taxorder_b %>%
    dplyr::group_by(bmisig) %>%
    dplyr::summarise(xpos = max(as.numeric(!!sym(feature_name))) + 0.5) %>%
    head(nrow(.)-1)

  write_tsv(usedf2, file = paste0(outdir, "loliplot_", "_", as.character(plim), "_fc", as.character(lfclim), "_", name, ".tsv"))
  write_tsv(linedf, file = paste0(outdir, "loliplot_", "_", as.character(plim), "_fc", as.character(lfclim), "_", name, "_linedf.tsv"))

  g0 <- ggplot(usedf2, aes(x=!!sym(feature_name), y=!!sym(lfc_name), col=Sig, fill=Sig))+
      facet_grid(as.formula(glue::glue("~ {contrast_name}"))) +
      geom_hline(yintercept = 0, linetype=2, col="lightgray")

  if(add_lines) g0 <- g0 + geom_vline(data=linedf, aes(xintercept=xpos), linetype=2, col="lightgray")

  g0 <- g0 + geom_segment(aes(x=!!sym(feature_name), xend = !!sym(feature_name), y=0, yend=!!sym(lfc_name))) +
      geom_point() +
      theme_bw() +
      coord_flip() +
      scale_color_manual(values = c("Down"="steelblue", "Up"="tomato", "NS"="darkgray")) +
      theme(axis.text.y= element_text(face=axis_face, size=10),
            axis.text.x= element_text( size=12),
            axis.title = element_text(size=12),
            strip.text = element_text(size=strip_fontsize))

  ggsave(filename = paste0(outdir, "loliplot_", "_", as.character(plim), "_fc", as.character(lfclim), "_", name, ".pdf"), g0,
         width = w, height = h)
  result <- list(
    usedf=usedf2,
    linedf=linedf,
    plot=g0,
    plotname=paste0(outdir, "loliplot_", "_", as.character(plim), "_fc", as.character(lfclim), "_", name, ".pdf"),
    dfname= paste0(outdir, "loliplot_", "_", as.character(plim), "_fc", as.character(lfclim), "_", name, ".tsv"),
    dfname_lines=paste0(outdir, "loliplot_", "_", as.character(plim), "_fc", as.character(lfclim), "_", name, "_linedf.tsv"),
    input=daa_df,
    args = list(vars2sort = vars2sort,
                vars2filter = vars2filter,
                name=name,
                plim = plim, plim_col = plim_col, lfclim = lfclim,
                w=w, h=h)
  )
  save(result, file = paste0(outdir, "loliplot_", "_", as.character(plim), "_fc", as.character(lfclim), "_", name, ".RData"))
  return(result)

}

makeReplacementsDF <- function(df, replace_list, totitlecase=c()){
  dfmod <- purrr::reduce(replace_list, ~ .x %>%
                           dplyr::mutate(!!sym(.y[3]) := gsub(.y[1], .y[2], !!sym(.y[3]), perl=TRUE)),
                         .init = df) %>%
    dplyr::mutate(across(all_of(totitlecase), tools::toTitleCase))

  return(dfmod)
}


replace_strings <- list(
  c("_clr$", ""),
  c("_g$", " (g)"),
  c("_kcal", " (Kcal)"),
  c("z_t1", "Z-score BMI T1"),
  c("z_waist_01", "Z-score waist T1"),
  c("z_t0", "Z-score BMI T0"),
  c("z_t1", "Z-score BMI T1"),
  c("inc_z_waist", "Change in Z-score waist"),
  c("inc_z_bmi", "Change in Z-score BMI"),
  c("_c1", " T0"),
  c("_00", " T0"),
  c("_01", " T1"),
  c("mg_p", "Body fat %"),
  c("nreads", "seq. depth"),
  c("sauces_con", "sauces, con"),
  c("juices_so", "juices, so"),
  c("fats_oils", "fats, oils"),
  c("sweets_pas", "sweets, pas"),
  c("^z_", "Z-score "),
  c("bmi", "BMI"),
  c("status_c2", "status T1"),
  #c("\\.", " "),
  c("_", " ")
)
replace_list <- c(
  list(
    c("\\.\\.", "-", "feature"),
    c("\\.", " ", "feature")
  ),
  purrr::map(replace_strings, \(x) c(x, "metadata")),
  purrr::map(replace_strings, \(x) c(x, "variable"))
)

plotdf <- funcres_all %>% filter(path_code %in% proc2use) %>%
  dplyr::mutate(var_full = paste(variable, value, sep=":")) %>%
  filter(! value %in% c("age_T0", "Girl") ) %>%
  makeReplacementsDF(replace_list = replace_list, totitlecase = c("metadata", "variable"))
head(plotdf)
write_tsv(plotdf, file = paste0(opt$out, "Maaslin2_allMerged_modnames1.tsv"))

#plotdf <- read_tsv(paste0(opt$out, "Maaslin2_allMerged_modnames1.tsv"))
## first, plot related to BMI

table(plotdf$metadata)

plotdf2 <- plotdf %>% dplyr::mutate(metadata = ifelse(dataset != "all", paste(metadata, " (Normal T0)", sep=""), metadata ))
write_tsv(plotdf2, file = paste0(opt$out, "Maaslin2_allMerged_modnames2.tsv"))

met2plot1 <- c("Z-Score BMI T0", "Z-Score Waist T0", "Z-Score BMI T1 (Normal T0)", "Z-Score Waist T1 (Normal T0)")
# use this with plotdf2

met2plot1_t1 <- c("Z-Score BMI T1", "Change in Z-Score BMI", "Z-Score Waist T1", "Change in Z-Score Waist")
df_t1 <- plotdf %>% filter(dataset != "all" & metadata %in% met2plot1_t1)

met2plot1 %in% plotdf2$metadata
met2plot1_t1 %in% df_t1$metadata

outdir <- paste0(opt$out, "loliplots/")
if(!dir.exists(outdir)) dir.create(outdir)

plots1 <- makeLoliplotFromMaaslin2(plotdf2,
                         vars2loliplot = met2plot1,
                         vars2sort = paste0(met2plot1, "__coef"),
                         vars2filter = c("Z-Score BMI T0"),
                         outdir = outdir,
                         name="Loliplot_BMIT0p05",
                         plim = 0.05, plim_col = 0.05, lfclim = 0,
                         strip_fontsize=9,
                         w=14, h=7,
                         contrast_name = "metadata",
                         pval_name =  "qval",
                         lfc_name = "coef",
                         feature_name = "feature",
                         axis_face = "plain",
                         add_lines = FALSE
)

plots1 <- makeLoliplotFromMaaslin2(df_t1,
                                   vars2loliplot = met2plot1_t1,
                                   vars2sort = paste0(met2plot1_t1, "__coef"),
                                   vars2filter =met2plot1_t1[1:2],
                                   outdir = outdir,
                                   name="Loliplot_BMIT1p05",
                                   plim = 0.05, plim_col = 0.05, lfclim = 0,
                                   strip_fontsize=9,
                                   w=12, h=6,
                                   contrast_name = "metadata",
                                   pval_name =  "qval",
                                   lfc_name = "coef",
                                   feature_name = "feature",
                                   axis_face = "plain",
                                   add_lines = FALSE
)


met2plot2 <- c(met2plot1[c(1, 3)], patnames)
plots1 <- makeLoliplotFromMaaslin2(plotdf2,
                                   vars2loliplot = met2plot2,
                                   vars2sort = paste0(met2plot2, "__coef"),
                                   vars2filter = met2plot2,
                                   outdir = outdir,
                                   name="Loliplot_BMIT0_FoodPatterns_p05",
                                   plim = 0.05, plim_col = 0.05, lfclim = 0,
                                   strip_fontsize=9,
                                   w=16, h=12,
                                   contrast_name = "metadata",
                                   pval_name =  "qval",
                                   lfc_name = "coef",
                                   feature_name = "feature",
                                   axis_face = "plain",
                                   add_lines = FALSE
)
plots1 <- makeLoliplotFromMaaslin2(plotdf2,
                                   vars2loliplot = met2plot2,
                                   vars2sort = paste0(met2plot2, "__coef"),
                                   vars2filter = met2plot2,
                                   outdir = outdir,
                                   name="Loliplot_BMIT0_FoodPatterns_p01",
                                   plim = 0.01, plim_col = 0.05, lfclim = 0,
                                   strip_fontsize=9,
                                   w=16, h=7,
                                   contrast_name = "metadata",
                                   pval_name =  "qval",
                                   lfc_name = "coef",
                                   feature_name = "feature",
                                   axis_face = "plain",
                                   add_lines = FALSE
)


food_mod_vars <- c("Energy (Kcal)",
                   "Carbohydrates (g)",
                   "Fiber (g)",
                   "Protein (g)",
                   "Total Fat (g)",
                   "Dairy Derivatives",
                   "Eggs",
                   "Meat",
                   "Fish",
                   "Vegetables",
                   "Tubers",
                   "Nuts",
                   "Oleaginous Fruits",
                   "Refined Cereals",
                   "Whole Grain Cereals",
                   "Legumes",
                   "Fats, Oils",
                   "Sweets, Pastries",
                   "Sugars and Sweets",
                   "Snacks Savory",
                   "Prepared Foods",
                   "Sauces, Condiments",
                   "Water",
                   "Juices, Softdrinks")
met2plot3 <- c(met2plot2, food_mod_vars)
plots1 <- makeLoliplotFromMaaslin2(plotdf2,
                                   vars2loliplot = met2plot3,
                                   vars2sort = paste0(met2plot2, "__coef"),
                                   vars2filter = met2plot2[1:2],
                                   outdir = outdir,
                                   name="Loliplot_FoodItems_p05",
                                   plim = 0.05, plim_col = 0.05, lfclim = 0,
                                   strip_fontsize=9,
                                   w=46, h=8,
                                   contrast_name = "metadata",
                                   pval_name =  "qval",
                                   lfc_name = "coef",
                                   feature_name = "feature",
                                   axis_face = "plain",
                                   add_lines = FALSE
)
plots1 <- makeLoliplotFromMaaslin2(plotdf2,
                                   vars2loliplot = food_mod_vars,
                                   vars2sort = paste0(food_mod_vars, "__coef"),
                                   vars2filter = food_mod_vars,
                                   outdir = outdir,
                                   name="Loliplot_FoodItems2_p05",
                                   plim = 0.05, plim_col = 0.05, lfclim = 0,
                                   strip_fontsize=9,
                                   w=46, h=16,
                                   contrast_name = "metadata",
                                   pval_name =  "qval",
                                   lfc_name = "coef",
                                   feature_name = "feature",
                                   axis_face = "plain",
                                   add_lines = FALSE
)


plots1 <- makeLoliplotFromMaaslin2(plotdf2,
                                   vars2loliplot = patnames,
                                   vars2sort = paste0(patnames, "__coef"),
                                   vars2filter = patnames,
                                   outdir = outdir,
                                   name="Loliplot_FoodPatternsOnly_p05",
                                   plim = 0.05, plim_col = 0.05, lfclim = 0,
                                   strip_fontsize=9,
                                   w=12, h=8,
                                   contrast_name = "metadata",
                                   pval_name =  "qval",
                                   lfc_name = "coef",
                                   feature_name = "feature",
                                   axis_face = "plain",
                                   add_lines = FALSE
)
