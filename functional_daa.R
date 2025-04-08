
source(opt$functional_functions)

opt <- restaurar(opt)
load(paste0(opt$out, "DeSEQ2_v8/DESEQ2_all.RData"))
opt$out <- paste0(opt$out, "Functional/")
if(!dir.exists(opt$out)) dir.create(opt$out)



samples2keep <- sample_data(all_phyloseq$MGBC_plus_NRA_raw) %>% data.frame %>% pull(sampleID)
metacyc_ab <-  readFunctionalMatrix(opt, "humann3_merged_abundances_CPM_named.tsv", sample_substring_index=2) #%>% dplyr::select(all_of(c("Pathway",samples2keep)))
#metacyc_rxn <- readFunctionalMatrix(opt, "humann3_merged_genetables_RXN_CPM.renamed.tsv") %>% dplyr::select(all_of(c("Pathway",samples2keep)))
#ko <- readFunctionalMatrix(opt, "humann3_merged_genetables_KO_CPM.renamed.tsv") %>% dplyr::select(all_of(c("Pathway",samples2keep)))
#go <- readFunctionalMatrix(opt, "humann3_merged_genetables_GO_CPM.renamed.tsv") %>% dplyr::select(all_of(c("Pathway",samples2keep)))
#cazy <- readFunctionalMatrix(opt, "humann3_merged_genetables_CAZY_CPM.renamed.tsv") %>% dplyr::select(all_of(c("Pathway",samples2keep)))



tabs2test <- list("MetaCyc"=metacyc_ab) #, "MetaCyc Reactions"=metacyc_rxn, "KEGG"=ko, "Gene Ontology"=go, "CAZy" = cazy, "GOmixerModules"=modules_nosp2 )
save(tabs2test, file = paste0(opt$out, "functional_tabs.RData"))
tabs2annot <- list("MetaCyc"=NULL) #, "MetaCyc Reactions"=metacyc_mapnames, "KEGG"=ko_mapnames, "Gene Ontology"=go_mapnames, "CAZy" = cazy_mapnames, "GOmixerModules"=modules_mapnames_nosp )

walk2(tabs2test, paste0(opt$out, "/functTabInput_",names(tabs2test), ".tsv"), \(df, name)if(!is.null(df))write_tsv(df, file=name))


#### Limma
#met2use <- sample_data(all_phyloseq$filt) %>% data.frame %>%  
#  dplyr::filter(status_c2 != "Insufficient gain") %>% 
#  dplyr::filter(status_c2 != "initially_overweight")

met2use <- sample_data(all_phyloseq$remove_tanda2) %>% data.frame 

opt$minsampleswithcount <- opt$minfreq*nrow(met2use)
filtered_samples <- lapply(tabs2test, \(x) x %>% dplyr::select(Pathway, all_of(met2use$sampleID)))
filtered_byproc <- lapply(filtered_samples, FUN=filterWholeProcessesAndFreq, opt)
walk2(filtered_byproc, paste0(opt$out, "/functTabInput_",names(tabs2test), "_filteredByProcess.tsv"), \(df, name)if(!is.null(df))write_tsv(df, file=name))

limmares_byproc <- lapply(filtered_byproc, FUN=limma4functional, met2use, "bmi_t0") #"status_c2" "edad_00"

limmares_byproc_annot <- mapply(limmares_byproc, tabs2annot, FUN=nameProcesses, SIMPLIFY = FALSE)
walk2(limmares_byproc_annot, paste0(opt$out, "/DAAlimma_process_",names(limmares_byproc_annot), ".tsv"), \(df, name)write_tsv(df%>% rownames_to_column("Process"), file=name))

restabs_limma <- lapply(limmares_byproc, getSummaryTablesDeseq, opt)
volcanos_limma <- lapply(names(limmares_byproc), 
                         FUN=function(name){
                           fname <- paste0("limma_volcano_SumProcesses_praw_", gsub(" ", "", name), sep="")
                           make_volcano(res = limmares_byproc[[name]], opt=opt, name=fname, pcol="pvalue")}
)
volcanos_limma_padj <- lapply(names(limmares_byproc), 
                              FUN=function(name){
                                fname <- paste0("limma_volcano_SumProcesses_padj_", gsub(" ", "", name), sep="")
                                make_volcano(res = limmares_byproc[[name]], opt=opt, name=fname, pcol="padj")}
)

plims <- c(0.1) #, 0.001, 0.001, 0.001, 0.05, 0.05)
names(plims) <- names(limmares_byproc)
barplots_limma <- lapply(names(limmares_byproc_annot), \(x){
  tab <- limmares_byproc_annot[[x]]
  include_longnames <- TRUE
  if(x == "CAZy"){
    include_longnames <- FALSE
    rownames(tab) <- paste(getCazyClass(rownames(tab)), rownames(tab), sep=" - ")
  }
  makeBarplotFunctional(tab, 
                        plim=plims[x], 
                        paste0("BarplotProcess_", x), opt$out, 
                        include_longnames = include_longnames, 
                        w=12, h=8)
})
# heatmaps_byproc <- lapply(1:length(limmares_byproc), FUN=function(i){
#             makeHeatmapFunctional(limmares_byproc[[i]], met2use, filtered_byproc[[i]],
#                         variable = "Condition",
#                         opt, 
#                         name = paste0("heatmap_byProcess_padj_", gsub(" ", "", names(limmares_byproc)[i]), sep=""), 
#                         logscale=FALSE, 
#                         ptype = "padj", w=20, h=14)
# })


### Functional with Age


