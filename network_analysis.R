library(tidyverse)
library(janitor)
library(DESeq2)
library(phyloseq)
library(minet)
library(igraph)
library(ggpubr)
library(car)

SEED <- 123
MODE = "LOCAL"

if(MODE == "IATA"){
  opt <- list()
}else{
  opt <- list(out ="/home/carmoma/Desktop/202311_DEPRESION/results_rstudio_8/networks1/",
              indir = "/home/carmoma/Desktop/202311_DEPRESION/results_rstudio_8/",
              phyloseq_list = "/home/carmoma/Desktop/202311_DEPRESION/results_rstudio_v2_4/phyloseq/phyloseq_all_list.RData",
              phyloseq_name = "remove_tanda2",
              r_functions="/home/carmoma/Desktop/202311_DEPRESION/depression_scripts/metagenomics_core_functions.R",
              r_functions_networks="/home/carmoma/Desktop/202311_DEPRESION/depression_scripts/functions_networks.R",
              metadata = "/home/carmoma/Desktop/202311_DEPRESION/metadatos_MC_AL12042023_CMcopy.xlsx",
              rewrite=TRUE,
              fc=1, 
              pval=0.05,
              adjust_pvals = TRUE
  )
}
if(! dir.exists(opt$out)){dir.create(opt$out)}

### LOAD DATA
source(opt$r_functions)
source(opt$r_functions_networks)
restaurar <- restauraropt_mk(opt)

load(opt$phyloseq_list)
phobj <- all_phyloseq[[opt$phyloseq_name]]
phobj <- updatePsWithLogs(phobj, c("Edad", "IMC"))
metadata <- sample_data(phobj) %>% data.frame

#load(paste0(opt$indir, "DESeq2_ControlVars/DeSEQ2/remove_tanda2_IMC_log/DESEQ2_all_results_remove_tanda2_IMC_log.R"))
#load(paste0(opt$indir, "DESeq2_ControlVarsMany/LFC_Comparison_AgeAndBMI_allCombos.RData"))
vstdf <- read_tsv(paste0(opt$indir, "DeSEQ2/remove_tanda2/remove_tanda2_vst_counts.tsv"))
normdf <- read_tsv(paste0(opt$indir, "DeSEQ2/remove_tanda2/remove_tanda2_norm_counts.tsv"))
#daatab <- read_tsv(paste0(opt$indir, "DeSEQ2/remove_tanda2/remove_tanda2_Condition_Depression_vs_Control_DAAshrinkNormal.tsv"))
#daataxa <- daatab %>% dplyr::filter(padj<0.05) %>% pull(taxon)

load(paste0("/home/carmoma/Desktop/202311_DEPRESION/results_rstudio_v2_4/", 
            "DESeq2_ControlVarsMany/LFC_Comparison_AgeAndBMI_allCombos.RData"))
daalist <- list(
  "D_vs_C" = dea2contrasts$firstContrast$resdf,
  "D_vs_C_adj_BMI" = dea2contrasts$contrastlist2$Condition_corrIMC$resdf,
  "BMI" = dea2contrasts$contrastlist2$BMI_alone$resdf,
  "BMI_adj_Depr" = dea2contrasts$contrastlist2$BMI_corrCond$resdf
)
taxalist <- map(daalist, \(x)x %>% dplyr::filter(padj < opt$pval) %>% pull(taxon)) %>% unlist %>% unique
taxalist2 <- dea2contrasts$contrastlist2$Condition_corrIMC$resdf %>% dplyr::filter(padj < opt$pval) %>% pull(taxon)
# https://github.com/ryanjw/co-occurrence
# https://rdrr.io/bioc/minet/man/minet.html  


# closeness centrality
# betweenness centrality
# Edge Degree
# keystone taxa -> highest combination of all 3 measures and average abn > 1%

net_estimator <- "spearman"
net_method <- "mrnet"

outdir <- paste0(opt$out, "/nets_4groups/")
if(!dir.exists(outdir)) dir.create(outdir)

metfilt <- metadata %>% dplyr::filter(!is.na(IMC)) %>% 
  dplyr::mutate(depr_and_ob = paste(as.character(Condition), 
                                    ifelse(IMC>25, "overweight", "normal"), sep=":")) %>% 
  group_by(depr_and_ob) %>% 
  group_split %>% 
  map( \(x){nn <- unique(x$depr_and_ob)[1]; y<-x %>% pull(sampleID) %>% list(.); names(y)<- nn; return(y)}) %>% 
  flatten()

nets <- map(names(metfilt), \(x){
  getGraphFromSamples(phobj, 
                      samples = metfilt[[x]], 
                      daataxa = taxalist, 
                      vstdf = vstdf, 
                      net_estimator = net_estimator, 
                      net_method = net_method, 
                      daatab = daalist$D_vs_C_adj_BMI,
                      daatab2 = NULL,
                      outdir=outdir,
                      filter_empty=FALSE, 
                      filt_quantile=0.95, 
                      name=paste0(gsub(":", "_", x), "_2cols"), 
                      w=12, h=12)
})
names(nets) <- names(metfilt)
save(nets, file = paste0(outdir, "/networks.RData"))

nets2 <- map(names(metfilt), \(x){
  getGraphFromSamples(phobj, 
                      samples = metfilt[[x]], 
                      daataxa = taxalist, 
                      vstdf = vstdf, 
                      net_estimator = net_estimator, 
                      net_method = net_method, 
                      daatab = daalist$D_vs_C_adj_BMI,
                      daatab2 = daalist$BMI_adj_Depr,
                      outdir=outdir,
                      filter_empty=FALSE, 
                      filt_quantile=0.95, 
                      name=paste0(gsub(":", "_", x), "_4cols"), 
                      w=12, h=12)
})
names(nets2) <- names(metfilt)
save(nets2, file = paste0(outdir, "/networks2.RData"))

nodeprops <- map(names(nets2), \(x){
  df <- getGraphProps(nets2[[x]])
  df$class <- x
  return(df)
}) %>% bind_rows() %>% 
  dplyr::mutate(
    class=factor(class, levels=c(
      "Control:normal", 
      "Depression:normal", 
      "Control:overweight", 
      "Depression:overweight"
    ))
  )
write_tsv(nodeprops, file = paste0(outdir, "/node_properties.tsv"))
bplot <- make_boxplot_nodeprops(nodeprops, outdir, "4groups", w=12, h=6, correct_pvals = TRUE)

map(names(nets), \(x){
  plotNetWithTopNames(nets[[x]], nodeprops, x, ntop=10, outdir, paste0(x, "_names"), add_labels = FALSE, w=12, h=12)
})

map(names(nets), \(x){
  plotNetWithTopNames(nets[[x]], nodeprops, x, ntop=10, outdir, paste0(x, "_nonames"), add_labels = FALSE, w=12, h=12)
})

list2recode <- list("#A1C6EA"="Control-normal", 
                    "#FD8B2F"="Depressed-normal",
                    "#00AA5A"="Control-overweight",
                    "#8E7BFF"="Depressed-overweight", 
                    "#A5ABBD"="NS")
nodeprops_b <- nodeprops %>% 
  dplyr::mutate(color=dplyr::recode(color, 
                                    !!!{list2recode})
  ) %>% 
  dplyr::mutate(color=factor(color, levels=c("Control-normal", 
                                             "Depressed-normal",
                                             "Control-overweight",
                                             "Depressed-overweight", 
                                             "NS")))
makeTopTaxaPlot(nodeprops_b, 10, outdir, "DeprAndOb", w = 8, h = 10, manual_scale=names(list2recode)[-4])

## Repeeat with 2 colors
nodeprops <- map(names(nets), \(x){
  df <- getGraphProps(nets[[x]])
  df$class <- x
  return(df)
}) %>% bind_rows() %>% 
  dplyr::mutate(
    class=factor(class, levels=c(
      "Control:normal", 
      "Depression:normal", 
      "Control:overweight", 
      "Depression:overweight"
    ))
  )
write_tsv(nodeprops, file = paste0(outdir, "/node_properties_2cols.tsv"))

map(names(nets), \(x){
  plotNetWithTopNames(nets[[x]], nodeprops, x, ntop=10, outdir, paste0(x, "_names2"), add_labels = FALSE, w=12, h=12)
})

map(names(nets), \(x){
  plotNetWithTopNames(nets[[x]], nodeprops, x, ntop=10, outdir, paste0(x, "_nonames2"), add_labels = FALSE, w=12, h=12)
})


list2recode <- list("#A1C6EA"="Control",
                    "#FD8B2F"="Depression", 
                    "#A5ABBD"="NS")
nodeprops_b <- nodeprops %>% 
  dplyr::mutate(color=dplyr::recode(color, 
                                    !!!{list2recode})
  ) %>% 
  dplyr::mutate(color=factor(color, levels=c("Control", 
                                             "Depression",
                                             "NS")))
makeTopTaxaPlot(nodeprops_b, 10, outdir, "DeprAndOb2cols", w = 8, h = 10, manual_scale=names(list2recode))

### ONLY DEPR AND CONTROL

outdir <- paste0(opt$out, "/nets_2groups/")
if(!dir.exists(outdir)) dir.create(outdir)

metfilt2 <- metadata %>% 
  group_by(Condition) %>% 
  group_split %>% 
  map( \(x){nn <- unique(x$Condition)[1]; y<-x %>% pull(sampleID) %>% list(.); names(y)<- nn; return(y)}) %>% 
  flatten()

nets_2g <- map(names(metfilt2), \(x){
  getGraphFromSamples(phobj, 
                      samples = metfilt2[[x]], 
                      daataxa = taxalist2, 
                      vstdf = vstdf, 
                      net_estimator = net_estimator, 
                      net_method = net_method, 
                      daatab = daalist$D_vs_C_adj_BMI,
                      daatab2 = NULL,
                      outdir=outdir,
                      filter_empty=FALSE, 
                      filt_quantile=0.95, 
                      name=paste0(gsub(":", "_", x), "_2groups"), 
                      w=12, h=12)
})
names(nets_2g) <- names(metfilt2)
save(nets_2g, file = paste0(outdir, "/networks.RData"))

nodeprops <- map(names(nets_2g), \(x){
  df <- getGraphProps(nets_2g[[x]])
  df$class <- x
  return(df)
}) %>% bind_rows() %>% 
  dplyr::mutate(
    class=factor(class, levels=c(
      "Control", 
      "Depression"
    ))
  )
write_tsv(nodeprops, file = paste0(outdir, "/node_properties.tsv"))
bplot <- make_boxplot_nodeprops(nodeprops, outdir, name, w=8, h=8,
                                correct_pvals = TRUE)



map(names(nets_2g), \(x){
  plotNetWithTopNames(nets_2g[[x]], nodeprops, x, ntop=10, outdir, paste0(x, "_names"))
})

map(names(nets_2g), \(x){
  plotNetWithTopNames(nets_2g[[x]], nodeprops, x, ntop=10, outdir, paste0(x, "_nonames"), add_labels = FALSE)
})

list2recode <- list("#A1C6EA"="Control",
                    "#FD8B2F"="Depression", 
                    "#A5ABBD"="NS")
nodeprops_b <- nodeprops %>% 
  dplyr::mutate(color=dplyr::recode(color, 
                                    !!!{list2recode})
  ) %>% 
  dplyr::mutate(color=factor(color, levels=c("Control", 
                                             "Depression",
                                             "NS")))

makeTopTaxaPlot(nodeprops_b, 10, outdir, "DeprVsCtrl", w = 8, h = 6)
