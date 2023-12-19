library(tidyverse)
library(DESeq2)
library(phyloseq)
library(janitor)
library(ggvenn)

## Mediation analysis based on tutorials:
## https://advstats.psychstat.org/book/mediation/index.php 
## https://rpubs.com/Momen/485122
## https://cran.r-project.org/web/packages/manymome/vignettes/med_lm.html
## https://stats.stackexchange.com/questions/19599/multiple-mediation-analysis-in-r

### Requirements:

### 1) Show that X is correlated with Y. Regress Y on X to estimate and test the path c
### 2) This step establishes that there is an effect that may be mediated.
### Show that X is correlated with M. Regress M on X to estimate and test path a
### 3) This step essentially involves treating the mediator as if it were an outcome variable.
### Show that M affects Y. Regress Y on both X and M to estimate and test path b
### 4) Note that it is not sufficient just to correlate the mediator with the outcome; the mediator and the outcome may be correlated because they are both caused by the input variable X. Thus, the input variable X must be controlled in establishing the effect of the mediator on the outcome.
### To establish that M completely mediates the X-Y relationship, the effect of X on Y controlling for M (path c′
### ) should be zero. The effects in both Steps 3 and 4 are estimated in the same equation.

### Only 2 and 3 are required. 

### How to meet them? 
### Our functional model is bacteria --> BMI --> depression AND bacteria --> depression. 
### 1) -> Differentially abundant bacteria (i.e, X bacteria are correlated with Y depression). We will select all bacteria in simple DAA
### 2) -> Mediator (BMI) is correlated to X (bacterial species). 
### 3) 

SEED <- 123
MODE = "LOCAL"

if(MODE == "IATA"){
  opt <- list()
}else{
  opt <- list(out ="/home/carmoma/Desktop/202311_DEPRESION/results_rstudio_v2_4/mediation_analysis2/",
              indir = "/home/carmoma/Desktop/202311_DEPRESION/results_rstudio_v2_4/",
              phyloseq_list = "/home/carmoma/Desktop/202311_DEPRESION/results_rstudio_v2_4/phyloseq/phyloseq_all_list.RData",
              phyloseq_name = "remove_tanda2",
              r_functions="/home/carmoma/Desktop/202311_DEPRESION/depression_scripts/metagenomics_core_functions.R",
              r_functions_mediation="/home/carmoma/Desktop/202311_DEPRESION/depression_scripts/mediation_functions.R",
              metadata = "/home/carmoma/Desktop/202311_DEPRESION/metadatos_MC_AL12042023_CMcopy.xlsx",
              rewrite=TRUE,
              fc=1, 
              pval=0.05
  )
}
if(! dir.exists(opt$out)){dir.create(opt$out)}

### LOAD DATA
source(opt$r_functions)
source(opt$r_functions_mediation)
load(opt$phyloseq_list)
phobj <- all_phyloseq[[opt$phyloseq_name]]
phobj <- updatePsWithLogs(phobj, c("Edad", "IMC"))
metadata <- sample_data(phobj) %>% data.frame

alphadiv <- read_tsv(paste0(opt$indir, "AlphaDiversity/remove_tanda2_AlphaDiv.tsv")) %>% 
  dplyr::select(sampleID, Observed, Chao1, Shannon, InvSimpson)

metadata2 <- merge(metadata, alphadiv, by="sampleID")

load(paste0(opt$indir, "DESeq2_ControlVarsMany/LFC_Comparison_AgeAndBMI_allCombos.RData"))
vstdf <- read_tsv(paste0(opt$indir, "DeSEQ2/remove_tanda2/remove_tanda2_vst_counts.tsv"))
normdf <- read_tsv(paste0(opt$indir, "DeSEQ2/remove_tanda2/remove_tanda2_norm_counts.tsv"))

### Do tests
#first do a small Venn diagram

vars2venn <- list(
  "D vs C" = dea2contrasts$firstContrast$resdf %>% dplyr::filter(padj < 0.05) %>% pull(taxon),
  #"D vs C adj. BMI" = summary_df %>% dplyr::filter(depr_adjimc_padj < plim) %>% pull(variable),
  "adj BMI" = dea2contrasts$contrastlist2$Condition_corrIMC$resdf %>% dplyr::filter(padj < 0.05) %>% pull(taxon),
  "adj Age" = dea2contrasts$contrastlist2$Condition_corrAge$resdf %>% dplyr::filter(padj < 0.05) %>% pull(taxon),
  "adj both" = dea2contrasts$contrastlist2$Condition_corr2$resdf %>% dplyr::filter(padj < 0.05) %>% pull(taxon)
)
gv <- ggvenn(
  vars2venn, columns = names(vars2venn),
  stroke_size = 0.5,
  stroke_color = C_NS,
  fill_color = c(C_CASE, C_CTRL2, C_CTRL, C_CASE2),show_elements = F
)
ggsave(filename = paste0(opt$out, "VennDiagram_CvsD_control.pdf"), gv, width = 8, height = 8)

#### Read data 

expr_df <- vstdf %>% column_to_rownames("gene") %>% 
  as.matrix() %>% t %>% data.frame %>% 
  rownames_to_column("sampleID")

df_all <- merge(metadata2, expr_df, by="sampleID")
df_all$Condition_bin <- ifelse(df_all$Condition=="Depression", TRUE, FALSE)
names(df_all) <- gsub("^X\\.", "", names(df_all)) %>% 
                 gsub("\\._", "_", .) %>% 
                 gsub("\\.$", "", .)

expr_df_norm <- normdf %>% column_to_rownames("gene") %>% 
  as.matrix() %>% t %>% data.frame %>% 
  rownames_to_column("sampleID")

df_all_norm <- merge(metadata2, expr_df_norm, by="sampleID")
df_all_norm$Condition_bin <- ifelse(df_all_norm$Condition=="Depression", TRUE, FALSE)
names(df_all_norm) <- gsub("^X\\.", "", names(df_all_norm)) %>% 
  gsub("\\._", "_", .) %>% 
  gsub("\\.$", "", .)


write_tsv(df_all, file=paste0(opt$out, "merged_data_vst.tsv"))
write_tsv(df_all_norm, file=paste0(opt$out, "merged_data_norm.tsv"))

### MEDIATION ANALYSIS WITH IMC

summary_df <- data.frame(variable = gsub("-", ".", vstdf$gene) %>% 
                                    gsub("\\[|\\]", "", .) %>% 
                                    gsub("\\._", "_", .) %>% 
                                    gsub("\\/|\\(|\\)", ".", .) %>% 
                                    gsub("\\.$", "", .) 
                         )

assertthat::assert_that(all(summary_df$variable %in% names(df_all)))

list2merge <- list(
  depr_only_padj = dea2contrasts$firstContrast$resdf,
  imc_only_padj = dea2contrasts$contrastlist2$BMI_alone$resdf,
  depr_adjimc_padj = dea2contrasts$contrastlist2$Condition_corrIMC$resdf,
  imc_adjdepr_padj = dea2contrasts$contrastlist2$BMI_corrCond$resdf
)

merged_pvals <- list2merge %>% 
  lapply(\(x){
  x$taxon <- x$taxon %>% 
    gsub("-", ".", .) %>% 
    gsub("\\[|\\]", "", .) %>% 
    gsub("\\._", "_", .) %>% 
    gsub("\\/|\\(|\\)", ".", .) %>% 
    gsub("\\.$", "", .) 
  x <- x[match(summary_df$variable, x$taxon), ]
  return(x$padj)
}) %>% bind_cols()
summary_df <- cbind(summary_df, merged_pvals)


mediator_name <- "IMC"
y_name <- "Condition_bin"
plim <- 0.05

vars2test <- summary_df %>% 
  dplyr::filter(depr_only_padj < plim &
                  #imc_only_padj < plim  & 
                  imc_adjdepr_padj < plim) %>% 
  pull(variable)

vars2venn <- list(
  "D vs C" = summary_df %>% dplyr::filter(depr_only_padj < plim) %>% pull(variable),
  #"D vs C adj. BMI" = summary_df %>% dplyr::filter(depr_adjimc_padj < plim) %>% pull(variable),
  "BMI" = summary_df %>% dplyr::filter(imc_only_padj < plim) %>% pull(variable),
  "BMI adj. Depr" = summary_df %>% dplyr::filter(imc_adjdepr_padj < plim) %>% pull(variable)
)

gv <- ggvenn(
  vars2venn, columns = names(vars2venn),
  stroke_size = 0.5,
  stroke_color = C_NS,
  fill_color = c(C_CASE, C_CASE2, C_CTRL2, C_CTRL),show_elements = F
)
ggsave(filename = paste0(opt$out, "VennDiagram_IMC.pdf"), gv)

medresults <- lapply(vars2test, \(x, df_all){
  df <- df_all %>% dplyr::select(all_of(c(x, y_name, mediator_name)))
  with_nas <-  df %>% apply(MAR=1, \(x)any(is.na(x)))
  df <- df[!with_nas, ]
  param_names <- c("b", "cp", "a", "a*b", "cp+a*b")
  res <- tryCatch({makeMediationSimple(df, x, y_name, mediator_name)}, error=\(x){
    
      xx <- data.frame(Estimate=rep(NA, 5), 
                       S.E. = rep(NA, 5), 
                       `Z-score`=rep(NA, 5), 
                       p.value=rep(NA, 5))
      rownames(xx) <- c(param_names)
      return(list(estimates=xx))   
  })
  res2 <- res$estimates %>% rownames_to_column("param") %>% 
    dplyr::select(param, Estimate, p.value) %>% 
    gather(key="var", value="value", Estimate, p.value) %>% 
    filter(param %in% param_names) %>% 
    unite("tmp", param, var, sep="_") %>% 
    spread(tmp, value) %>% 
    mutate(Xvar = x, Yvar = y_name, Mediator=mediator_name, fullmodel=list(res$estimates))
  return(res2)
  
}, df_all) %>% bind_rows()
#all(medresults$Xvar %in% summary_df$variable)
medresults_merged <- merge(medresults %>% dplyr::select(-fullmodel), summary_df, by.x="Xvar", by.y="variable")
write_tsv(medresults_merged, file=paste0(opt$out, "mediation_analysis_IMC.tsv"))
save(medresults_merged, file=paste0(opt$out, "mediation_analysis_IMC.RData"))

## Transform to plot
def_effects <- c('a', 'b', 'cp','a*b', 'cp+a*b') 
res_trans <- map(1:nrow(medresults), \(i){
  res_i <- medresults$fullmodel[[i]]
  oldnames2 <- rownames(res_i) %>% 
    sapply(\(x)strsplit(x, "\\*|\\+",perl = TRUE)[[1]][1]) 
  resdf_i <- res_i %>% 
    rownames_to_column("param") %>% 
    mutate(x_labels = ifelse(param %in% def_effects, 
                             medresults$Xvar[i],
                             gsub("V\\[|\\]", "", param)),
           param = ifelse(param %in% def_effects, 
                          gsub("(b|a|cp)", paste0("\\1", as.character(i)), param, perl=T), param) 
           )

}) %>% bind_rows()
bs_only <- res_trans %>% dplyr::filter(grepl("^b[0-9]+$", param))
(g_bs <- ggplot(bs_only,aes(x=Estimate))+geom_histogram()+mytheme+ggtitle("Histogram of b (IMC->Depr)"))
ggsave(filename = paste0(opt$out, "histogram_bs_IMC_separate.pdf"), g_bs)

res_final <- rbind(
  bs_only %>% 
    mutate_if(is.character, \(x)"b") %>% 
    group_by(param, x_labels) %>% 
    dplyr::summarise_all(mean) %>% 
    dplyr::select(all_of(names(res_trans))),
  res_trans %>% dplyr::filter(!grepl("^b[0-9]+$", param)) %>% 
    dplyr::mutate(param=gsub("b[0-9]+$", "b", param))
) 
write_tsv(res_final, file=paste0(opt$out, "mediation_analysis_IMC_separate_reformat.tsv"))
write_tsv(bs_only, file=paste0(opt$out, "mediation_analysis_IMC_separate_bs.tsv"))

name<-"analysis_IMC_separateModel_vjust"
plotres <- plotIMCMediationSimple(res_final, vars2test, opt$out, name, plim_plot = 0.05, 
                       use_color_scale = FALSE, w=14, h=10)

## Make also boxplot
bacnames <- plotres$bacorder$x_labels
plots <- makePlotBySpecies(bacnames, df_all, opt$out, "IMC_separate_BySpeciesPearson", quantvar="IMC_log", 
                           quantvar_name = "log(IMC)", corrmethod = "pearson", w=8, h=10)
plots <- makePlotBySpecies(bacnames, df_all, opt$out, "IMC_separate_BySpeciesSpearman", quantvar="IMC_log", 
                           quantvar_name = "log(IMC)", corrmethod = "spearman", w=8, h=10)
plots <- makePlotBySpecies(bacnames, df_all, opt$out, "IMC_separate_BySpeciesKendall", quantvar="IMC_log", 
                           quantvar_name = "log(IMC)", corrmethod = "kendall", w=8, h=10)
# 
# plots <- makePlotBySpecies(bacnames, df_all_norm, opt$out, "IMCseparate_NormCountsBySpeciesPearson", quantvar="IMC_log", 
#                            quantvar_name = "log(IMC)", corrmethod = "pearson", w=8, h=10)
# plots <- makePlotBySpecies(bacnames, df_all_norm, opt$out, "IMCseparate_NormCountsBySpeciesSpearman", quantvar="IMC_log", 
#                            quantvar_name = "log(IMC)", corrmethod = "spearman", w=8, h=10)
# plots <- makePlotBySpecies(bacnames, df_all_norm, opt$out, "IMCseparate_NormCountsBySpeciesKendall", quantvar="IMC_log", 
#                            quantvar_name = "log(IMC)", corrmethod = "kendall", w=8, h=10)

###########################################################################
### MEDIATION WITH IMC, MERGED

summary_df <- data.frame(variable = gsub("-", ".", vstdf$gene) %>% 
                           gsub("\\[|\\]", "", .) %>% 
                           gsub("\\._", "_", .) %>% 
                           gsub("\\/|\\(|\\)", ".", .) %>% 
                           gsub("\\.$", "", .) 
)

assertthat::assert_that(all(summary_df$variable %in% names(df_all)))

list2merge <- list(
  depr_only_padj = dea2contrasts$firstContrast$resdf,
  imc_only_padj = dea2contrasts$contrastlist2$BMI_alone$resdf,
  depr_adjimc_padj = dea2contrasts$contrastlist2$Condition_corrIMC$resdf,
  imc_adjdepr_padj = dea2contrasts$contrastlist2$BMI_corrCond$resdf
)

merged_pvals <- list2merge %>% 
  lapply(\(x){
    x$taxon <- x$taxon %>% 
      gsub("-", ".", .) %>% 
      gsub("\\[|\\]", "", .) %>% 
      gsub("\\._", "_", .) %>% 
      gsub("\\/|\\(|\\)", ".", .) %>% 
      gsub("\\.$", "", .) 
    x <- x[match(summary_df$variable, x$taxon), ]
    return(x$padj)
  }) %>% bind_cols()
summary_df <- cbind(summary_df, merged_pvals)


mediator_name <- "IMC"
y_name <- "Condition_bin"
plim <- 0.05

vars2test <- summary_df %>% 
  dplyr::filter(depr_only_padj<plim & 
                  imc_only_padj < plim & 
                  imc_adjdepr_padj < plim) %>% 
  pull(variable)
length(vars2test)
#vars2test <- "Actinomyces_naeslundii"
df <- df_all %>% dplyr::select(all_of(c(vars2test, y_name, mediator_name)))
with_nas <-  df %>% apply(MAR=1, \(x)any(is.na(x)))
df <- df[!with_nas, ]
res <- makeMediationSimple_mergedX(df, vars2test, y_name, mediator_name)
write_tsv(res, file=paste0(opt$out, "mediation_analysis_IMC_mergedModel.tsv"))

### PLOT 

name<-"analysis_IMC_mergedModel_vjust"
plotres2 <- plotIMCMediationSimple(res, vars2test, opt$out, name, plim_plot = 0.05, use_color_scale = FALSE)


## Make also boxplot
bacnames <- plotres2$bacorder$x_labels
plots <- makePlotBySpecies(bacnames, df_all, opt$out, "IMC_merged_BySpeciesPearson", quantvar="IMC_log", 
                           quantvar_name = "log(IMC)", corrmethod = "pearson", w=8, h=10)
plots <- makePlotBySpecies(bacnames, df_all, opt$out, "IMC_merged_BySpeciesSpearman", quantvar="IMC_log", 
                           quantvar_name = "log(IMC)", corrmethod = "spearman", w=8, h=10)
plots <- makePlotBySpecies(bacnames, df_all, opt$out, "IMC_merged_BySpeciesKendall", quantvar="IMC_log", 
                           quantvar_name = "log(IMC)", corrmethod = "kendall", w=8, h=10)
# 

###################################################################
### MEDIATION ANALYSIS WITH AGE

summary_df <- data.frame(variable = gsub("-", ".", vstdf$gene) %>% 
                           gsub("\\[|\\]", "", .) %>% 
                           gsub("\\._", "_", .) %>% 
                           gsub("\\/|\\(|\\)", ".", .) %>% 
                           gsub("\\.$", "", .) 
)

assertthat::assert_that(all(summary_df$variable %in% names(df_all)))

list2merge <- list(
  depr_only_padj = dea2contrasts$firstContrast$resdf,
  age_only_padj = dea2contrasts$contrastlist2$Age_alone$resdf,
  depr_adjage_padj = dea2contrasts$contrastlist2$Condition_corrAge$resdf,
  age_adjdepr_padj = dea2contrasts$contrastlist2$Age_corrCond$resdf
)

merged_pvals <- list2merge %>% 
  lapply(\(x){
    x$taxon <- x$taxon %>% 
      gsub("-", ".", .) %>% 
      gsub("\\[|\\]", "", .) %>% 
      gsub("\\._", "_", .) %>% 
      gsub("\\/|\\(|\\)", ".", .) %>% 
      gsub("\\.$", "", .) 
    x <- x[match(summary_df$variable, x$taxon), ]
    return(x$padj)
  }) %>% bind_cols()
summary_df <- cbind(summary_df, merged_pvals)


x_name <- "Edad"
y_name <- "Condition_bin"
plim <- 0.05

vars2test <- summary_df %>% 
  dplyr::filter(depr_only_padj < plim & 
                  #age_only_padj < plim &
                  age_adjdepr_padj < plim) %>% 
  pull(variable)

vars2venn <- list(
  "D vs C" = summary_df %>% dplyr::filter(depr_only_padj < plim) %>% pull(variable),
  "Age" = summary_df %>% dplyr::filter(age_only_padj < plim) %>% pull(variable),
  "Age adj. Depr" = summary_df %>% dplyr::filter(age_adjdepr_padj < plim) %>% pull(variable)
)

gv <- ggvenn(
  vars2venn, columns = names(vars2venn),
  stroke_size = 0.5,
  stroke_color = C_NS,
  fill_color = c(C_CASE, C_CASE2, C_CTRL2),show_elements = F
)
ggsave(filename = paste0(opt$out, "VennDiagram_Age.pdf"), gv)


medresultsAge <- lapply(vars2test, \(mediator_name, df_all){
  df <- df_all %>% dplyr::select(all_of(c(x_name, y_name, mediator_name)))
  with_nas <-  df %>% apply(MAR=1, \(x)any(is.na(x)))
  df <- df[!with_nas, ]
  param_names <- c("b", "cp", "a", "a*b", "cp+a*b")
  res <- tryCatch({makeMediationSimple(df, x_name, y_name, mediator_name)}, error=\(x){
    
    xx <- data.frame(Estimate=rep(NA, 5), 
                     S.E. = rep(NA, 5), 
                     `Z-score`=rep(NA, 5), 
                     p.value=rep(NA, 5))
    rownames(xx) <- c(param_names)
    return(list(estimates=xx))   
  })
  res2 <- res$estimates %>% 
    rownames_to_column("param") %>% 
    dplyr::select(param, Estimate, p.value) %>% 
    gather(key="var", value="value", Estimate, p.value) %>% 
    dplyr::filter(param %in% param_names) %>% 
    unite("tmp", param, var, sep="_") %>% 
    spread(tmp, value) %>% 
    dplyr::mutate(Xvar = x_name, Yvar = y_name, Mediator=mediator_name, fullmodel=list(res$estimates))
  return(res2)
  
}, df_all) %>% bind_rows() 
#all(medresults$Xvar %in% summary_df$variable)

medresultsAge_merged <- merge(medresultsAge, summary_df, by.x="Mediator", by.y="variable")
write_tsv(medresultsAge_merged, file=paste0(opt$out, "mediation_analysis_Edad.tsv"))

## Transform to plot ##falta
def_effects <- c('a', 'b', 'cp','a*b', 'cp+a*b') 
res_trans <- map(1:nrow(medresultsAge), \(i){
  res_i <- medresultsAge$fullmodel[[i]]
  #oldnames2 <- rownames(res_i) %>% 
  #  sapply(\(x)strsplit(x, "\\*|\\+",perl = TRUE)[[1]][1]) 
  resdf_i <- res_i %>% 
    rownames_to_column("param") %>% 
    mutate(x_labels = ifelse(param %in% def_effects, 
                             medresultsAge$Mediator[i],
                             gsub("V\\[|\\]", "", param)),
           param = ifelse(param %in% def_effects, 
                          gsub("(b|a|cp)", paste0("\\1", as.character(i)), param, perl=T), param) 
    )
  
}) %>% bind_rows()
cps_only <- res_trans %>% dplyr::filter(grepl("^cp[0-9]+$", param))
(g_bs <- ggplot(cps_only,aes(x=Estimate))+geom_histogram()+mytheme+ggtitle("Histogram of b (IMC->Depr)"))
ggsave(filename = paste0(opt$out, "histogram_bs_IMC_separate.pdf"), g_bs)

res_final <- rbind(
  cps_only %>% 
    mutate(param= "cp", x_labels="cp") %>% 
    group_by(param, x_labels) %>% 
    dplyr::summarise_if(is.numeric, mean) %>% 
    dplyr::select(all_of(names(res_trans))),
  res_trans %>% dplyr::filter(!grepl("^cp[0-9]+$", param)) %>% 
    dplyr::mutate(param=gsub("cp[0-9]+", "cp", param))
) 
write_tsv(res_final, file=paste0(opt$out, "mediation_analysis_Edad_separate_reformat.tsv"))
write_tsv(cps_only, file=paste0(opt$out, "mediation_analysis_Edad_separate_cps.tsv"))


name<-"analysis_Age_separateModel_vjust"
plotres3 <- plotAgeMediationSimple(res_final, vars2test, opt$out, name, plim_plot = 0.05, use_color_scale = FALSE, 
                                   w=10, h=6)

## Make also boxplot
bacnames <- plotres3$bacorder$x_labels
plots <- makePlotBySpecies(bacnames, df_all, opt$out, "Age_separate_BySpeciesPearson", quantvar="Edad_log", 
                           quantvar_name = "log(Age)", corrmethod = "pearson", w=8, h=10)
plots <- makePlotBySpecies(bacnames, df_all, opt$out, "Age_separate_BySpeciesSpearman", quantvar="Edad_log", 
                           quantvar_name = "log(Age)", corrmethod = "spearman", w=8, h=10)
plots <- makePlotBySpecies(bacnames, df_all, opt$out, "Age_separate_BySpeciesKendall", quantvar="Edad_log", 
                           quantvar_name = "log(Age)", corrmethod = "kendall", w=8, h=10)
# 

##################################################################################33
### MEDIATION ANALYSIS WITH AGE, MERGED

summary_df <- data.frame(variable = gsub("-", ".", vstdf$gene) %>% 
                           gsub("\\[|\\]", "", .) %>% 
                           gsub("\\._", "_", .) %>% 
                           gsub("\\/|\\(|\\)", ".", .) %>% 
                           gsub("\\.$", "", .) 
)

assertthat::assert_that(all(summary_df$variable %in% names(df_all)))

list2merge <- list(
  depr_only_padj = dea2contrasts$firstContrast$resdf,
  age_only_padj = dea2contrasts$contrastlist2$Age_alone$resdf,
  depr_adjage_padj = dea2contrasts$contrastlist2$Condition_corrAge$resdf,
  age_adjdepr_padj = dea2contrasts$contrastlist2$Age_corrCond$resdf
)

merged_pvals <- list2merge %>% 
  lapply(\(x){
    x$taxon <- x$taxon %>% 
      gsub("-", ".", .) %>% 
      gsub("\\[|\\]", "", .) %>% 
      gsub("\\._", "_", .) %>% 
      gsub("\\/|\\(|\\)", ".", .) %>% 
      gsub("\\.$", "", .) 
    x <- x[match(summary_df$variable, x$taxon), ]
    return(x$padj)
  }) %>% bind_cols()
summary_df <- cbind(summary_df, merged_pvals)


x_name <- "Edad"
y_name <- "Condition_bin"
plim <- 0.05

vars2test <- summary_df %>% 
  dplyr::filter(depr_only_padj < plim  & 
                  age_only_padj < plim & 
                  age_adjdepr_padj) %>% 
  pull(variable)

df <- df_all %>% dplyr::select(all_of(c(x_name, y_name, vars2test)))
with_nas <-  df %>% apply(MAR=1, \(x)any(is.na(x)))
df <- df[!with_nas, ]
res <- makeMediationSimple_mergedMediat(df, x_name, y_name, vars2test)
write_tsv(res, file=paste0(opt$out, "mediation_analysis_Edad_mergedModel.tsv"))

## Plot network
name<-"analysis_Age_mergedModel_vjust"
plotres4 <- plotAgeMediationSimple(res, vars2test, opt$out, name, plim_plot = 0.05, use_color_scale = FALSE, 
                                   w=10, h=6)

## Make also boxplot
bacnames <- plotres4$bacorder$x_labels
plots <- makePlotBySpecies(bacnames, df_all, opt$out, "Age_merged_BySpeciesPearson", quantvar="Edad_log", 
                           quantvar_name = "log(Age)", corrmethod = "pearson", w=8, h=6)
plots <- makePlotBySpecies(bacnames, df_all, opt$out, "Age_merged_BySpeciesSpearman", quantvar="Edad_log", 
                           quantvar_name = "log(Age)", corrmethod = "spearman", w=8, h=6)
plots <- makePlotBySpecies(bacnames, df_all, opt$out, "Age_merged_BySpeciesKendall", quantvar="Edad_log", 
                           quantvar_name = "log(Age)", corrmethod = "kendall", w=8, h=6)
# 

### MEDIATION ANALYSIS WITH AGE AND IMC

summary_df <- data.frame(variable = gsub("-", ".", vstdf$gene) %>% 
                           gsub("\\[|\\]", "", .) %>% 
                           gsub("\\._", "_", .) %>% 
                           gsub("\\/|\\(|\\)", ".", .) %>% 
                           gsub("\\.$", "", .) 
)

assertthat::assert_that(all(summary_df$variable %in% names(df_all)))

list2merge <- list(
  depr_only_padj = dea2contrasts$firstContrast$resdf,
  age_only_padj = dea2contrasts$contrastlist2$Age_alone$resdf,
  imc_only_padj = dea2contrasts$contrastlist2$BMI_alone$resdf,
  
  depr_adjage_padj = dea2contrasts$contrastlist2$Condition_corrAge$resdf,
  age_adjdepr_padj = dea2contrasts$contrastlist2$Age_corrCond$resdf,
  
  depr_adjimc_padj = dea2contrasts$contrastlist2$Condition_corrIMC$resdf,
  imc_adjdepr_padj = dea2contrasts$contrastlist2$BMI_corrCond$resdf,
  
  depr_adj2 = dea2contrasts$contrastlist2$Condition_corr2$resdf,
  imc_adj2 = dea2contrasts$contrastlist2$BMI_corr2$resdf,
  age_adj2 = dea2contrasts$contrastlist2$Age_corr2$resdf
)

merged_pvals <- list2merge %>% 
  lapply(\(x){
    x$taxon <- x$taxon %>% 
      gsub("-", ".", .) %>% 
      gsub("\\[|\\]", "", .) %>% 
      gsub("\\._", "_", .) %>% 
      gsub("\\/|\\(|\\)", ".", .) %>% 
      gsub("\\.$", "", .) 
    x <- x[match(summary_df$variable, x$taxon), ]
    return(x$padj)
  }) %>% bind_cols()
summary_df <- cbind(summary_df, merged_pvals)

x_name <- "Edad"
y_name <- "Condition_bin"
mediator_name1 <- "IMC"
plim <- 0.05

# signvars <- summary_df %>% dplyr::select(-variable) %>% 
#   as.matrix() %>% apply(MAR=1, \(x)any(x<=plim)) 
# vars2test <- summary_df %>% dplyr::filter(signvars) %>% 
#   pull(variable)

vars2test <- summary_df %>% dplyr::filter(depr_adj2<plim & 
                                     ((imc_adjdepr_padj < plim ) |  #& imc_only_padj < plim
                                        (age_adjdepr_padj < plim ))) %>%  #& age_only_padj < plim
     pull(variable)


vars2venn <- list(
  "D vs C" = summary_df %>% dplyr::filter(depr_only_padj < plim) %>% pull(variable),
  "BMI adj." = summary_df %>% dplyr::filter(imc_adjdepr_padj < plim) %>% pull(variable),
  "Age adj." = summary_df %>% dplyr::filter(age_adjdepr_padj < plim) %>% pull(variable)
)

gv <- ggvenn(
  vars2venn, columns = names(vars2venn),
  stroke_size = 0.5,
  stroke_color = C_NS,
  fill_color = c(C_CASE, C_CASE2, C_CTRL2),show_elements = F
)
ggsave(filename = paste0(opt$out, "VennDiagram_AgeAndBMI.pdf"), gv)


param_names <- c('a', 'b', 'cp', 'd', 'e', 'fp', 'a*b', 'cp+a*b', 'd*b', 'fp+d*b', 'e*cp', 'fp+e*cp', 'd*b+fp+e*cp+e*a*b') 
medresultsBoth <- lapply(vars2test, \(mediator_name2, df_all){
  df <- df_all %>% dplyr::select(all_of(c(x_name, y_name, mediator_name1, mediator_name2)))
  with_nas <-  df %>% apply(MAR=1, \(x)any(is.na(x)))
  df <- df[!with_nas, ]
  
  res <- tryCatch({makeMediationComplex(df, xname=x_name, 
                                        yname=y_name, 
                                        medname=mediator_name1, 
                                        medname2=mediator_name2)}, 
                  error=\(x){
    
    xx <- data.frame(Estimate=rep(NA, length(param_names)), 
                     S.E. = rep(NA, length(param_names)), 
                     `Z-score`=rep(NA, length(param_names)), 
                     p.value=rep(NA, length(param_names)))
    rownames(xx) <- c(param_names)
    return(list(estimates=xx))   
  })
  res2 <- res$estimates %>% 
    rownames_to_column("param") %>% 
    dplyr::select(param, Estimate, p.value) %>% 
    gather(key="var", value="value", Estimate, p.value) %>% 
    dplyr::filter(param %in% param_names) %>% 
    unite("tmp", param, var, sep="_") %>% 
    spread(tmp, value) %>% 
    dplyr::mutate(Xvar = x_name, 
                  Yvar = y_name, 
                  Mediator=mediator_name1, 
                  Mediator2=mediator_name2,
                  fullmodel=list(res$estimates))
  return(res2)
  
}, df_all) %>% bind_rows() 

#all(medresults$Xvar %in% summary_df$variable)
medresultsBoth_merged <- merge(medresultsBoth, summary_df, by.x="Mediator2", by.y="variable")
write_tsv(medresultsBoth_merged, file=paste0(opt$out, "mediation_analysis_EdadAndIMC.tsv"))

## Transform to plot ##falta
def_effects <- c('a', 'b', 'cp', 'd', 'fp', 'e', 'a*b', 'cp+a*b', 
                 'cp+a*b', 'fp+d*b', 'e*cp', 'fp+e*cp', 'd*b+fp+e*cp+e*a*b') 
res_trans <- map(1:nrow(medresultsBoth), \(i){
  res_i <- medresultsBoth$fullmodel[[i]]
  #oldnames2 <- rownames(res_i) %>% 
  #  sapply(\(x)strsplit(x, "\\*|\\+",perl = TRUE)[[1]][1]) 
  resdf_i <- res_i %>% 
    rownames_to_column("param") %>% 
    mutate(x_labels = ifelse(param %in% def_effects, 
                             medresultsBoth$Mediator2[i],
                             gsub("V\\[|\\]", "", param)),
           param = ifelse(param %in% def_effects, 
                          gsub("(b|a|cp|d|e|fp)", paste0("\\1", as.character(i)), param, perl=T), param) 
    )
  
}) %>% bind_rows()
edbs_only <- res_trans %>% dplyr::filter(grepl("^(d|b|fp)[0-9]+$", param)) %>% 
  dplyr::mutate(x_labels=gsub("[0-9]+", "", param))
(edbs_hist <- ggplot(edbs_only,aes(x=Estimate))+
    facet_wrap(~x_labels, scales="free")+geom_histogram()+
    mytheme+ggtitle("Histogram of params"))
ggsave(filename = paste0(opt$out, "histogram_edbfs_AgeAndIMC_separate.pdf"), edbs_hist)

res_final <- rbind(
  edbs_only %>% 
    mutate(param = x_labels) %>% 
    group_by(param, x_labels) %>% 
    dplyr::summarise_if(is.numeric, mean) %>% 
    dplyr::select(all_of(names(res_trans))),
  res_trans %>% dplyr::filter(!grepl("^(d|b|fp)[0-9]+$", param)) %>% 
    dplyr::mutate(param=gsub("(d|b|fp)[0-9]+", "\\1", param))
) 
write_tsv(res_final, file=paste0(opt$out, "mediation_analysis_EdadAndIMC_separate_reformat.tsv"))
write_tsv(cps_only, file=paste0(opt$out, "mediation_analysis_EdadAndIMC_separate_cps.tsv"))


name<-"analysis_AgeAndIMC_separateModel_vjust"
plotres3 <- plotAgeAndIMCMediationComplex(res_final, vars2test, opt$out, name, plim_plot = 0.05, use_color_scale = FALSE, 
                                   w=14, h=10)

## Make also boxplot
bacnames <- plotres3$bacorder$x_labels
plots <- makePlotBySpecies(bacnames, df_all, opt$out, "AgeAndIMC_plotAge_separate_BySpeciesPearson", quantvar="Edad_log", 
                           quantvar_name = "log(Age)", corrmethod = "pearson", w=8, h=10)
plots <- makePlotBySpecies(bacnames, df_all, opt$out, "AgeAndIMC_plotAge_separate_BySpeciesSpearman", quantvar="Edad_log", 
                           quantvar_name = "log(Age)", corrmethod = "spearman", w=8, h=10)
plots <- makePlotBySpecies(bacnames, df_all, opt$out, "AgeAndIMC_plotAge_separate_BySpeciesKendall", quantvar="Edad_log", 
                           quantvar_name = "log(Age)", corrmethod = "kendall", w=8, h=10)
plots <- makePlotBySpecies(bacnames, df_all, opt$out, "AgeAndIMC_plotIMC_separate_BySpeciesPearson", quantvar="IMC_log", 
                           quantvar_name = "log(BMI)", corrmethod = "pearson", w=8, h=10)
plots <- makePlotBySpecies(bacnames, df_all, opt$out, "AgeAndIMC_plotIMC_separate_BySpeciesSpearman", quantvar="IMC_log", 
                           quantvar_name = "log(BMI)", corrmethod = "spearman", w=8, h=10)
plots <- makePlotBySpecies(bacnames, df_all, opt$out, "AgeAndIMC_plotIMC_separate_BySpeciesKendall", quantvar="IMC_log", 
                           quantvar_name = "log(BMI)", corrmethod = "kendall", w=8, h=10)
## 

### MEDIATION ANALYSIS WITH AGE AND IMC, MERGED

summary_df <- data.frame(variable = gsub("-", ".", vstdf$gene) %>% 
                           gsub("\\[|\\]", "", .) %>% 
                           gsub("\\._", "_", .) %>% 
                           gsub("\\/|\\(|\\)", ".", .) %>% 
                           gsub("\\.$", "", .) 
)

assertthat::assert_that(all(summary_df$variable %in% names(df_all)))

list2merge <- list(
  depr_only_padj = dea2contrasts$firstContrast$resdf,
  age_only_padj = dea2contrasts$contrastlist2$Age_alone$resdf,
  imc_only_padj = dea2contrasts$contrastlist2$BMI_alone$resdf,
  
  depr_adjage_padj = dea2contrasts$contrastlist2$Condition_corrAge$resdf,
  age_adjdepr_padj = dea2contrasts$contrastlist2$Age_corrCond$resdf,
  
  depr_adjimc_padj = dea2contrasts$contrastlist2$Condition_corrIMC$resdf,
  imc_adjdepr_padj = dea2contrasts$contrastlist2$BMI_corrCond$resdf,
  
  depr_adj2 = dea2contrasts$contrastlist2$Condition_corr2$resdf,
  imc_adj2 = dea2contrasts$contrastlist2$BMI_corr2$resdf,
  age_adj2 = dea2contrasts$contrastlist2$Age_corr2$resdf
)

merged_pvals <- list2merge %>% 
  lapply(\(x){
    x$taxon <- x$taxon %>% 
      gsub("-", ".", .) %>% 
      gsub("\\[|\\]", "", .) %>% 
      gsub("\\._", "_", .) %>% 
      gsub("\\/|\\(|\\)", ".", .) %>% 
      gsub("\\.$", "", .) 
    x <- x[match(summary_df$variable, x$taxon), ]
    return(x$padj)
  }) %>% bind_cols()
summary_df <- cbind(summary_df, merged_pvals)



x_name <- "Edad"
y_name <- "Condition_bin"
mediator_name1 <- "IMC"
plim <- 0.05
# 
# signvars <- summary_df %>% dplyr::select(-variable) %>% 
#   as.matrix() %>% apply(MAR=1, \(x)any(x<=plim)) 

vars2test <- summary_df %>% dplyr::filter(depr_adj2<plim & 
                                            ((imc_adjdepr_padj < plim ) | #& imc_only_padj < plim
                                               (age_adjdepr_padj < plim ))) %>%  #& age_only_padj < plim
  pull(variable)

df <- df_all %>% dplyr::select(all_of(c(x_name, y_name, mediator_name1, vars2test)))
with_nas <-  df %>% apply(MAR=1, \(x)any(is.na(x)))
df <- df[!with_nas, ]

res <- makeMediationComplex_mergedMEdiat2(df, xname=x_name, 
                            yname=y_name, 
                            medname=mediator_name1, 
                            mednames2=vars2test)

#all(medresults$Xvar %in% summary_df$variable)
write_tsv(res, file=paste0(opt$out, "mediation_analysis_EdadAndIMC_mergedModel.tsv"))

## PLOT
name<-"analysis_AgeAndIMC_mergedModel_vjust"
plotres4 <- plotAgeAndIMCMediationComplex(res, vars2test, opt$out, name, plim_plot = 0.05, use_color_scale = FALSE, 
                                          w=14, h=10)

## Make also boxplot
bacnames <- plotres4$bacorder$x_labels
plots <- makePlotBySpecies(bacnames, df_all, opt$out, "AgeAndIMC_plotAge_merged_BySpeciesPearson", quantvar="Edad_log", 
                           quantvar_name = "log(Age)", corrmethod = "pearson", w=8, h=10)
plots <- makePlotBySpecies(bacnames, df_all, opt$out, "AgeAndIMC_plotAge_merged_BySpeciesSpearman", quantvar="Edad_log", 
                           quantvar_name = "log(Age)", corrmethod = "spearman", w=8, h=10)
plots <- makePlotBySpecies(bacnames, df_all, opt$out, "AgeAndIMC_plotAge_merged_BySpeciesKendall", quantvar="Edad_log", 
                           quantvar_name = "log(Age)", corrmethod = "kendall", w=8, h=10)
plots <- makePlotBySpecies(bacnames, df_all, opt$out, "AgeAndIMC_plotIMC_merged_BySpeciesPearson", quantvar="IMC_log", 
                           quantvar_name = "log(BMI)", corrmethod = "pearson", w=8, h=10)
plots <- makePlotBySpecies(bacnames, df_all, opt$out, "AgeAndIMC_plotIMC_merged_BySpeciesSpearman", quantvar="IMC_log", 
                           quantvar_name = "log(BMI)", corrmethod = "spearman", w=8, h=10)
plots <- makePlotBySpecies(bacnames, df_all, opt$out, "AgeAndIMC_plotIMC_merged_BySpeciesKendall", quantvar="IMC_log", 
                           quantvar_name = "log(BMI)", corrmethod = "spearman", w=8, h=10)
