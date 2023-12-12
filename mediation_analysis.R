library(tidyverse)
library(bmem)
library(DESeq2)
library(phyloseq)

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
  opt <- list(out ="/home/carmoma/Desktop/202311_DEPRESION/results_rstudio_v2_4/mediation_analysis",
              indir = "/home/carmoma/Desktop/202311_DEPRESION/results_rstudio_v2_4/",
              phyloseq_list = "/home/carmoma/Desktop/202311_DEPRESION/results_rstudio_v2_4/phyloseq/phyloseq_all_list.RData",
              phyloseq_name = "remove_tanda2",
              r_functions="/home/carmoma/Desktop/202311_DEPRESION/depression_scripts/metagenomics_core_functions.R",
              metadata = "/home/carmoma/Desktop/202311_DEPRESION/metadatos_MC_AL12042023_CMcopy.xlsx",
              rewrite=TRUE,
              fc=1, 
              pval=0.05
  )
}
if(! dir.exists(opt$out)){dir.create(opt$out)}

### LOAD DATA
source(opt$r_functions)
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

makeMediationSimple <- function(df, xname, yname, medname){
  library(bmem)
  library(sem)
  
  eq1 = paste0(yname, " = b*", medname, " + cp*", xname)
  eq2 = paste0(medname, " = a*", xname)
  eqs <- paste(eq1, eq2, sep="\n", collapse="\n")
  
  iris.model<-specifyEquations(exog.variances=T, text=eqs)

  effects<-c('a*b', 'cp+a*b') 
  nlsy.res<-bmem.sobel(df, iris.model, effects) 
  return(nlsy.res)

}

makeMediationSimple_mergedX <- function(df, xnames, yname, medname){
  library(bmem)
  library(sem)
  
  a_params <- paste("a", 1:length(xnames), sep="")
  cp_params <- paste("cp", 1:length(xnames), sep="")
  
  eq1 = paste(yname, " = b*", medname, " + ", cp_params, "*", xnames, sep="") %>% 
    paste(collapse = "\n")
  eq2 = paste0(medname, " = ", a_params, "*", xnames) %>% 
    paste(collapse = "\n")
  eqs <- paste(eq1, eq2, sep="\n", collapse="\n")
  
  iris.model<-specifyEquations(exog.variances=T, text=eqs)
  
  effects<-c(paste(a_params, "*b", sep="" ), paste(cp_params, '+', a_params,'*b', sep="")) 
  nlsy.res<-bmem.sobel(df, iris.model, effects) 
  
  oldnames2 <- rownames(nlsy.res$estimates) %>% 
    sapply(\(x)strsplit(x, "\\*|\\+",perl = TRUE)[[1]][1]) 
  # newnames <- ifelse(oldnames %in% a_params, 
  #        paste(oldnames, xnames[match(oldnames, a_params)], sep="_"),
  #        ifelse(oldnames %in% cp_params, 
  #               paste(oldnames, xnames[match(oldnames, cp_params)], sep="_"), 
  #               oldnames))
  resdf <- nlsy.res$estimate %>% 
    rownames_to_column("param")
  resdf$x_labels <- ifelse(oldnames2 %in% a_params, 
                     xnames[match(oldnames2, a_params)],
                     ifelse(oldnames2 %in% cp_params, 
                            xnames[match(oldnames2, cp_params)], 
                            gsub("V\\[|\\]", "", oldnames2)))
  
  #rownames(nlsy.res$estimates) <- newnames
  return(resdf)
}

makeMediationSimple_mergedMediat <- function(df, x_name, yname, mednames){
  library(bmem)
  library(sem)
  
  a_params <- paste("a", 1:length(mednames), sep="")
  b_params <- paste("b", 1:length(mednames), sep="")
  
  eq1 = paste(yname, " = ", b_params, "*", mednames, " + ", "cp*", x_name, sep="") %>% 
    paste(collapse = "\n")
  eq2 = paste0(mednames, " = ", a_params, "*", x_name) %>% 
    paste(collapse = "\n")
  eqs <- paste(eq1, eq2, sep="\n", collapse="\n")
  
  iris.model<-specifyEquations(exog.variances=T, text=eqs)
  
  effects<-c(paste(a_params, "*", b_params, sep="" ), paste('cp+', a_params,'*', b_params, sep="")) 
  nlsy.res<-bmem.sobel(df, iris.model, effects) 
  
  oldnames2 <- rownames(nlsy.res$estimates) %>% 
    sapply(\(x){y<-strsplit(x, "\\*|\\+",perl = TRUE)[[1]];y[length(y)]}) 
  # newnames <- ifelse(oldnames %in% a_params, 
  #        paste(oldnames, xnames[match(oldnames, a_params)], sep="_"),
  #        ifelse(oldnames %in% cp_params, 
  #               paste(oldnames, xnames[match(oldnames, cp_params)], sep="_"), 
  #               oldnames))
  resdf <- nlsy.res$estimate %>% 
    rownames_to_column("param")
  resdf$x_labels <- ifelse(oldnames2 %in% a_params, 
                           xnames[match(oldnames2, a_params)],
                           ifelse(oldnames2 %in% b_params, 
                                  xnames[match(oldnames2, b_params)], 
                                  gsub("V\\[|\\]", "", oldnames2)))
  
  #rownames(nlsy.res$estimates) <- newnames
  return(resdf)
}

makeMediationComplex <- function(df, xname="Edad", yname="Condition_bin", 
         medname="IMC", medname2="Faecalibacterium_prausnitzii"){
  library(bmem)
  library(sem)
  
  eq1 = paste0(yname, " = b*", medname, " + cp*", medname2, " + fp*", xname)
  eq2 = paste0(medname, " = a*", medname2, " + d*", xname)
  eq3 = paste0(medname2, " = e*", xname)
  eqs <- paste(eq1, eq2, eq3, sep="\n", collapse="\n")
  
  iris.model<-specifyEquations(exog.variances=T, text=eqs)
  
  effects<-c('a*b', 'cp+a*b', 'd*b', 'fp+d*b', 'e*cp', 'fp+e*cp', 'd*b+fp+e*cp+a*b') 
  nlsy.res<-bmem.sobel(df, iris.model, effects) 
  return(nlsy.res)
  
}

makeMediationComplex_mergedMEdiat2 <- function(df, xname="Edad", 
                                               yname="Condition_bin", 
                                               medname="IMC", 
                                               mednames2=c("Faecalibacterium_prausnitzii")){
  library(bmem)
  library(sem)
  
  a_params <- paste("a", 1:length(mednames2), sep="")
  cp_params <- paste("cp", 1:length(mednames2), sep="")
  e_params <- paste("e", 1:length(mednames2), sep="")
  
  eq1 = paste(yname, " = b*", medname, " + ", cp_params, "*", mednames2, " + fp*", xname, sep="") %>% 
    paste(collapse = "\n")
  eq2 = paste0(medname, " = ", a_params, "*", mednames2, " + d*", xname) %>% 
    paste(collapse = "\n")
  eq3 = paste0(mednames2, " = ", e_params, "*", xname)%>% 
    paste(collapse = "\n")
  eqs <- paste(eq1, eq2, eq3, sep="\n", collapse="\n")
  
  iris.model<-specifyEquations(exog.variances=T, text=eqs)
  
  effects<-c(paste(a_params, "*b", sep="" ), 
             paste(cp_params, '+', a_params,'*b', sep=""),
             'fp+d*b',
             paste(e_params, '*', cp_params, sep=""),
             paste('fp + ', e_params, '*', cp_params, sep=""),
             paste('d*b + fp + ', e_params, '*', cp_params, ' + ', a_params, "*b", sep="")
             ) 
  nlsy.res<-bmem.sobel(df, iris.model, effects) 
  
  resdf <- nlsy.res$estimate %>% 
    rownames_to_column("param")
  resdf$x_labels <- resdf$param %>% 
    sapply(\(x){
      y <-strsplit(x, "\\*|\\+",perl = TRUE)[[1]] %>% gsub(" ", "", .) %>% sort 
      
      if(any(stringr::str_starts(y, "a"))){
        pthis <- y[stringr::str_starts(y, "a")]
        return(paste(mednames2[a_params==pthis], collapse=":"))
      }else if(any(stringr::str_starts(y, "cp"))){
        pthis <- y[stringr::str_starts(y, "cp")]
        return(paste(mednames2[cp_params==pthis], collapse=":"))
      }else if(any(stringr::str_starts(y, "e"))){
        pthis <- y[stringr::str_starts(y, "e")]
        return(paste(mednames2[e_params==pthis], collapse=":"))
      }else{
        return(gsub("V\\[|\\]", "", x))
      }
      
    }) 
  return(resdf)
}

#### Read data 

expr_df <- vstdf %>% column_to_rownames("gene") %>% 
  as.matrix() %>% t %>% data.frame %>% 
  rownames_to_column("sampleID")

df_all <- merge(metadata2, expr_df, by="sampleID")
df_all$Condition_bin <- ifelse(df_all$Condition=="Depression", TRUE, FALSE)
names(df_all) <- gsub("^X\\.", "", names(df_all)) %>% 
                 gsub("\\._", "_", .) %>% 
                 gsub("\\.$", "", .)


write_tsv(df_all, file=paste0(opt$out, "merged_data.tsv"))

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
  dplyr::filter(depr_only_padj<plim | imc_only_padj < plim | depr_adjimc_padj < plim) %>% 
  pull(variable)



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
    mutate(Xvar = x, Yvar = y_name, Mediator=mediator_name)
  return(res2)
  
}, df_all) %>% bind_rows()
#all(medresults$Xvar %in% summary_df$variable)
medresults_merged <- merge(medresults, summary_df, by.x="Xvar", by.y="variable")
write_tsv(medresults_merged, file=paste0(opt$out, "mediation_analysis_IMC.tsv"))

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
  dplyr::filter(depr_only_padj<plim & imc_only_padj < plim & imc_adjdepr_padj < plim) %>% 
  pull(variable)
length(vars2test)

df <- df_all %>% dplyr::select(all_of(c(vars2test, y_name, mediator_name)))
with_nas <-  df %>% apply(MAR=1, \(x)any(is.na(x)))
df <- df[!with_nas, ]
res <- makeMediationSimple_mergedX(df, vars2test, y_name, mediator_name)
write_tsv(res, file=paste0(opt$out, "mediation_analysis_IMC_mergedModel.tsv"))


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
  dplyr::filter(depr_only_padj < plim | age_only_padj < plim | depr_adjage_padj < plim) %>% 
  pull(variable)

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
    dplyr::mutate(Xvar = x_name, Yvar = y_name, Mediator=mediator_name)
  return(res2)
  
}, df_all) %>% bind_rows() 
#all(medresults$Xvar %in% summary_df$variable)

medresultsAge_merged <- merge(medresultsAge, summary_df, by.x="Mediator", by.y="variable")
write_tsv(medresultsAge_merged, file=paste0(opt$out, "mediation_analysis_Edad.tsv"))

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
  dplyr::filter(depr_only_padj < plim  & age_adjdepr_padj < plim) %>% 
  pull(variable)

df <- df_all %>% dplyr::select(all_of(c(x_name, y_name, vars2test)))
with_nas <-  df %>% apply(MAR=1, \(x)any(is.na(x)))
df <- df[!with_nas, ]
res <- makeMediationSimple_mergedMediat(df, x_name, y_name, vars2test)
write_tsv(medresultsAge_merged, file=paste0(opt$out, "mediation_analysis_Edad_mergedModel.tsv"))

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

summary_df %>% filter(depr_adj2<plim & 
                        imc_adj2<plim & 
                        age_adj2<plim) ## just curious

x_name <- "Edad"
y_name <- "Condition_bin"
mediator_name1 <- "IMC"
plim <- 0.05

signvars <- summary_df %>% dplyr::select(-variable) %>% 
  as.matrix() %>% apply(MAR=1, \(x)any(x<=plim)) 
vars2test <- summary_df %>% dplyr::filter(signvars) %>% 
  pull(variable)

medresultsBoth <- lapply(vars2test, \(mediator_name2, df_all){
  df <- df_all %>% dplyr::select(all_of(c(x_name, y_name, mediator_name1, mediator_name2)))
  with_nas <-  df %>% apply(MAR=1, \(x)any(is.na(x)))
  df <- df[!with_nas, ]
  param_names <- c('a', 'b', 'cp', 'd', 'e', 'fp', 'a*b', 'cp+a*b', 'd*b', 'fp+d*b', 'e*cp', 'fp+e*cp', 'd*b+fp+e*cp+a*b') 
  
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
                  Mediator2=mediator_name2)
  return(res2)
  
}, df_all) %>% bind_rows() 
#all(medresults$Xvar %in% summary_df$variable)
medresultsBoth_merged <- merge(medresultsBoth, summary_df, by.x="Mediator2", by.y="variable")
write_tsv(medresultsBoth_merged, file=paste0(opt$out, "mediation_analysis_EdadAndIMC.tsv"))



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

vars2test <-summary_df %>% filter(depr_only_padj<plim & 
                                    imc_adjdepr_padj<plim ) %>% pull(variable)

df <- df_all %>% dplyr::select(all_of(c(x_name, y_name, mediator_name1, vars2test)))
with_nas <-  df %>% apply(MAR=1, \(x)any(is.na(x)))
df <- df[!with_nas, ]

res <- makeMediationComplex_mergedMEdiat2(df, xname=x_name, 
                            yname=y_name, 
                            medname=mediator_name1, 
                            mednames2=vars2test)

#all(medresults$Xvar %in% summary_df$variable)
write_tsv(res, file=paste0(opt$out, "mediation_analysis_EdadAndIMC_mergedModel.tsv"))


