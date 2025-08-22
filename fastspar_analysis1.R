library(tidyverse)
library(pheatmap)

PLIM <- 0.01
put_back_names <- function(tab, otus){
  names(tab)[1] <- "taxon"
  assertthat::assert_that(all( as.character(tab$taxon) == names(tab)[2:ncol(tab)]))
  
  tab$taxon <- otus$taxonomy[match(as.character(tab$taxon), as.character(otus$ID))]
  names(tab)[2:ncol(tab)] <- tab$taxon
  return(tab)
}
  
meta <- read_tsv("/home/carlos/Documentos/CORALS/results_rstudio/results_Feb25_gam1/prepare_fastspar/METADATA_noTanda2.tsv")
otus <- read_tsv("/home/carlos/Documentos/CORALS/results_rstudio/results_Feb25_gam1/prepare_fastspar/OTUS_DF_filt_noTanda2.tsv")

fpall_corr <- read_tsv("/home/carlos/Documentos/CORALS/results_rstudio/results_Feb25_gam1/fastspar_results/CORALS_full1/median_correlation_splitAbn_OTUS_DF_filt_noTanda2_full.tsv")
fpall_corr <- put_back_names(fpall_corr, otus)
fpall_cov <- read_tsv("/home/carlos/Documentos/CORALS/results_rstudio/results_Feb25_gam1/fastspar_results/CORALS_full1/median_covariance_splitAbn_OTUS_DF_filt_noTanda2_full.tsv")
fpall_cov <- put_back_names(fpall_cov, otus)
fpall_p <- read_tsv("/home/carlos/Documentos/CORALS/results_rstudio/results_Feb25_gam1/fastspar_results/CORALS_full1/pvalues_splitAbn_OTUS_DF_filt_noTanda2_full.tsv")
fpall_p <- put_back_names(fpall_p, otus)

agesplit<- list.files("/home/carlos/Documentos/CORALS/results_rstudio/results_Feb25_gam1/fastspar_results/CORALS_ageclass1/", pattern = "tsv", full.names = T) %>% 
  map(read_tsv) %>% 
  map(put_back_names, otus=otus)
names(agesplit) <- list.files("/home/carlos/Documentos/CORALS/results_rstudio/results_Feb25_gam1/fastspar_results/CORALS_ageclass1/", pattern = "tsv", full.names = F) 
agelist <- names(agesplit) %>% map( \(x) {y <- strsplit(x, "_")[[1]]; return(y[length(y)])}) %>% gsub(".tsv", "", .)
typelist <- names(agesplit) %>% map( \(x) {y <- strsplit(x, "_")[[1]]; return(ifelse(y[1]=="pvalues", y[1], y[2]))}) %>% gsub(".tsv", "", .)


mat <- fpall_corr %>% column_to_rownames("taxon") %>% as.matrix
pmat <-  fpall_p %>% column_to_rownames("taxon") %>% as.matrix
assertthat::assert_that(all(rownames(pmat) == rownames(mat) ))
mat2 <- mat 
mat2[pmat>PLIM ] <- NA
hm <- pheatmap(mat)

mat2 <- mat2[hm$tree_row$order, hm$tree_col$order]
pheatmap(mat2, cluster_cols=F, cluster_rows = F)

## correlations of correlations and pvalues


cors <- map(unique(agelist), \(x){
  tab <-agesplit[[which(agelist == x & typelist == "correlation")]] %>% column_to_rownames("taxon") %>% as.matrix
  cor(mat %>% as.vector, tab %>% as.vector)
  
}) %>% unlist
names(cors) <- unique(agelist)
#3         4         5         6       gt7 
#0.7095323 0.8812139 0.8991866 0.8800905 0.7420784

allxall <- expand.grid(unique(agelist), unique(agelist)) %>% 
  mutate(cors = map2(Var1, Var2, \(x, y){
  tab <-agesplit[[which(agelist == x & typelist == "correlation")]] %>% column_to_rownames("taxon") %>% as.matrix
  tab2 <-agesplit[[which(agelist == y & typelist == "correlation")]] %>% column_to_rownames("taxon") %>% as.matrix
  cor(tab %>% as.vector, tab2 %>% as.vector)
  
  }) %>% unlist
  ) %>% spread(key=Var2, value=cors) %>% 
  column_to_rownames("Var1") %>% 
  as.matrix()

pheatmap(allxall)

matvec <- mat %>% as.vector
pmatvec <- pmat %>% as.vector
coincid <- map(unique(agelist), \(x){
  tab <-agesplit[[which(agelist == x & typelist == "correlation")]] %>% column_to_rownames("taxon") %>% 
    as.matrix %>% as.vector
  ptab <-agesplit[[which(agelist == x & typelist == "pvalues")]] %>% column_to_rownames("taxon") %>% 
    as.matrix %>% as.vector
  cc <- ifelse(ptab <= PLIM & matvec <= PLIM, 
               ifelse(matvec > 0 & tab > 0, "POS_POS", 
                      ifelse(matvec < 0 & tab < 0, "NEG_NEG", "SIG_Disc")), 
               ifelse(ptab > PLIM & matvec > PLIM, "NS_NS", 
                      ifelse(matvec > 0 & tab > 0, "NSdisc_POS_POS", 
                             ifelse(matvec < 0 & tab < 0, "NSdisc_NEG_NEG", "NSdisc_Disc"))
               )
               ) %>% 
    table %>% prop.table
  cc <- 100*cc
  yy <- as.vector(cc)
  names(yy) <- names(cc)
  return(yy)
}) %>% bind_rows() %>% 
  mutate(AgeGroup = unique(agelist))

write_tsv(coincid, file = paste0("/home/carlos/Documentos/CORALS/results_rstudio/results_Feb25_gam1/fastspar_results/", "Coincidence_AgeGroup_vs_All.tsv"))
