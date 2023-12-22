library(tidyverse)
library(janitor)
library(readxl)

SEED <- 123
MODE = "LOCAL"

if(MODE == "IATA"){
  opt <- list()
}else{
  opt <- list(out ="/home/carmoma/Desktop/202311_DEPRESION/results_rstudio_v2_4/mediation_analysis4/",
              indir = "/home/carmoma/Desktop/202311_DEPRESION/results_rstudio_v2_4/",
              phyloseq_list = "/home/carmoma/Desktop/202311_DEPRESION/results_rstudio_v2_4/phyloseq/phyloseq_all_list.RData",
              phyloseq_name = "remove_tanda2",
              r_functions="/home/carmoma/Desktop/202311_DEPRESION/depression_scripts/metagenomics_core_functions.R",
              r_functions_mediation="/home/carmoma/Desktop/202311_DEPRESION/depression_scripts/mediation_functions.R",
              metadata = "/home/carmoma/Desktop/202311_DEPRESION/metadatos_MC_AL12042023_CMcopy.xlsx",
              rewrite=TRUE,
              fc=1, 
              pval=0.05,
              adjust_pvals = TRUE
  )
}

hmdadfile <- "/home/carmoma/Desktop/202311_DEPRESION/microbiome_databases/HMDAD/data_download.txt"
hmdad <- read_tsv(hmdadfile) %>% dplyr::filter(Position == "Gastrointestinal tract")

load(paste0(opt$indir, "DESeq2_ControlVarsMany/LFC_Comparison_AgeAndBMI_allCombos.RData"))

calcScore <- function(name, ev, hmdad){
  if(! name %in% hmdad$Microbe) return(0)
  n <- hmdad %>% dplyr::filter(Microbe == name & Evidence == ev) %>% pull(PMID) %>% unique %>% length
  return(n)
}


daatab <- dea2contrasts$firstContrast$resdf %>% dplyr::mutate(taxon=gsub("_", " ", taxon), 
                                                        taxon = gsub("[\\[\\]]", "", taxon)) %>% 
  dplyr::mutate(goodScore = map_vec(taxon, calcScore, "Decrease", hmdad),
                badScore = map_vec(taxon, calcScore, "Increase", hmdad),
                psig = padj < 0.05,
                direction = ifelse(padj > 0.05, "NS", ifelse(log2FoldChangeShrink > 0, "Up", "Down")))

t1 <- daatab %>% dplyr::filter(psig) %>% group_by(direction) %>%  dplyr::summarise(meanGood = mean(goodScore, na.rm = T),
                                                                                   meanBad = mean(badScore, na.rm = T))
daafilt <- daatab %>% dplyr::filter(psig)

conttab_bad <- c( 
  length(which(daafilt$direction == "Up" & daafilt$badScore>0)), 
  length(which(daafilt$direction == "Up" & daafilt$badScore==0)),
  length(which(daafilt$direction == "Down" & daafilt$badScore>0)), 
  length(which(daafilt$direction == "Down" & daafilt$badScore==0))
) %>% matrix(nrow=2, byrow=T)
chisq.test(conttab_bad)
fisher.test(conttab_bad)


###
positive_strings <- c("Microflora promote disease's progression",
                      "Positive correlation", "Positive correlation due to antibiotic", "Unassociated", "Associated")
amadis <- read_excel("/home/carmoma/Desktop/202311_DEPRESION/microbiome_databases/Amadis/GIFTED.xlsx") %>% 
  dplyr::filter(Organ %in% c("Liver and Intestines", "Intestines", "Large Intestines", "Brain")) %>% 
  dplyr::mutate(Evidence = ifelse(amadis$`Association between disease and microflora` %in% positive_strings, "Increase", "Decrease"),
                Microbe = Flora,
                PMID = PubmedID)

table(daatab$taxon %in% amadis$Flora)

daatab <- daatab %>% 
  dplyr::mutate(goodScore_amadis = map_vec(taxon, calcScore, "Decrease", amadis),
                badScore_amadis = map_vec(taxon, calcScore, "Increase", amadis),
                psig = padj < 0.05,
                direction = ifelse(padj > 0.05, "NS", ifelse(log2FoldChangeShrink > 0, "Up", "Down")))

t2 <- daatab %>% dplyr::filter(psig) %>% group_by(direction) %>%  dplyr::summarise(meanGood = mean(goodScore_amadis, na.rm = T),
                                                                                   meanBad = mean(badScore_amadis, na.rm = T))
daafilt <- daatab %>% dplyr::filter(psig)

conttab_bad2 <- c( 
  length(which(daafilt$direction == "Up" & daafilt$badScore_amadis>0)), 
  length(which(daafilt$direction == "Up" & daafilt$badScore_amadis==0)),
  length(which(daafilt$direction == "Down" & daafilt$badScore_amadis>0)), 
  length(which(daafilt$direction == "Down" & daafilt$badScore_amadis==0))
) %>% matrix(nrow=2, byrow=T)
conttab_good2 <- c( 
  length(which(daafilt$direction == "Up" & daafilt$goodScore_amadis>0)), 
  length(which(daafilt$direction == "Up" & daafilt$goodScore_amadis==0)),
  length(which(daafilt$direction == "Down" & daafilt$goodScore_amadis>0)), 
  length(which(daafilt$direction == "Down" & daafilt$goodScore_amadis==0))
) %>% matrix(nrow=2, byrow=T)
chisq.test(conttab_bad2)
fisher.test(conttab_bad2)
chisq.test(conttab_good2)
fisher.test(conttab_good2)

####

daatab <- daatab %>% mutate(goodScore2 = goodScore + goodScore_amadis, badScore2 = badScore + badScore_amadis)
daafilt <- daatab %>% dplyr::filter(psig)


conttab_bad2 <- c( 
  length(which(daafilt$direction == "Up" & daafilt$badScore2>0)), 
  length(which(daafilt$direction == "Up" & daafilt$badScore2==0)),
  length(which(daafilt$direction == "Down" & daafilt$badScore2>0)), 
  length(which(daafilt$direction == "Down" & daafilt$badScore2==0))
) %>% matrix(nrow=2, byrow=T)
conttab_good2 <- c( 
  length(which(daafilt$direction == "Up" & daafilt$goodScore2>0)), 
  length(which(daafilt$direction == "Up" & daafilt$goodScore2==0)),
  length(which(daafilt$direction == "Down" & daafilt$goodScore2>0)), 
  length(which(daafilt$direction == "Down" & daafilt$goodScore2==0))
) %>% matrix(nrow=2, byrow=T)
chisq.test(conttab_bad2)
fisher.test(conttab_bad2)
chisq.test(conttab_good2)
fisher.test(conttab_good2)

t.test(badScore2 ~ direction, data=daafilt)
wilcox.test(badScore2 ~ direction, data=daafilt)
