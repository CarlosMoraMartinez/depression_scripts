library(ANCOMBC)
library(tidyverse)
library(wesanderson)
library(phyloseq)
# use this! use without log

pal <- wes_palette("AsteroidCity1", 2)
outdir <- "/home/carlos/Escritorio/202311_DEPRESION/ReviewJune2025/Results_rstudio/results2/DAA_ANCOMBC2/"
if(!dir.exists(outdir)) dir.create(outdir)

#opt <- restaurar(opt)
#load("/home/carlos/Escritorio/202311_DEPRESION/ReviewJune2025/Results_rstudio/results2/phyloseq_original//phyloseq_all_list.RData")
load("/home/carlos/Escritorio/202311_DEPRESION/ReviewJune2025/Results_rstudio/results2/DESeq2_AgeSexInteraction/Integrate1/phyloseq_used_remove_tanda2.RData")


anres_cond <- ancombc2(phobj, fix_formula = "Condition", p_adj_method = "BH")
save(anres_cond, file = paste0(outdir, "ancombc2_res_Condition.RData"))
write_tsv(anres_cond$res, file = paste0(outdir, "ancombc2_res_Condition.tsv"))

anres_bmi <- ancombc2(phobj, fix_formula = "BMI_log", p_adj_method = "BH")
save(anres_bmi, file = paste0(outdir, "ancombc2_res_BMI.RData"))
write_tsv(anres_bmi$res, file = paste0(outdir, "ancombc2_res_BMI.tsv"))

#interaction doesn't work
anres_cond_adj1 <- ancombc2(phobj, fix_formula = "Condition+Sex*Age_log", p_adj_method = "BH")
anres_bmi_adj1 <- ancombc2(phobj, fix_formula = "BMI_log + Sex * Age_log", p_adj_method = "BH")
anres_adj_all <- ancombc2(phobj, fix_formula = "Condition + BMI_log + Sex * Age_log", p_adj_method = "BH")


anres_adj_all <- ancombc2(phobj, fix_formula = "Condition + BMI_log + Sex + Age_log", p_adj_method = "BH")
save(anres_adj_all, file = paste0(outdir, "ancombc2_res_adjBMISexAge.RData"))
write_tsv(anres_adj_all$res, file = paste0(outdir, "ancombc2_res_adjBMISexAge.tsv"))

adjres <- anres_adj_all$res
fromMed <- c(
  "Blautia_hansenii",
  "Actinomyces_naeslundii",
  "Veillonella_nakazawae",
  "Ligilactobacillus_ruminis",
  "[Clostridium]_innocuum",
  "Anaerostipes_caccae",
  "Enterocloster_clostridioformis",
  "Streptococcus_parasanguinis",
  "Streptococcus_salivarius",
  "Streptococcus_sp._LPB0220",
  "Streptococcus_sp._FDAARGOS_192",
  "Streptococcus_lactarius",
  "Streptococcus_sp._HSISM1",
  "Streptococcus_gordonii",
  "Streptococcus_vestibularis",
  "Lactobacillus_gasseri"
)
fromMed[!fromMed %in% adjres$taxon]

fromMedx <- adjres %>% filter(taxon %in% fromMed) %>% 
  select(taxon, p_ConditionDepression, p_BMI_log)

write_tsv(fromMedx, file = paste0(outdir, "ancombc2_res_adjBMISexAge_filt16SpeciesFromMediation.tsv"))

anres_adj_all$res %>% filter(grepl("prausnit", taxon))
