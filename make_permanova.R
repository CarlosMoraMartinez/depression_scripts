
# Note:
# z_t0, z_t1, status_t0, status_t1 --> calculated by us
# z_bmi_00, z_bmi_01, Category_T0, Category_T1, z_waist_00, z_waist_01 --> calculated by collaborators

#all_phyloseq:
load("/home/carlos/Documentos/CORALS/results_rstudio/results_Abril25_2/foodPCA/phyloseq_list_foodPCA_withNMF.RData")

for(ph in names(all_phyloseq)){
 sample_data(all_phyloseq[[ph]])$inc_z_waist <- sample_data(all_phyloseq[[ph]])$z_waist_01 - sample_data(all_phyloseq[[ph]])$z_waist_00
 sample_data(all_phyloseq[[ph]])$inc_z_bmiCOR <- sample_data(all_phyloseq[[ph]])$z_bmi_01 - sample_data(all_phyloseq[[ph]])$z_bmi_00
 sample_data(all_phyloseq[[ph]])$inc_z_bmi <- sample_data(all_phyloseq[[ph]])$z_t1 - sample_data(all_phyloseq[[ph]])$z_t0
}
save(all_phyloseq, file = "/home/carlos/Documentos/CORALS/results_rstudio/results_Agosto25_1//foodPCA/phyloseq_list_foodPCA_withNMF_withInc.RData")

load("/home/carlos/Documentos/CORALS/results_rstudio/results_Agosto25_1//foodPCA/phyloseq_list_foodPCA_withNMF_withInc.RData")

exclude_vars <- c("sampleID", "sampleID2", "tanda")
dtypes <- c("bray")
ph2use <- names(all_phyloseq)[5:4]

food_variables <- names(sample_data(all_phyloseq$remove_tanda2 ))[5:30]
outdir <- paste0(opt$out, "PERMANOVA_dietPatterns/")
if(!dir.exists(outdir)) dir.create(outdir)
permaresults <- list()
for(i in ph2use){
  cat("Doing PERMANOVA of: ", i)
  phobj <- all_phyloseq[[i]]
  permaresults[[i]] <- lapply(dtypes, FUN=function(dd, phobj, exclude_vars, SEED){
    oname <- paste0(outdir, "permanova_results2_", i, "_", dd, ".tsv")
    makePermanova(phobj,dist_method = dd,
                  seed = SEED,
                  exclude_vars = exclude_vars,
                  outname = oname)
  }, phobj, exclude_vars, SEED)
  names(permaresults[[i]]) <- dtypes
}

save(permaresults, file = paste0(outdir, "permanova_results_bray_250814_2.RData"))
write_tsv(permaresults$remove_tanda2_rarefied_min$bray,
          file = paste0(outdir, "permanova_results_remove_tanda2_rarefied_min_bray_250814_2.tsv"))
write_tsv(permaresults$remove_tanda2$bray,
          file = paste0(outdir, "permanova_results_remove_tanda2_bray_250814_2.tsv"))
#xx <- read_tsv(paste0(outdir, "permanova_results_remove_tanda2_rarefied_min_bray.tsv"))
#yy <- read_tsv(paste0(outdir, "permanova_results2_remove_tanda2_rarefied_min_bray.tsv"))
#
#xx <- rbind(xx, yy)
#write_tsv(xx, file=paste0(outdir, "permanova_results2_remove_tanda2_rarefied_min_bray_full.tsv"))

# xx <- read_tsv(paste0(outdir, "permanova_results2_remove_tanda2_rarefied_min_bray_full.tsv")) # from April Results
xx <- read_tsv(paste0(outdir, "permanova_results_remove_tanda2_rarefied_min_bray_250814_2.tsv"))
##


##### only Normal at T0

permaresults <- list()
include_vars <- c("z_t0", "z_t1", "inc_z_bmi",
                  "z_waist_00", "z_waist_01", "inc_z_waist",
                  "status_c2")
exclude_vars <- names(sample_data(all_phyloseq$remove_tanda2_rarefied_min))
assertthat::assert_that(all(include_vars %in% exclude_vars))
exclude_vars <- exclude_vars[! exclude_vars %in% include_vars]

permaresults_normalT0 <- list()
for(i in ph2use){
  cat("Doing PERMANOVA of: ", i)
  phobj <- all_phyloseq[[i]]
  normalsamples <- sample_data(phobj) %>% data.frame() %>%
    filter(status_c1 == "normal")
  phobj_normt0 <- phyloseq::prune_samples(normalsamples$sampleID, phobj)

  permaresults_normalT0[[i]] <- lapply(dtypes, FUN=function(dd, phobj_normt0, exclude_vars, SEED){
    oname <- paste0(outdir, "permanova_results2_onlyNormalT0", i, "_", dd, ".tsv")
    makePermanova(phobj_normt0,
                  dist_method = dd,
                  seed = SEED,
                  exclude_vars = exclude_vars,
                  outname = oname)
  }, phobj_normt0, exclude_vars, SEED)
  names(permaresults_normalT0[[i]]) <- dtypes
}

save(permaresults_normalT0, file = paste0(outdir, "permanova_results_bray_onlyNormalT0_250814_2.RData"))
write_tsv(permaresults_normalT0$remove_tanda2_rarefied_min$bray,
          file = paste0(outdir, "permanova_results_remove_tanda2_rarefied_min_bray_normalT0_250814_2.tsv"))
write_tsv(permaresults_normalT0$remove_tanda2$bray,
          file = paste0(outdir, "permanova_results_remove_tanda2_bray_normalT0_250814_2.tsv"))

### ADONIS with multiple variables

#permaformulas <- c(
#  "braydist ~ edad_00 + Sex + hospital + nreads_filt + Category_T0 + educ_m_discrete",
#  "braydist ~ edad_00 + Sex + hospital + nreads_filt + Category_T0 + z_cintura_00 + educ_m_discrete",
#  "braydist ~ edad_00 + Sex + hospital + nreads_filt + z_imc_00 + z_cintura_00 + educ_m_discrete",
#  "braydist ~ edad_00 + Sex + hospital + nreads_filt + Category_T0 + z_cintura_00 + educ_m_discrete + lacteos_00 + ffq_h_carb_00 + grasas_00 + refresc_00 + ffq_energia_00 + legum_00 + frutas_00 + ffq_prot_00 + cereref_00",
#  "braydist ~ edad_00 + Sex + hospital + nreads_filt + z_imc_00 + z_cintura_00 + educ_m_discrete +  lacteos_00 + ffq_h_carb_00 + grasas_00 + refresc_00 + ffq_energia_00 + legum_00 + frutas_00 + ffq_prot_00 + cereref_00",
#  "braydist ~ PC2 + PC4 + PC11 + PC6 + PC1 + PC9 + PC12 + PC21",
#  "braydist ~ PC2 + PC4",
#  "braydist ~ PC2 + PC4 + PC1+ PC11",
#  "braydist ~ edad_00 + PC2 + PC4 + PC11 + PC6 + PC1 + PC9 + PC12 + PC21",
#  "braydist ~ edad_00 + PC2 + PC4",
#  "braydist ~ edad_00 + PC2 + PC4 + PC1+ PC11",
#  "braydist ~ edad_00 + Sex + hospital + nreads_filt + Category_T0 + educ_m_discrete + PC2 + PC4 + PC11 + PC6 + PC1 + PC9 + PC12 + PC21",
#  "braydist ~ edad_00 + Sex + hospital + nreads_filt + z_imc_00 + z_cintura_00+ educ_m_discrete + PC2 + PC4",
#  "braydist ~ edad_00 + Sex + hospital + nreads_filt + Category_T0 + educ_m_discrete + PC2 + PC4 + PC1+ PC11",
#  "braydist ~ edad_00 + Sex + hospital + nreads_filt + z_imc_00 + z_cintura_00 + educ_m_discrete + PC2 + PC4 + PC1+ PC11",
#)

outdir <- paste0(opt$out, "PERMANOVA_dietPatterns/")
if(!dir.exists(outdir)) dir.create(outdir)

s_meta <- sample_data(all_phyloseq$remove_tanda2_rarefied_min) %>% data.frame
vars2test1 <- c(#"status_c1",
               # "status_c2",
                "Sex",
               # "mg_p_00",
               "hospital",
               "mother_educ",
                "age_T0",
               "z_t0",
               #"z_bmi_00",
               #"z_bmi_01",
               "nreads",
               "z_waist_00"
               #"z_waist_01",
               #"inc_z_waist",s
               #"inc_z_bmi"
               ) # "exercise_cat", # too many NAs "age_class1",
quant_vars_diet_clr <- c(names(s_meta)[5], names(s_meta)[grep("_clr$", names(s_meta))])
pcvars_clr <- names(s_meta)[grep("^clr_PC", names(s_meta))]
patnames <- c("Preprocessed", "Mediterranean", "Western")

#xx <- read_tsv(paste0(outdir, "permanova_results_remove_tanda2_rarefied_min_bray.tsv"))
#xx <- read_tsv(paste0(outdir, "permanova_results2_remove_tanda2_rarefied_min_bray_full.tsv"))
xx <- read_tsv(paste0(outdir, "permanova_results_remove_tanda2_rarefied_min_bray_250814_2.tsv"))

sig_food_vars <- xx %>% filter(variable %in% quant_vars_diet_clr) %>% arrange(padj) %>% filter(padj<0.05) %>% pull(variable)
sig_pop_vars <- xx %>% filter(variable %in% vars2test1) %>% arrange(padj) %>% filter(padj<0.05) %>% pull(variable)
sig_PC_clrs <-  xx %>% filter(grepl("^clr_PC", variable, perl=T)) %>% arrange(padj) %>% filter(padj<0.05) %>% pull(variable)


permaformulas_food <- paste0("braydist ~ ", paste(quant_vars_diet_clr, sep=" + ", collapse=" + "))
permaformulas_pop <- paste0("braydist ~ ", paste(vars2test1, sep=" + ", collapse=" + "))
permaformulas_PCclr <- paste0("braydist ~ ", paste(pcvars_clr, sep=" + ", collapse=" + "))
permaformulas_pop_food <- paste0("braydist ~ ", paste(c(vars2test1, quant_vars_diet_clr), sep=" + ", collapse=" + "))
permaformulas_pop_PCclr <- paste0("braydist ~ ", paste(c(vars2test1, pcvars_clr), sep=" + ", collapse=" + "))

permaformulas_patterns <- paste0("braydist ~ ", paste(patnames, sep=" + ", collapse=" + "))
permaformulas_pop_patterns <- paste0("braydist ~ ", paste(c(vars2test1, patnames), sep=" + ", collapse=" + "))

permaformulas_food_sig <- paste0("braydist ~ ", paste(sig_food_vars, sep=" + ", collapse=" + "))
permaformulas_pop_sig <- paste0("braydist ~ ", paste(sig_pop_vars, sep=" + ", collapse=" + "))
permaformulas_PCclr_sig <- paste0("braydist ~ ", paste(sig_PC_clrs, sep=" + ", collapse=" + "))
permaformulas_pop_food_sig <- paste0("braydist ~ ", paste(c(sig_pop_vars, sig_food_vars), sep=" + ", collapse=" + "))
permaformulas_pop_PCclr_sig <- paste0("braydist ~ ", paste(c(sig_pop_vars, sig_PC_clrs), sep=" + ", collapse=" + "))

make_formulas_excluding1 <- function(varvec){
  res <- map_vec(varvec, \(x){ paste0("braydist ~ ", paste(varvec[varvec != x], sep=" + ", collapse=" + ")) })
  names(res) <- varvec
  return(res)
}
allbut1_all_food <- make_formulas_excluding1(c(vars2test1, quant_vars_diet_clr))
names(allbut1_all_food) <- paste("AllBut1", names(allbut1_all_food), sep="_")
allbut1_sig_food <- make_formulas_excluding1(c(sig_pop_vars, sig_food_vars))
names(allbut1_sig_food) <- paste("AllBut1Sig", names(allbut1_sig_food), sep="_")
allbut1_all_pcs <- make_formulas_excluding1(c(vars2test1, pcvars_clr))
names(allbut1_all_pcs) <- paste("AllBut1PCs", names(allbut1_all_pcs), sep="_")
allbut1_sig_pcs <- make_formulas_excluding1(c(sig_pop_vars, sig_PC_clrs))
names(allbut1_sig_pcs) <- paste("AllBut1SigPCs", names(allbut1_sig_pcs), sep="_")

permaformulas <- c(
  allfood=permaformulas_food,
  allpop=permaformulas_pop,
  PCclr=permaformulas_PCclr,
  PopFood=permaformulas_pop_food,
  PopPCclr=permaformulas_pop_PCclr,

  FoodSig=permaformulas_food_sig,
  PopSig=permaformulas_pop_sig,
  PCclSig=permaformulas_PCclr_sig,
  PopFoodSig=permaformulas_pop_food_sig,
  PopPCclrSig=permaformulas_pop_PCclr_sig, #,
  Patterns=permaformulas_patterns,
  PopPatterns=permaformulas_pop_patterns

  #allbut1_all_food,
  #allbut1_sig_food,
  #allbut1_all_pcs,
  #allbut1_sig_pcs
)


ph2use <- c("remove_tanda2_rarefied_min", "remove_tanda2" )
permaresults_mult <- list()
permaresults_mult_byVar <- list()
for(i in ph2use){
  phobj <- all_phyloseq[[i]]
  #phobj <- updatePsWithLogs(phobj, c("Edad", "BMI"))

  cat("Doing PERMANOVA of: ", i, " with by=margin\n")
  permaresults_mult_byVar[[i]] <- lapply(dtypes, FUN=function(dd, phobj, exclude_vars, SEED){
    oname <- paste0(outdir, "permanova_resultsPatternsMult_", i, "_", dd, ".tsv")
    makePermanovaFormulas(phobj,
                          permaformulas,
                          dist_method = dd,
                          seed = SEED,
                          outname = oname, parallel=12, by="margin")
  }, phobj, exclude_vars, SEED)
  names(permaresults_mult_byVar[[i]]) <- dtypes

  cat("Doing PERMANOVA of: ", i, "\n")
  permaresults_mult[[i]] <- lapply(dtypes, FUN=function(dd, phobj, exclude_vars, SEED){
    oname <- paste0(outdir, "permanova_resultsPatternsMult_", i, "_", dd, ".tsv")
    makePermanovaFormulas(phobj,
                          permaformulas,
                          dist_method = dd,
                          seed = SEED,
                          outname = oname, parallel=12)
  }, phobj, exclude_vars, SEED)
  names(permaresults_mult[[i]]) <- dtypes


}
permaresults_mult$remove_tanda2_rarefied_min$bray$res$model_name <- names(permaformulas)
names(permaresults_mult_byVar$remove_tanda2_rarefied_min$bray$modelos) <- names(permaformulas)
names(permaresults_mult$remove_tanda2_rarefied_min$bray$modelos) <- names(permaformulas)


save(permaresults_mult, file = paste0(outdir, "PERMANOVA_MULT_250814.RData"))
save(permaresults_mult_byVar, file = paste0(outdir, "PERMANOVA_MULT_margin1_250814.RData" ))

load(paste0(outdir, "PERMANOVA_MULT_250814.RData"))
load(paste0(outdir, "PERMANOVA_MULT_margin1_250814.RData"))
load(paste0(outdir, "permanova_results_bray_250814_2.RData"))

## Code for merging with previous permanovas
#reserva1 <- permaresults_mult
#reserva2 <- permaresults_mult_byVar
#load("/home/carlos/Documentos/CORALS/results_rstudio/results_Abril25_2/PERMANOVA/PERMANOVA_MULT.RData")
#load("/home/carlos/Documentos/CORALS/results_rstudio/results_Abril25_2/PERMANOVA/PERMANOVA_MULT_margin1.RData")
#
#
#for(i in names(reserva1$remove_tanda2_rarefied_min$bray$modelos)){
#  permaresults_mult$remove_tanda2_rarefied_min$bray$modelos[[i]] <- reserva1$remove_tanda2_rarefied_min$bray$modelos[[i]]
#}
#for(i in names(reserva2$remove_tanda2_rarefied_min$bray$modelos)){
#  permaresults_mult_byVar$remove_tanda2_rarefied_min$bray$modelos[[i]] <- reserva2$remove_tanda2_rarefied_min$bray$modelos[[i]]
#}
#permaresults_mult_byVar$remove_tanda2_rarefied_min$bray$res <- rbind(permaresults_mult_byVar$remove_tanda2_rarefied_min$bray$res,
#                                                                     reserva2$remove_tanda2_rarefied_min$bray$res)
#
#permaresults_mult$remove_tanda2_rarefied_min$bray$res <- rbind(permaresults_mult$remove_tanda2_rarefied_min$bray$res,
#                                                               reserva1$remove_tanda2_rarefied_min$bray$res)
#
#save(permaresults_mult, file = paste0(outdir, "PERMANOVA_MULT.RData"))
#save(permaresults_mult_byVar, file = paste0(outdir, "PERMANOVA_MULT_margin1.RData"))


mm <- permaresults_mult_byVar$remove_tanda2_rarefied_min$bray$modelos
mm_full <- permaresults_mult$remove_tanda2_rarefied_min$bray$modelos

vars2test1 <- c("status_c1",
                "Sex",
                "hospital",
                "mother_educ",
                "age_T0",
                "z_bmi_00",
                "nreads",
                "mg_p_00",
                "status_c2",
                "z_waist_00")
vars2test1 <-  c("status_c1", "Sex",
                 "z_t0", "z_t1", "inc_z_bmi",
                 "z_waist_00", "z_waist_01", "inc_z_waist",
                 "z_waist_01_T1", "inc_z_waist_T1",
                 "z_t1_T1", "inc_z_bmi_T1",
                 "exercise_cat", "age_T0", #"mg_p_00",
                 "hospital",  "mother_educ",
                 "status_c2_T1", "nreads")

patnames <- c("Preprocessed",  "Mediterranean", "Western")
quant_vars_diet <- names(s_meta)[5:30]
diet_aggr_vars <- quant_vars_diet[1:5]
dem_all_vars <- c("Demographic all", "Demographic + Diet", "Diet all", "All patterns", "Demographic + patterns")

type_levels <- rev(c("Aggregated",
                     "Demographic",
                     "Demographic T1", "T1 (normal BMI at T0)",
                     "Diet pattern",
                     "Macronutrients",
                      "Food groups"))

xx <- read_tsv(paste0(outdir, "permanova_results_remove_tanda2_rarefied_min_bray_250814_2.tsv"))
#or:
#xx <- permaresults$remove_tanda2_rarefied_min$bray
xx_normT0 <- read_tsv(paste0(outdir, "permanova_results_remove_tanda2_rarefied_min_bray_normalT0_250814_2.tsv"))
# or:
#xx_normT0 <- permaresults_normalT0$remove_tanda2_rarefied_min$bray
xx_normT0 <- xx_normT0 %>%
  dplyr::mutate(variable = paste0(variable, "_T1")) %>%
  dplyr::mutate(model_name = variable)

betadf <- rbind(xx %>% mutate(model_name = variable),
                xx_normT0,
                permaresults_mult$remove_tanda2_rarefied_min$bray$res) %>%
  filter(model_name %in% c(quant_vars_diet_clr, vars2test1, patnames,
                           "exercise_cat", "PCclr","allpop", "PopPCclr",
                           "Patterns", "PopPatterns")) %>%
  filter(model_name != "nreads") %>%
  dplyr::mutate(model_name = ifelse(model_name == "allpop", "Demographic all", model_name)) %>%
  dplyr::mutate(model_name = ifelse(model_name == "Patterns", "All patterns", model_name)) %>%
  dplyr::mutate(model_name = ifelse(model_name == "PopPatterns", "Demographic + patterns", model_name)) %>%
  dplyr::mutate(model_name = ifelse(model_name == "PopPCclr", "Demographic + Diet", model_name)) %>%
  dplyr::mutate(model_name = ifelse(model_name == "PCclr", "Diet all", model_name))%>%
  dplyr::mutate(model_name = gsub("_clr$", "", model_name),
                type = ifelse(model_name %in% quant_vars_diet,
                              ifelse(model_name %in% diet_aggr_vars, "Macronutrients", "Food groups"),
                              "Demographic"),
                type = ifelse( grepl("[Pp]attern",model_name, perl=T) | (model_name %in% patnames), "Diet pattern", type ) #grepl("[Pp]attern",model_name, perl=T) |
                ) %>%
  dplyr::mutate(type = ifelse(model_name %in% dem_all_vars, "Aggregated", type)) %>%
  dplyr::mutate(
    model_name = gsub("_g$", " (g)", model_name, perl=T),
    model_name = gsub("_kcal", " (Kcal)", model_name),

    model_name = gsub("z_t1_T1", "Z-score BMI T1 (normal T0)", model_name, perl=T),
    model_name = gsub("z_waist_01_T1", "Z-score waist T1 (normal T0)", model_name),

    model_name = gsub("z_t0", "Z-score BMI T0", model_name, perl=T),
    model_name = gsub("z_t1", "Z-score BMI T1", model_name, perl=T),

    model_name = gsub("inc_z_waist_T1", "Change in Z-score waist (normal T0)", model_name),
    model_name = gsub("inc_z_waist", "Change in Z-score waist", model_name, perl=T),
    model_name = gsub("inc_z_bmi_T1", "Change in Z-score BMI (normal T0)", model_name, perl=T),
    model_name = gsub("inc_z_bmi", "Change in Z-score BMI", model_name, perl=T),

    model_name = gsub("_c1", " T0", model_name),
    model_name = gsub("_00", " T0", model_name),
    model_name = gsub("_01", " T1", model_name),
    model_name = gsub("mg_p", "Body fat %", model_name),
    model_name = gsub("nreads", "seq. depth", model_name),
    model_name = gsub("sauces_con", "sauces, con", model_name),
    model_name = gsub("juices_so", "juices, so", model_name),
    model_name = gsub("fats_oils", "fats, oils", model_name),
    model_name = gsub("sweets_pas", "sweets, pas", model_name),
    model_name = gsub("^z_", "Z-score ", model_name, perl=T),
    model_name = gsub("bmi", "BMI", model_name, perl=T),
    model_name = gsub("status_c2_T1", "status T1 (normal T0)", model_name, perl=T),
    model_name = gsub("_", " ", model_name)
  ) %>%
  dplyr::mutate(plog = -log10(P),
                 padj_log = -log10(padj)) %>%
  dplyr::mutate(type = ifelse( grepl("T1$|Change", model_name, perl=T), "Demographic T1", type),
                type = ifelse(grepl("\\(normal T0\\)", model_name, perl=T), "T1 (normal BMI at T0)", type)) %>%
  dplyr::mutate(type = factor(type, levels=type_levels)) %>%
  dplyr::mutate(model_name = tools::toTitleCase(model_name)) %>%
  group_by(type) %>%
  arrange(R2_var) %>%
  dplyr::mutate(model_name = factor(model_name, levels=model_name)) %>%
  dplyr::mutate(Sig = map_vec(padj_log, \(x) paste(rep("*", min(c(as.integer(x), 3)) ), collapse=""))) # %>%
  #dplyr::mutate(Sig = ifelse(padj <= 0.05, Sig, ""))


linedf <- betadf %>% group_by(type) %>%
  dplyr::summarise(xpos = max(as.numeric(model_name)) + 0.5) %>%
  head(nrow(.) - 1)

write_tsv(betadf, file = paste0(outdir, "BETADF_TABLE_USED_FOR_PLOT.tsv"))
write_tsv(linedf, file = paste0(outdir, "BETADF_TABLE_USED_FOR_PLOT_HelperLinedf.tsv"))
#colors <- ggsci::pal_lancet()(length(levels(betadf$type)))
#colors <- colors[c(1, length(colors):2)]
TSIZE <- 6
(g0 <- ggplot(betadf, aes(x=model_name, y=R2_var, col=type, fill=type))+
  geom_segment(aes(x = model_name, xend = model_name, y=0, yend=R2_var)) +
  #geom_hline(yintercept = -log10(0.05), linetype=2, col="tomato") +
  #geom_hline(yintercept = 0, linetype=1, col="black") +
  #geom_vline(aes(xintercept = Variable), col="gray", linetype=2, linewidth=0.1) +
  geom_point() +
  theme_classic() +
  ylab(expression(Adonis2 ~ R^2)) +
  ggsci::scale_color_lancet() +
  ggsci::scale_fill_lancet() +
  geom_vline(data=linedf, aes(xintercept = xpos),
             linetype=2, col="gray") +
    geom_text( #valor salario
      col = "black",
      aes(label=Sig),
      inherit.aes = TRUE,
      size = TSIZE,
      angle = 0,
      #face="bold",
      hjust = -0.8,
      vjust = 0.8
    ) +
    ylim(0, 0.1)+
    coord_flip() +
    xlab("Variable")

)
ggsave( paste0(outdir, "PERMANOVA_rsquared_merged_plot2_patterns_PatsIntoAggr_250814_padj0.1.pdf"), g0,
        width = 7, height = 7)

write_tsv(betadf, paste0(outdir, "PERMANOVA_rsquared_merged_plot_patterns_250814_padj0.1.tsv"))
