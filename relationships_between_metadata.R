

load("/home/carlos/Documentos/CORALS/results_rstudio/results_Abril25_2/foodPCA/phyloseq_list_foodPCA_withNMF.RData")
outdir <- paste0(opt$out, "/test_differences_metadata/")
if(!dir.exists(outdir)) dir.create(outdir)

s_meta <- sample_data(all_phyloseq$remove_tanda2_rarefied_min) %>% data.frame %>%
  dplyr::mutate(energy_kcal_log = log(energy_kcal+1))
this_metadata <- s_meta


s_meta_tmp <- s_meta %>% filter(!is.na(Mediterranean) & ! is.na(Preprocessed) & ! is.na(Western))
cor(s_meta_tmp$Mediterranean, s_meta_tmp$Preprocessed)
cor(s_meta_tmp$Mediterranean, s_meta_tmp$Western)
cor(s_meta_tmp$Western, s_meta_tmp$Preprocessed)
plot(s_meta_tmp$Mediterranean, s_meta_tmp$Preprocessed)
plot(s_meta_tmp$Mediterranean, s_meta_tmp$Western)
plot(s_meta_tmp$Western, s_meta_tmp$Preprocessed)



patnames <- c("Preprocessed", "Mediterranean", "Western")
foodvars_clr <- c("energy_kcal_log", names(s_meta)[grep("_clr$", names(s_meta), perl=T)])
foodvars_macro <- foodvars_clr[1:5]

othervars_cat <- c("Sex", "Category_T0",  "Category_T1", "mother_educ",
                   "exercise_cat", "hospital")

othervars_quant <- c("edad_00",  "z_bmi_00", "z_bmi_01",
                     "nreads", "mg_p_00", "z_waist_00")

allvars_num <- c(patnames, foodvars_clr, othervars_quant)
allvars <- allvars_num

(hist_clr <- s_meta %>% select(sampleID, all_of(patnames), all_of(foodvars_clr)) %>% #all_of(foodvars_clr)
  gather(var, value, -sampleID) %>%
  dplyr::mutate(var = factor(var, levels = c(patnames, foodvars_clr))) %>%
  ggplot(aes(x=value)) +
  facet_wrap(~ var, scales="free") +
  geom_histogram(stat = "bin", bins=30, fill="steelblue4")+
  theme_bw()
)
ggsave(paste0(outdir, "histogram_food_clr.pdf"), hist_clr, width = 12, height = 10)


allmods <- map(allvars_num, \(x){
  varsthis <- allvars_num[allvars_num != x]
  makeLinearModelsSingleVariable(this_metadata, varsthis[1],
                                 varsthis[2:length(varsthis)],
                                 c(x),
                                 combos=1,
                                 outdir = outdir,
                                 name = paste0("models_", x) )


})
names(allmods) <- allvars_num

allmods_tab <- map(allmods, \(x) x$single_anovas) %>%
  bind_rows() %>%
  dplyr::mutate(testvar = map_vec(model, \(x) strsplit(x, " ~ ")[[1]][2]))
write_tsv(file = paste0(outdir, "all_models_quant_vars.tsv") , allmods_tab)


alldif <- map(allvars, \(x){
  varsthis <- allvars[allvars != x]
  testDiversityDifferences(this_metadata,
                           c(x),
                           varsthis,
                           outdir,
                           paste0("models_", x))


})
names(alldif) <- allvars
alldif_df <- alldif %>% bind_rows()
write_tsv(file = paste0(outdir, "all_models_categ_vars.tsv") , alldif_df)

################33


library(outliers)

foodvars_all <- c(names(s_meta)[5:30],
              "energy_kcal_log",
              names(s_meta)[grep("_clr$", names(s_meta), perl=T)])

P_GRUBBS <- 0.0001
outliers_food <- map(foodvars_all, \(vv){
  aux <- s_meta %>%
    filter(!is.na(!!sym(vv)))
  outs <- data.frame()
  while(grubbs.test(aux[, vv])$p.value < P_GRUBBS & nrow(aux) > 2){
   outs <- rbind(outs,
             aux %>% filter(!!sym(vv) == max(!!sym(vv) )) %>%
               select(sampleID, all_of(vv))
        )
   aux <- aux %>% filter(!!sym(vv) < max(!!sym(vv) ))
  }
  return(outs)
})
names(outliers_food) <- foodvars_all
sapply(outliers_food, nrow)
save(outliers_food, file = paste0(outdir, "outliers_food_Grubbs_pe4"))




allmods_noOuts <- map(allvars_num, \(x){
  varsthis <- allvars[allvars != x]
  if(gsub("_clr$", "", x, perl=T) %in% names(outliers_food) ){
    out_this <- outliers_food[[gsub("_clr$", "", x, perl=T)]]
    aux <- this_metadata %>% filter(!sampleID %in% out_this$sampleID)
  }else{
    aux <- this_metadata
  }

  makeLinearModelsSingleVariable(aux, varsthis[1],
                                 varsthis[2:length(varsthis)],
                                 c(x),
                                 combos=1,
                                 outdir = outdir,
                                 name = paste0("models_noOutliersGrubbsPe4_", x) )


})
names(allmods_noOuts) <- allvars_num

allmods_tab_noOuts <- map(allmods_noOuts, \(x) x$single_anovas) %>%
  bind_rows() %>%
  dplyr::mutate(testvar = map_vec(model, \(x) strsplit(x, " ~ ")[[1]][2]))
write_tsv(file = paste0(outdir, "all_models_quant_vars_noOuts.tsv") , allmods_tab_noOuts)


alldif_noOuts <- map(allvars, \(x){
  varsthis <- allvars[allvars != x]
  if(gsub("_clr$", "", x, perl=T) %in% names(outliers_food) ){
    out_this <- outliers_food[[gsub("_clr$", "", x, perl=T)]]
    aux <- this_metadata %>% filter(!sampleID %in% out_this$sampleID)
  }else{
    aux <- this_metadata
  }
  testDiversityDifferences(aux,
                           c(x),
                           varsthis,
                           outdir,
                           paste0("models_", x))


})
names(alldif_noOuts) <- allvars
alldif_df_noOuts <- alldif_noOuts %>% bind_rows()
write_tsv(file = paste0(outdir, "all_models_categ_vars_noOuts.tsv") , alldif_df_noOuts)

# plot

cat_pop <- c("Sex", "status_c1",  "status_c2", "mother_educ",
                   "exercise_cat", "hospital")

quant_pop <-  c("edad_00",  "z_bmi_00", "z_bmi_01",
                 "mg_p_00", "z_waist_00")
quant_food <- c(patnames, foodvars_clr)
macro_food <- foodvars_clr[1:5]

plotdf_q <- s_meta %>% select(sampleID,
                              all_of(quant_pop),
                              all_of(quant_food)) %>%
  gather("pop_var", "pop_val", all_of(quant_pop)) %>%
  gather("diet_var", "diet_val", all_of(quant_food)) %>%
  dplyr::mutate(
    diet_var_orig = diet_var,
    pop_var_orig = pop_var
  ) %>%
  dplyr::mutate(diet_type = ifelse(diet_var %in% patnames, "Diet pattern", "Food groups"),
                diet_type = ifelse(diet_var %in% macro_food, "Macronutrients", diet_type),
                diet_type = factor(diet_type, levels=c("Macronutrients", "Diet pattern", "Food groups")),
  ) %>% dplyr::mutate(
    diet_var = gsub("_clr$", "", diet_var, perl=T),
    diet_var = gsub("_g$", " (g)", diet_var, perl=T),
    diet_var = gsub("_kcal", " (Kcal)", diet_var),
    pop_var = gsub("_c1", " T0", pop_var),
    pop_var = gsub("_00", " T0", pop_var),
    pop_var = gsub("mg_p", "Body fat %", pop_var),
    pop_var = gsub("nreads", "seq. depth", pop_var),
    diet_var = gsub("sauces_con", "sauces, con", diet_var),
    diet_var = gsub("juices_so", "juices, so", diet_var),
    diet_var = gsub("fats_oils", "fats, oils", diet_var),
    diet_var = gsub("sweets_pas", "sweets, pas", diet_var),
    pop_var = gsub("^z_", "Z-score ", pop_var, perl=T),
    pop_var = gsub("bmi", "BMI", pop_var, perl=T),
    pop_var = gsub("status_c2_T1", "status T1 (norm. w. T0)", pop_var, perl=T),
    pop_var = gsub("_", " ", pop_var),
    diet_var = gsub("_", " ", diet_var)
  )

assertthat::assert_that( all(plotdf_q$pop_var_orig %in% allmods_tab_noOuts$Index ))
assertthat::assert_that(all(plotdf_q$diet_var_orig %in% allmods_tab_noOuts$testvar ))

REG_PLIM <- 0.01
plotdf_q <- plotdf_q %>%
  dplyr::mutate(x_y = paste(diet_var_orig, pop_var_orig, sep="_")) %>%
  merge(allmods_tab_noOuts %>%
          dplyr::mutate(x_y = paste(testvar, Index, sep="_")) %>%
          select(-Index, -testvar),
        by="x_y", all.x=T, all.y=F) %>%
  group_by(pop_var_orig) %>%
  dplyr::mutate(padj = p.adjust(`Pr(>F)`, method = "BH")) %>%
  dplyr::mutate(Sig = ifelse(`Pr(>F)` < REG_PLIM, paste0("p<", as.character(REG_PLIM)), "NS"),
                Sig_padj = ifelse(padj < REG_PLIM, paste0("p<", as.character(REG_PLIM)), "ns"))


(gg <- ggplot(plotdf_q, aes(x=diet_val, y=pop_val, col=Sig, fill=Sig )) +
  facet_grid(diet_var ~ pop_var) +
    geom_smooth(col="darkgray", alpha=0.2) +
  geom_point(alpha=0.2, size=0.05) +
  theme_classic()
)
ggsave(filename = paste0(outdir, "relation_between_BMI_and_other_variables.pdf"), gg, width = 10, height = 22)
