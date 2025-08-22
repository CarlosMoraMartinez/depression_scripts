# Alpha 4 each

## Cualitativas

#load(paste0(opt$out, "foodPCA/phyloseq_list_foodPCA_withNMF.RData"))
load("/home/carlos/Documentos/CORALS/results_rstudio/results_Agosto25_1//foodPCA/phyloseq_list_foodPCA_withNMF_withInc.RData")

s_meta <- sample_data(all_phyloseq$remove_tanda2_rarefied_min) %>% data.frame
alpha_indices <- c("Observed", "Chao1", "Shannon", "InvSimpson")
vars2test <- c("status_c2", "status_c1","Sex", "Category_T0", "Category_T1",
               "hospital", "age_class1", "mother_educ", "exercise_cat",
               "cat_peso_01", "cat_peso_00", "pattern")

quant_vars <- c("age_T0", "z_bmi_00", "z_bmi_01", "nreads", "mg_p_00", "z_waist_00", "z_waist_01",
                "z_t0", "z_t1",
                "inc_z_waist", "inc_z_bmi", "inc_z_bmiCOR",
                "Preprocessed", "Mediterranean", "Western")
patnames <- c("Preprocessed", "Mediterranean", "Western")
quant_vars_diet <- names(s_meta)[5:30]
quant_vars_diet_clr <- names(s_meta)[grep("_clr$", names(s_meta))]
quant_vars_diet_PCs <- names(s_meta)[grep("^PC", names(s_meta))]
quant_vars_diet_PCs_clr <-names(s_meta)[grep("clr_PC", names(s_meta))]

vars2log <- c()

if(length(vars2log) > 0){
  quant_vars_ext <- c(quant_vars, paste(vars2log, "_log", sep=""))
}else{
  quant_vars_ext <- quant_vars
}
interestvar <- "age_T0" # "status_c2"
extravars <- c(quant_vars, vars2test)
extravars <- extravars[extravars != interestvar]
extravars <- extravars[extravars != "age_class1"]
extravars <- extravars[extravars != "edad_00_meses"]

outdir <- paste0(opt$out, "/AlphaDiversity/")
if(!dir.exists(outdir)) dir.create(outdir)

extravars2 <- c("Sex", "z_bmi_00",  "z_bmi_01", "mother_educ", "z_waist_00",  "z_waist_01",
                "inc_z_waist", "inc_z_bmi",
                "exercise_cat") #"edad_00_meses_log"

phseq_to_use <- c("remove_tanda2", "remove_tanda2_rarefied_min") #names(all_phyloseq)
#load(allphyloseqlist_fname)

for(phname in phseq_to_use){
  cat("Alpha diversity in ", phname, "\n")
  phobj <- all_phyloseq[[phname]]
  phobj <- updatePsWithLogs(phobj, vars2log)
  sample_data(phobj)$tanda[is.na(sample_data(phobj)$tanda)] <- 2

  #sample_data(phobj)$tanda <- as.character(sample_data(phobj)$tanda)
  #sample_data(phobj)$educ_m_discrete <- ifelse(sample_data(phobj)$edu_m_00 <4, "1-3",
  #                                             ifelse(sample_data(phobj)$edu_m_00 <7, "4-6", "7-10"))

  sample_data(phobj)$status_c2 <- factor(sample_data(phobj)$status_c2,
                                         levels = c("Insufficient gain", "Normal", "Excessive gain"))
  divtab <- calculateAlphaDiversityTable(phobj, outdir, alpha_indices, paste0(phname, "_AlphaDiv") )

  models1 <- makeLinearModelsSingleVariable(divtab, interestvar,
                                            extravars,
                                            alpha_indices,
                                            combos=1,
                                            outdir = outdir, name = paste0(phname, "_AlphaDiv_linMod1var") )

  models2 <- makeLinearModelsSingleVariable(divtab, interestvar,
                                            extravars2,
                                            alpha_indices,
                                            combos=1:4,
                                            outdir = outdir, name = paste0(phname, "_AlphaDiv_linModManyVars") )


  alphadif <- testDiversityDifferences(divtab, alpha_indices, vars2test, outdir, "AlphaDiv_rawdata")


  divplots <- getAlphaDiversityCustomPlot(phobj, vars2test, quant_vars_ext,
                                opt,
                                indices= alpha_indices,
                                correct_pvalues = T, correct_pvalues_indices = F,
                                name = paste0(phname, "_AlphaDivCusP"), w = 10, h = 4)
  divplots <- getAlphaDiversityCustomPlot(phobj, vars2test, quant_vars_ext,
                                opt,
                                indices= alpha_indices,
                                correct_pvalues = T, correct_pvalues_indices = T,
                                name = paste0(phname, "_AlphaDivCusPAdjInd"), w = 10, h = 4)
  divplots <- getAlphaDiversityCustomPlot(phobj, vars2test, quant_vars_ext,
                                          opt,
                                          indices= alpha_indices,
                                          correct_pvalues = T, correct_pvalues_indices = F,
                                          name = paste0(phname, "_AlphaDivCusPraw"), w = 10, h = 3) #3 better for regressions

  # Models for status_c2
  new_divtab <- divtab %>% filter(status_c1 == "normal")
  new_phobj <- subset_samples(phobj, status_c1 == "normal")

  models1_onlyOw <- makeLinearModelsSingleVariable(new_divtab, "status_c2",
                                            extravars,
                                            alpha_indices,
                                            combos=1,
                                            outdir = outdir, name = paste0(phname, "_AlphaDiv_linMod1var_onlyT0Normal") )

  models2_onlyOw <- makeLinearModelsSingleVariable(new_divtab, "status_c2",
                                            extravars2,
                                            alpha_indices,
                                            combos=1:4,
                                            outdir = outdir, name = paste0(phname, "_AlphaDiv_linModManyVars_onlyT0Normal") )


  alphadif_onlyOw <- testDiversityDifferences(new_divtab, alpha_indices, vars2test, outdir, "AlphaDiv_rawdata_onlyT0Normal")


  divplots_onlyOw <- getAlphaDiversityCustomPlot(new_phobj, vars2test, quant_vars_ext,
                                          opt,
                                          indices= alpha_indices,
                                          correct_pvalues = T, correct_pvalues_indices = F,
                                          name = paste0(phname, "_AlphaDivCusP_onlyT0Normal"), w = 10, h = 4)
  divplots_onlyOw <- getAlphaDiversityCustomPlot(new_phobj, vars2test, quant_vars_ext,
                                          opt,
                                          indices= alpha_indices,
                                          correct_pvalues = T, correct_pvalues_indices = T,
                                          name = paste0(phname, "_AlphaDivCusPAdjInd_onlyT0Normal"), w = 10, h = 4)
  divplots_onlyOw <- getAlphaDiversityCustomPlot(new_phobj, vars2test, quant_vars_ext,
                                          opt,
                                          indices= alpha_indices,
                                          correct_pvalues = T, correct_pvalues_indices = F,
                                          name = paste0(phname, "_AlphaDivCusPraw_onlyT0Normal"), w = 10, h = 3) #3 better for regressions

  #Models for diet
  models1_diet <- makeLinearModelsSingleVariable(divtab, "status_c2",
                                                   quant_vars_diet,
                                                   alpha_indices,
                                                   combos=1,
                                                   outdir = outdir, name = paste0(phname, "_AlphaDiv_linMod1var_Diet") )
  models1_diet2 <- makeLinearModelsSingleVariable(divtab, "age_T0",
                                                 quant_vars_diet,
                                                 alpha_indices,
                                                 combos=1,
                                                 outdir = outdir, name = paste0(phname, "_AlphaDiv_linMod1var_DietAgeAdj") )
  models1_dietclr <- makeLinearModelsSingleVariable(divtab, "status_c2",
                                                    quant_vars_diet_clr,
                                                 alpha_indices,
                                                 combos=1,
                                                 outdir = outdir, name = paste0(phname, "_AlphaDiv_linMod1var_DietCLR") )
  models1_diet2clr <- makeLinearModelsSingleVariable(divtab, "age_T0",
                                                     quant_vars_diet_clr,
                                                  alpha_indices,
                                                  combos=1,
                                                  outdir = outdir, name = paste0(phname, "_AlphaDiv_linMod1var_DietAgeAdjCLR") )
  models1_dietPCs <- makeLinearModelsSingleVariable(divtab, "status_c2",
                                                 quant_vars_diet_PCs,
                                                 alpha_indices,
                                                 combos=1,
                                                 outdir = outdir, name = paste0(phname, "_AlphaDiv_linMod1var_DietPCs") )
  models1_dietPCs_clr <- makeLinearModelsSingleVariable(divtab, "status_c2",
                                                    quant_vars_diet_PCs_clr,
                                                    alpha_indices,
                                                    combos=1,
                                                    outdir = outdir, name = paste0(phname, "_AlphaDiv_linMod1var_DietPCsCLR") )
  divplots_diet <- getAlphaDiversityCustomPlot(phobj, vars2test[1], quant_vars_diet,
                                                 opt,
                                                 indices= alpha_indices,
                                                 correct_pvalues = T, correct_pvalues_indices = F,
                                                 name = paste0(phname, "_AlphaDivCusPraw_Diet"), w = 10, h = 3)

  divplots_diet_clr <- getAlphaDiversityCustomPlot(phobj, vars2test[1], quant_vars_diet_clr,
                                               opt,
                                               indices= alpha_indices,
                                               correct_pvalues = T, correct_pvalues_indices = F,
                                               name = paste0(phname, "_AlphaDivCusPraw_DietCLR"), w = 10, h = 3)
  divplots_dietPCs <- getAlphaDiversityCustomPlot(phobj, vars2test[1], quant_vars_diet_PCs,
                                               opt,
                                               indices= alpha_indices,
                                               correct_pvalues = T, correct_pvalues_indices = F,
                                               name = paste0(phname, "_AlphaDivCusPraw_DietPCs"), w = 10, h = 3)
  divplots_dietPCs_clr <- getAlphaDiversityCustomPlot(phobj, vars2test[1], quant_vars_diet_PCs_clr,
                                                  opt,
                                                  indices= alpha_indices,
                                                  correct_pvalues = T, correct_pvalues_indices = F,
                                                  name = paste0(phname, "_AlphaDivCusPraw_DietPCsCLR"), w = 10, h = 3)
}

# Plot alpha div significance

names_models <- names(models1_dietclr$models)[sapply(paste0("\\~ ",  quant_vars_diet_clr, "$"),
                                                     grep,
                                                     names(models1_dietclr$models), perl=T) %>%
  as.vector]

food_model_coefs <- map(models1_dietclr$models[names_models], \(x) {data.frame(Intercept=x$coefficients[1],
                                                           Slope = x$coefficients[2],
                                                           Variable = names(x$coefficients)[2] ) }) %>%
  bind_rows %>%
  mutate(model = names_models)

names_models2 <- names(models1_diet$models)[sapply(paste0("\\~ ", quant_vars_diet, "$"),
                                                     grep,
                                                     names(models1_diet$models), perl=T) %>%
                                                as.vector]
food_model_coefs2 <- map(models1_diet$models[names_models2], \(x) {data.frame(Intercept=x$coefficients[1],
                                                                                Slope = x$coefficients[2],
                                                                                Variable = names(x$coefficients)[2] ) }) %>%
  bind_rows %>%
  mutate(model = names_models2)

food_model_coefs <- food_model_coefs %>%
  rbind(food_model_coefs2)

all_models <- rbind(models1_dietclr$single_anovas %>% filter(sapply(model, \(x)any(sapply(quant_vars_diet_clr, grepl, x)) )),
                    models1$single_anovas ,
                    models1_diet$single_anovas,
                    models1_onlyOw$single_anovas %>%
                      filter(grepl("status_c2|z_waist_01|z_bmi_01|inc_z_waist|inc_z_bmi|z_t0|z_t1", model, perl=T)) %>%
                      #dplyr::mutate(model=gsub("status_c2", "status_c2_T1", model)) %>%
                      dplyr::mutate(model=paste0(model, "_T1")) %>%
                      distinct
                    ) %>%
  dplyr::mutate(Variable = map_vec(model, \(x)strsplit(x, " ~ ")[[1]][2]))


model_df <- merge(all_models, food_model_coefs, by=c("Variable", "model"), all.x=TRUE, all.y=TRUE)
write_tsv(model_df, paste0(opt$out, "AlphaDiversity/alpha_div_models_merged_250814.tsv"))

## brief tests
ggplot(divtab, aes(x=inc_z_bmi, y = Shannon, col=status_c2)) +
  geom_point() +
  geom_smooth(method="lm") + facet_grid(~ status_c1)

ggplot(divtab, aes(x=bmi_t0, y = inc_z_bmi, col=Shannon)) +
  geom_point() +
  geom_smooth()

ggplot(divtab, aes(x=z_bmi_00, y = inc_z_bmi, col=Shannon)) +
  geom_point() +
  geom_smooth()

ggplot(divtab, aes(x=z_bmi_00, y = Shannon , col=inc_z_bmi)) +
  geom_point() +
  geom_smooth()

ggplot(divtab, aes(x=z_bmi_00, y = Shannon , col=status_c1)) +
  geom_point() +
  geom_smooth(method="lm")

ggplot(divtab, aes(x=age_T0, y = Shannon , col=hospital)) +
  geom_point() +
  geom_smooth(method="lm")
ggplot(divtab, aes(x=age_T0, y = Observed , col=hospital)) +
  geom_point() +
  geom_smooth(method="lm")
########
vars2plot <- c(quant_vars_diet_clr,
               "energy_kcal", patnames,
               c("status_c1", "Sex",
                 "z_t0", "z_t1", "inc_z_bmi",
                 "z_waist_00", "z_waist_01", "inc_z_waist",
                 "z_waist_01_T1", "inc_z_waist_T1",
                 "z_t1_T1", "inc_z_bmi_T1",
                 "exercise_cat", "age_T0", #"mg_p_00",
                 "hospital",  "mother_educ",
                "status_c2_T1", "nreads") #"age_class1",
         )

diet_aggr_vars <- quant_vars_diet[1:5]

models2plot <- model_df %>%
  filter(Variable %in% vars2plot) %>%
  filter(Variable != "nreads") %>%
  dplyr::mutate(padj = p.adjust(`Pr(>F)`, method = "BH"),
                plog = -log10(`Pr(>F)`),
                padj_log = -log10(padj)) %>%
  group_by(Index) %>%
  dplyr::mutate(padj_groupInd = p.adjust(`Pr(>F)`, method = "BH"),
                padj_log_groupInd = -log10(padj_groupInd)) %>%
  ungroup() %>%
  dplyr::mutate(Variable = gsub("_clr$", "", Variable),
                type = ifelse(Variable %in% quant_vars_diet,
                              ifelse(Variable %in% diet_aggr_vars, "Macronutrients", "Food groups"),
                              "Demographic"),
                type = ifelse(Variable %in% patnames, "Diet pattern", type)) %>%
  dplyr::mutate(
                Variable = gsub("_g$", " (g)", Variable, perl=T),
                Variable = gsub("_kcal", " (Kcal)", Variable),
                Variable = gsub("z_t1_T1", "Z-score BMI T1 (normal T0)", Variable, perl=T),
                Variable = gsub("z_waist_01_T1", "Z-score waist T1 (normal T0)", Variable),

                Variable = gsub("z_t0", "Z-score BMI T0", Variable, perl=T),
                Variable = gsub("z_t1", "Z-score BMI T1", Variable, perl=T),

                Variable = gsub("inc_z_waist_T1", "Change in Z-score waist (normal T0)", Variable),
                Variable = gsub("inc_z_waist", "Change in Z-score waist", Variable, perl=T),
                Variable = gsub("inc_z_bmi_T1", "Change in Z-score BMI (normal T0)", Variable, perl=T),
                Variable = gsub("inc_z_bmi", "Change in Z-score BMI", Variable, perl=T),

                Variable = gsub("_c1", " T0", Variable),
                Variable = gsub("_00", " T0", Variable),
                Variable = gsub("_01", " T1", Variable),
                Variable = gsub("mg_p", "Body fat %", Variable),
                Variable = gsub("nreads", "seq. depth", Variable),
                Variable = gsub("sauces_con", "sauces, con", Variable),
                Variable = gsub("juices_so", "juices, so", Variable),
                Variable = gsub("fats_oils", "fats, oils", Variable),
                Variable = gsub("sweets_pas", "sweets, pas", Variable),
                Variable = gsub("^z_", "Z-score ", Variable, perl=T),
                Variable = gsub("bmi", "BMI", Variable, perl=T),
                Variable = gsub("status_c2_T1", "status T1 (normal T0)", Variable, perl=T),
                Variable = gsub("_", " ", Variable)
                ) %>%
  dplyr::mutate(type = ifelse( grepl("T1$|Change", Variable, perl=T), "Demographic T1", type),
                type = ifelse(grepl("\\(normal T0\\)", Variable, perl=T), "T1 (normal BMI at T0)", type)) %>%

  dplyr::mutate(type = factor(type, levels=rev(c("Demographic", "Demographic T1", "T1 (normal BMI at T0)",
                                                 "Diet pattern", "Macronutrients", "Food groups")))) %>%
  dplyr::mutate(Variable = tools::toTitleCase(Variable)) %>%
  group_by(type) %>%
  arrange(type, plog) %>%
  mutate(Variable = factor(Variable, levels=Variable[Index=="Observed"])) %>%
  mutate(Index = factor(Index, levels=alpha_indices))

write_tsv(models2plot, paste0(opt$out, "AlphaDiversity/alpha_div_models_merged_2plot3_250814.tsv"))

#models2plot %>% select(Variable, type) %>% distinct() %>% view

linedf <- models2plot %>% group_by(type) %>%
  dplyr::summarise(xpos = max(as.numeric(Variable))+0.5) %>%
  head(nrow(.)-1)

colors <- ggsci::pal_lancet()(7)
colors <- colors[-7]
assertthat::assert_that(length(colors) == length(unique(models2plot$type)))

(g0 <- ggplot(models2plot, aes(x=Variable, y=plog, col=type, fill=type))+
  facet_grid(~ Index) +
  geom_hline(yintercept = -log10(0.05), linetype=2, col="tomato") +
  geom_hline(yintercept = 0, linetype=1, col="black") +
  geom_vline(data = linedf, inherit.aes=F, aes(xintercept = xpos), linetype=2, col="gray") +
  geom_vline(aes(xintercept = Variable), col="gray", linetype=2, linewidth=0.1) +
  geom_segment(aes(x = Variable, xend = Variable, y=0, yend=plog)) +
  geom_point() +
  theme_classic() +
  ylab("-10log(p value)") +
  coord_flip() +
  #ggsci::scale_color_lancet() +
  #ggsci::scale_fill_lancet()
  scale_color_manual(values=colors) +
  scale_fill_manual(values=colors)
)
ggsave( paste0(opt$out, "AlphaDiversity/alpha_div_models_merged_2plot_withDietPat_250814.pdf"), g0,
        width = 9, height = 6.5) # for old version: idth = 6.5, height = 6

(g0 <- ggplot(models2plot, aes(x=Variable, y=padj_log, col=type, fill=type))+
    facet_grid(~ Index) +
    geom_hline(yintercept = -log10(0.05), linetype=2, col="tomato") +
    geom_hline(yintercept = 0, linetype=1, col="black") +
    geom_vline(data = linedf, inherit.aes=F, aes(xintercept = xpos), linetype=2, col="gray") +
    geom_vline(aes(xintercept = Variable), col="gray", linetype=2, linewidth=0.1) +
    geom_segment(aes(x = Variable, xend = Variable, y=0, yend=padj_log)) +
    geom_point() +
    theme_classic() +
    ylab("-10log(p value)") +
    coord_flip() +
    #ggsci::scale_color_lancet() +
    #ggsci::scale_fill_lancet()
    scale_color_manual(values=colors) +
    scale_fill_manual(values=colors)
)
ggsave( paste0(opt$out, "AlphaDiversity/alpha_div_models_merged_2plot_withDietPat_padj3_250814.pdf"), g0,
        width = 9, height = 6.5)

(g0 <- ggplot(models2plot, aes(x=Variable, y=padj_log_groupInd, col=type, fill=type))+
    facet_grid(~ Index) +

    geom_hline(yintercept = -log10(0.05), linetype=2, col="tomato") +
    geom_hline(yintercept = 0, linetype=1, col="black") +
    geom_vline(data = linedf, inherit.aes=F, aes(xintercept = xpos), linetype=2, col="gray") +
    geom_vline(aes(xintercept = Variable), col="gray", linetype=2, linewidth=0.1) +
    geom_segment(aes(x = Variable, xend = Variable, y=0, yend=padj_log)) +
    geom_point() +
    theme_classic() +
    ylab("-10log(p value)") +
    coord_flip() +
    #ggsci::scale_color_lancet() +
    #ggsci::scale_fill_lancet()
    scale_color_manual(values=colors) +
    scale_fill_manual(values=colors) +
    scale_y_continuous(breaks = scales::pretty_breaks(n=5))
)
ggsave( paste0(opt$out, "AlphaDiversity/alpha_div_models_merged_2plot_withDietPat_padj3_normByIndex_250814.pdf"), g0,
        width = 9, height = 6.5)

only_sig <- models2plot %>% filter(`Pr(>F)` <= 0.05)
write_tsv(only_sig, paste0(opt$out, "AlphaDiversity/alpha_div_models_merged_2plot_onlySig3_250814.tsv"))

only_sig <- models2plot %>% filter(padj <= 0.05)
write_tsv(only_sig, paste0(opt$out, "AlphaDiversity/alpha_div_models_merged_2plot_onlySig3_padj_250814.tsv"))

only_sig <- models2plot %>% filter(padj_groupInd <= 0.05)
write_tsv(only_sig, paste0(opt$out, "AlphaDiversity/alpha_div_models_merged_2plot_onlySig3_padj_groupByIndex_250814.tsv"))

#####################################
# Beta 4 each
outdir <- paste0(opt$out, "/BetaDiversity/")
if(!dir.exists(outdir)) dir.create(outdir)

phseq_to_use <- names(all_phyloseq)[4:5]  #[c(9,10,2,7)] # [c(2,3,7,9)]

dists <- c("bray") # "jaccard"
METHODS <- c("PCoA", "NMDS") #, "NMDS"
vars2pcoa <- c(vars2test, quant_vars_ext)
ccaplots <- list()
for(phname in phseq_to_use){
  for(method in METHODS){
    for(dist in dists){
      name <- paste0(phname, "_", dist, "_", method)
      cat("Beta diversity for ", name, "\n")
      if(method != "NMDS"){
        extradims_use <- 2:3
        w <- 12
        h <- 4
      }else{
        extradims_use <- c(2)
        w <- 6
        h <- 4
      }
      logvars <- vars2pcoa[grepl("_log$", vars2pcoa, perl=T)]
      origvars <- gsub("_log", "", logvars)
      phobj <- updatePsWithLogs(all_phyloseq[[phname]], origvars)

      ccaplots[[name]] <- makeAllPCoAs(phobj, outdir,
                                       method = method,
                                       name = name,
                                       dist_type = dist,
                                       dist_name = dist,
                                       vars2plot = vars2pcoa,
                                       extradims = extradims_use,
                                       create_pdfs = T, w=w, h=h) #w=16, h=12
    }}}

# Composition 4 each

outdir <- paste0(opt$out, "/DescriptiveAbundances/")
if(!dir.exists(outdir)) dir.create(outdir)
tops <- c(5, 10)
interestvar <- "status_c2"

for(phname in phseq_to_use){
  cat("Doing Abundance Plots for: ", phname, "\n")
  abund_plots <- plotAbundanceFullPipeline(all_phyloseq[[phname]], interestvar, outdir, phname, c("Control", "Depression"), tops)
}
