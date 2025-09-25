library(tidyverse)
library(phyloseq)

makeReplacementsDF <- function(df, replace_list, totitlecase=c()){
  dfmod <- purrr::reduce(replace_list, ~ .x %>%
                           dplyr::mutate(!!sym(.y[3]) := gsub(.y[1], .y[2], !!sym(.y[3]), perl=TRUE)),
                         .init = df)%>%
    dplyr::mutate(across(all_of(totitlecase), tools::toTitleCase))
  
  return(dfmod)
}


# opt <- list(out="/home/carlos/Documentos/CORALS/results_rstudio/results_Agosto25_1/")

load("/home/carlos/Documentos/CORALS/results_rstudio/results_Agosto25_1//foodPCA/phyloseq_list_foodPCA_withNMF_withInc.RData")
outdir <- paste0(opt$out, "/test_differences_metadata3/")
if(!dir.exists(outdir)) dir.create(outdir)

s_meta <- sample_data(all_phyloseq$remove_tanda2_rarefied_min) %>% data.frame %>%
  dplyr::mutate(energy_kcal_log = log(energy_kcal+1))
this_metadata <- s_meta



vars2test1 <-  c("status_c1", "Sex",
                 "z_t0", "z_t1", "inc_z_bmi",
                 "z_waist_00", "z_waist_01", "inc_z_waist",
                 #"z_waist_01_T1", "inc_z_waist_T1",
                 #"z_t1_T1", "inc_z_bmi_T1",
                 "exercise_cat", "age_T0", #"mg_p_00",
                 "hospital",  "mother_educ",
                 #"status_c2_T1", 
                 "nreads")

patnames <- c("Preprocessed",  "Mediterranean", "Western")
quant_vars_diet <- names(s_meta)[5:30]
diet_aggr_vars <- quant_vars_diet[1:5]

table(vars2test1 %in% names(s_meta))
vars2test1[!vars2test1 %in% names(s_meta)]

alphadif <- testDiversityDifferences(s_meta, c(patnames, quant_vars_diet), 
                                     groupvars = vars2test1, 
                                     outdir = outdir, name = "TestDifferencesFood1")

alphadif_tab2 <- alphadif %>% dplyr::mutate(
  dataset = "T0",
  pval  = ifelse(comparison == "all" & is.na(wilcox_test), anova_p, wilcox_test),
  padj = ifelse(comparison == "all" & is.na(wilcox_test), anova_p, wilcox_corrected),
  test = ifelse(comparison == "all" & is.na(wilcox_test), "ANOVA", "Wilcox"), 
  type =  ifelse(comparison == "all" & is.na(wilcox_test), "categorical_all", "categorical")
) %>% 
  dplyr::select(dataset, type, test, variable, groups, comparison, pval, padj)

models <- makeLinearModelsSingleVariable(s_meta, vars2test1[1],
                                         vars2test1[-1],
                                         c(patnames, quant_vars_diet),
                                         combos=1,
                                         outdir = outdir, 
                                         name = "linmodels1")
models_tab <- models$single_anovas %>% 
  dplyr::mutate(variable = Index,
                groups = gsub(".+ ~ ", "", model),
                comparison = "all", 
                type = "numerical",
                test= "ANOVA", 
                dataset = "T0") %>% 
  dplyr::select(-Index) %>% 
  dplyr::mutate(pval = `Pr(>F)`, 
                padj = p.adjust(`Pr(>F)`, method="BH")) %>% 
  dplyr::select(dataset, type, test, variable, groups, comparison, pval, padj)

alphadif_normT0 <- testDiversityDifferences(s_meta %>% filter(status_c1 == "normal"), 
                                            c(patnames, quant_vars_diet), 
                                     groupvars = c("status_c2", vars2test1[-1]), 
                                     outdir = outdir, name = "TestDifferencesFood_normT0")

alphadif_normT0_tab2 <- alphadif_normT0 %>% dplyr::mutate(
  dataset = "T1_normT0",
  pval  = ifelse(comparison == "all" & is.na(wilcox_test), anova_p, wilcox_test),
  padj = ifelse(comparison == "all" & is.na(wilcox_test), anova_p, wilcox_corrected),
  test = ifelse(comparison == "all" & is.na(wilcox_test), "ANOVA", "Wilcox"), 
  type =  ifelse(comparison == "all" & is.na(wilcox_test), "categorical_all", "categorical")
) %>% 
  dplyr::select(dataset, type, test, variable, groups, comparison, pval, padj)

models_normT0 <- makeLinearModelsSingleVariable(s_meta %>% filter(status_c1 == "normal"), 
                                         "status_c2",
                                         vars2test1[-1],
                                         c(patnames, quant_vars_diet),
                                         combos=1,
                                         outdir = outdir, 
                                         name = "linmodels1_normT0")

models_tab_T0 <- models_normT0$single_anovas %>% 
  dplyr::mutate(variable = Index,
                groups = gsub(".+ ~ ", "", model),
                comparison = "all", 
                type = "numerical",
                test= "ANOVA", 
                dataset = "T1_normT0") %>% 
  dplyr::select(-Index) %>% 
  dplyr::mutate(pval = `Pr(>F)`, 
                padj = p.adjust(`Pr(>F)`, method="BH")) %>% 
  dplyr::select(dataset, type, test, variable, groups, comparison, pval, padj)


mergedtabs <- rbind(alphadif_tab2, models_tab, alphadif_normT0_tab2, models_tab_T0)
mergedtabs <- mergedtabs %>% dplyr::mutate(
  #padj_global = p.adjust(pval, method = "BH"), # make this after removing uninteresting tests
  log10p = -log10(pval),
  log10padj = -log10(padj),
  #log10padj_global = -log10(padj_global),
  Sig = ifelse(padj<=0.01, "**", ifelse(padj<="0.05", "*", ""))
  #Sig_global = ifelse(padj_global<=0.01, "**", ifelse(padj_global<="0.05", "*", "")),
)
write_tsv(mergedtabs, paste0(outdir, "merged_statistical_differences.tsv"))


type_levels <- rev(c("Aggregated",
                     "Demographic",
                     "Demographic T1", "T1 (normal BMI at T0)",
                     "Diet pattern",
                     "Macronutrients",
                     "Food groups"))
replace_strings <- list(
  c("_clr$", ""),
  c("_g$", " (g)"),
  c("_kcal", " (Kcal)"),
  c("z_t1", "Z-score BMI T1"),
  c("z_waist_01", "Z-score waist T1"),
  c("z_t0", "Z-score BMI T0"),
  c("z_t1", "Z-score BMI T1"),
  c("inc_z_waist", "Change in Z-score waist"),
  c("inc_z_bmi", "Change in Z-score BMI"),
  c("_c1", " T0"),
  c("_00", " T0"),
  c("_01", " T1"),
  c("mg_p", "Body fat %"),
  c("nreads", "seq. depth"),
  c("sauces_con", "sauces, con"),
  c("juices_so", "juices, so"),
  c("fats_oils", "fats, oils"),
  c("sweets_pas", "sweets, pas"),
  c("^z_", "Z-score "),
  c("bmi", "BMI"),
  c("status_c2", "status T1"),
  c("\\.", " "),
  c("_", " ")
)

replace_list <- c(
  purrr::map(replace_strings, \(x) c(x, "variable")),
  purrr::map(replace_strings, \(x) c(x, "groups")), 
  list(c("4-6-1-3", "4-6 vs 1-3", "comparison"), 
       c("7-10-1-3", "7-10 vs 1-3", "comparison"),
       c("7-10-4-6", "7-10 vs 4-6", "comparison"),
      c("h-", "h vs ", "comparison")
  )
)

merged_mod <- makeReplacementsDF(mergedtabs, replace_list, c("variable", "groups"))
write_tsv(mergedtabs, paste0(outdir, "merged_statistical_differences_mod1.tsv"))

test_at_t1 <- c("Z-Score BMI T1", "Z-Score Waist T1", "Status T1", "Change in Z-Score BMI", "Change in Z-Score Waist")
merged2 <- merged_mod %>% 
  dplyr::filter(variable %in% patnames) %>% 
  dplyr::filter( (dataset == "T0" & !grepl("T1", groups) ) | (dataset == "T1_normT0" & groups %in% test_at_t1)) %>% 
  dplyr::filter(! (groups == "Hospital" & type != "categorical_all")) %>% 
  dplyr::filter(! (groups == "Exercise Cat" & type == "numerical")) %>% # reoetidos porque aparecen tanto en tests cat como en numerical
  dplyr::filter(! (groups == "Status T0" & type == "numerical")) %>% 
  dplyr::filter(! (groups == "Mother Educ" & type == "numerical")) %>% 
  dplyr::filter(! (groups == "Sex" & type == "numerical")) %>% 
  dplyr::filter(! (groups == "Status T1" & type == "numerical")) %>% 
  dplyr::filter(!grepl("NaN", comparison)) %>% 
  dplyr::mutate(groups_mod = ifelse(dataset == "T0", groups, paste(groups, "(Normal T0)"))) %>% 
  
  dplyr::mutate(
    comp2show = ifelse(comparison == "all", groups_mod, paste(groups_mod, comparison, sep=":")), 
    padj_global = p.adjust(pval, method = "BH"),
    log10padj_global = -log10(padj_global), 
    Sig_global = ifelse(padj_global<=0.01, "p<0.01", ifelse(padj_global<="0.05", "p<0.05", "NS")),
    Sig_global = factor(Sig_global, levels = c("NS", "p<0.05", "p<0.01"))
  )

groups_order <- c("Age T0", "Sex", "Hospital", "Mother Educ", "Exercise Cat", "Seq  Depth", 
                  "Status T0", "Status T1 (Normal T0)", 
                  "Z-Score BMI T0", "Z-Score BMI T1 (Normal T0)", "Z-Score Waist T0", "Z-Score Waist T1 (Normal T0)", 
                  "Change in Z-Score BMI", "Change in Z-Score BMI (Normal T0)", "Change in Z-Score Waist",
                  "Change in Z-Score Waist (Normal T0)")

merged3 <- merged2 %>% 
  dplyr::mutate(groups_mod = factor(groups_mod, levels=groups_order)) %>% 
  arrange(desc(groups_mod), desc(comparison) )

complevs <- merged3 %>% filter(variable == "Preprocessed") %>% 
  pull(comp2show)

merged4 <- merged3 %>% dplyr::mutate(comp2show=factor(comp2show, levels=complevs))

## faltaria ordenar correctamente las variables
g0 <- ggplot(merged4, aes(x=comp2show, y=log10padj_global, col=Sig_global)) + 
  facet_grid(. ~ variable ) +
  geom_hline(yintercept = -log10(0.05), linetype=2, col="gray") +
  scale_color_manual(values=c("gray", "firebrick4", "firebrick1")) +
  #geom_line() +
  geom_point() +
  coord_flip() + 
  ylab("-log10(adj. p-value)") + 
  xlab("Comparison") +
  theme_bw()
ggsave(filename = paste0(outdir, "dotplot1.pdf"), plot = g0, width = 8, height = 5)
write_tsv(merged2, paste0(outdir, "merged_statistical_differences_mod2.tsv"))
write_tsv(merged4, paste0(outdir, "merged_statistical_differences_mod4.tsv"))


library(ggpmisc)
library(broom)
use_vars_all <- c("z_t0") #, "z_waist_00"

use_vars_t1 <- c("z_t1", "inc_z_bmi") # "z_waist_01", "inc_z_waist"

var_levorder <- c("Z-Score BMI T0",
                  #"Z-Score Waist T0",
                  "Z-Score BMI T1 (Normal T0)",
                  #"Z-Score Waist T1 (Normal T0)",
                  "Change in Z-Score BMI (Normal T0)"
                  #"Change in Z-Score Waist  (Normal T0)"
                  )

df_bmi_all <- s_meta %>% select(sampleID, all_of(use_vars_all), all_of(patnames), status_c1, status_c2) %>% 
  gather("variable", "value", all_of(use_vars_all)) %>% 
  gather("pattern", "pattern_value", all_of(patnames)) 

df_bmi_t1 <- s_meta %>% select(sampleID, all_of(use_vars_t1), all_of(patnames), status_c1, status_c2) %>% 
  filter(status_c1 == "normal") %>% 
  gather("variable", "value", all_of(use_vars_t1)) %>% 
  gather("pattern", "pattern_value", all_of(patnames)) %>% 
  dplyr::mutate(variable = paste(variable, "(Normal T0)"))

df_bmi <- rbind(df_bmi_all, df_bmi_t1)

replace_list_df <- purrr::map(replace_strings, \(x) c(x, "variable"))

df_bmi <- makeReplacementsDF(df_bmi, replace_list_df, c("variable"))

regrs <- df_bmi %>% 
  group_by(variable, pattern) %>% 
  group_map( \(x, y) {
    
    m <- lm(pattern_value ~ value, data = x) 
    m %>% broom::glance() %>% cbind(y) %>% 
      dplyr::mutate(Coef = m$coefficients[2])
    
    }) %>% 
  bind_rows() %>% 
  dplyr::mutate(padj = p.adjust(p.value, method = "BH")) %>% 
  dplyr::mutate(sig = ifelse(padj<= 0.05, ifelse(Coef < 0, "Neg", "Pos"), "NS")) %>% 
  dplyr::mutate(sig = factor(sig, levels=c("Neg", "Pos", "NS"))) %>% 
  dplyr::mutate(labs = ifelse(padj < 0.001, "p<0.001", paste0("p=", as.character(round(padj, 3)))))

aux <- merged4 %>% filter(comp2show %in% regrs$variable) %>% 
  dplyr::select(pattern = variable, 
                variable = comp2show, 
                orig_log10padj_global = log10padj_global, 
                orig_padj_global= padj_global, orig_pval = pval)
nrow(regrs)
regrs <- merge(regrs, aux, by=c("pattern", "variable"), all.x=TRUE, all.y=TRUE)%>% 
  dplyr::mutate(sig_orig = ifelse(orig_padj_global<= 0.05, ifelse(Coef < 0, "Neg", "Pos"), "NS")) %>% 
  dplyr::mutate(sig_orig = factor(sig_orig, levels=c("Neg", "Pos", "NS"))) %>% 
  dplyr::mutate(labs_orig = ifelse(orig_padj_global < 0.001, "p<0.001", paste0("p=", as.character(round(orig_padj_global, 3)))))
nrow(regrs)
regrs <-regrs %>% dplyr::mutate(variable = factor(variable, levels=var_levorder))

df_bmi_final <- df_bmi %>% merge(regrs, by = c("variable", "pattern")) %>% 
  dplyr::mutate(variable = factor(variable, levels=var_levorder))

(g0 <- ggplot(df_bmi_final, aes(x=value, y=pattern_value, col=sig_orig)) + 
  facet_grid(pattern ~ variable) + 
  geom_point(size=0.1, alpha=0.5, col="gray") + 
  # stat_poly_eq(use_label(c("p")), #c("eq", "R2", "f", "p", "n")
  #              method="lm", small.p=T, small.r=F, label.y=0.99, label.x=0.8)+
  geom_smooth(method="lm") + 
  geom_text(data = regrs, aes(x = -Inf, y = Inf, label = labs_orig, col=sig_orig),
            hjust = -0.1,  
            vjust = 1.5,
            inherit.aes = FALSE, size=4) +
  scale_color_manual(values=c("tomato", "steelblue", "darkgray")) +
  theme_classic() + 
  ylab("Diet Pattern Score") +
  xlab("Variable Value") +
  theme(strip.text = element_text(size = 12)) +
    theme(axis.title = element_text(size = 14))
)
ggsave(paste0(outdir, "regressions_dietpatterns_bmi_alldata.pdf"), g0, width = 10.5, height = 6)
ggsave(paste0(outdir, "regressions_dietpatterns_bmi_alldata.png"), g0, width = 10.5, height = 6)
write_tsv(x = regr, file = paste0(outdir, "regressions_dietpatterns_bmi_alldata_REGR.tsv"))
write_tsv(x = df_bmi_final, file = paste0(outdir, "regressions_dietpatterns_bmi_alldata_DF.tsv"))

## remove outliers
qlim <- 0.01

df_bmi_filt <- df_bmi %>% group_by(variable, pattern) %>% 
  group_split() %>% 
  map(\(x){
    qs_pat <- quantile(x$pattern_value, c(qlim, 1-qlim), na.rm=TRUE)
    qs_val <- quantile(x$value, c(qlim, 1-qlim), na.rm=TRUE)
    
    x %>% 
      dplyr::filter(pattern_value > qs_pat[1] & pattern_value < qs_pat[2]) %>% 
      dplyr::filter(value > qs_val[1] & value < qs_val[2])
    
  }) %>% bind_rows()

regrs_filt <- df_bmi_filt %>% 
  group_by(variable, pattern) %>% 
  group_map( \(x, y) {
    
    m <- lm(pattern_value ~ value, data = x) 
    m %>% broom::glance() %>% cbind(y) %>% 
      dplyr::mutate(Coef = m$coefficients[2])
    
  }) %>% 
  bind_rows() %>% 
  dplyr::mutate(padj = p.adjust(p.value, method = "BH")) %>% 
  dplyr::mutate(sig = ifelse(padj<= 0.05, ifelse(Coef < 0, "Neg", "Pos"), "NS")) %>% 
  dplyr::mutate(sig = factor(sig, levels=c("Neg", "Pos", "NS"))) %>% 
  dplyr::mutate(labs = ifelse(padj < 0.001, "p<0.001", paste0("p=", as.character(round(padj, 3))))) %>% 
  dplyr::mutate(variable=factor(variable, levels=var_levorder))

nrow(regrs_filt)
regrs_filt <- merge(regrs_filt, aux, by=c("pattern", "variable"), all.x=TRUE, all.y=TRUE)%>% 
  dplyr::mutate(sig_orig = ifelse(orig_padj_global<= 0.05, ifelse(Coef < 0, "Neg", "Pos"), "NS")) %>% 
  dplyr::mutate(sig_orig = factor(sig_orig, levels=c("Neg", "Pos", "NS"))) %>% 
  dplyr::mutate(labs_orig = ifelse(orig_padj_global < 0.001, "p<0.001", paste0("p=", as.character(round(orig_padj_global, 3)))))
nrow(regrs_filt)

df_bmi_final_filt <- df_bmi_filt %>% merge(regrs_filt, by = c("variable", "pattern")) %>% 
  dplyr::mutate(variable = factor(variable, levels=var_levorder))

g1 <- ggplot(df_bmi_final_filt, aes(x=value, y=pattern_value, col=sig_orig)) + 
  facet_grid(pattern ~ variable) + 
  geom_point(size=0.1, alpha=0.5, col="gray") + 
  # stat_poly_eq(use_label(c("p")), #c("eq", "R2", "f", "p", "n")
  #              method="lm", small.p=T, small.r=F, label.y=0.99, label.x=0.8)+
  geom_smooth(method="lm") + 
  geom_text(data = regrs_filt, aes(x = -Inf, y = Inf, label = labs_orig, col=sig_orig),
            hjust = -0.1,  
            vjust = 1.5,
            inherit.aes = FALSE, size=5) +
  scale_color_manual(values=c("tomato", "steelblue", "darkgray")) +
  theme_classic() + 
  ylab("Diet Pattern Score") +
  xlab("Variable Value") +
  theme(strip.text = element_text(size = 12)) +
  theme(axis.title = element_text(size = 14))

ggsave(paste0(outdir, "regressions_dietpatterns_bmi_filt", as.character(qlim), ".pdf"), g1, width = 11, height = 6)
write_tsv(regrs_filt, paste0(outdir, "regressions_dietpatterns_bmi_filt", as.character(qlim), "_REGR.tsv"))
write_tsv(df_bmi_final_filt, paste0(outdir, "regressions_dietpatterns_bmi_filt", as.character(qlim), "_DF.tsv"))


### Now boxplots
library(ggpubr)

cols2 <- c( "#EE0000FF", "#3B4992FF", "#008B45FF")

df_educ_all <- s_meta %>% select(sampleID, mother_educ, all_of(patnames), status_c1, status_c2) %>% 
  filter(!is.na(mother_educ)) %>% 
  gather("pattern", "pattern_value", all_of(patnames)) %>% 
  dplyr::mutate(pattern = factor(pattern, levels = c("Mediterranean", "Preprocessed", "Western")))

combs <- combn(unique(df_educ_all$mother_educ), 2, simplify = FALSE)

gh <- ggplot(df_educ_all, aes(x=mother_educ, y=pattern_value, col=pattern)) + 
  facet_grid(. ~ pattern ) + 
  #geom_point(size=0.2, alpha=0.5, position=position_jitterdodge()) +
  geom_violin(alpha=0) +
  geom_boxplot(width=0.2, col="black", fill="gray", outliers = FALSE, alpha=0.5) +
  
    theme_classic() + 
    ylab("Diet Pattern Score") +
    xlab("Mother Education (years)") +
    theme(strip.text = element_text(size = 12)) +
    theme(axis.title = element_text(size = 14)) +
  theme(axis.text = element_text(size = 14)) +
    scale_color_manual(values=cols2)+
    scale_fill_manual(values=cols2) +
  stat_compare_means(
    vjust=0.5,
    hide.ns = TRUE,
    tip.length = 0.01,
    method = "wilcox.test",
    comparisons = combs,  # <-- replace with your groups
    label = "p.signif"   # or "p.format" for numeric p-values
  )
ggsave(paste0(outdir, "boxplot_dietpatterns_educ", ".pdf"), gh, width = 8, height = 3)
write_tsv(df_educ_all, paste0(outdir, "boxplot_dietpatterns_educ", ".tsv"))


df_wx_all <- s_meta %>% dplyr::select(sampleID, exercise_cat, all_of(patnames), status_c1, status_c2) %>% 
  filter(!is.na(exercise_cat)) %>% 
  gather("pattern", "pattern_value", all_of(patnames)) %>% 
  filter(!is.na(pattern_value)) %>% 
  dplyr::mutate(pattern = factor(pattern, levels = c("Mediterranean", "Preprocessed", "Western"))) # %>% 
  #dplyr::mutate(exercise_cat=as.character(exercise_cat))

combs <- combn(levels(df_wx_all$exercise_cat), 2, simplify = FALSE)

gh <- ggplot(df_wx_all, aes(x=exercise_cat, y=pattern_value, col=pattern)) + 
  facet_grid(. ~ pattern ) + 
  #geom_point(size=0.2, alpha=0.5, position=position_jitterdodge()) +
  geom_violin(alpha=0) +
  geom_boxplot(width=0.2, col="black", fill="gray", outliers = FALSE, alpha=0.5) +
  
  theme_classic() + 
  ylab("Diet Pattern Score") +
  xlab("Weekly Exercise") +
  theme(strip.text = element_text(size = 12)) +
  theme(axis.title = element_text(size = 14)) +
  theme(axis.text = element_text(size = 14)) +
  scale_color_manual(values=cols2)+
  scale_fill_manual(values=cols2) +
  ylim(0, 3) +
  stat_compare_means(
    vjust=0.5,
    hide.ns = TRUE,
    tip.length = 0.01,
    method = "wilcox.test",
    comparisons = combs,  # <-- replace with your groups
    label = "p.signif"   # or "p.format" for numeric p-values
  )
ggsave(paste0(outdir, "boxplot_dietpatterns_exercise", ".pdf"), gh, width = 8, height = 4)
write_tsv(df_educ_all, paste0(outdir, "boxplot_dietpatterns_exercise", ".tsv"))




(gl <- ggplot(s_meta %>% dplyr::mutate(exercise=af_extraesc_m_00/60), aes(x=edu_m_00, y=exercise)) + 
  geom_point(size=0.2, alpha=0.5) +
  geom_smooth(method="lm") +
   stat_poly_eq(use_label(c("p")), #c("eq", "R2", "f", "p", "n")
                method="lm", small.p=T, small.r=F, label.y=0.99, label.x=0.9)+
  theme_classic() + 
  ylab("Weekly Exercise (hours)") +
  xlab("Mother Education (years)") +
  theme(strip.text = element_text(size = 12)) +
  theme(axis.title = element_text(size = 14)) +
  theme(axis.text = element_text(size = 14)) +
  scale_color_manual(values=cols2)+
  scale_fill_manual(values=cols2) +
    ylim(0, 8)
)
ggsave(paste0(outdir, "regr_mothereduc_exercise", ".pdf"), gl, width = 5, height = 3)

## plot with Hospital
df_hosp_all <- s_meta %>% select(sampleID, hospital, all_of(patnames), status_c1, status_c2) %>% 
  filter(!is.na(hospital)) %>% 
  gather("pattern", "pattern_value", all_of(patnames)) %>% 
  dplyr::mutate(pattern = factor(pattern, levels = c("Mediterranean", "Preprocessed", "Western"))) %>% 
  dplyr::mutate(hospital = gsub("Santiago_", "", hospital))

combs <- combn(unique(df_hosp_all$hospital), 2, simplify = FALSE)

(gh <- ggplot(df_hosp_all, aes(x=hospital, y=pattern_value, col=pattern)) + 
  facet_grid(. ~ pattern ) + 
  #geom_point(size=0.2, alpha=0.5, position=position_jitterdodge()) +
  geom_violin(alpha=0) +
  geom_boxplot(width=0.2, col="black", fill="gray", outliers = FALSE, alpha=0.5) +
  
  theme_classic() + 
  ylab("Diet Pattern Score") +
  xlab("Hospital") +
  theme(strip.text = element_text(size = 12)) +
  theme(axis.title = element_text(size = 14)) +
  theme(axis.text = element_text(size = 14)) +
  theme(axis.text.x = element_text(size = 14, angle=45, vjust = 1, hjust=1)) +
  scale_color_manual(values=cols2)+
  scale_fill_manual(values=cols2) #+
   #stat_compare_means(
   #  vjust=0.5,
   #  hide.ns = TRUE,
   #  tip.length = 0.01,
   #  method = "wilcox.test",
   #  comparisons = combs,  # <-- replace with your groups
   #  label = "p.signif"   # or "p.format" for numeric p-values
   #)
)
ggsave(paste0(outdir, "boxplot_dietpatterns_hospital", ".pdf"), gh, width = 10, height = 3)
write_tsv(df_hosp_all, paste0(outdir, "boxplot_dietpatterns_hospital", ".tsv"))

## plot with Age
df_age_all <- s_meta %>% select(sampleID, age_T0, all_of(patnames), status_c1, status_c2) %>% 
  filter(!is.na(age_T0)) %>% 
  gather("pattern", "pattern_value", all_of(patnames)) %>% 
  dplyr::mutate(pattern = factor(pattern, levels = c("Mediterranean", "Preprocessed", "Western"))) 


(gh <- ggplot(df_age_all, aes(x=age_T0, y=pattern_value, col=pattern)) + 
    facet_grid(. ~ pattern ) + 
    stat_poly_eq(use_label(c("p")), #c("eq", "R2", "f", "p", "n")
                 method="lm", small.p=T, small.r=F, label.y=0.99, label.x=0.9) +
    geom_smooth(method="lm", alpha=0.5) +
    geom_point(size=0.2, alpha=0.5) +
    theme_classic() + 
    ylab("Diet Pattern Score") +
    xlab("Age at T0") +
    theme(strip.text = element_text(size = 12)) +
    theme(axis.title = element_text(size = 14)) +
    theme(axis.text = element_text(size = 14)) +
    scale_color_manual(values=cols2)+
    scale_fill_manual(values=cols2) #+
  #stat_compare_means(
  #  vjust=0.5,
  #  hide.ns = TRUE,
  #  tip.length = 0.01,
  #  method = "wilcox.test",
  #  comparisons = combs,  # <-- replace with your groups
  #  label = "p.signif"   # or "p.format" for numeric p-values
  #)
)
ggsave(paste0(outdir, "regr_dietpatterns_age", ".pdf"), gh, width = 10, height = 3)
ggsave(paste0(outdir, "regr_dietpatterns_age", ".png"), gh, width = 10, height = 3)
write_tsv(df_age_all, paste0(outdir, "boxplot_dietpatterns_age", ".tsv"))


## plot again diet patterns

load("/home/carlos/Documentos/CORALS/results_rstudio/results_Abril25_2/foodPCA/NMF_food.RData")


coef <- coef(nmf_res)        
apply(coef, 2, which.max) %>% sort

pattern_contrib <- apply(coef, 2, which.max)
patnames <-   c(  "Preprocessed", "Western", "Mediterranean") #c("Mediterranean", "Sugars", "Preprocessed")  
rownames(coef) <- patnames
coefdf <- coef %>% t %>% data.frame %>% rownames_to_column("food_group") %>% 
  dplyr::mutate(Main_pattern = patnames[pattern_contrib]) %>% 
  dplyr::mutate(
    food_group = gsub("_g$", " (g)", food_group, perl=T),
    food_group = gsub("_kcal", " (Kcal)", food_group),
    food_group = gsub("_prop", "", food_group),
    food_group = gsub("_c1", " T0", food_group),
    food_group = gsub("_00", " T0", food_group),
    food_group = gsub("nreads", "seq. depth", food_group),
    food_group = gsub("sauces_con", "sauces, con", food_group),
    food_group = gsub("juices_so", "juices, so", food_group),
    food_group = gsub("fats_oils", "fats, oils", food_group),
    food_group = gsub("sweets_pas", "sweets, pas", food_group),
    food_group = gsub("^z_", "Z-score ", food_group, perl=T),
    food_group = gsub("bmi", "BMI", food_group, perl=T),
    food_group = gsub("status_c2", "status T1 (norm. w. T0)", food_group, perl=T),
    food_group = gsub("_", " ", food_group)
  ) %>% 
  dplyr::mutate(food_group = tools::toTitleCase(food_group)) %>% 
  dplyr::mutate(main_score = apply(coef, 2, max)) %>%  
  group_by(Main_pattern) %>% 
  arrange(desc(Main_pattern), main_score) %>% 
  dplyr::mutate(food_group = factor(food_group, levels=food_group)) %>% 
  gather("Food pattern", "Score", all_of(patnames))

linedf <- coefdf %>% 
  group_by(Main_pattern) %>% 
  dplyr::summarise(pos = max(as.numeric(food_group))+0.5) %>% 
  head(2)

(g0 <-ggplot(coefdf, aes(x=food_group, y=Score, col = `Main_pattern` ))+ #, fill = `Main_pattern`
    facet_grid(~ `Food pattern`) +
    geom_col(fill="gray", alpha=0.5) + # 
    theme_bw() +
    geom_vline(xintercept = linedf$pos, linetype=2, col="gray") +
    coord_flip() +
    xlab("food group") +
    scale_color_manual(values=cols2) +
    theme_classic() +
    theme(
      axis.text = element_text(size = 14),
      strip.text = element_text(size = 14) 
    )
)
ggsave(filename = paste0(outdir, "food_patterns_Mod.pdf"), g0, width = 8, height = 5)

#### Now make a few alpha div plots

## get statistics 
alpha_indices <- c("Observed", "Chao1", "Shannon", "InvSimpson")
alphastats <- read_tsv("/home/carlos/Documentos/CORALS/results_rstudio/results_Agosto25_1/AlphaDiversity/alpha_div_models_merged_2plot3_250814.tsv")
divtab <- read_tsv("/home/carlos/Documentos/CORALS/results_rstudio/results_Agosto25_1/AlphaDiversity/remove_tanda2_rarefied_min_AlphaDiv.tsv")

divlong <- divtab %>% gather("alpha_index", "alpha_value", all_of(alpha_indices) ) %>% 
  gather("pattern", "pattern_value", all_of(patnames)) %>% 
  dplyr::mutate(pattern = factor(pattern, levels = c("Mediterranean", "Preprocessed", "Western"))) %>% 
  dplyr::mutate(alpha_index = factor(alpha_index, levels=alpha_indices))

regr_txt <- alphastats %>% filter(Variable %in% patnames) %>% 
  dplyr::mutate(str2show = ifelse( `Pr(>F)` < 0.05, paste("p=", round(`Pr(>F)`, 3), ", p. adj=", round(padj_groupInd, 2), sep=""), 
                                   paste("p=", round(`Pr(>F)`, 2), sep="") 
                                   )) %>% 
  select(pattern=Variable, alpha_index=Index, `Pr(>F)`, padj_groupInd, str2show) %>% 
  dplyr::mutate(pattern = factor(pattern, levels=patnames)) %>% 
  dplyr::mutate(alpha_index = factor(alpha_index, levels=alpha_indices)) %>% 
  dplyr::mutate(textcol = ifelse(padj_groupInd <= 0.05, as.character(pattern), "NS")) 

(g0 <- ggplot(divlong, aes(x=pattern_value, y=alpha_value, col=pattern)) + 
    facet_grid(alpha_index ~ pattern , scales="free") + 
    geom_point(size=0.2, alpha=0.5, col="gray") + 
    # stat_poly_eq(use_label(c("p")), #c("eq", "R2", "f", "p", "n")
    #              method="lm", small.p=T, small.r=F, label.y=0.99, label.x=0.8)+
    geom_smooth(method="lm") + 
    scale_color_manual(values=cols2) +
    geom_text(data = regr_txt, aes(x = -Inf, y = Inf, label = str2show),
              col= "#3B4992FF",
              hjust = -0.1,  
              vjust = 1.5,
              inherit.aes = TRUE, size=4) +
    scale_color_manual(values=c(cols2, "darkgray")) +
    theme_classic() + 
    ylab("Alpha Index Value") +
    xlab("Diet Pattern Score") +
    theme(strip.text = element_text(size = 14)) +
    theme(axis.title = element_text(size = 14)) +
    theme(axis.text = element_text(size = 12))
)
ggsave(paste0(outdir, "alpha_diversity_dietpatterns_2.pdf"), g0, width = 9, height = 7)
ggsave(paste0(outdir, "alpha_diversity_dietpatterns_2.png"), g0, width = 9, height = 7)
write_tsv(divlong, paste0(outdir, "alpha_diversity_dietpatterns.tsv"))
write_tsv(regr_txt, paste0(outdir, "alpha_diversity_dietpatterns_REGR.tsv"))

## now BMI
library(ggpmisc)
library(broom)
use_vars_all <- c("z_t0") #, "z_waist_00"

use_vars_t1 <- c("z_t1", "inc_z_bmi") # "z_waist_01", "inc_z_waist"

var_levorder <- c("Z-Score BMI T0",
                  #"Z-Score Waist T0",
                  "Z-Score BMI T1 (Normal T0)",
                  #"Z-Score Waist T1 (Normal T0)",
                  "Change in Z-Score BMI (Normal T0)"
                  #"Change in Z-Score Waist  (Normal T0)"
)

df_bmi_all_alpha <- divtab %>% select(sampleID, all_of(use_vars_all), all_of(alpha_indices), status_c1, status_c2) %>% 
  gather("variable", "value", all_of(use_vars_all)) %>% 
  gather("alpha_index", "alpha_value", all_of(alpha_indices)) 

df_bmi_t1_alpha <- divtab %>% select(sampleID, all_of(use_vars_t1), all_of(alpha_indices), status_c1, status_c2) %>% 
  filter(status_c1 == "normal") %>% 
  gather("variable", "value", all_of(use_vars_t1)) %>% 
  gather("alpha_index", "alpha_value", all_of(alpha_indices)) %>% 
  dplyr::mutate(variable = paste(variable, "(Normal T0)"))

df_bmi_alpha <- rbind(df_bmi_all_alpha, df_bmi_t1_alpha) %>% 
  dplyr::mutate(colby = ifelse(grepl("Normal T0", variable), status_c2, status_c1)) %>% 
  dplyr::mutate(alpha_index = factor(alpha_index, levels = alpha_indices))

replace_list_df <- purrr::map(replace_strings, \(x) c(x, "variable"))

df_bmi_alpha <- makeReplacementsDF(df_bmi_alpha, replace_list_df, c("variable"))
newnames <- c("Z-Score BMI T0", "Z-Score BMI T1 (Normal T0)", "Change in Z-Score BMI (Normal T0)")   

regr_txt <- alphastats %>% filter(Variable %in% newnames) %>% 
  dplyr::mutate(str2show = ifelse( `Pr(>F)` < 0.05, paste("p=", round(`Pr(>F)`, 3), ", p. adj=", round(padj_groupInd, 2), sep=""), 
                                   paste("p=", round(`Pr(>F)`, 2), sep="") 
  )) %>% 
  select(pattern=Variable, alpha_index=Index, `Pr(>F)`, padj_groupInd, str2show) %>% 
  dplyr::mutate(pattern = factor(pattern, levels=newnames)) %>% 
  dplyr::mutate(alpha_index = factor(alpha_index, levels=alpha_indices)) %>% 
  dplyr::mutate(textcol = ifelse(padj_groupInd <= 0.05, "Sig.", "NS"))

(g0 <- ggplot(df_bmi_alpha, aes(x=value, y=alpha_value)) + 
    facet_grid(alpha_index ~ pattern, scales = "free") + 
    geom_point(size=0.2, alpha=0.5, col="gray") + # , aes(col=colby)
    # stat_poly_eq(use_label(c("p")), #c("eq", "R2", "f", "p", "n")
    #              method="lm", small.p=T, small.r=F, label.y=0.99, label.x=0.8)+
    geom_smooth(method="lm") + 
    geom_text(data = regr_txt, aes(x = -Inf, y = Inf, label = str2show, colour =textcol),
              hjust = -0.1,  
              vjust = 1.5,
              inherit.aes = TRUE, size=4) +
    ggsci::scale_color_aaas() +
    theme_classic() + 
    ylab("Alpha Index Value") +
    xlab("Variable Value") +
    theme(strip.text = element_text(size = 10)) +
    theme(axis.title = element_text(size = 14))
)
ggsave(paste0(outdir, "regressions_alphadiv_bmi_alldata.pdf"), g0, width = 9, height = 6)
ggsave(paste0(outdir, "regressions_alphadiv_bmi_alldata.png"), g0, width = 9, height = 6)
write_tsv(x = regr_txt, file = paste0(outdir, "regressions_alphadiv_bmi_alldata_REGR.tsv"))
write_tsv(x = df_bmi_alpha, file = paste0(outdir, "regressions_alphadiv_bmi_alldata_DF.tsv"))

ggplot(divtab, aes(x=Observed, y=Shannon, col=status_c1))+
  geom_point(size=0.5) +
  stat_ellipse() +
  theme_classic( )
  
## Age
divlong <- divtab %>% gather("alpha_index", "alpha_value", all_of(alpha_indices) ) %>% 
  dplyr::mutate(alpha_index = factor(alpha_index, levels=alpha_indices))

regr_txt <- alphastats %>% filter(Variable == "Age T0") %>% 
  dplyr::mutate(str2show = ifelse( `Pr(>F)` < 0.05, paste( ifelse(`Pr(>F)` < 0.001, "p<0.001", 
                                                                  paste("p=", round(`Pr(>F)`, 3), sep="")), 
                                                           ", p. adj=", round(padj_groupInd, 2), sep=""), 
                                   paste("p=", round(`Pr(>F)`, 2), sep="") 
  )) %>% 
  select(pattern=Variable, alpha_index=Index, `Pr(>F)`, padj_groupInd, str2show) %>% 
  dplyr::mutate(alpha_index = factor(alpha_index, levels=alpha_indices)) %>% 
  dplyr::mutate(textcol = ifelse(padj_groupInd <= 0.05, "Sig.", "NS"))

(g0 <- ggplot(divlong, aes(x=age_T0, y=alpha_value)) + 
    facet_wrap(. ~ alpha_index, scales="free", ncol=4) + 
    geom_point(size=0.2, alpha=0.5, col="gray") + 
    # stat_poly_eq(use_label(c("p")), #c("eq", "R2", "f", "p", "n")
    #              method="lm", small.p=T, small.r=F, label.y=0.99, label.x=0.8)+
    geom_smooth(method="lm") + 
    geom_text(data = regr_txt, aes(x = -Inf, y = Inf, label = str2show, colour =textcol),
              hjust = -0.1,  
              vjust = 1.5,
              inherit.aes = TRUE, size=4) +
    ggsci::scale_color_aaas() +
    theme_classic() + 
    ylab("Alpha Index Value") +
    xlab("Age at T0") +
    theme(strip.text = element_text(size = 14)) +
    theme(axis.title = element_text(size = 14)) +
    theme(axis.text = element_text(size = 12))
)
ggsave(paste0(outdir, "alpha_diversity_age.pdf"), g0, width = 10, height = 2.5)
ggsave(paste0(outdir, "alpha_diversity_age.png"), g0, width = 10, height = 2.5)
write_tsv(divlong, paste0(outdir, "alpha_diversity_age.tsv"))
write_tsv(regr_txt, paste0(outdir, "alpha_diversity_age_REGR.tsv"))
