

opt <- restaurar(opt)
#load(paste0(opt$out, "foodPCA/phyloseq_list_foodPCA_withNMF.RData"))
load("/home/carlos/Documentos/CORALS/results_rstudio/results_Agosto25_1//foodPCA/phyloseq_list_foodPCA_withNMF_withInc.RData")

deseqname <- "DESeq_FoodVars_Patterns_250814/"
deseqname_normal <- "DESeq_FoodVars_Patterns_onlyNormal_250814/"

s_meta <- sample_data(all_phyloseq$remove_tanda2) %>% data.frame

food_variables<-c(names(s_meta)[5], names(s_meta)[grep("_clr$", names(s_meta))]) # names(sample_data(all_phyloseq$remove_tanda2))[5:30]
patnames <- c("Preprocessed", "Mediterranean", "Western")
quant_vars <- c("age_T0", "af_extraesc_m_00",
                "z_bmi_00", "z_bmi_01", "inc_z_bmiCOR", # calc by collabs
                "z_t0", "z_t1","inc_z_bmi", # calc by us
                "mg_p_00", "z_waist_00", "z_waist_01","inc_z_waist",
                "nreads"
                )
cat_vars <-  c("pattern", "Sex", "hospital", "age_class1", "mother_educ", "exercise_cat",
              "Category_T0", "Category_T1",  # calc by collabs
              "status_c2", "status_c1" # calc by us
                )
vars2deseq <- c(food_variables, patnames, cat_vars, quant_vars)



daa_all <- list()
daa_all_onlyNormalT0 <- list()

opt$mincount <- 1
phseq_to_use <-  "remove_tanda2" #names(all_phyloseq)#c("remove_tanda2", "rmbatch_tanda", "filt")
opt <- restaurar(opt)
for(phname in phseq_to_use){
  daa_all[[phname]] <- list()
  daa_all_onlyNormalT0[[phname]] <- list()
  for(var in vars2deseq[39:51]){
    cat("Doing DESeq2 Analysys for: ", phname, ", var=",var, "(", which(var==vars2deseq), " of ", length(vars2deseq), "), all data\n")
    phobj <- all_phyloseq[[phname]]
    samples <- sample_data(phobj)$sampleID[! is.na(sample_data(phobj)[, var ])]
    phobj_filt <- phyloseq::prune_samples(samples, phobj)

    daa_all[[phname]][[var]] <-deseq_full_pipeline(phobj_filt,
                                                   phname, var, opt, deseqname=deseqname)

    if(var == "status_c1") next
    cat("Doing DESeq2 Analysys for: ", phname, ", var=",var, "(", which(var==vars2deseq), " of ", length(vars2deseq), "), only normal at T0 \n")
    samples_normal <- sample_data(phobj_filt) %>% data.frame() %>%
      filter(status_c1 == "normal") %>% pull(sampleID)
    phobj_filt_normal <- phyloseq::prune_samples(samples_normal, phobj_filt)

    daa_all_onlyNormalT0[[phname]][[var]] <-deseq_full_pipeline(phobj_filt_normal,
                                                                phname, var, opt, deseqname=deseqname_normal)


  }
}
save(daa_all, file = paste0(opt$out, deseqname, "/DESEQ2_all_food_variables_Patterns_250814.RData"))
save(daa_all_onlyNormalT0, file = paste0(opt$out, deseqname, "/DESEQ2_all_food_variables_Patterns_onlyNormalT0_250814.RData"))
#load(paste0(opt$out, deseqname, "/DESEQ2_all_food_variables_CLR.RData"))

rm(daa_all)
rm(daa_all_onlyNormalT0)
#########################################################################################
## Including confounders

conf_sets <- list(
  SexHosAge = c("Sex",  "hospital", "age_T0"),
  SexHosAgeEdu = c("Sex",  "hospital", "age_T0", "mother_educ")
  #SexHosAgeKcal = c("Sex",  "hospital", "age_T0", "energy_kcal_log"),
  #SexHosAgeEduKcal = c("Sex",  "hospital", "age_T0", "edu_m_00", "energy_kcal_log")
)
vars2test <- c(
  patnames,
  "z_bmi_00", "z_bmi_01",
  "mg_p_00", "z_waist_00",
  "z_t0", "z_t1","inc_z_bmi", # calc by us
  "z_waist_00", "z_waist_01","inc_z_waist",
  food_variables
  #names(s_meta)[grep("_clr$", names(s_meta))],
)
vars2log <- c("energy_kcal")
vars2test <- c("edu_m_00")

opt <- restaurar(opt)
deseqname <- "DESeq_Patterns_withCovars_250818/"
#deseqname_normal <- "DESeq_Patterns_withCovars_onlyNormal_250814/"
#samples2use <- sample_data(all_phyloseq$remove_tanda2) %>% data.frame %>%
#  filter(status_c1 == "normal") %>%
#  filter(!is.na(status_c2)) %>%
#  filter(!is.nan(status_c2)) %>%
#  filter(status_c2 != "NaN") %>%
#  pull(sampleID)

daa_all_cov <- list()

opt$mincount <- 1
phseq_to_use <-  "remove_tanda2" #names(all_phyloseq)#c("remove_tanda2", "rmbatch_tanda", "filt")
for(phname in phseq_to_use){
  daa_all_cov[[phname]] <- list()
  phobj <- all_phyloseq[[phname]]
  sample_data(phobj)$mother_educ <- gsub("-", "_", as.character(sample_data(phobj)$mother_educ))

  for(v2log in vars2log) sample_data(phobj)[, paste0(v2log, "_log")] <- log(sample_data(phobj)[, v2log]+1)

  for(confn in names(conf_sets)[1]){
    conf <- conf_sets[[confn]]
    opt$out
    vname <-paste(phname, "CovsOnly", confn, sep="_", collapse="_")
    cat("Doing DESeq2 Analysys for: ", vname, "(", 0, " of ", length(vars2test), ")", "\n")
    start <- Sys.time()
    daa_all_cov[[phname]][[confn]]  <-deseq_full_pipeline(phobj,
                                                          vname,
                                                          conf,
                                                          opt,
                                                          deseqname=deseqname,
                                                          plot_all=0,
                                                          return_all = TRUE,
                                                          get_all_shrinkages = FALSE,
                                                          get_all_contrasts=TRUE)
    print( Sys.time() - start )
    for(var in vars2test){
      vars2deseq <- c(var, conf)
      vname <-paste(phname, var, confn, sep="_", collapse="_")

      cat("Doing DESeq2 Analysys for: ", vname, "(", which(var == vars2test), " of ", length(vars2test), ")", "\n")
      phobj_filt <- phobj
      for(v in vars2deseq){
        samples_notnas <- sample_data(phobj_filt) %>% data.frame %>%
          dplyr::filter(!is.na(!!sym(v))) %>%
          pull(sampleID)
        phobj_filt <- phyloseq::prune_samples(samples_notnas, phobj_filt)
      }

      start <- Sys.time()
      daa_all_cov[[phname]][[vname]]  <-deseq_full_pipeline(phobj_filt,
                                                           vname,
                                                           vars2deseq,
                                                           opt,
                                                           deseqname=deseqname,
                                                           plot_all=0,
                                                           return_all = TRUE,
                                                           get_all_shrinkages = FALSE,
                                                           get_all_contrasts=TRUE)
      print( Sys.time() - start )
    }
  }
}
#save(daa_all_cov, file = paste0(opt$out, deseqname, "/DESEQ2_FoodPatterns_covars_onlySexHosAge_onlyEducVars.RData"))
save(daa_all_cov, file = paste0(opt$out, deseqname, "/DESEQ2_FoodPatterns_covars_onlySexHosAge.RData"))
load( paste0(opt$out, deseqname, "/DESEQ2_FoodPatterns_covars_full.RData"))



#########################################################################################
## Including confounders only Normal at T0

conf_sets <- list(
  SexHosAge = c("Sex",  "hospital", "edad_00"),
  SexHosAgeEdu = c("Sex",  "hospital", "edad_00", "edu_m_00"),
  SexHosAgeEduCat = c("Sex",  "hospital", "edad_00", "mother_educ"),
  SexHosAgeKcal = c("Sex",  "hospital", "edad_00", "energy_kcal_log"),
  SexHosAgeEduKcal = c("Sex",  "hospital", "edad_00", "edu_m_00", "energy_kcal_log")

)
vars2test <- c(
  patnames,
  "status_c2",
  "z_bmi_00", "z_bmi_01",
  "mg_p_00", "z_waist_00",
  "z_t0", "z_t1","inc_z_bmi", # calc by us
  "z_waist_00", "z_waist_01","inc_z_waist"
  #names(s_meta)[grep("_clr$", names(s_meta))],
)
patnames <- c("Preprocessed", "Mediterranean", "Western")
vars2log <- c("energy_kcal")

opt <- restaurar(opt)
deseqname <- "DESeq_Patterns_withCovars_onlyNormalT0_250818/"
samples2use <- sample_data(all_phyloseq$remove_tanda2) %>% data.frame %>%
  filter(status_c1 == "normal") %>%
  filter(!is.na(status_c2)) %>%
  filter(!is.nan(status_c2)) %>%
  filter(status_c2 != "NaN") %>%
  pull(sampleID)

daa_all_cov <- list()

opt$mincount <- 1
phseq_to_use <-  "remove_tanda2" #names(all_phyloseq)#c("remove_tanda2", "rmbatch_tanda", "filt")
for(phname in phseq_to_use){
  daa_all_cov[[phname]] <- list()
  phobj <- all_phyloseq[[phname]]
  for(v2log in vars2log) sample_data(phobj)[, paste0(v2log, "_log")] <- log(sample_data(phobj)[, v2log]+1)

  for(confn in names(conf_sets)[1]){
    conf <- conf_sets[[confn]]
    cat(conf)
    for(var in vars2test){
      vars2deseq <- c(var, conf)
      vname <-paste(phname, var, confn, sep="_", collapse="_")

      cat("Doing DESeq2 Analysys for: ", vname, "\n")
      phobj_filt <- phyloseq::prune_samples(samples2use, phobj)
      for(v in vars2deseq){
        samples_notnas <- sample_data(phobj_filt) %>% data.frame %>%
          dplyr::filter(!is.na(!!sym(v))) %>%
          pull(sampleID)
        phobj_filt <- phyloseq::prune_samples(samples_notnas, phobj_filt)
      }
      start <- Sys.time()
      daa_all_cov[[phname]][[vname]] <-deseq_full_pipeline(phobj_filt,
                                                           vname,
                                                           vars2deseq,
                                                           opt,
                                                           deseqname=deseqname,
                                                           plot_all=0,
                                                           return_all = TRUE,
                                                           get_all_shrinkages = FALSE,
                                                           get_all_contrasts=TRUE)
     print( Sys.time() - start )
    }
  }
}
save(daa_all_cov, file = paste0(opt$out, deseqname, "/DESEQ2_FoodPatterns_covars_onlyNormalT0.RData"))
load( paste0(opt$out, deseqname, "/DESEQ2_FoodPatterns_covars_onlyNormalT0.RData"))


##33333333333#3###################################3
## Only normal at T0 (old way)

opt <- restaurar(opt)
deseqname <- "DESeq_normalT0/"
samples2use <- sample_data(all_phyloseq$remove_tanda2) %>% data.frame %>%
  filter(status_c1 == "normal") %>%
  filter(!is.na(status_c2)) %>%
  filter(!is.nan(status_c2)) %>%
  filter(status_c2 != "NaN") %>%
  pull(sampleID)

vars2deseq <- c("status_c2")

daa_all_normalt0 <- list()

opt$mincount <- 1
phseq_to_use <-  "remove_tanda2" #names(all_phyloseq)#c("remove_tanda2", "rmbatch_tanda", "filt")
for(phname in phseq_to_use){
  daa_all_normalt0[[phname]] <- list()
  for(var in vars2deseq){
    cat("Doing DESeq2 Analysys for: ", phname, "\n")
    phobj <- all_phyloseq[[phname]]
    phobj_filt <- phyloseq::prune_samples(samples2use, phobj)

    daa_all_normalt0[[phname]][[var]] <-deseq_full_pipeline(phobj_filt, phname, var, opt, deseqname=deseqname)
  }
}
#save(daa_all_normalt0, file = paste0(opt$out, deseqname, "/DESEQ2_normalT0.RData"))
load(paste0(opt$out, deseqname, "/DESEQ2_normalT0.RData"))

## Only normal at T0, with covariates

opt <- restaurar(opt)
deseqname <- "DESeq_normalT0_withCovars/"
samples2use <- sample_data(all_phyloseq$remove_tanda2) %>% data.frame %>%
  filter(status_c1 == "normal") %>%
  filter(!is.na(status_c2)) %>%
  filter(!is.nan(status_c2)) %>%
  filter(status_c2 != "NaN") %>%
  pull(sampleID)

vars2deseq <- c("status_c2", "edad_00", "Sex",  "mother_educ", "hospital")
vname <- paste0(toupper(substr(vars2deseq, 1, 1)), tolower(substr(vars2deseq, 2, 2)), collapse = "")

daa_all_normalt0_cov <- list()

opt$mincount <- 1
phseq_to_use <-  "remove_tanda2" #names(all_phyloseq)#c("remove_tanda2", "rmbatch_tanda", "filt")
for(phname in phseq_to_use){
  daa_all_normalt0_cov[[phname]] <- list()
    cat("Doing DESeq2 Analysys for: ", phname, "\n")
    phobj <- all_phyloseq[[phname]]
    phobj_filt <- phyloseq::prune_samples(samples2use, phobj)
    for(v in vars2deseq){
      samples_notnas <- sample_data(phobj_filt) %>% data.frame %>%
        dplyr::filter(!is.na(!!sym(v))) %>%
        pull(sampleID)
      phobj_filt <- phyloseq::prune_samples(samples_notnas, phobj_filt)
    }

    daa_all_normalt0_cov[[phname]][[vname]] <-deseq_full_pipeline(phobj_filt, phname, vars2deseq, opt, deseqname=deseqname)
}
#save(daa_all_normalt0_cov, file = paste0(opt$out, deseqname, "/DESEQ2_normalT0_covars.RData"))
load( paste0(opt$out, deseqname, "/DESEQ2_normalT0_covars.RData"))

#########################################################################3
# Make Heatmap with LFC all

dfmerged <- map2(daa_all$remove_tanda2, names(daa_all$remove_tanda2),
                 \(x, xname){
                   map2(x$all_contrasts, names(x$all_contrasts), \(y, yname){
                     y$resdf %>% dplyr::mutate(Comparison = yname)
                   }) %>% bind_rows() %>%
                     dplyr::mutate(Variable = xname)
                 }) %>% bind_rows %>%
  select(Variable, Comparison, everything()) %>%
  dplyr::filter(Variable != "status_c2") %>%
  dplyr::filter(Variable != "hospital")

mos <- dfmerged %>% filter(padj < 0.000001 & abs(log2FoldChangeShrink)>3) %>% pull(taxon) %>% unique

mergedfilt <- dfmerged %>% filter(taxon %in% mos) %>%
  mutate(gsub("_", " ", taxon))
taxord <- mergedfilt %>%
  group_by(taxon) %>%
  dplyr::summarise(maxl = max(log2FoldChangeShrink)) %>%
  arrange(maxl)

mergedfilt <- mergedfilt %>%
  dplyr::mutate(taxon=factor(taxon, levels=taxord$taxon)) %>%
  dplyr::mutate(Sig = ifelse(padj < 0.01,
                             ifelse(log2FoldChangeShrink < 0, "Down", "Up"),
                             "NS")) %>%
  dplyr::mutate(Sig = factor(Sig, levels = c("Down", "Up", "NS")))

ggplot(mergedfilt, aes(x=taxon, y=log2FoldChangeShrink, col=Sig, fill=Sig))+
  facet_grid(Comparison ~ .) +
  geom_point() +
  scale_color_manual(values = c("steelblue", "firebrick4", "gray")) +
  scale_fill_manual(values = c("steelblue", "firebrick4", "gray")) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, vjust=1, hjust=1),
        strip.text.y = element_text(angle = 0, vjust=1, hjust=1))


################3 Integrate

for(nn in names(reserva$remove_tanda2)){
  daa_all_cov$remove_tanda2[[nn]] <- reserva$remove_tanda2[[nn]]
}

PLIM=0.0001
PLIM_COL=0.05
LFCLIM=1

deseqname <- "DESeq_Patterns_withCovars/"
load(paste0(opt$out, deseqname, "/DESEQ2_FoodPatterns_covars_full2.RData"))

corr_used_short <- "SexHosAge"

for(corr_used_short in names(conf_sets)[2:5]){
corr_used <- paste0("_", corr_used_short, "$")

uselist <-daa_all_cov$remove_tanda2[grepl(corr_used, names(daa_all_cov$remove_tanda2), perl=T)]
usedf <- map(vars2test, \(vv){
  uu <- uselist[grepl(paste0("_", vv, corr_used), names(uselist), perl=T)]
  uu[[1]]$all_contrasts[[vv]]$resdf %>%
    dplyr::mutate(Variable = vv, Corrdataset = corr_used_short)
}) %>% bind_rows() %>%
  dplyr::mutate(
    Variable = gsub("_00", " T0", Variable),
    Variable = gsub("_01", " T1", Variable),
    Variable = gsub("mg_p", "Body fat %", Variable),
    Variable = gsub("^z_", "Z-score ", Variable, perl=T),
    Variable = gsub("bmi", "BMI", Variable, perl=T),
    Variable = gsub("_", " ", Variable)
  ) %>%
  dplyr::mutate(
    Variable = factor(Variable, levels=c(
      "Mediterranean",
      "Preprocessed",
      "Western",
      "Z-score BMI T0",
      "Z-score BMI T1",
      "Z-score waist T0",
      "Body fat % T0"
    ))
  )%>%
  dplyr::mutate(
    Sig = ifelse(padj <= PLIM_COL, ifelse(log2FoldChangeShrink < 0 , "Down", "Up"), "NS")
  )%>%
  dplyr::mutate(
    taxon_old = taxon,
    taxon = gsub("_", " ", taxon),
    taxon = gsub("[\\[\\]]", "", taxon, perl=TRUE)
  )

table(usedf$Variable)


tax2use <- usedf %>%
  filter(padj < PLIM & abs(log2FoldChangeShrink) > LFCLIM) %>%
  pull(taxon) %>% unique
tax2use %>% length


taxorder <- usedf %>%
  filter(taxon %in% tax2use) %>%
  select(taxon, log2FoldChangeShrink) %>%
  group_by(taxon) %>%
  dplyr::summarise(meanlfc=mean(log2FoldChangeShrink)) %>%
  arrange(meanlfc)

taxorder2 <- usedf %>%
  filter(taxon %in% tax2use) %>%
  select(taxon, Variable, log2FoldChangeShrink) %>%
  spread(Variable, log2FoldChangeShrink) %>%
  dplyr::mutate(tmp = Mediterranean -1*Preprocessed -1*Western) %>%
  arrange(Mediterranean) #tmp

taxorder3 <- usedf %>%
  filter(taxon %in% tax2use) %>%
  select(taxon, Variable, log2FoldChangeShrink, padj) %>%
  gather("vart", "valt", log2FoldChangeShrink, padj) %>%
  unite("vart2", Variable, vart, sep="__") %>%
  spread(vart2, valt) %>%
  dplyr::mutate(bmisig = ifelse(`Z-score BMI T0__padj` < PLIM_COL, "Sig", "NS")) %>%
  dplyr::mutate(bmisig = factor(bmisig, levels=c("NS", "Sig"))) %>%
  group_by(bmisig) %>%
  arrange(bmisig,Mediterranean__log2FoldChangeShrink) %>%
  dplyr::mutate(taxon=factor(taxon, levels=taxon))

TAXORDER <- taxorder3$taxon

usedf2 <- usedf %>%
  filter(taxon %in% tax2use) %>%
  dplyr::mutate(taxon=factor(taxon, levels=TAXORDER))

linedf <- taxorder3 %>%
  dplyr::summarise(xpos = max(as.numeric(taxon)) + 0.5) %>%
  head(nrow(.)-1)

(g0 <- ggplot(usedf2, aes(x=taxon, y=log2FoldChangeShrink, col=Sig, fill=Sig))+
  facet_grid(~ Variable) +
  geom_hline(yintercept = 0, linetype=2, col="lightgray") +
  geom_vline(data=linedf, aes(xintercept=linedf$xpos), linetype=2, col="lightgray")+
  geom_segment(aes(x=taxon, xend = taxon, y=0, yend=log2FoldChangeShrink)) +
  geom_point() +
  theme_bw() +
  coord_flip() +
  scale_color_manual(values = c("Down"="steelblue", "Up"="tomato", "NS"="darkgray")) +
  theme(axis.text.y= element_text(face="italic", size=10),
        axis.text.x= element_text( size=12),
        axis.title = element_text(size=12),
        strip.text = element_text(size=10))
)

ggsave(paste0(opt$out, deseqname, "/DESEQ2_FoodPatterns_covars_Lolliplot_", corr_used_short, ".pdf"), g0,
       height = 7, width = 12)

tax2use_oldname_pe4 <- usedf %>%
  filter(Variable %in% patnames) %>%
  filter(padj < PLIM & abs(log2FoldChangeShrink) > LFCLIM) %>%
  pull(taxon_old) %>% unique

tax2use_oldname_pe2 <- usedf %>%
    filter(Variable %in% patnames) %>%
  filter(padj < 0.01 & abs(log2FoldChangeShrink) > LFCLIM) %>%
    pull(taxon_old) %>% unique

df2pca <- daa_all_cov$remove_tanda2[grepl(corr_used, names(daa_all_cov$remove_tanda2), perl=T)]
df2pca <- df2pca[[1]]$vst_counts_df

pp1 <- makeAllPCAs(all_phyloseq$remove_tanda2, df2pca, tax2use_oldname_pe4,
            c("pattern", patnames),
            opt, paste0(corr_used_short, "_DiffTaxaDietPat_p1e4"))

pp2 <- makeAllPCAs(all_phyloseq$remove_tanda2, df2pca, tax2use_oldname_pe4,
            c("pattern", patnames),
            opt, paste0(corr_used_short, "_DiffTaxaDietPat_p1e2"))
}
