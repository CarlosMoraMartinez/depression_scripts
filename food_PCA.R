
## Add Diet PCA to metadata
source(opt$predictive_functions)
load("/home/carlos/Documentos/CORALS/results_rstudio/results_Abril25_2//phyloseq/phyloseq_all_list.RData")

outdir <- paste0(opt$out, "foodPCA_rmOutliers_Median/")
if(!dir.exists(outdir)) dir.create(outdir)

var2predict <- c("Category_T0", "status_c1")

for(phname in names(all_phyloseq)){
  phobj <- all_phyloseq[[phname]]

  #sample_data(phobj)$eggs[sample_data(phobj)$eggs > 60] <- 60
  this_metadata <- sample_data(phobj) %>% data.frame
  
  food_variables<- names(this_metadata)[5:30]
  
  this_metadata <- replaceAllOutliers(food_variables, this_metadata)
  
  food_mat <- this_metadata %>% select( all_of(food_variables)) %>% 
    #column_to_rownames("sampleID") %>% 
    as.matrix 
  
  no_food_samples <- food_mat %>% rowSums()
  no_food_samples <- no_food_samples[is.na(no_food_samples)] %>% names
  
  this_metadata <- this_metadata %>% filter(! sampleID %in% no_food_samples)
  
  
  food_PCA <- make_meta_PCA(this_metadata, food_variables, 
                            var2predict,
                            outdir,
                            make_log=TRUE,
                            name=paste0(phname, "foodVars_CorrOutl")
                            )
  pca_meta <- food_PCA$pca$x %>% as.data.frame 
  
  for(pc in names(pca_meta)){
    sample_data(all_phyloseq[[phname]])[, pc] <- pca_meta[sample_data(phobj)$sampleID, pc] 
  }

}
#### Replace outliers in food variables

for(phname in names(all_phyloseq)){
  phobj <- all_phyloseq[[phname]]
  
  #sample_data(phobj)$eggs[sample_data(phobj)$eggs > 60] <- 60
  this_metadata <- sample_data(phobj) %>% data.frame
  food_variables<- names(this_metadata)[5:30]
  this_metadata <- replaceAllOutliers(food_variables, this_metadata, p_lim = 0.01, funct_name = "median")
  sample_data(all_phyloseq[[phname]]) <- sample_data(this_metadata)
  
}

save(all_phyloseq, file = paste0(outdir, "/phyloseq_all_list_CorrOutl_median.RData"))
#### Transform food data to compositional

library(compositions)
library(NMF)

this_metadata <- sample_data(all_phyloseq$remove_tanda2_rarefied_min) %>% data.frame
total_name <- names(this_metadata)[5]
macro_names <- names(this_metadata)[6:9]
specific_names <- names(this_metadata)[c(10:28, 30)] #not water

for(phname in names(all_phyloseq)){
  phobj <- all_phyloseq[[phname]]
  this_metadata <- sample_data(phobj) %>% data.frame %>% 
    rowwise() %>% 
    dplyr::mutate(total_g_macro = sum(c_across(all_of(macro_names)))) %>% 
    dplyr::mutate(total_g_foods = sum(c_across(all_of(specific_names)))) %>% 
    ungroup() %>%
    dplyr::mutate(across(
      all_of(macro_names),
      ~ .x / total_g_macro,
      .names = "{.col}_prop"
    )) %>%
    dplyr::mutate(across(
      all_of(specific_names),
      ~ .x / total_g_foods,
      .names = "{.col}_prop"
    ))


  macro_prop_names <- paste0(macro_names, "_prop")
  spec_prop_names <- paste0(specific_names, "_prop")
  macro_prop_matrix <- this_metadata %>%
    select(all_of(macro_prop_names)) %>%
    as.matrix()
  spec_prop_matrix <- this_metadata %>%
    select(all_of(spec_prop_names)) %>%
    as.matrix()

  macro_prop_matrix[macro_prop_matrix == 0] <- 1e-6
  spec_prop_matrix[spec_prop_matrix == 0] <- 1e-6
  
  # Apply CLR transformation
  macro_clr <- compositions::clr(macro_prop_matrix)
  spec_clr <- compositions::clr(spec_prop_matrix)
  
  colnames(macro_clr) <- paste0(macro_names, "_clr")
  colnames(spec_clr) <- paste0(specific_names, "_clr")
  
  this_metadata <- bind_cols(this_metadata, as.data.frame(macro_clr))
  this_metadata <- bind_cols(this_metadata, as.data.frame(spec_clr))
  this_metadata <- as.data.frame(this_metadata)
  
  ## Now make nonnegative mattrix factorization
  #nnmat <- this_metadata %>% select(all_of(c(spec_prop_names)))
  #nas <- apply(nnmat, 1, \(x)any(is.na(x)))
  #nnmat <- nnmat[!nas, ]
  #
  #all_nmf <- list()
  #for(rank in 2:10){
  #  all_nmf[[rank]] <- nmf(nnmat, rank = rank, method = "brunet", nrun = 100)
  #}
  #coph <- map(all_nmf, \(x)summary(x)) %>% 
  #  bind_rows %>% 
  #  filter(!is.na(rank))
  #par(mfrow=c(1,3))
  #plot(coph$rank, coph$residuals)
  #plot(coph$rank, coph$cophenetic)
  #plot(coph$rank, coph$silhouette.basis)
#
  #nmf_res <- all_nmf[[4]] #nmf(nnmat, rank = 2, method = "brunet", nrun = 100)
  #
  #basis <- basis(nmf_res)         # "components" (food patterns)
  #coef <- coef(nmf_res)           # "weights" per subject
  #
  ## Assign dominant pattern per subject (optional)
  #dominant_pattern <- apply(coef, 2, which.max)
  #diet_df$pattern <- factor(dominant_pattern)
  #
  rownames(this_metadata) <- this_metadata$sampleID
  sample_data(phobj) <- sample_data(this_metadata)
  all_phyloseq[[phname]] <- phobj
  
  
  
}


### PCA with clr
phname <- "remove_tanda2_rarefied_min"
this_metadata <- sample_data(all_phyloseq$remove_tanda2_rarefied_min) %>% data.frame
food_variables<- names(this_metadata)[c(grep("_clr", names(this_metadata)))]

food_mat <- this_metadata %>% select( all_of(food_variables)) %>% 
  #column_to_rownames("sampleID") %>% 
  as.matrix 

no_food_samples <- food_mat %>% rowSums()
no_food_samples <- no_food_samples[is.na(no_food_samples)] %>% names

this_metadata <- this_metadata %>% filter(! sampleID %in% no_food_samples)

var2predict <- c("status_c1", "status_c2", "Category_T0", "Category_T1", "age_class1", "mother_educ", "Sex", "hospital")
food_PCA <- make_meta_PCA(this_metadata, food_variables, 
                          var2predict,
                          outdir,
                          make_log=FALSE,
                          name=paste0(phname, "_foodVarsCLR")
)
pca_meta <- food_PCA$pca$x %>% as.data.frame 


for(phname in names(all_phyloseq)){
  for(pc in names(pca_meta)){
    sample_data(all_phyloseq[[phname]])[, paste0("clr_", pc)] <- pca_meta[sample_data(all_phyloseq[[phname]])$sampleID, pc] 
  }
}

fname <- paste0(outdir, "phyloseq_list_foodPCA_CorrOutl_median.RData")
save(all_phyloseq, file=fname)
load(fname)

#modify exercise
for(phname in names(all_phyloseq)){
    metdf <- sample_data(all_phyloseq[[phname]]) %>% data.frame %>% 
      dplyr::mutate(exercise_cat = (af_extraesc_m_00/60) %>% floor,
             exercise_cat = ifelse(exercise_cat > 4, 4, exercise_cat ),
             exercise_cat = ifelse(exercise_cat == 0, "<1h", 
                                   ifelse(exercise_cat == 4, ">4h", paste0(exercise_cat, "h"))),
             exercise_cat = factor(exercise_cat, levels=c("<1h", "2h", "3h", ">4h"))
             ) 
    sample_data(all_phyloseq[[phname]])[, "exercise_cat"] <-metdf$exercise_cat

}

fname <- paste0(outdir, "phyloseq_list_foodPCA_CorrOutl_median.RData")
save(all_phyloseq, file=fname)
load(fname)

## Now make nonnegative mattrix factorization

library(NMF)
this_metadata <- sample_data(all_phyloseq$remove_tanda2) %>% data.frame
spec_prop_names <- names(this_metadata)[grepl("prop", names(this_metadata)) & ! grepl("_g_", names(this_metadata)) ]
nnmat <- this_metadata %>% select(all_of(c(spec_prop_names)))
nas <- apply(nnmat, 1, \(x)any(is.na(x)))
nnmat <- nnmat[!nas, ] %>% as.matrix 

#all_nmf <- list()
#for(rank in 2:10){
#  all_nmf[[rank]] <- nmf(nnmat, rank = rank, method = "brunet", nrun = 100)
#}
#coph <- map(all_nmf, \(x)summary(x)) %>% 
#  bind_rows %>% 
#  filter(!is.na(rank))
#par(mfrow=c(1,3))
#plot(coph$rank, coph$residuals)
#plot(coph$rank, coph$cophenetic)
#plot(coph$rank, coph$silhouette.basis)
#
#nmf_res <- all_nmf[[4]] #nmf(nnmat, rank = 2, method = "brunet", nrun = 100)

nmf_res <-  nmf(nnmat, rank = 3, method = "brunet", nrun = 1000)
basis <- basis(nmf_res)         # "components" (food patterns)
coef <- coef(nmf_res)           # "weights" per subject
save(nmf_res, file = paste0(outdir, "NMF_3factors_food_ModOutliersMedian.RData"))
apply(coef, 2, which.max) %>% sort
# Assign dominant pattern per subject (optional)
patnames <- c("Mediterranean", "Preprocessed", "Western") # Western

dominant_pattern <- apply(basis, 1, which.max)
this_metadata$pattern<- patnames[factor(dominant_pattern)[this_metadata$sampleID]]

nas <- this_metadata$sampleID[!this_metadata$sampleID %in% rownames(basis)]
nasmat <- matrix(nrow=length(nas), ncol=ncol(basis))
rownames(nasmat) <- nas
basis2 <- rbind(basis, nasmat)
for(i in 1:ncol(basis)){
  this_metadata[[patnames[i]]] <- basis2[this_metadata$sampleID, i]
}

sample_data(all_phyloseq$remove_tanda2) <- sample_data(this_metadata)

this_metadata <-sample_data(all_phyloseq$remove_tanda2_rarefied_min) %>% data.frame
this_metadata$pattern<-  patnames[factor(dominant_pattern)[this_metadata$sampleID]]
for(i in 1:ncol(basis)){
  this_metadata[[patnames[i]]] <- basis2[this_metadata$sampleID, i]
}
sample_data(all_phyloseq$remove_tanda2_rarefied_min) <- sample_data(this_metadata)


fname <- paste0(outdir, "phyloseq_list_foodPCA_withNMF_AdjOutliers.RData")
save(all_phyloseq, file=fname)

pattern_contrib <- apply(coef, 2, which.max)
patnames <-   c("Mediterranean", "Preprocessed", "Western") #c("Mediterranean", "Sugars", "Preprocessed")  
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
    dplyr::mutate(main_score = apply(coef, 2, max)) %>%  
    group_by(Main_pattern) %>% 
    arrange(main_score) %>% 
    dplyr::mutate(food_group = factor(food_group, levels=food_group)) %>% 
  gather("Food pattern", "Score", all_of(patnames))

linedf <- coefdf %>% 
  group_by(Main_pattern) %>% 
  dplyr::summarise(pos = max(as.numeric(food_group))+0.5) %>% 
  head(2)

(g0 <-ggplot(coefdf, aes(x=food_group, y=Score, col = `Main_pattern`, fill = `Main_pattern` ))+
  facet_grid(~ `Food pattern`) +
  geom_col(col="black") +
  theme_bw() +
  geom_vline(xintercept = linedf$pos, linetype=2, col="gray") +
  coord_flip() +
  xlab("food group") +
  ggsci::scale_fill_lancet() +
  theme(
    axis.text = element_text(size = 14),
    strip.text = element_text(size = 14) 
  )
)
ggsave(filename = paste0(outdir, "food_patterns1.pdf"), g0, width = 8, height = 6)

metlong <- this_metadata %>% 
  gather("Pattern_score", "Pattern_score_value", all_of(patnames)) 


library(ggsci)


(g1 <- ggplot(metlong, aes(x = age_T0, y=Pattern_score_value, 
                           col=Pattern_score, fill=Pattern_score)) +
    facet_grid(~ Pattern_score) +
    geom_point(, size=0.05) + 
    geom_smooth() +
    theme_bw() + 
    ggsci::scale_fill_lancet() +
    ggsci::scale_color_lancet() +
    theme(
      axis.text = element_text(size = 14),
      axis.title = element_text(size = 14),
      strip.text = element_text(size = 14) 
    ) + 
    ylab("Score")
  
  )

ggsave(filename = paste0(outdir, "food_patterns_vs_age.pdf"), g1, width = 8, height = 3)

othervars_cat <- c("Sex", "Category_T0",  "Category_T1", "mother_educ", 
                    "exercise_cat", "status_c2", "hospital")

othervars_quant <- c("edad_00",  "z_bmi_00", "z_bmi_01", 
               "nreads", "mg_p_00", "z_waist_00")

othervars <- c(othervars_cat, othervars_quant)

library(ggpmisc)
models1 <- makeLinearModelsSingleVariable(this_metadata, "age_T0", 
                                          othervars, 
                                          patnames, 
                                          combos=1,
                                          outdir = outdir, 
                                          name = paste0("FoodNMF") )

alphadif <- testDiversityDifferences(this_metadata, 
                                     patnames, 
                                     othervars, 
                                     outdir, "this_metadata")


qplots <- list()

for(var in othervars_quant){
  
  (qplots[[var]] <- ggplot(metlong, aes(x = !!sym(var), 
                             y=Pattern_score_value, 
                             col=Pattern_score, 
                             fill=Pattern_score)) +
     facet_grid(~ Pattern_score) +
     geom_point(, size=0.05) + 
     geom_smooth(method = "lm") +
     stat_poly_eq(use_label(c("R2", "p")), #c("eq", "R2", "f", "p", "n")
                  method="lm", small.p=T, small.r=F, label.y=0.99) +
     theme_bw() + 
     ggsci::scale_fill_lancet() +
     ggsci::scale_color_lancet() +
     theme(
       axis.text = element_text(size = 14),
       axis.title = element_text(size = 14),
       strip.text = element_text(size = 14) 
     ) + 
     ylab("Score")
   
  )
  
}

pdf(paste0(outdir, "food_patterns_vs_DemoVars.pdf"), width = 8, height = 3)
for(gg in qplots){
  print(gg)
}

metlong <- metlong %>% 
  dplyr::mutate(status_c2 = ifelse(status_c2 == "NaN", NA, status_c2)) %>% 
  dplyr::mutate(status_c2 = factor(status_c2, levels=c("Insufficient gain", "Normal", "Excessive gain")))
  
qplots <- list()
for(var in othervars_cat){
  
  aux <- metlong %>% filter(!is.na(!!sym(var))) %>% 
    filter(!!sym(var) != "NaN")
  
  (qplots[[var]] <- ggplot(aux, 
                           aes(x = !!sym(var), 
                                        y=Pattern_score_value, 
                                        col=Pattern_score)) +
     facet_grid(~ Pattern_score) +
     geom_boxplot(fill="white") +
     geom_jitter( size=0.05, alpha=0.5) +
     theme_bw() + 
     ggsci::scale_fill_lancet() +
     ggsci::scale_color_lancet() +
     theme(
       axis.text = element_text(size = 14),
       axis.text.x = element_text(size = 14, angle=45, 
                                 vjust=1, hjust=1),
       axis.title = element_text(size = 14),
       strip.text = element_text(size = 14) 
     ) + 
     ylab("Score")
   
  )
  
}

pdf(paste0(outdir, "food_patterns_vs_DemoVarsCat.pdf"), width = 8, height = 4)
for(gg in qplots){
  print(gg)
}
dev.off()



qplots <- list()
for(var in othervars_cat){
  
  aux <- metlong %>% filter(!is.na(!!sym(var))) %>% 
    filter(!!sym(var) != "NaN") %>% 
    filter(status_c1 == "normal")
  
  (qplots[[var]] <- ggplot(aux, 
                           aes(x = !!sym(var), 
                               y=Pattern_score_value, 
                               col=Pattern_score)) +
      facet_grid(~ Pattern_score) +
      geom_boxplot(fill="white") +
      geom_jitter( size=0.05, alpha=0.5) +
      theme_bw() + 
      ggsci::scale_fill_lancet() +
      ggsci::scale_color_lancet() +
      theme(
        axis.text = element_text(size = 14),
        axis.text.x = element_text(size = 14, angle=45, 
                                   vjust=1, hjust=1),
        axis.title = element_text(size = 14),
        strip.text = element_text(size = 14) 
      ) + 
      ylab("Score")
    
  )
  
}

pdf(paste0(outdir, "food_patterns_vs_DemoVarsCat_onlyNormalT0.pdf"), width = 8, height = 4)
for(gg in qplots){
  print(gg)
}
dev.off()

aux <- this_metadata%>% 
  filter(status_c2!= "NaN") %>% 
  filter(status_c1 == "normal")
models1 <- makeLinearModelsSingleVariable(aux, "status_c2", 
                                          othervars, 
                                          patnames, 
                                          combos=1,
                                          outdir = outdir, 
                                          name = paste0("FoodNMF_onlyNormalT0") )

alphadif <- testDiversityDifferences(aux, 
                                     patnames, 
                                     othervars, 
                                     outdir, "FoodNMF_onlyNormalT0")

