library(tidyverse)
library(phyloseq)
library(mixOmics)
library(DESeq2)

getVstFromPhyloseq <- function(phobj, rows_as_samples=TRUE, outdir="./", name=""){
  dds <- phyloseq::phyloseq_to_deseq2(phobj, ~ 1)  # no design needed, just for VST
  dds <- DESeq2::estimateSizeFactors(dds)
  dds <- DESeq2::estimateDispersions(dds)
  vst_counts <- assay(varianceStabilizingTransformation(dds, blind = TRUE))

  if(rows_as_samples){
    vstdf_write <- vst_counts %>% t %>% as.data.frame %>%
      rownames_to_column("sampleID")
  }else{
    vstdf_write <- vst_counts %>% as.data.frame %>%
      rownames_to_column("sampleID")
  }
  if(outdir != ""){
    write_tsv(vstdf_write, paste0(outdir, "VST_counts_", name, ".tsv"))
  }
  return(vstdf_write)
}

makeRCC_MixOmics <- function(microbiome_mat, diet_mat, tax2plot=c(), meta2plot = c(),
                             outdir="./", name="all"){

  assertthat::assert_that(all(tax2plot %in% colnames(microbiome_mat)))
  #Some NAs are already encoded as 0s
  diet_mat[is.na(diet_mat)] <- 0
  not_nas <- rowSums(diet_mat != 0) > 0
  not_nas <- names(not_nas)[not_nas > 0]

  microbiome_mat_filt <- microbiome_mat[not_nas, ]
  if(length(tax2plot)>0){
    microbiome_mat_filt <- microbiome_mat[, tax2plot]
  }

  diet_mat_filt <- diet_mat[not_nas, ]
  if(length(meta2plot)>0){
    diet_mat_filt <- diet_mat_filt[, meta2plot]
  }

  microbiome_mat_filt_scale <- scale(microbiome_mat_filt)
  with_nas <- apply(microbiome_mat_filt_scale, MAR=2, \(x)length(which(is.na(x)))) %>% sort
  with_nas <- with_nas[with_nas > 0]
  microbiome_mat_filt_scale <- microbiome_mat_filt_scale[, ! colnames(microbiome_mat_filt_scale) %in% names(with_nas)]

  diet_mat_filt_scale <- scale(diet_mat_filt)

  rcc_res <- mixOmics::rcc(microbiome_mat_filt_scale, diet_mat_filt_scale, method = "shrinkage")
  save(rcc_res, file = paste0(outdir, "mixomics_RCC_", name, ".RData"))

  return(rcc_res)
}

makeReplacementsDF <- function(df, replace_list, totitlecase=c()){
  dfmod <- purrr::reduce(replace_list, ~ .x %>%
                           dplyr::mutate(!!sym(.y[3]) := gsub(.y[1], .y[2], !!sym(.y[3]), perl=TRUE)),
                         .init = df)%>%
    dplyr::mutate(across(all_of(totitlecase), tools::toTitleCase))

  return(dfmod)
}

prepare_RCC_plot <- function(rcc_res,
                             types_ordered,
                             replace_strings,
                             all_var_types,
                             face_bold = c(),
                             plotn=20,
                             plottype= "corr",
                             outdir="./", name="all"){

  df1 <- rcc_res$loadings$X  %>% as.data.frame %>% rownames_to_column("variable") %>%
    dplyr::mutate(class = "Microbiome")
  df2 <- rcc_res$loadings$Y  %>% as.data.frame %>% rownames_to_column("variable")%>%
    dplyr::mutate(class = "Variable")

  # Correlation  with canonical variates
  corY <- cor(rcc_res$Y, rcc_res$variates$X) %>%
    as.data.frame %>% rownames_to_column("variable")
  names(corY) <- c("variable", "corV1", "corV2")
  corX <- cor(rcc_res$X, rcc_res$variates$Y)%>%
    as.data.frame %>% rownames_to_column("variable")
  names(corX) <- c("variable", "corV1", "corV2")

  df_cor <- rbind(corX, corY)

  arrange_by_dist <- ifelse(plottype == "corr", "dist_to_center_cor", "dist_to_center")

  df3 <- rbind(df1, df2) %>%
    merge(df_cor, by="variable") %>%
    dplyr::mutate(label_var = ifelse(class=="Variable", variable, "")) %>%
    dplyr::mutate(dist_to_center = sqrt(as.numeric(scale(V1))^2 + as.numeric(scale(V2))^2)) %>%
    dplyr::mutate(dist_to_center_cor = sqrt(as.numeric(scale(corV1))^2 + as.numeric(scale(corV2))^2)) %>%
    dplyr::group_by(class) %>%
    dplyr::arrange(desc(!!sym(arrange_by_dist))) %>%
    dplyr::mutate(dist_order = 1:n(),
                  label_mic = ifelse(dist_order <= plotn & class != "Variable", variable, "")) %>%
    dplyr::ungroup()


  # List containing replacement in the label_mic column and in the label_var column
  replace_list <- c(
    list(
      c("_", " ", "label_mic"),
      c("[\\[\\]]", "", "label_mic")),
    purrr::map(replace_strings, \(x) c(x, "label_var"))
  )

  df3b <- makeReplacementsDF(df3, replace_list, c("label_var"))

  df3b <- df3b %>%
    dplyr::mutate(type = ifelse(class == "Variable",
                                all_var_types[label_var],
                                "Microbiome")) %>%
    dplyr::mutate(type = factor(type, levels=c( types_ordered))) %>%
    dplyr::mutate(label_full = ifelse(class == "Microbiome" & label_mic== "", "",
                                      ifelse(class=="Variable",
                                             ifelse(variable %in% face_bold, paste0("bold('", label_var, "')"), paste0("plain('", label_var, "')")),
                                             paste0("italic('", label_mic, "')")
                                      )
    )
    )

  write_tsv(df3b, file = paste0(outdir, name, "_RCC_table2plot.tsv"))
  return(df3b)
}

rotate_points <- function(mat, vnames, element_to_right){
  assertthat::assert_that(element_to_right %in% vnames, msg = "rotate_points: element_to_right to calculate angle not in vnames")
  rownames(mat) <- vnames

  rotation <- pi - atan2(mat[element_to_right, 2], mat[element_to_right, 1])
  
  R <- matrix(c(cos(rotation), -sin(rotation),
                sin(rotation),  cos(rotation)),
              nrow = 2)
  rotated <- mat %*% R
  return(rotated)
}


plotRCC <- function(df3, cols=c(), outdir="./", name="plot",
                    alpha_not_shown=0.2,
                    shape_special = c(),
                    w=8, h=8,
                    topn=20,
                    element_to_right = NULL,
                    plottype= "corr"){
  # plottype: plot correlation to loadings or loadings
  if(length(cols) == 0){
    cols <- ggsci::pal_aaas()(length(unique(df3$type))-1)
    cols <- c(cols, "gray20")  ## Microbiome always gray
  }

  df3 <- df3 %>% dplyr::mutate(
    #Size = ifelse(dist_order <= topn, 0.5, 0.1),
    Alpha = ifelse(type != "Microbiome" | dist_order <= topn, 1, alpha_not_shown),
    PointShape = ifelse(variable %in% names(shape_special), shape_special[variable], 16)
  )

  if(plottype == "corr"){
    v1 <- "corV1"
    v2 <- "corV2"
    name <- paste0(name, "_corr2CC")
  }else{
    v1 <- "V1"
    v2 <- "V2"
  }
  if(! is.null(element_to_right )){
    rotated <- rotate_points(as.matrix(df3[, c(v1, v2)]), df3$variable, element_to_right)
    newnames <- paste0(c(v1, v2), "_rot")
    df3 <- df3 %>% dplyr::mutate(
      !!sym(newnames[1]) := rotated[, 1],
      !!sym(newnames[2]) := rotated[, 2]
    )
    v1 = newnames[1]
    v2 = newnames[2]
  }
  
  g1 <- ggplot(df3, aes(x=!!sym(v1), y=!!sym(v2), col=type)) +
    geom_vline(xintercept = 0, linetype=2, col="lightgray") +
    geom_hline(yintercept = 0, linetype=2, col="lightgray") +
    geom_point(alpha=df3$Alpha, shape=df3$PointShape) + #size=Size, aes(alpha=Alpha)
    ggrepel::geom_text_repel(data=df3 %>% filter(label_full != "") , 
                             aes(label = label_full), parse=TRUE, 
                             max.overlaps = Inf,
                             max.time = 2,
                             max.iter = 5000) +
    theme_classic(base_size = 14) +
    scale_color_manual(values=cols) +
    xlab("Canonical Variate 1")  + 
    ylab("Canonical Variate 2")  + 
    scale_fill_manual(values=cols)

  ggsave(filename = paste0(outdir_cca, "RCC_", name, ".pdf"), g1, width = w, height = h)
  return(g1)
}

RCC_mixOmixs_fullPipeline <- function(phobj, vars2select,
                                      types_ordered=c(), replace_strings=c(),
                                      all_var_types = c(),
                                      tax2plot=c(), vstdf = NULL, cols=c(),
                                      face_bold = c(),
                                      shape_special = c(),
                                      alpha_not_shown=0.2, plotn=20,
                                      element_to_right = NULL,
                                      outdir="./", plottype= "corr",
                                      name="all", w=8, h=8){
  if(is.null(vstdf)){
    vstdf <- getVstFromPhyloseq(phobj, outdir = outdir, name=name)
  }
  microbiome_mat <- vstdf %>% column_to_rownames("sampleID") %>% as.matrix

  diet_mat <- sample_data(phobj) %>% data.frame %>%
    dplyr::select(all_of(vars2select)) %>% as.matrix()

  rcc_res <- makeRCC_MixOmics(microbiome_mat, diet_mat, tax2plot, c(), outdir=outdir, name=name)

  df3 <- prepare_RCC_plot(rcc_res,
                          types_ordered=types_ordered,
                          replace_strings=replace_strings,
                          all_var_types=all_var_types,
                          plotn=plotn, 
                          face_bold=face_bold,
                          plottype= plottype, outdir=outdir_cca, name=name)

  g1 <- plotRCC(df3, cols=cols, outdir=outdir_cca, name=name,
                plottype = plottype, topn=plotn,
                alpha_not_shown=alpha_not_shown, 
                shape_special = shape_special, 
                element_to_right= element_to_right,
                w=w, h=h)

  return(list(
    microbiome_mat = microbiome_mat,
    diet_mat = diet_mat,
    rcc=rcc_res,
    df=df3,
    plot=g1,
    name=name,
    plottype=plottype
  ))
}


#################################
OUT_BASE <- "/home/carlos/Documentos/CORALS/results_rstudio/results_Agosto25_1/"
outdir <-paste0(OUT_BASE, "DAA_integrate_food_etc_250924/")
if(!dir.exists(outdir)) dir.create(outdir)
outdir_cca <- paste0(outdir, "CCA_4/")
if(!dir.exists(outdir_cca)) dir.create(outdir_cca)

 all_daa <- read_tsv(paste0(OUT_BASE, "DAA_integrate_food_etc_250818/all_DAA_long.tsv"))
 all_daa_mod <- read_tsv(paste0(outdir, "all_DAA_long_NamesModified.tsv"))

load("/home/carlos/Documentos/CORALS/results_rstudio/results_Agosto25_1//foodPCA/phyloseq_list_foodPCA_withNMF_withInc.RData")
phobj <- all_phyloseq$remove_tanda2
s_meta <- sample_data(phobj) %>% data.frame()
if(! file.exists( paste0(outdir, "VST_counts.tsv"))){
  vstdf_write <- getVstFromPhyloseq(phobj, outdir = outdir)
}else{
  vstdf_write <-read_tsv(paste0(outdir, "VST_counts.tsv"))
}
################################################################

# Original variable names
patnames <- c("Preprocessed",  "Mediterranean", "Western")
popvarscca <- c("z_t0", "z_waist_00", "inc_z_bmi", "inc_z_waist")
popvarscca_t1 <- c("z_t1", "z_waist_01", "inc_z_bmi", "inc_z_waist")
quant_vars_diet <- names(s_meta)[grepl("_clr$", names(s_meta))]
diet_aggr_vars <- quant_vars_diet[1:4]
food_vars <- quant_vars_diet[5:length(quant_vars_diet)]

vars2select <- c(food_vars, popvarscca)
vars2select_withPats <-  c(patnames, food_vars, popvarscca)

vars2select_t1 <- c(food_vars, popvarscca_t1)
vars2select_withPats_t1 <-  c(patnames, food_vars, popvarscca_t1)

# Modified variable names
types_ordered <- c("Demographic","Mediterranean", "Preprocessed", "Western", "Microbiome")


food_western <- c("Dairy","Refined Cereals", "Sweets, Pastries", "Meat", "Tubers", "Sugars and Sweets", "Sauces, Condiments")
food_prep <- c("Juices, Softdrinks", "Dairy Derivatives" , "Prepared Foods", "Snacks Savory")
food_med <- c("Fruits", "Vegetables", "Fish", "Fats, Oils", "Eggs" , "Legumes", "Whole Grain Cereals", "Nuts", "Oleaginous Fruits")
food_macro <- c("Energy (Kcal)", "Carbohydrates (g)", "Fiber (g)", "Protein (g)", "Total Fat (g)")
patnames <- c("Preprocessed", "Mediterranean", "Western")

popvarscca_newnames <- c("Z-Score BMI T0", "Z-Score Waist T0", "Change in Z-Score BMI", "Change in Z-Score Waist")
popvarscca_newnames_t1 <- c("Z-Score BMI T1", "Z-Score Waist T1", "Change in Z-Score BMI", "Change in Z-Score Waist")

# Vector with the 'type' of each variable, used to color in plots
all_var_types <- c(rep("Mediterranean", length(food_med)),
                   rep("Preprocessed", length(food_prep)),
                   rep("Western", length(food_western)),
                   rep("Demographic", length(popvarscca_newnames))
)
names(all_var_types) <- c(food_med, food_prep, food_western, popvarscca_newnames)

all_var_types_t1 <- c(rep("Mediterranean", length(food_med)),
                   rep("Preprocessed", length(food_prep)),
                   rep("Western", length(food_western)),
                   rep("Demographic", length(popvarscca_newnames_t1))
)
names(all_var_types_t1) <- c(food_med, food_prep, food_western, popvarscca_newnames_t1)


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

topn <- 20

##########################################################

#tax2plot <- all_daa %>% filter(Contrast %in% colnames(diet_mat)) %>%
#  filter(padj <= 0.001) %>% pull(taxon) %>% unique
#length(tax2plot)
#tax2plot <- c()

#
#res1 <- RCC_mixOmixs_fullPipeline(phobj =  phobj,
#                          vars2select = vars2select,
#                          types_ordered = types_ordered,
#                          replace_strings = replace_strings,
#                          all_var_types = all_var_types,
#                          tax2plot = c(),
#                          vstdf = vstdf_write,
#                          cols = c(),
#                          outdir = outdir_cca,
#                          plotn = 15,
#                          plottype = "corr",
#                          name = "allMicroFoodVarsPlusBMI_Corr2")
#
#res1$plot
#
#
#all_var_types_pats <- c(Mediterranean="Mediterranean", Western="Western", Preprocessed="Preprocessed",
#                        all_var_types[21:24])
#res1 <- RCC_mixOmixs_fullPipeline(phobj =  phobj,
#                                  vars2select = c(patnames, popvarscca),
#                                  types_ordered = types_ordered,
#                                  replace_strings = replace_strings,
#                                  all_var_types = all_var_types_pats,
#                                  tax2plot = c(),
#                                  vstdf = vstdf_write,
#                                  element_to_right = "Mediterranean",
#                                  cols = c(),
#                                  outdir = outdir_cca,
#                                  plotn = 15,
#                                  plottype = "corr",
#                                  name = "allMicroFoodPatternsPlusBMI_Corr2")
#
#res1$plot


## Food + patterns + BMI at T0
shapes_special <- c(rep(8, length(patnames)), rep(23, length(popvarscca)))
names(shapes_special) <- c(patnames, popvarscca)

#cols <- ggsci::pal_lancet()(7)
cols1 <- c("#00468BFF", "#ED0000FF", "#42B540FF", "#0099B4FF", "#925E9FFF", "#FDAF91FF", "#AD002AFF")
#cols2 <- ggsci::pal_aaas()(3)
cols2 <- c("#3B4992FF", "#EE0000FF", "#008B45FF")

cols2use <- c("Demographic"=cols1[4], 
              "Mediterranean" = cols2[2],
              "Preprocessed" = cols2[1],
              "Western" = cols2[3],
              "Microbiome"="gray20")
names(cols2use) <- NULL
all_var_types_pats <- c(Mediterranean="Mediterranean", Western="Western", Preprocessed="Preprocessed",
                        all_var_types[21:24])
res1 <- RCC_mixOmixs_fullPipeline(phobj =  phobj,
                                  vars2select = vars2select_withPats,
                                  types_ordered = types_ordered,
                                  replace_strings = replace_strings,
                                  all_var_types = c(all_var_types, all_var_types_pats[1:3]),
                                  tax2plot = c(),
                                  vstdf = vstdf_write,
                                  face_bold = c(patnames, popvarscca), 
                                  shape_special = shapes_special,
                                  cols = cols2use,
                                  plotn = 20,
                                  element_to_right = "Mediterranean",
                                  outdir = outdir_cca,
                                  plottype = "corr", 
                                  name = "allMicroFoodAndPatsAndBMI_Corr5_top20", 
                                  h=8, w=10)

res1$plot 

tax2use_cca <- all_daa %>% filter(Contrast %in% c("patnames", "inc_z_bmi", "z_t0", "z_waist_00")) %>% 
  dplyr::filter(padj <= 0.01) %>% pull(taxon) %>% unique
vstdf_write_filt <- vstdf_write[, c("sampleID", tax2use_cca)] 

res1_filt <- RCC_mixOmixs_fullPipeline(phobj =  phobj,
                                  vars2select = vars2select_withPats,
                                  types_ordered = types_ordered,
                                  replace_strings = replace_strings,
                                  all_var_types = c(all_var_types, all_var_types_pats[1:3]),
                                  tax2plot = c(),
                                  vstdf = vstdf_write_filt,
                                  face_bold = c(patnames, popvarscca), 
                                  shape_special = shapes_special,
                                  cols = cols2use,
                                  plotn = 20,
                                  element_to_right = "Mediterranean",
                                  outdir = outdir_cca,
                                  alpha_not_shown = 0.3,
                                  plottype = "corr", 
                                  name = "SigMicroFoodAndPatsAndBMI_Corr5_top20", 
                                  h=8, w=10)

#### only T1

shapes_special_t1 <- c(rep(8, length(patnames)), rep(23, length(popvarscca_t1)))
names(shapes_special_t1) <- c(patnames, popvarscca_t1)

s2use <- s_meta %>% filter(status_c1=="normal") %>% pull(sampleID)
phobj_filt <- phyloseq::prune_samples(s2use,phobj)

res2 <- RCC_mixOmixs_fullPipeline(phobj =  phobj_filt,
                                  vars2select = vars2select_withPats_t1,
                                  types_ordered = types_ordered,
                                  replace_strings = replace_strings,
                                  all_var_types = c(all_var_types_t1, all_var_types_pats[1:3]),
                                  tax2plot = c(),
                                  vstdf = vstdf_write,
                                  face_bold = c(patnames, popvarscca_t1), 
                                  shape_special = shapes_special_t1,
                                  element_to_right = "Mediterranean",
                                  cols = cols2use,
                                  plotn = 20,
                                  outdir = outdir_cca,
                                  plottype = "corr", 
                                  name = "allMicroFoodAndPatsAndBMI_T1onlyNormalT0_Corr5_top20", 
                                  h=8, w=10)

## filtering taxa
vstdf_write_filt_t1 <- vstdf_write_filt %>% filter(sampleID %in% s2use)

res2_filt <- RCC_mixOmixs_fullPipeline(phobj =  phobj_filt,
                                  vars2select = vars2select_withPats_t1,
                                  types_ordered = types_ordered,
                                  replace_strings = replace_strings,
                                  all_var_types = c(all_var_types_t1, all_var_types_pats[1:3]),
                                  tax2plot = c(),
                                  vstdf = vstdf_write_filt_t1,
                                  face_bold = c(patnames, popvarscca_t1), 
                                  shape_special = shapes_special_t1,
                                  element_to_right = "Mediterranean",
                                  cols = cols2use,
                                  plotn = 20,
                                  outdir = outdir_cca,
                                  plottype = "corr", 
                                  name = "SigMicroFoodAndPatsAndBMI_T1onlyNormalT0_Corr5_top20", 
                                  h=8, w=10)


### Finally, combine CC with DAA

library(pheatmap)

#### heatmaps were saved to: heatmaps_most_distant_to_center_CCA3

ccaresult_t0 <- read_tsv(paste0(outdir_cca, "allMicroFoodAndPatsAndBMI_Corr5_top20_RCC_table2plot.tsv"))
ccaresult_t1normt0 <- read_tsv(paste0(outdir_cca, "allMicroFoodAndPatsAndBMI_T1onlyNormalT0_Corr5_top20_RCC_table2plot.tsv"))
all_daa <- read_tsv(paste0(outdir, "all_DAA_long.tsv"))
all_daa_mod <- read_tsv(paste0(outdir, "all_DAA_long_NamesModified.tsv"))


## Now plot


PLIM <- 0.01
LFCLIM <- 1
PLIM_PLOT <- 0.05
MIN_COMP_LIM <- "byBMI" # 10 # usually a number to filter the number of comparisons in which each taxa is signif
TOPN <- 20
#tax2plot <- all_daa %>%
#  filter(!is.na(padj) & padj <= PLIM) %>%
#  filter(!is.na(log2FoldChangeShrink) & abs(log2FoldChangeShrink) >= LFCLIM) %>%
#  pull(taxon) %>% unique
#length(tax2plot)


tax2plot <-  ccaresult_t0 %>%
  filter(type=="Microbiome") %>% 
  filter(dist_order <= TOPN) %>% 
  dplyr::mutate(taxon=variable)

nrow(tax2plot)

all_daa_mod1 <- all_daa_mod %>%
  dplyr::mutate(LFC = ifelse(padj <= PLIM_PLOT, log2FoldChangeShrink, log2FoldChangeShrink)) %>%
  dplyr::filter(taxon %in% tax2plot$taxon) %>%
  dplyr::select(taxon, Contrast, LFC)

mat1 <- all_daa_mod1 %>%
  #dplyr::mutate(Contrast = ifelse(Samples == "Normal_T0", paste(Contrast, Samples, sep="__") , Contrast)) %>%
  #select(-Samples) %>%
  pivot_wider(names_from = Contrast, values_from = LFC) %>%
  column_to_rownames("taxon") %>%
  as.matrix

mat1_sc <- apply(mat1, MAR=2, scale)
rownames(mat1_sc) <- rownames(mat1)

pheatmap(mat1_sc %>% t,
         cluster_rows = TRUE, cluster_cols = TRUE,
         filename = paste0(outdir, "heatmap1_CCA_maxDistTop", as.character(TOPN), "comps.pdf"),
         height = 6, width = 12,
         show_colnames = FALSE)

mat2 <- mat1_sc
qs <- quantile(mat1, c(0.025, 0.975), na.rm = TRUE)
mat2[mat2 < qs[1]] <- qs[1]
mat2[mat2 > qs[2]] <- qs[2]

mat2[is.na(mat2)] <- 0

hm <- pheatmap(mat2) #just to cluster

pheatmap(mat2 %>% t,
         cluster_rows = TRUE, cluster_cols = TRUE,
         filename = paste0(outdir, "heatmap1_CCA_maxDistTop", as.character(TOPN), "comps2.pdf"),
         height = 6, width = 12,
         show_colnames = FALSE)


## hm with NS in gray
mat_p <- all_daa_mod %>%
  dplyr::filter(taxon %in% tax2plot$taxon) %>%
  dplyr::select(taxon, Contrast, padj) %>%
  #dplyr::mutate(Contrast = ifelse(Samples == "Normal_T0", paste(Contrast, Samples, sep="__") , Contrast)) %>%
  #select(-Samples) %>%
  pivot_wider(names_from = Contrast, values_from = padj) %>%
  column_to_rownames("taxon") %>%
  as.matrix

mat_p <- mat_p[rownames(mat2), colnames(mat2)]
mat3 <- mat2
mat3[mat_p > PLIM_PLOT] <- NA

mat3 <- mat3[hm$tree_row$order, hm$tree_col$order]
pheatmap(t(mat3),
         cluster_rows = FALSE, cluster_cols = FALSE,
         filename = paste0(outdir, "heatmap1_CCA_maxDistTop", as.character(TOPN), "comps3.pdf"),
         height = 6, width = 12,
         show_colnames = FALSE)


#####

new_food_names <- all_daa_mod %>% filter(type == "Food groups") %>%
  pull(Contrast) %>% unique

matfood <- mat2[, new_food_names]
pheatmap(matfood %>% t,
         cluster_rows = TRUE, cluster_cols = TRUE,
         filename = paste0(outdir, "heatmap1_CCA_maxDistTop", as.character(TOPN), "_foodOnly_comps4.pdf"),
         #clustering_method = "ward.D2",
         height = 5, width = 12,
         show_colnames = FALSE)


allvars_names <- all_daa_mod$Contrast %>% unique
food_western <- c("Dairy","Refined Cereals", "Sweets, Pastries", "Meat", "Tubers", "Sugars and Sweets", "Sauces, Condiments")
food_prep <- c("Juices, Softdrinks", "Dairy Derivatives" , "Prepared Foods", "Snacks Savory")
food_med <- c("Fruits", "Vegetables", "Fish", "Fats, Oils", "Eggs" , "Legumes", "Whole Grain Cereals", "Nuts", "Oleaginous Fruits")
food_macro <- c("Energy (Kcal)", "Carbohydrates (g)", "Fiber (g)", "Protein (g)", "Total Fat (g)")
patnames <- c("Preprocessed", "Mediterranean", "Western")
popvars <- allvars_names[!allvars_names %in% c(new_food_names, patnames, food_macro)]
popvars_all <- popvars[!grepl("\\(Normal T0\\)", popvars)] %>% sort
popvars_all <- c(popvars_all[grepl("Age", popvars_all)],
                 popvars_all[grepl("Sex", popvars_all)],
                 popvars_all[!grepl("Sex|Age|Change", popvars_all, perl=TRUE)],
                 popvars_all[grepl("Change", popvars_all, perl=TRUE)]
)

popvars_all_t0 <- popvars_all[!grepl("T1|Change", popvars_all, perl=TRUE)]
popvars_all_t1 <- popvars_all[grep("T1|Change", popvars_all, perl=TRUE)]

popvars_nt0 <- popvars[grepl("\\(Normal T0\\)", popvars)] %>% sort
popvars_nt0 <- c(
  popvars_nt0[!grepl("Status|Change", popvars_nt0, perl=TRUE)],
  popvars_nt0[grepl("Change", popvars_nt0, perl=TRUE)],
  popvars_nt0[grepl("Status", popvars_nt0, perl=TRUE) & grepl("vs Normal", popvars_nt0, perl=TRUE)],
  popvars_nt0[grepl("Status", popvars_nt0, perl=TRUE) & !grepl("vs Normal", popvars_nt0, perl=TRUE)]
)

varlist <- list(popvars_all_t0,
                popvars_all_t1, popvars_nt0,
                patnames, food_macro,
                food_western, food_prep, food_med)
names(varlist) <- c("Demographic", "Demographic T1", "T1 (normal BMI at T0)",
                    "Diet pattern", "Macronutrients",
                    "Food groups - Western", "Food groups - Preprocessed", "Food groups - Mediterranean")
varorder <- unlist(varlist)

## Now order bacteria
bacorder_ind <- mat2
assertthat::assert_that(all(colnames(bacorder_ind) == colnames(mat_p)))
assertthat::assert_that(all(rownames(bacorder_ind) == rownames(mat_p)))
bacorder_ind[mat_p >= PLIM_PLOT] <- 0
bacorder_ind[bacorder_ind > 0] <- 1
bacorder_ind[bacorder_ind < 0] <- -1

bacorder_ind[, "Mediterranean"] <- -1*bacorder_ind[, "Mediterranean"]
bacorder_ind[, food_med] <- -1*bacorder_ind[, food_med]

vars2score <- c(patnames,
                popvars_all_t0[!grepl("Sex|Age", popvars_all_t0, perl=TRUE)],
                popvars_all_t1,
                popvars_nt0[!grepl("Insufficient", popvars_nt0)],
                food_western, food_prep, food_med
)
bacorder_ind <- bacorder_ind[, vars2score] %>% rowSums() %>% sort

mat_ord <- mat2[names(bacorder_ind), varorder]

labels_col <- gsub("_", " ", rownames(mat_ord))
labels_col <- lapply(labels_col, \(x){
  bquote(italic(.(x)))
}) %>% as.expression()

anncol <- map2(varlist, names(varlist), \(x, xn) data.frame(var=x, type=xn)) %>%
  bind_rows %>% column_to_rownames("var")

colors <- ggsci::pal_lancet()(7)
colors <- rev(colors[-7])

library(colorspace)
library(scico)
shades <- colorspace::lighten(colors[6], amount = seq(0.4, 0, length.out = 3))
#"#6C8BCCFF" "#4367A9FF" "#00468BFF"

color_list_cols <-  list(type=c(colors[1:5], shades))
names(color_list_cols$type) <- names(varlist)

mat_p2 <- mat_p[rownames(mat_ord), colnames(mat_ord)]
mat_pchar <- matrix(character(prod(dim(mat_p2))), ncol=ncol(mat_p2))
mat_pchar[mat_p2 <= PLIM_PLOT] <- "*"
colnames(mat_pchar) <- colnames(mat_p2)
rownames(mat_pchar) <- rownames(mat_p2)

# scico palettes
# “acton”, “bam”, “bamako”, “bamO”, “batlow”, “batlowK”, “batlowW”,
#“berlin”, “bilbao”, “broc”, “brocO”, “buda”, “bukavu”, “cork”, “corkO”, “davos”, “devon”, “fes”, “glasgow”, “grayC”, “hawaii”, “imola”, “lajolla”, “lapaz”, “lipari”, “lisbon”, “managua”, “navia”, “nuuk”, “oleron”, “oslo”, “roma”, “romaO”, “tofino”, “tokyo”, “turku”, “vanimo”, “vik”, “vikO”

#outdir <- paste0(outdir, "/all_colors2/")
#if(!dir.exists(outdir)) dir.create(outdir)

#allcols <- c("acton", "bam", "bamako", "bamO", "batlow", "batlowK", "batlowW", "berlin",
#"bilbao", "broc", "brocO", "buda", "bukavu", "cork", "corkO", "davos", "devon", "fes", "glasgow",
#"grayC", "hawaii", "imola", "lajolla", "lapaz", "lipari", "lisbon", "managua", "navia", "nuuk", "oleron", "oslo", "roma",
#"romaO", "tofino", "tokyo", "turku", "vanimo", "vik", "vikO")

#for(pname in allcols){
pname <-"vik"
col_fun <- scico(100, palette = pname)

#cluster_rows = TRUE --> IGNORE ORDERING
pheatmap(mat_ord,
         cluster_rows = TRUE, cluster_cols = FALSE,
         filename = paste0(outdir, "heatmap1_CCA_maxDistTop", as.character(TOPN), "_comps5.pdf"),
         #clustering_method = "ward.D2",
         height = 6, width = 18, # height = 14 MIN_COMP_LIM = 4 (80 y algo taxa)
         gaps_col = sapply(varlist, length) %>% cumsum(),
         annotation_col = anncol,
         show_colnames = TRUE,
         angle_col = 45,           # rotate column names
         fontsize_col = 10,
         fontsize_row = 10,
         labels_row = labels_col,
         #display_numbers = mat_pchar,
         number_color = "white",
         color = col_fun,
         annotation_colors = color_list_cols)

pheatmap(mat_ord,
         cluster_rows = TRUE, cluster_cols = FALSE,
         filename = paste0(outdir, "heatmap1_CCA_maxDistTop", as.character(TOPN), "_comps5b.pdf"),
         #clustering_method = "ward.D2",
         height = 6, width = 18, # height = 14 MIN_COMP_LIM = 4 (80 y algo taxa)
         gaps_col = sapply(varlist, length) %>% cumsum(),
         annotation_col = anncol,
         show_colnames = TRUE,
         angle_col = 45,           # rotate column names
         fontsize_col = 10,
         fontsize_row = 10,
         labels_row = labels_col,
         display_numbers = mat_pchar,
         number_color = "white",
         fontsize_number=12,
         color = col_fun,
         annotation_colors = color_list_cols)

mat_ord2 <- mat_ord
mat_ord2[mat_p2 >= 0.1] <- 0
pheatmap(mat_ord2,
         cluster_rows = TRUE, cluster_cols = FALSE,
         filename = paste0(outdir, "heatmap1_CCA_maxDistTop", as.character(TOPN), "_comps6.pdf"),
         #clustering_method = "ward.D2",
         height = 6, width = 16,  # height = 14 MIN_COMP_LIM = 4 (80 y algo taxa)
         gaps_col = sapply(varlist, length) %>% cumsum(),
         annotation_col = anncol,
         show_colnames = TRUE,
         angle_col = 45,           # rotate column names
         fontsize_col = 10,
         fontsize_row = 10,
         labels_row = labels_col,
         #display_numbers = mat_pchar,
         number_color = "white",
         color = col_fun,
         annotation_colors = color_list_cols)
#}

########### Until here in August 2025.
########### September 24, 2025: make CCA only with bacteria that are DAA


