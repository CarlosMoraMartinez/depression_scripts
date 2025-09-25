library(tidyverse)
library(phyloseq)


OUT_BASE <- "/home/carlos/Documentos/CORALS/results_rstudio/results_Agosto25_1/"
outdir <-paste0(OUT_BASE, "DAA_integrate_food_etc_250818/")
if(!dir.exists(outdir)) dir.create(outdir)

# phyloseq objects
load("/home/carlos/Documentos/CORALS/results_rstudio/results_Agosto25_1//foodPCA/phyloseq_list_foodPCA_withNMF_withInc.RData")
s_meta <- all_phyloseq$remove_tanda2 %>% sample_data %>% data.frame
rm(all_phyloseq)
write_tsv(s_meta, paste0(outdir, "metadata_used.tsv"))
# DAA adjusting for Sex, Age and Hospital
load(paste0(OUT_BASE, "DESeq_Patterns_withCovars_250818/DESEQ2_FoodPatterns_covars_onlySexHosAge.RData"))

all_daa <- rbind(
  daa_all_cov$remove_tanda2$SexHosAge$all_contrasts$Sex_Girl_vs_Boy$resdf %>% dplyr::mutate(Contrast = "Sex_Girl_vs_Boy"),
  daa_all_cov$remove_tanda2$SexHosAge$all_contrasts$age_T0$resdf %>% dplyr::mutate(Contrast = "Age")
)
names2bind <- names(daa_all_cov$remove_tanda2)
names2bind <- names2bind[names2bind != "SexHosAge"]

for(nn in names2bind){
  cat(nn, ": ", daa_all_cov$remove_tanda2[[nn]]$all_contrasts[[1]]$contrast_name, "\n")
  aux <- daa_all_cov$remove_tanda2[[nn]]$all_contrasts[[1]]$resdf %>%
    dplyr::mutate(Contrast = daa_all_cov$remove_tanda2[[nn]]$all_contrasts[[1]]$contrast_name)
  all_daa <- rbind(all_daa, aux)
}

write_tsv(all_daa, file = paste0(outdir, "all_DAA_long.tsv"))
rm(daa_all_cov)
# DAA only in normal at T0; adjusting for Sex, Age and Hospital
load(paste0(OUT_BASE, "DESeq_Patterns_withCovars_onlyNormalT0_250818/DESEQ2_FoodPatterns_covars_onlyNormalT0.RData"))

all_daa$Samples <- "All"

names2bind <- names(daa_all_cov$remove_tanda2)
names2bind <- names2bind[names2bind != "remove_tanda2_status_c2_SexHosAge"]
for(nn in names2bind){
  cat(nn, ": ", daa_all_cov$remove_tanda2[[nn]]$all_contrasts[[1]]$contrast_name, "\n")
  aux <- daa_all_cov$remove_tanda2[[nn]]$all_contrasts[[1]]$resdf %>%
    dplyr::mutate(Contrast = daa_all_cov$remove_tanda2[[nn]]$all_contrasts[[1]]$contrast_name,
                  Samples = "Normal_T0")
  all_daa <- rbind(all_daa, aux)
}

status_c2 <- rbind(
  daa_all_cov$remove_tanda2$remove_tanda2_status_c2_SexHosAge$all_contrasts$status_c2_Insufficient.gain_vs_Excessive.gain$resdf %>%
    dplyr::mutate(log2FoldChange = -1*log2FoldChange,
                  log2FoldChangeShrink = -1*log2FoldChangeShrink,
                  Contrast = "status_c2_Excessive.gain_vs_Insufficient.gain"
                  ),
  daa_all_cov$remove_tanda2$remove_tanda2_status_c2_SexHosAge$all_contrasts$status_c2_Normal_vs_Excessive.gain$resdf %>%
    dplyr::mutate(log2FoldChange = -1*log2FoldChange,
                  log2FoldChangeShrink = -1*log2FoldChangeShrink,
                  Contrast = "status_c2_Excessive.gain_vs_Normal"
    ),
  daa_all_cov$remove_tanda2$remove_tanda2_status_c2_SexHosAge$all_contrasts$status_c2_Normal_vs_Insufficient.gain$resdf %>%
    dplyr::mutate(Contrast = "status_c2_Normal_vs_Insufficient.gain")
) %>% dplyr::mutate(Samples = "Normal_T0")

all_daa <- rbind(all_daa, status_c2)

write_tsv(all_daa, file = paste0(outdir, "all_DAA_long.tsv"))
all_daa <- read_tsv(paste0(outdir, "all_DAA_long.tsv"))

rm(daa_all_cov)
rm(aux, status_c2)

######################################3
## modify Contrast names

type_levels <- rev(c("Demographic",
                     "Demographic T1", "T1 (normal BMI at T0)",
                     "Diet pattern",
                     "Macronutrients",
                     "Food groups"))

patnames <- c("Preprocessed",  "Mediterranean", "Western")
quant_vars_diet <- names(s_meta)[5:30]
diet_aggr_vars <- quant_vars_diet[1:5]


vars2exclude <- c("z_bmi_01", "z_bmi_00")
vars2exclude_normal_t0 <- c("z_bmi_00", "z_waist_00", "z_t0", "mg_p_00", patnames)


all_daa_mod <- all_daa %>%
  filter(! Contrast %in% vars2exclude) %>%
  filter( ! (Contrast %in% vars2exclude_normal_t0 & Samples == "Normal_T0")) %>%
  dplyr::mutate(Contrast = gsub("_clr$", "", Contrast),
              type = ifelse(Contrast %in% quant_vars_diet,
                            ifelse(Contrast %in% diet_aggr_vars, "Macronutrients", "Food groups"),
                            "Demographic"),
              type = ifelse( grepl("[Pp]attern",Contrast, perl=T) | (Contrast %in% patnames), "Diet pattern", type ) #grepl("[Pp]attern",Contrast, perl=T) |
) %>%
  dplyr::mutate(
    Contrast = gsub("_g$", " (g)", Contrast, perl=T),
    Contrast = gsub("_kcal", " (Kcal)", Contrast),

    Contrast = ifelse(Samples == "Normal_T0", paste(Contrast, " (normal T0)", sep="") , Contrast),

    Contrast = gsub("z_t1", "Z-score BMI T1", Contrast, perl=T),
    Contrast = gsub("z_waist_01", "Z-score waist T1", Contrast),

    Contrast = gsub("z_t0", "Z-score BMI T0", Contrast, perl=T),
    Contrast = gsub("z_t1", "Z-score BMI T1", Contrast, perl=T),

    Contrast = gsub("inc_z_waist", "Change in Z-score waist", Contrast, perl=T),
    Contrast = gsub("inc_z_bmi", "Change in Z-score BMI", Contrast, perl=T),

    Contrast = gsub("_c1", " T0", Contrast),
    Contrast = gsub("_00", " T0", Contrast),
    Contrast = gsub("_01", " T1", Contrast),
    Contrast = gsub("mg_p", "Body fat %", Contrast),
    Contrast = gsub("nreads", "seq. depth", Contrast),
    Contrast = gsub("sauces_con", "sauces, con", Contrast),
    Contrast = gsub("juices_so", "juices, so", Contrast),
    Contrast = gsub("fats_oils", "fats, oils", Contrast),
    Contrast = gsub("sweets_pas", "sweets, pas", Contrast),
    Contrast = gsub("^z_", "Z-score ", Contrast, perl=T),
    Contrast = gsub("bmi", "BMI", Contrast, perl=T),
    Contrast = gsub("status_c2", "status T1", Contrast, perl=T),
    Contrast = gsub("\\.", " ", Contrast, perl=T),
    Contrast = gsub("_", " ", Contrast)
  ) %>%
  dplyr::mutate(plog = -log10(pvalue),
                padj_log = -log10(padj)) %>%
  dplyr::mutate(type = ifelse( grepl("T1$|Change", Contrast, perl=T), "Demographic T1", type),
                type = ifelse(grepl("\\(normal T0\\)", Contrast, perl=T), "T1 (normal BMI at T0)", type)) %>%
  dplyr::mutate(type = factor(type, levels=type_levels)) %>%
  dplyr::mutate(Contrast = tools::toTitleCase(Contrast))

write_tsv(all_daa_mod, file = paste0(outdir, "all_DAA_long_NamesModified.tsv"))
all_daa_mod <- read_tsv(paste0(outdir, "all_DAA_long_NamesModified.tsv"))
unique(all_daa_mod$Contrast)


## Now plot
library(pheatmap)

PLIM <- 0.01
LFCLIM <- 1
PLIM_PLOT <- 0.05
MIN_COMP_LIM <- "byBMI" # 10 # usually a number to filter the number of comparisons in which each taxa is signif
#tax2plot <- all_daa %>%
#  filter(!is.na(padj) & padj <= PLIM) %>%
#  filter(!is.na(log2FoldChangeShrink) & abs(log2FoldChangeShrink) >= LFCLIM) %>%
#  pull(taxon) %>% unique
#length(tax2plot)

tax_numcomp <- all_daa_mod %>%
  dplyr::mutate(sig_p05 = !is.na(padj) & padj <= 0.05,
                sig_p01 = !is.na(padj) & padj <= 0.01,
                sig_p001 = !is.na(padj) & padj <= 0.001,
                sig_p05_and_lfc = sig_p05 & !is.na(log2FoldChangeShrink) & abs(log2FoldChangeShrink) >= LFCLIM,
                sig_p01_and_lfc = sig_p01 & !is.na(log2FoldChangeShrink) & abs(log2FoldChangeShrink) >= LFCLIM,
                sig_p001_and_lfc = sig_p001 & !is.na(log2FoldChangeShrink) & abs(log2FoldChangeShrink) >= LFCLIM) %>%
  dplyr::mutate(Contrast = ifelse(Samples == "Normal_T0", paste(Contrast, Samples, sep="__") , Contrast)) %>%
  group_by(taxon) %>%
  dplyr::summarise_if(is.logical, sum) %>%
  arrange(desc(sig_p05))

write_tsv(tax_numcomp, file = paste0(outdir, "number_of_comparisons_significant_per_taxa_p01_LFC1.tsv"))


tax2plot <- tax_numcomp %>% filter(sig_p001 >= MIN_COMP_LIM)
nrow(tax2plot)

reserva_outdir <- outdir

#Instead, select those that are significant with BMI

tax2plot <-  all_daa_mod %>% filter(Contrast == "Z-Score BMI T0" | grepl("Change", Contrast) | Contrast == "Z-Score BMI T1") %>%
  filter(padj <= 0.001 & !is.na(padj)) %>%
  group_by(taxon) %>%
  dplyr::summarise(n_sig = n())
name2hms <- "SigBMIT0OrT1All"

#Instead, select those that are significant with BMI only at T0

tax2plot <-  all_daa_mod %>% filter(Contrast == "Z-Score BMI T0") %>%
  filter(padj <= 0.01 & !is.na(padj)) %>%
  group_by(taxon) %>%
  dplyr::summarise(n_sig = n())
name2hms <- "SigOnlyBMIT0"


#Instead, select those that are significant with BMI at T1
tax2plot <-  all_daa_mod %>% filter(Contrast == "Z-Score BMI T1 (Normal T0)" | Contrast == "Change in Z-Score BMI (Normal T0)") %>%
  filter(padj <= 0.01 & !is.na(padj)) %>%
  group_by(taxon) %>%
  dplyr::summarise(n_sig = n())
name2hms <- "SigT1NormalT0"

#Instead, select those that are significant with diet patterns
patnames <- c("Preprocessed", "Mediterranean", "Western" )
tax2plot <-  all_daa_mod %>% filter(Contrast %in% patnames) %>%
  filter(padj <= 0.01 & !is.na(padj)) %>%
  group_by(taxon) %>%
  dplyr::summarise(n_sig = n())
name2hms <- "SigPatternsOnly"

#Instead, select those that are significant with diet patterns AND at least one BMI variable
patnames <- c("Preprocessed", "Mediterranean", "Western" )
tax2plot_a <-  all_daa_mod %>% filter(Contrast %in% patnames) %>%
  filter(padj <= 0.01 & !is.na(padj)) %>% 
  group_by(taxon) %>%
  dplyr::summarise(n_sig = n())

tax2plot_b <-  all_daa_mod %>% filter(Contrast == "Z-Score BMI T1 (Normal T0)" | Contrast == "Z-Score BMI T0") %>%
  filter(padj <= 0.01 & !is.na(padj)) %>% 
  group_by(taxon) %>%
  dplyr::summarise(n_sig = n())

tax2plot <- tax2plot_a %>% filter(taxon %in% tax2plot_b$taxon)

name2hms <- "SigOnePatternAndOneBMI"

####
nrow(tax2plot)
outdir <- paste0(reserva_outdir, name2hms, "/")
if(!dir.exists(outdir)) dir.create(outdir)

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
         filename = paste0(outdir, "heatmap1_scaledFC_", as.character(MIN_COMP_LIM), "comps.pdf"),
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
         filename = paste0(outdir, "heatmap2_trimmed0975_", as.character(MIN_COMP_LIM), "comps.pdf"),
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
         filename = paste0(outdir, "heatmap3_NSgray_", as.character(MIN_COMP_LIM), "comps.pdf"),
         height = 6, width = 12,
         show_colnames = FALSE)


## Sig per contrast:
sig_per_contrast <- all_daa_mod %>%
  dplyr::mutate(sig_p05 = !is.na(padj) & padj <= 0.05,
                sig_p01 = !is.na(padj) & padj <= 0.01,
                sig_p001 = !is.na(padj) & padj <= 0.001,
                sig_p05_and_lfc = sig_p05 & !is.na(log2FoldChangeShrink) & abs(log2FoldChangeShrink) >= LFCLIM,
                sig_p01_and_lfc = sig_p01 & !is.na(log2FoldChangeShrink) & abs(log2FoldChangeShrink) >= LFCLIM,
                sig_p001_and_lfc = sig_p001 & !is.na(log2FoldChangeShrink) & abs(log2FoldChangeShrink) >= LFCLIM) %>%
  dplyr::mutate(Contrast = ifelse(Samples == "Normal_T0", paste(Contrast, Samples, sep="__") , Contrast)) %>%
  group_by(Contrast) %>%
  dplyr::summarise_if(is.logical, sum) %>%
  arrange(desc(sig_p05))

write_tsv(sig_per_contrast, file = paste0(outdir, "number_of_taxa_significant_per_per_contrast.tsv"))

## Now cluster
## first find optimal number of clusters for both bacteria and variables
#MAX_CLUSTS<-30
#wcss <- numeric(MAX_CLUSTS)
#wcss_vars <- numeric(MAX_CLUSTS)
#for (k in 1:MAX_CLUSTS) {
#  km <- kmeans(mat2, centers = k, nstart = 100)
#  km_vars <- kmeans(t(mat2), centers = k, nstart = 100)
#  wcss[k] <- km$tot.withinss
#  wcss_vars[k] <- km_vars$tot.withinss
#}
#
## Plot elbow method
#par(mfrow=c(1, 2))
#plot(1:MAX_CLUSTS, wcss, type = "b", pch = 19, frame = FALSE,
#     main = "Clustering taxa",
#     xlab = "Number of clusters K",
#     ylab = "Total within-clusters sum of squares")
#plot(1:MAX_CLUSTS, wcss_vars, type = "b", pch = 19, frame = FALSE,
#     main = "Clustering variables",
#     xlab = "Number of clusters K",
#     ylab = "Total within-clusters sum of squares")
#
#
#library(cluster)
#
#sil_width <- numeric(MAX_CLUSTS-1)
#sil_width_vars <- numeric(MAX_CLUSTS-1)
#for (k in 2:MAX_CLUSTS) {
#  km <- kmeans(mat2, centers = k, nstart = 100)
#  km_vars <- kmeans(t(mat2), centers = k, nstart = 100)
#
#  ss <- silhouette(km$cluster, dist(mat2))
#  ss_vars <- silhouette(km_vars$cluster, dist(t(mat2)))
#  sil_width[k-1] <- mean(ss[, 3])
#  sil_width_vars[k-1] <- mean(ss_vars[, 3])
#}
#
#par(mfrow=c(1, 2))
#plot(2:MAX_CLUSTS , sil_width, type = "b", pch = 19, frame = FALSE,
#     xlab = "Number of clusters K",
#     ylab = "Average silhouette width")
#plot(2:MAX_CLUSTS , sil_width_vars, type = "b", pch = 19, frame = FALSE,
#     xlab = "Number of clusters K",
#     ylab = "Average silhouette width")
#best_k <- which.max(sil_width) + 1
#best_k
## not very good results
#
#library(factoextra)
#
#set.seed(123)
#gap_stat <- clusGap(mat2, FUN = kmeans, nstart = 25, K.max = 30, B = 50)
#fviz_gap_stat(gap_stat)
#
##################
#
#NCLUS <- 5
#clusts <- kmeans(mat2, centers = NCLUS)
#clusts_vars <- kmeans(mat2 %>% t, centers = NCLUS)
#
#annrow <- data.frame(taxon= names(clusts$cluster[rownames(mat2)]),
#                     cluster= paste0("Cluster ", as.character(clusts$cluster[rownames(mat2)]))) %>%
#  column_to_rownames("taxon")
#
#anncol <- data.frame(var= names(clusts_vars$cluster[colnames(mat2)]),
#                     cluster= paste0("Cluster ", as.character(clusts_vars$cluster[colnames(mat2)]))) %>%
#  column_to_rownames("var")
#
#bacpca <- prcomp(mat2)
#pcadf <- bacpca$x %>% as.data.frame %>%
#  rownames_to_column("taxon") %>%
#  dplyr::mutate(cluster=annrow[taxon, "cluster"])
#
#
#plots <- map(paste0("PC", 2:5), \(PC){
#  ggplot(pcadf, aes(x=PC1, y=!!sym(PC), col=cluster, fill=cluster)) +
#    geom_point() +
#    stat_ellipse() +
#    theme_minimal() +
#    #ggsci::scale_color_lancet() +
#    xlab(G4Micro::getPropVar(bacpca, "PC1")) +
#    ylab(G4Micro::getPropVar(bacpca, PC))
#
#})
#pdf(paste0(outdir, "PCA_and_kmeans_clustering_bacteria.pdf"), width = 10, height = 7)
#cowplot::plot_grid(plotlist = plots, ncol=2)
#dev.off()
#
### try another one
#
#library(mclust)
#mc <- Mclust(mat2)
#summary(mc)
#plot(mc)
#
#library(dbscan)
#db <- dbscan(mat2, eps = 0.5, minPts = 5)
#plot(db, data = mat2)
#
#####
#
#pheatmap(mat2 %>% t,
#         cluster_rows = TRUE, cluster_cols = TRUE,
#         annotation_row = anncol,
#         annotation_col = annrow,
#         filename = paste0(outdir, "heatmap2_annotKmeansK7.pdf"),
#         height = 7, width = 12,
#         show_colnames = FALSE)
#
#pheatmap(mat2 %>% t,
#         cluster_rows = TRUE, cluster_cols = TRUE,
#         annotation_row = anncol,
#         annotation_col = annrow,
#         filename = paste0(outdir, "heatmap2_annotKmeansK7_wardD2.pdf"),
#         clustering_method = "ward.D2",
#         height = 7, width = 12,
#         show_colnames = FALSE)

# Ok, cluster only Food variables

new_food_names <- all_daa_mod %>% filter(type == "Food groups") %>%
  pull(Contrast) %>% unique

matfood <- mat2[, new_food_names]
pheatmap(matfood %>% t,
         cluster_rows = TRUE, cluster_cols = TRUE,
         filename = paste0(outdir, "heatmap4_foodOnly_", as.character(MIN_COMP_LIM), "comps.pdf"),
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
         filename = paste0(outdir, "heatmap5_preOrder_", as.character(MIN_COMP_LIM), "comps_all_", pname, "_clust.pdf"),
         #clustering_method = "ward.D2",
         height = 8, width = 18, # height = 14 MIN_COMP_LIM = 4 (80 y algo taxa)
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
         filename = paste0(outdir, "heatmap5_preOrder_", as.character(MIN_COMP_LIM), "comps_asterisks_", pname, "_clust.pdf"),
         #clustering_method = "ward.D2",
         height = 8, width = 18, # height = 14 MIN_COMP_LIM = 4 (80 y algo taxa)
         gaps_col = sapply(varlist, length) %>% cumsum(),
         annotation_col = anncol,
         show_colnames = TRUE,
         angle_col = 45,           # rotate column names
         fontsize_col = 10,
         fontsize_row = 10,
         labels_row = labels_col,
         display_numbers = mat_pchar,
         fontsize_number=12,
         number_color = "white",
         color = col_fun,
         annotation_colors = color_list_cols)

mat_ord2 <- mat_ord
mat_ord2[mat_p2 >= 0.1] <- 0
pheatmap(mat_ord2,
         cluster_rows = TRUE, cluster_cols = FALSE,
         filename = paste0(outdir, "heatmap5_preOrder_", as.character(MIN_COMP_LIM), "_SigOnly_comps_asterisks_", pname, "_clust.pdf"),
         #clustering_method = "ward.D2",
         height = 7, width = 16,  # height = 14 MIN_COMP_LIM = 4 (80 y algo taxa)
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

################################################################33
######## END OF HEATMAPS

##################Lolipop plots
clean_names <- function(tax){
  gsub("_", " ", tax) %>% 
    gsub("[\\[\\]]", "", .)
}

makeLoliplot <- function(daa_df, 
                         vars2loliplot = c(), # Variables which LFC include in plot (all in column named 'Contrast')
                         vars2sort = c(), # Variables to sort taxa, ordered
                         vars2filter = c(), # Variables to select significant taxa 
                         outdir = "./",
                         name="test",
                         plim = 0.01, plim_col = 0.05, lfclim = 0,
                         strip_fontsize=9,
                         w=12, h=12
){
  
  if(length(vars2sort) == 0) vars2sort <- vars2loliplot
  if(length(vars2filter) == 0) vars2filter <- vars2loliplot
  
  usedf <- daa_df %>% 
    dplyr::filter(Contrast %in% vars2loliplot) %>% 
    dplyr::mutate(
      Sig = ifelse(!is.na(padj) & padj <= plim_col, ifelse(log2FoldChangeShrink < 0 , "Down", "Up"), "NS"), 
      Contrast = factor(Contrast, levels=vars2loliplot)
    )
  
  tax2use <- daa_df %>% 
    dplyr::filter(Contrast %in% vars2filter) %>% 
    dplyr::filter(padj < plim & abs(log2FoldChangeShrink) > lfclim) %>%
    dplyr::pull(taxon) %>% unique
  tax2use %>% length
  
  taxorder <- usedf %>%
    dplyr::filter(taxon %in% tax2use) %>%
    dplyr::select(taxon, Contrast, log2FoldChangeShrink, padj) %>%
    tidyr::gather("vart", "valt", log2FoldChangeShrink, padj) %>%
    unite("vart2", Contrast, vart, sep="__") %>%
    tidyr::spread(vart2, valt) %>% 
    dplyr::mutate(bmisig = "NS")
  
  taxorder_b <- purrr::reduce(vars2sort, ~ .y %>% dplyr::mutate(bmisig := ifelse(!!sym(.x) < plim_col, 
                                                                                 .x, bmisig )), 
                              .init= taxorder, .dir="backward") 
  
  newlevs <- if("NS" %in% taxorder_b$bmisig){c("NS", rev(vars2sort))}else{ c("NS", rev(vars2sort))}
  
  taxorder_b <- taxorder_b %>% 
    dplyr::mutate(bmisig = factor(bmisig, levels = newlevs)) %>% 
    dplyr::group_by(bmisig) %>% 
    dplyr::arrange(bmisig) %>% 
    group_split() 
  
  sortlevs <- purrr::map_vec(taxorder_b, ~ unique(.x[["bmisig"]])) %>% as.character
  sortlevs[sortlevs=="NS"] <- sortlevs[length(sortlevs)]
  taxorder_b <- purrr::map2(taxorder_b, sortlevs, ~ .x %>% 
                    dplyr::arrange(.data[[gsub("__padj", "__log2FoldChangeShrink", .y)]])) %>% 
    #purrr::map( \(x) {
    #  grname <- x %>% pull(bmisig) %>% unique()
    #  if(grname == "NS") return(x)
    #  x %>% 
    #    dplyr::arrange(x, .data[[gsub("__padj", "__log2FoldChangeShrink", grname)]])
    #  }) %>% #vars2sort[1] -> avoid NS to fail
    bind_rows() %>% 
    #dplyr::arrange(bmisig,`Z-Score BMI T0__log2FoldChangeShrink`) %>%
    dplyr::mutate(taxon=clean_names(taxon)) %>% 
    dplyr::mutate(taxon=factor(taxon, levels=taxon))
  
  TAXORDER <- taxorder_b$taxon
  
  usedf2 <- usedf %>%
    filter(taxon %in% tax2use) %>%
    dplyr::mutate(taxon=clean_names(taxon)) %>% 
    dplyr::mutate(taxon=factor(taxon, levels=TAXORDER))
  
  linedf <- taxorder_b %>%
    dplyr::group_by(bmisig) %>% 
    dplyr::summarise(xpos = max(as.numeric(taxon)) + 0.5) %>%
    head(nrow(.)-1)
  
  write_tsv(usedf2, file = paste0(outdir, "loliplot_", "_", as.character(plim), "_fc", as.character(lfclim), "_", name, ".tsv"))
  write_tsv(linedf, file = paste0(outdir, "loliplot_", "_", as.character(plim), "_fc", as.character(lfclim), "_", name, "_linedf.tsv"))
  
  (g0 <- ggplot(usedf2, aes(x=taxon, y=log2FoldChangeShrink, col=Sig, fill=Sig))+
      facet_grid(~ Contrast) +
      geom_hline(yintercept = 0, linetype=2, col="lightgray") +
      geom_vline(data=linedf, aes(xintercept=xpos), linetype=2, col="lightgray")+
      geom_segment(aes(x=taxon, xend = taxon, y=0, yend=log2FoldChangeShrink)) +
      geom_point() +
      theme_bw() +
      coord_flip() +
      scale_color_manual(values = c("Down"="steelblue", "Up"="tomato", "NS"="darkgray")) +
      theme(axis.text.y= element_text(face="italic", size=10),
            axis.text.x= element_text( size=12),
            axis.title = element_text(size=12),
            strip.text = element_text(size=strip_fontsize))
  )
  ggsave(filename = paste0(outdir, "loliplot_", "_", as.character(plim), "_fc", as.character(lfclim), "_", name, ".pdf"), g0, 
         width = w, height = h)
  result <- list(
    usedf=usedf2, 
    linedf=linedf, 
    plot=g0,
    plotname=paste0(outdir, "loliplot_", "_", as.character(plim), "_fc", as.character(lfclim), "_", name, ".pdf"),
    dfname= paste0(outdir, "loliplot_", "_", as.character(plim), "_fc", as.character(lfclim), "_", name, ".tsv"),
    dfname_lines=paste0(outdir, "loliplot_", "_", as.character(plim), "_fc", as.character(lfclim), "_", name, "_linedf.tsv"),
    input=daa_df, 
    args = list(vars2sort = vars2sort,
                vars2filter = vars2filter,
                name=name,
                plim = plim, plim_col = plim_col, lfclim = lfclim,
                w=w, h=h)
  )
  save(result, file = paste0(outdir, "loliplot_", "_", as.character(plim), "_fc", as.character(lfclim), "_", name, ".RData"))
  return(result)
  
}

all_daa <- read_tsv(paste0(outdir, "all_DAA_long.tsv"))
all_daa_mod <- read_tsv(paste0(outdir, "all_DAA_long_NamesModified.tsv"))

vars2loliplot <- c("Z-Score BMI T0",
                  "Z-Score Waist T0",
                 # "Change in Z-Score BMI",
                 # "Change in Z-Score Waist",
                  "Z-Score BMI T1 (Normal T0)",
                  "Z-Score Waist T1 (Normal T0)"
                  #"Change in Z-Score BMI (Normal T0)",
                 # "Change in Z-Score Waist (Normal T0)"
                 #"Status T1 Excessive Gain vs Normal (Normal T0)"
                  #"Body Fat % T0",
                 )

vars2filtersig <- vars2loliplot[1]
vars2sort <- paste0(vars2loliplot, "__padj")


lol1 <- makeLoliplot(all_daa_mod, vars2loliplot = vars2loliplot, 
                         vars2sort = vars2sort,
                         vars2filter = vars2filtersig,
                         outdir = outdir, 
                         name="BMI_T0_4",
                         plim = 0.01, plim_col = 0.05, lfclim = 0,
                         w=12, h=8
)

lol1 <- makeLoliplot(all_daa_mod, vars2loliplot = vars2loliplot, 
                     vars2sort = vars2sort,
                     vars2filter = vars2filtersig,
                     outdir = outdir, 
                     name="BMI_T0_2",
                     plim = 0.001, plim_col = 0.05, lfclim = 0.5,
                     w=12, h=6
)


### focus on T1

vars2loliplot_t1 <- c("Change in Z-Score BMI (Normal T0)",
                   "Change in Z-Score Waist (Normal T0)",
                   "Z-Score BMI T1 (Normal T0)",
                   "Z-Score Waist T1 (Normal T0)",
                   "Z-Score BMI T0",
                   "Z-Score Waist T0"
                   #"Change in Z-Score BMI (Normal T0)",
                   # "Change in Z-Score Waist (Normal T0)"
                   #"Status T1 Excessive Gain vs Normal (Normal T0)"
                   #"Body Fat % T0",
)

vars2filtersig_t1 <- c("Change in Z-Score BMI (Normal T0)",
                    "Change in Z-Score Waist (Normal T0)",
                    "Z-Score BMI T1 (Normal T0)",
                    "Z-Score Waist T1 (Normal T0)"
)
vars2sort_t1 <- paste0(vars2loliplot_t1, "__padj")


lol2 <- makeLoliplot(all_daa_mod, vars2loliplot = vars2loliplot_t1, 
                     vars2sort = vars2sort_t1,
                     vars2filter = vars2filtersig_t1,
                     outdir = outdir, 
                     name="BMI_T1_1",
                     plim = 0.01, plim_col = 0.05, lfclim = 0,
                     w=18, h=12
)

vars2loliplot_t1_b <- c("Change in Z-Score BMI (Normal T0)",
                        "Z-Score BMI T1 (Normal T0)",
                      "Change in Z-Score Waist (Normal T0)",
                      "Z-Score Waist T1 (Normal T0)"
)
vars2filtersig_t1_b <-vars2loliplot_t1_b[c(1, 2)]
vars2sort_t1_b <- paste0(vars2loliplot_t1_b, "__padj")

daa_filtered <- all_daa_mod %>% filter(Contrast %in% vars2loliplot_t1_b) %>% 
  dplyr::mutate(Contrast = gsub(" \\(Normal T0\\)", "", Contrast))

lol3 <- makeLoliplot(daa_filtered, 
                     vars2loliplot =  gsub(" \\(Normal T0\\)", "", vars2loliplot_t1_b) , 
                     vars2sort = gsub(" \\(Normal T0\\)", "", vars2sort_t1_b) ,
                     vars2filter = gsub(" \\(Normal T0\\)", "", vars2filtersig_t1_b) ,
                     outdir = outdir, 
                     name="BMI_T1_3_onlyNormalT0_",
                     plim = 0.01, plim_col = 0.05, lfclim = 0,
                     strip_fontsize=11,
                     w=12, h=7
)
################################################################################################################ 
###### Now plot LFC of different variables  against LFC of other variabels (each species is a point)

all_daa <- read_tsv(paste0(outdir, "all_DAA_long.tsv"))
all_daa_mod <- read_tsv(paste0(outdir, "all_DAA_long_NamesModified.tsv"))

xx <- all_daa_mod %>% filter(Contrast == "Mediterranean")
yy <- all_daa_mod %>% filter(Contrast == "Z-Score BMI T0")
zz <- merge(xx, yy, by="taxon")
zz$sig_BMI <- ifelse(zz$padj.y <= 0.01, "Sig. BMI", "NS BMI")
zz$sig_any <- ifelse(zz$padj.y <= 0.01 | zz$padj.x <= 0.01, "Sig. Any", "NS")
ggplot(zz, aes(x=log2FoldChangeShrink.y, y=log2FoldChangeShrink.x, col=sig_BMI, fill=sig_BMI)) +
  geom_point() +
  geom_smooth(method="lm") +
  theme_bw()

## not very ogood result

## Try canonical correlations
library(phyloseq)      # microbiome data handling
library(DESeq2)        # for vst normalization
library(vegan)         # ecological ordination, including CCA/RDA
library(CCA)           # canonical correlation analysis
library(mixOmics)

type_levels <- rev(c("Demographic",
                     "Demographic T1", "T1 (normal BMI at T0)",
                     "Diet pattern",
                     "Macronutrients",
                     "Food groups"))
patnames <- c("Preprocessed",  "Mediterranean", "Western")
food_western <- c("Dairy","Refined Cereals", "Sweets, Pastries", "Meat", "Tubers", "Sugars and Sweets", "Sauces, Condiments")
food_prep <- c("Juices, Softdrinks", "Dairy Derivatives" , "Prepared Foods", "Snacks Savory")
food_med <- c("Fruits", "Vegetables", "Fish", "Fats, Oils", "Eggs" , "Legumes", "Whole Grain Cereals", "Nuts", "Oleaginous Fruits")
food_macro <- c("Energy (Kcal)", "Carbohydrates (g)", "Fiber (g)", "Protein (g)", "Total Fat (g)")
patnames <- c("Preprocessed", "Mediterranean", "Western")

popvarscca <- c("z_t0", "z_waist_00", "inc_z_bmi", "inc_z_waist")
popvarscca_newnames <- c("Z-Score BMI T0", "Z-Score Waist T0", "Change in Z-Score BMI", "Change in Z-Score Waist")

all_var_types <- c(rep("Mediterranean", length(food_med)),
                    rep("Preprocessed", length(food_prep)),
                    rep("Western", length(food_western)),
                   rep("Demographic", length(popvarscca_newnames))
)
names(all_var_types) <- c(food_med, food_prep, food_western, popvarscca_newnames)


outdir_cca <- paste0(outdir, "CCA/")
if(!dir.exists(outdir_cca)) dir.create(outdir_cca)

all_daa <- read_tsv(paste0(outdir, "all_DAA_long.tsv"))
all_daa_mod <- read_tsv(paste0(outdir, "all_DAA_long_NamesModified.tsv"))

load("/home/carlos/Documentos/CORALS/results_rstudio/results_Agosto25_1//foodPCA/phyloseq_list_foodPCA_withNMF_withInc.RData")
phobj <- all_phyloseq$remove_tanda2
dds <- phyloseq_to_deseq2(phobj, ~ 1)  # no design needed, just for VST
s_meta <- sample_data(phobj) %>% data.frame()

quant_vars_diet <- names(s_meta)[grepl("_clr$", names(s_meta))]
diet_aggr_vars <- quant_vars_diet[1:4]
food_vars <- quant_vars_diet[5:length(quant_vars_diet)]

vars2select <- c(food_vars, popvarscca)

# Run variance stabilizing transformation
dds <- estimateSizeFactors(dds)
dds <- estimateDispersions(dds)
vst_counts <- assay(varianceStabilizingTransformation(dds, blind = TRUE))
vstdf_write <- vst_counts %>% t %>% as.data.frame %>%
  rownames_to_column("sampleID")
write_tsv(vstdf_write, paste0(outdir, "VST_counts.tsv"))

microbiome_mat <- vstdf_write %>% column_to_rownames("sampleID") %>% as.matrix
dim(microbiome_mat)

diet_mat <- sample_data(phobj) %>% data.frame %>%
  #dplyr::select( ends_with("clr")) %>% # # already sampleID in rownames
  dplyr::select(all_of(vars2select)) %>%
  #column_to_rownames("sampleID") %>%
  as.matrix()
#not_nas <- rowSums(diet_mat != 0) # use when only diet is being used

not_nas <- (rowSums(diet_mat[, food_vars] != 0) > 0) & (rowSums( is.na(diet_mat[, popvarscca])) == 0)
not_nas <- names(not_nas)[not_nas > 0]

tax2plot <- all_daa %>% filter(Contrast %in% colnames(diet_mat)) %>%
  filter(padj <= 0.001) %>% pull(taxon) %>% unique
length(tax2plot)
tax2plot <- c()

assertthat::assert_that(all(tax2plot %in% colnames(microbiome_mat)))

microbiome_mat_filt <- microbiome_mat[not_nas, ]
if(length(tax2plot)>0){
  microbiome_mat_filt <- microbiome_mat[, tax2plot]
}
diet_mat_filt <- diet_mat[not_nas, ]

microbiome_mat_filt_scale <- scale(microbiome_mat_filt)
with_nas <- apply(microbiome_mat_filt_scale, MAR=2, \(x)length(which(is.na(x)))) %>% sort
with_nas <- with_nas[with_nas > 0]
microbiome_mat_filt_scale <- microbiome_mat_filt_scale[, ! colnames(microbiome_mat_filt_scale) %in% names(with_nas)]

diet_mat_filt_scale <- scale(diet_mat_filt)

#cca_res <- cancor(microbiome_mat_filt_scale, diet_mat_filt_scale)
rcc_res <- rcc(microbiome_mat_filt_scale, diet_mat_filt_scale, method = "shrinkage")

# Inspect canonical correlations
#rcc_res$cor

#---------------------------------------------------------
# 4. Plot results
#---------------------------------------------------------

## Correlation circle plot (variables contribution)
#plotVar(rcc_res, comp = 1:2)
#
## Sample plot (canonical variates)
#plotIndiv(rcc_res, comp = c(1,2), group = NULL, legend = TRUE)
#
## Clustered image map of correlations
#cim(rcc_res, comp = 1:2)
#
## Network view of associations
#network(rcc_res, comp = 1:2, cutoff = 0.1)

save(rcc_res, file = paste0(outdir, "mixomics_RCC_allMicro_DietAndWeightVars.RData"))


df1 <- rcc_df <- rcc_res$loadings$X  %>% as.data.frame %>% rownames_to_column("variable") %>%
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

topn <- 20

df3 <- rbind(df1, df2) %>%
  merge(df_cor, by="variable") %>%
  dplyr::mutate(label_var = ifelse(class=="Variable", variable, "")) %>%
  dplyr::mutate(dist_to_center = sqrt(scale(V1)^2 + scale(V2)^2)) %>%
  group_by(class) %>%
  dplyr::arrange(desc(dist_to_center)) %>%
  dplyr::mutate(dist_order = 1:n(),
                label_mic = ifelse(dist_order <= topn & class != "Variable", variable, "")) %>%
  ungroup() %>%
  dplyr::mutate(
    label_mic = gsub("_", " ", label_mic),
    label_mic = gsub("[\\[\\]]", "", label_mic, perl=TRUE),
  ) %>%
  dplyr::mutate(
    label_var = gsub("_clr$", "", label_var),
    label_var = gsub("_g$", " (g)", label_var, perl=T),
    label_var = gsub("_kcal", " (Kcal)", label_var),

    #label_var = ifelse(Samples == "Normal_T0", paste(label_var, " (normal T0)", sep="") , label_var),

    label_var = gsub("z_t1", "Z-score BMI T1", label_var, perl=T),
    label_var = gsub("z_waist_01", "Z-score waist T1", label_var),

    label_var = gsub("z_t0", "Z-score BMI T0", label_var, perl=T),
    label_var = gsub("z_t1", "Z-score BMI T1", label_var, perl=T),

    label_var = gsub("inc_z_waist", "Change in Z-score waist", label_var, perl=T),
    label_var = gsub("inc_z_bmi", "Change in Z-score BMI", label_var, perl=T),

    label_var = gsub("_c1", " T0", label_var),
    label_var = gsub("_00", " T0", label_var),
    label_var = gsub("_01", " T1", label_var),
    label_var = gsub("mg_p", "Body fat %", label_var),
    label_var = gsub("nreads", "seq. depth", label_var),
    label_var = gsub("sauces_con", "sauces, con", label_var),
    label_var = gsub("juices_so", "juices, so", label_var),
    label_var = gsub("fats_oils", "fats, oils", label_var),
    label_var = gsub("sweets_pas", "sweets, pas", label_var),
    label_var = gsub("^z_", "Z-score ", label_var, perl=T),
    label_var = gsub("bmi", "BMI", label_var, perl=T),
    label_var = gsub("status_c2", "status T1", label_var, perl=T),
    label_var = gsub("\\.", " ", label_var, perl=T),
    label_var = gsub("_", " ", label_var)
  ) %>% dplyr::mutate(label_var = tools::toTitleCase(label_var)) %>%
  dplyr::mutate(type = ifelse(class == "Variable",
                              all_var_types[label_var],
                              "Microbiome")) %>%
  dplyr::mutate(type = factor(type, levels=c( "Demographic","Mediterranean", "Preprocessed", "Western", "Microbiome"))) %>%
  dplyr::mutate(label_full = ifelse(class == "Microbiome" & label_mic== "", "",
                  ifelse(class=="Variable", paste0("plain('", label_var, "')"), paste0("italic('", label_mic, "')"))))

cols <- ggsci::pal_aaas()(4)
cols <- c(cols, "gray20")

df3 <- df3 %>% dplyr::mutate(
  #Size = ifelse(dist_order <= topn, 0.5, 0.1),
  Alpha = ifelse(dist_order <= topn, 1, 0.2)
)

(g0 <- ggplot(df3, aes(x=V1, y=V2, col=type)) +
    geom_vline(xintercept = 0, linetype=2, col="lightgray") +
    geom_hline(yintercept = 0, linetype=2, col="lightgray") +
  geom_point(aes(alpha=Alpha)) + #size=Size,
 ggrepel::geom_text_repel(data=df3 %>% filter(label_full != "") , aes(label = label_full), parse=TRUE) +
  #ggrepel::geom_text_repel(data=df3 %>% filter(class != "Variable" & dist_order <= topn),
  #                           aes(label = paste0("italic('", label_mic, "')")),
  #                         parse = TRUE) +
  theme_classic() +
  scale_color_manual(values=cols) +
  scale_fill_manual(values=cols)

)
  ggsave(filename = paste0(outdir_cca, "CCA_alltaxa_foodgroupsAndWeight.pdf"), g0, width = 8, height = 8)

  (g0 <- ggplot(df3, aes(x=corV1, y=corV2, col=type)) +
      geom_vline(xintercept = 0, linetype=2, col="lightgray") +
      geom_hline(yintercept = 0, linetype=2, col="lightgray") +
      geom_point(aes(alpha=Alpha)) + #size=Size,
      ggrepel::geom_text_repel(data=df3 %>% filter(label_full != "") , aes(label = label_full), parse=TRUE) +
      #ggrepel::geom_text_repel(data=df3 %>% filter(class != "Variable" & dist_order <= topn),
      #                           aes(label = paste0("italic('", label_mic, "')")),
      #                         parse = TRUE) +
      theme_classic() +
      scale_color_manual(values=cols) +
      scale_fill_manual(values=cols)

  )
  ggsave(filename = paste0(outdir_cca, "CCA_alltaxa_foodgroupsAndWeight_cor2CC.pdf"), g0, width = 8, height = 8)

# Other tests; not working very well
spls_res <- spls(microbiome_mat, diet_mat, ncomp = 2, keepX = 24, keepY = 10)
plotVar(spls_res, comp = c(1,2))




