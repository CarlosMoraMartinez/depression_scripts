library(tidyverse)
library(wesanderson)
library(DESeq2)
# use this! use without log

pal <- wes_palette("AsteroidCity1", 2)

opt <- restaurar(opt)
load("/home/carlos/Escritorio/202311_DEPRESION/ReviewJune2025/Results_rstudio/results2/phyloseq_original//phyloseq_all_list.RData")

phseq_to_correct <- names(all_phyloseq)[4] # "remove_tanda2"
interestvar <- "Condition"


DANAME <- "DESeq2_AgeSexInteraction_log/"
opt$reserva_0 <- opt$out
opt$out <- paste0(opt$out, DANAME)
opt$mincount <- 1
opt$minsampleswithcount <- opt$minfreq * nrow(sample_data(all_phyloseq$remove_tanda2))

if(!dir.exists(opt$out)) dir.create(opt$out)

s_meta <- sample_data(all_phyloseq$remove_tanda2) %>% data.frame

  phname <- phseq_to_correct
  cat("Doing DESeq2 Analysys with correction for: ", phname, "\n")
  phobj <- all_phyloseq[[phname]]
  phobj <- updatePsWithLogs(phobj, c("Edad", "BMI"))
  sample_data(phobj)$Age_log <- sample_data(phobj)$Edad_log
  sample_data(phobj)$Age <- sample_data(phobj)$Edad
  sample_data(phobj)$Sex <- factor(ifelse(sample_data(phobj)$Sexo == "Hombre", "Male", "Female"), levels=c("Male", "Female"))

  dds <- phyloseq_to_deseq2(phobj, design =~ Condition + Sex * Age_log)
  dds <- dds[rowSums(counts(dds) >= opt$mincount) >= opt$minsampleswithcount,]
  dds <- DESeq(dds, test = "LRT", reduced = ~ Condition + Sex + Age_log )
  results_interaction <- results(dds)

  resLFC <- lfcShrink(dds, coef= "SexFemale.Age_log", type="ashr", lfcThreshold = log2(opt$fc)) %>%
    data.frame
  dds_nobmi <- dds

save(dds, file=paste0(opt$out, "/DESEQ2_Condition_TestSexAgeInteraction.RData"))
write_tsv(resLFC, file = paste0(opt$out, "/DESEQ2_Condition_TestSexAgeInteraction.tsv"))


samples_with_nas <- sample_data(phobj) %>% data.frame
samples_with_nas <- is.na(samples_with_nas$BMI)
samples <- sample_data(phobj)$sampleID[!samples_with_nas]
phobj_filt <- phyloseq::prune_samples(samples, phobj)

dds <- phyloseq_to_deseq2(phobj_filt, design =~ Condition + BMI_log + Sex * Age_log)
dds <- dds[rowSums(counts(dds) >= opt$mincount) >= opt$minsampleswithcount,]

dds <- DESeq(dds, test = "LRT", reduced = ~ Condition + BMI_log + Sex + Age_log)
results_interaction_BMI <- results(dds)

resLFC_BMI <- lfcShrink(dds, coef= "SexFemale.Age_log", type="ashr", lfcThreshold = log2(opt$fc)) %>%
  data.frame


save(dds, file=paste0(opt$out, "/DESEQ2_ConditionBMI_TestSexAgeInteraction.RData"))
write_tsv(resLFC_BMI, file = paste0(opt$out, "/DESEQ2_ConditionBMI_TestSexAgeInteraction.tsv"))

raw_counts <- counts(dds)

PLIM <- 0.05
LFC_LIM <- 1
s_meta <- sample_data(phobj) %>% data.frame

norm_counts <- counts(dds, normalized = T)
norm_counts_df <- defWriteMatAsDF(norm_counts, opt, "ConditionBMI_TestSexAgeInteraction_norm_counts.tsv")

s_meta$sampleID[!(s_meta$sampleID %in% colnames(norm_counts_df))]

vstds <- varianceStabilizingTransformation(raw_counts, blind=F)
vst_counts_df <- defWriteMatAsDF(vstds, opt, "ConditionBMI_TestSexAgeInteraction_vst_counts.tsv" )

vst_counts_long <- vst_counts_df %>%
  gather("sample", "vst_abn", -gene)


norm_counts_long <- norm_counts_df %>%
  gather("sample", "norm_abn", -gene) %>%
  merge(vst_counts_long, by=c("gene","sample" )) %>%
  merge(s_meta, by.x="sample", by.y="sampleID", all.y=F, all.x=T)

write_tsv(x = norm_counts_long, paste0(opt$out, "vst_and_norm_counts_Long_withMetadata.tsv"))

norm_counts_long_filt <- norm_counts_long %>%
  filter(gene %in% (resLFC_BMI %>% filter(padj <  PLIM) %>% rownames)) %>%
  dplyr::mutate(gene = gsub("_", " ", gene),
                gene = gsub("[\\[\\]]", "", gene, perl=T))


(gage <- ggplot(norm_counts_long_filt, aes(
  x=Edad,
  y= norm_abn,
  group=Sex,
  col=Sex,
  fill=Sex,
  shape=Condition)) +
    facet_wrap(~ gene, scales="free") +
   geom_smooth(alpha=0.5, method="lm") +
    geom_point() +
    theme_bw()+
    xlab("Age") +
    ylab("VST abundance") +
    scale_color_manual(values=pal) +
    scale_fill_manual(values=pal) +
    theme(strip.text = element_text(face="italic"))
  )
ggsave(filename = paste0(opt$out, "/InteractionTaxa_abundance_line.pdf"), gage, width = 12, height = 8)


(gage <- ggplot(norm_counts_long_filt, aes(
  x=Edad,
  y= vst_abn,
  group=Sex,
  col=Sex,
  fill=Sex,
  shape=Condition)) +
    facet_wrap(~ gene, scales="free") +
    geom_smooth(alpha=0.5, method="lm") +
    geom_point() +
    theme_bw()+
    xlab("Age") +
    ylab("VST abundance") +
    scale_color_manual(values=pal) +
    scale_fill_manual(values=pal) +
    theme(strip.text = element_text(face="italic"))
)
ggsave(filename = paste0(opt$out, "/InteractionTaxa_abundance_line_VST.pdf"), gage, width = 11, height = 7)

(gage <- ggplot(norm_counts_long_filt, aes(
  x=Edad,
  y= vst_abn,
  group=Sex,
  col=Sex,
  fill=Sex,
  shape=Condition)) +
    facet_wrap(~ gene, scales="free") +
    geom_smooth(alpha=0.5) +
    geom_point() +
    theme_bw()+
    xlab("Age") +
    ylab("VST abundance") +
    scale_color_manual(values=pal) +
    scale_fill_manual(values=pal) +
    theme(strip.text = element_text(face="italic"))
)
ggsave(filename = paste0(opt$out, "/InteractionTaxa_abundance_line_VST_GAM.pdf"), gage, width = 11, height = 7)

## F prausnitzii alone

norm_counts_long_filt2 <- norm_counts_long %>%
  filter(grepl("prausnit", gene)) %>%
  dplyr::mutate(gene = gsub("_", " ", gene),
                gene = gsub("[\\[\\]]", "", gene, perl=T))


(gage <- ggplot(norm_counts_long_filt2, aes(
  x=Edad,
  y= norm_abn,
  group=Sex,
  col=Sex,
  fill=Sex,
  shape=Condition)) +
    facet_wrap(~ gene, scales="free") +
    geom_smooth(alpha=0.5, method="lm") +
    geom_point() +
    theme_bw()+
    xlab("Age") +
    ylab("VST abundance") +
    scale_color_manual(values=pal) +
    scale_fill_manual(values=pal) +
    theme(strip.text = element_text(face="italic"))
)
ggsave(filename = paste0(opt$out, "/InteractionTaxa_abundance_line_Fprausni.pdf"), gage, width = 12, height = 8)


(gage <- ggplot(norm_counts_long_filt2, aes(
  x=Edad,
  y= vst_abn,
  group=Sex,
  col=Sex,
  fill=Sex,
  shape=Condition)) +
    facet_wrap(~ gene, scales="free") +
    geom_smooth(alpha=0.5, method="lm") +
    geom_point() +
    theme_bw()+
    xlab("Age") +
    ylab("VST abundance") +
    scale_color_manual(values=pal) +
    scale_fill_manual(values=pal) +
    theme(strip.text = element_text(face="italic"))
)
ggsave(filename = paste0(opt$out, "/InteractionTaxa_abundance_line_VST_Fprausni.pdf"), gage, width = 5, height = 3)

(gage <- ggplot(norm_counts_long_filt2, aes(
  x=Edad,
  y= vst_abn,
  group=Sex,
  col=Sex,
  fill=Sex,
  shape=Condition)) +
    facet_wrap(~ gene, scales="free") +
    geom_smooth(alpha=0.5) +
    geom_point() +
    theme_bw()+
    xlab("Age") +
    ylab("VST abundance") +
    scale_color_manual(values=pal) +
    scale_fill_manual(values=pal) +
    theme(strip.text = element_text(face="italic"))
)
ggsave(filename = paste0(opt$out, "/InteractionTaxa_abundance_line_VST_GAM_Fprausni.pdf"), gage, width = 11, height = 7)




library(outliers)

norm_counts_long_filt2 <- norm_counts_long_filt %>%
  dplyr::mutate(
    norm_abn_olrm = norm_abn,
    sampleID = sample  ) %>%
  group_by(gene) %>%
  group_modify(.f = ~ {
    aux <- data.frame(.x)
    ols <- detectOutliers(aux, "norm_abn", p_lim = 0.05)
    .x$norm_abn_olrm[.x$sampleID %in% ols$sampleID ] <- NA
   return(.x)
  })

write_tsv(x = norm_counts_long, paste0(opt$out, "vst_and_norm_counts_Long_withMetadata_outliers_p05_Identified.tsv"))

(gage <- ggplot(norm_counts_long_filt2, aes(
  x=Edad,
  y= norm_abn_olrm,
  group=Sex,
  col=Sex,
  fill=Sex,
  shape=Condition)) +
    facet_wrap(~ gene, scales="free") +
    geom_smooth(alpha=0.5, method="lm") +
    geom_point() +
    theme_bw()+
    xlab("Age") +
    ylab("VST abundance") +
    scale_color_manual(values=pal) +
    scale_fill_manual(values=pal) +
    theme(strip.text = element_text(face="italic"))
)
ggsave(filename = paste0(opt$out, "/InteractionTaxa_abundance_line_rmOutliers.pdf"), gage, width = 11, height = 7)


(gage <- ggplot(norm_counts_long_filt2, aes(
  x=Edad,
  y= norm_abn_olrm,
  group=Sex,
  col=Sex,
  fill=Sex,
  shape=Condition)) +
    facet_wrap(~ gene, scales="free") +
    geom_smooth(alpha=0.5) +
    geom_point() +
    theme_bw()+
    xlab("Age") +
    ylab("VST abundance") +
    scale_color_manual(values=pal) +
    scale_fill_manual(values=pal) +
    theme(strip.text = element_text(face="italic"),
          panel.grid = element_blank(),                     # remove all gridlines
          panel.background = element_rect(fill = "white"),  # set panel bg to white
          panel.border = element_rect(color = "black"))
)
ggsave(filename = paste0(opt$out, "/InteractionTaxa_abundance_line_rmOutliers_GAM.pdf"), gage, width = 12, height = 8)


## Now do the analysis

opt$out <- paste0(opt$out, "DESeq2/")
if(!dir.exists(opt$out)) dir.create(opt$out)

designs <- list(
  Cond="~ Condition",
  BMI="~ BMI_log",
  CondSexIAge = "~ Condition + Sex * Age_log",
  BMISexIAge ="~ BMI_log + Sex * Age_log",
  CondBMISexIAge ="~ Condition + BMI_log + Sex * Age_log"
)

dds_all <- list()
for(ff in names(designs)){
  cat(paste0("Making DESeq2 for: ", designs[[ff]]), "\n")
  opt <- restaurar(opt)
  opt$out <- paste0(opt$out,  DANAME, ff, "/")
  opt$mincount <- 1
  opt$minsampleswithcount <- opt$minfreq * nrow(sample_data(all_phyloseq$remove_tanda2))

  if(!dir.exists(opt$out)) dir.create(opt$out)
  dds_all[[ff]] <- getDeseqResults(phobj_filt, opt = opt, name = ff, variables=NULL, formula=designs[[ff]])
  dds_all[[ff]]$formula <- designs[[ff]]
}
opt <- restaurar(opt)
save(dds_all, file = paste0(opt$out, DANAME, "/all_DESeq2.RData"))

####

#Sólo he guardado resdf
firstContrast <- dds_all$Cond
contrastlist2 <- list(
  dds_all$CondSexIAge$all_contrasts$Condition_Depression_vs_Control$resdf,
  dds_all$CondBMISexIAge$all_contrasts$Condition_Depression_vs_Control$resdf,
  dds_all$BMI$resdf,
  dds_all$BMISexIAge$all_contrasts$BMI_log$resdf,
  dds_all$CondBMISexIAge$all_contrasts$BMI_log$resdf

) %>%
  lapply(\(x)return(list(
      resdf= x
  ))) %>%
  lapply(\(x){
      x$resdf$log2FoldChangeShrink <- x$resdf$log2FoldChange
      x$resdf$lfcSE_Shrink <- x$resdf$lfcSE
      return(x)
  })

name2remove2 <- "xxx"
names(contrastlist2) <- c(
  "D_vs_C_adj.Age*Sex",
  "D_vs_C_adj.Age*Sex+BMI",
  "BMI_alone",
  "BMI_adj.Age*Sex",
  "BMI_adj.Age*Sex+Depr")

outdir <- paste0(opt$out, DANAME, "/Integrate1/")
if(!dir.exists(outdir) ) dir.create(outdir)

mainContrastName <- "Depression vs Control"
contrastNamesOrdered2 <- c("Depression vs Control",  gsub("_", " ", names(contrastlist2)))
compareLFCContrats(contrastlist2, firstContrast,
                   contrastNamesOrdered2, mainContrastName,
                   plim_select= 0.05, plim_plot=0.05,
                   name2remove = name2remove2,
                   resdfname="resdf",
                   outdir = outdir,
                   name="LFC_Comparison_AgeISex_BMI_allCombos_p05", w=12, h=12, scale_mode = "free")


compareLFCContrats(contrastlist2, firstContrast,
                   contrastNamesOrdered2, mainContrastName,
                   plim_select= 0.001, plim_plot=0.05,
                   name2remove = name2remove2,
                   resdfname="resdf",
                   outdir = outdir, name="LFC_Comparison_AgeISex_BMI_allCombos_pem3", w=12, h=12,
                   scale_mode = "free")
dea2contrasts <- list(firstContrast = firstContrast, contrastlist2=contrastlist2)
save(dea2contrasts, file = paste0(outdir, "/LFC_Comparison_AgeISex_BMI_allCombos.RData"))

phobj <- phobj_filt
save(phobj, file = paste0(outdir, "/phyloseq_used_remove_tanda2.RData"))

opt <- restaurar(opt)


## Check Age alone
ddsAge <- phyloseq_to_deseq2(phobj, design =~Age_log)
ddsAge <- ddsAge[rowSums(counts(ddsAge) >= opt$mincount) >= opt$minsampleswithcount,]
ddsAge <- DESeq(ddsAge)
resLFC_age <- results(ddsAge)

resLFC_age <- lfcShrink(ddsAge, coef= "Age_log", type="ashr", lfcThreshold = log2(opt$fc)) %>%
  data.frame

resLFC_age %>% rownames_to_column("taxon") %>% filter(grepl("prausn", taxon)) # padj=0.0825

write_tsv(resLFC_age, file = paste0(opt$out, "DESeqDAA_Age_log_alone.tsv"))

##
ddsAge <- phyloseq_to_deseq2(phobj, design =~Age)
ddsAge <- ddsAge[rowSums(counts(ddsAge) >= opt$mincount) >= opt$minsampleswithcount,]
ddsAge <- DESeq(ddsAge)
resLFC_age <- results(ddsAge)

resLFC_age <- lfcShrink(ddsAge, coef= "Age", type="ashr", lfcThreshold = log2(opt$fc)) %>%
  data.frame

resLFC_age %>% rownames_to_column("taxon") %>% filter(grepl("prausn", taxon)) # padj=0.099

write_tsv(resLFC_age, file = paste0(opt$out, "DESeqDAA_Age_alone.tsv"))


## Age only in women

samples <- sample_data(phobj)$sampleID[sample_data(phobj)$Sex == "Female"]
phobj_women <- phyloseq::prune_samples(samples, phobj)

ddsAge <- phyloseq_to_deseq2(phobj_women, design =~Age)
ddsAge <- ddsAge[rowSums(counts(ddsAge) >= opt$mincount) >= opt$minsampleswithcount,]
ddsAge <- DESeq(ddsAge)
resLFC_age <- results(ddsAge)

resLFC_age <- lfcShrink(ddsAge, coef= "Age", type="ashr", lfcThreshold = log2(opt$fc)) %>%
  data.frame

resLFC_age %>% rownames_to_column("taxon") %>% filter(grepl("prausn", taxon)) # padj=0.18

write_tsv(resLFC_age, file = paste0(opt$out, "DESeqDAA_Age_alone_WomenOnly.tsv"))

##

samples <- sample_data(phobj)$sampleID[sample_data(phobj)$Sex == "Female"]
phobj_women <- phyloseq::prune_samples(samples, phobj)

ddsAge <- phyloseq_to_deseq2(phobj_women, design =~Age_log)
ddsAge <- ddsAge[rowSums(counts(ddsAge) >= opt$mincount) >= opt$minsampleswithcount,]
ddsAge <- DESeq(ddsAge)
resLFC_age <- results(ddsAge)

resLFC_age <- lfcShrink(ddsAge, coef= "Age_log", type="ashr", lfcThreshold = log2(opt$fc)) %>%
  data.frame

resLFC_age %>% rownames_to_column("taxon") %>% filter(grepl("prausn", taxon)) # padj=0.138

write_tsv(resLFC_age, file = paste0(opt$out, "DESeqDAA_Age_alone_log_WomenOnly.tsv"))


