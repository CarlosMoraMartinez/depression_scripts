library(tidyverse)

outdir <- paste0("/home/carlos/Escritorio/202311_DEPRESION/ReviewJune2025/Results_rstudio/results2/", "read_number_plots/")
if(!dir.exists(outdir)) dir.create(outdir)

multiqc <- read_csv("/home/carlos/Escritorio/202311_DEPRESION/ReviewJune2025/Data/num_reads/multiqc_table_mod.csv") %>%
  dplyr::mutate(Classified  = 100 - pct_Unclassified) %>%
  dplyr::mutate(sampleID = gsub("G4M", "", Sample_Name))%>%
  dplyr::mutate(sampleID = gsub("^0", "", sampleID, perl=T))

raw_reads <- read_tsv("/home/carlos/Escritorio/202311_DEPRESION/ReviewJune2025/Data/num_reads/raw.count.txt") %>%
  mutate_if(is.numeric, ~ .x/2) %>%
  dplyr::mutate(sampleID = gsub("G4M", "", sampleID))%>%
  dplyr::mutate(sampleID = gsub("^0", "", sampleID, perl=T))
trim_reads <- read_tsv("/home/carlos/Escritorio/202311_DEPRESION/ReviewJune2025/Data/num_reads/trim.count.txt") %>%
  mutate_if(is.numeric, ~ .x/2)%>%
  dplyr::mutate(sampleID = gsub("G4M", "", sampleID))%>%
  dplyr::mutate(sampleID = gsub("^0", "", sampleID, perl=T))
filt_reads <- read_tsv("/home/carlos/Escritorio/202311_DEPRESION/ReviewJune2025/Data/num_reads/filt.count.txt") %>%
  mutate_if(is.numeric, ~ .x/2) %>%
  dplyr::mutate(sampleID = gsub("G4M", "", sampleID)) %>%
  dplyr::mutate(sampleID = gsub("^0", "", sampleID, perl=T))

all(s_meta$sampleID %in% raw_reads$sampleID)
all(s_meta$sampleID %in% trim_reads$sampleID)
all(s_meta$sampleID %in% filt_reads$sampleID)
all(s_meta$sampleID %in% multiqc$sampleID)

s_meta_read_data <- s_meta %>%
  filter(Tanda == 1) %>%
  merge(raw_reads, by="sampleID") %>%
  merge(trim_reads, by="sampleID") %>%
  merge(filt_reads, by="sampleID") %>%
  merge(multiqc, by="sampleID") %>%
  dplyr::mutate(
    pct_trimmed_calc = (raw_reads - trimmed_reads)/raw_reads,
    pct_filt_calc = (trimmed_reads - filt_reads)/raw_reads,
    pct_filt_calc_rel2trim = (trimmed_reads - filt_reads)/trimmed_reads,
    num_classified = Classified*filt_reads/100
  ) %>%
  dplyr::mutate(group = paste(Condition, ob_o_sobrepeso, sep=":"))

plot(round(s_meta_read_data$trimmed_reads, 1),
     2*s_meta_read_data$M_Seqs_trimmed_r1)
round(s_meta_read_data$trimmed_reads/1e6, 1) - 2*s_meta_read_data$M_Seqs_trimmed_r1

hist(s_meta_read_data$num_classified - 2*s_meta_read_data$reads)

s_meta_long <- s_meta_read_data %>%
  select(sampleID, group, Condition, ob_o_sobrepeso, raw_reads, trimmed_reads, filt_reads) %>%
  tidyr::gather(key="read_type", value="num_reads", raw_reads, trimmed_reads, filt_reads) %>%
  dplyr::mutate(`million reads` = num_reads/1e6) %>%
  filter(! is.na(ob_o_sobrepeso)) %>%
  dplyr::mutate(read_type = ifelse(
    read_type == "raw_reads", "Raw Reads",
    ifelse(read_type == "filt_reads", "Host filtered", "Trimmed Reads")
  )) %>%
  dplyr::mutate(read_type = factor(read_type, levels=c("Raw Reads", "Trimmed Reads", "Host filtered")))



signif_levels=c("***"=0.001, "**"=0.01, "*"=0.05, "ns"=1.1)

comp <- combn(unique(s_meta_long$group), 2, simplify = F) %>% lapply(as.character)
num_comparisons <- 3*length(comp)
signif_levels_bonferroni <- c(signif_levels[1:3]/num_comparisons, signif_levels[4])

greads <- ggplot(s_meta_long, aes(x=group, y=`million reads`, col=group)) +
  facet_grid( ~ read_type)+
  geom_boxplot(width=0.5) +
  geom_jitter() +
  ggsignif::stat_signif(test="wilcox.test", na.rm=T, comparisons = comp,
                        step_increase=0.06,
                        tip_length = 0.01,
                        map_signif_level=signif_levels_bonferroni,
                        vjust=0.3,
                        color = "black"
  ) +
  theme_bw() +
  theme(legend.position = "none") +
  theme(axis.text.x = element_blank(),
        strip.text = element_text(size=12),
        axis.title.x = element_text(size=12),
        axis.title.y = element_text(size=12),
        axis.text.y = element_text(size=12)) #element_text(size = 11, colour = "black", angle = 45, vjust=1, hjust=1)



df_filt <- s_meta_read_data %>%
  filter(!is.na(ob_o_sobrepeso))

comp <- combn(unique(df_filt$group), 2, simplify = F) %>% lapply(as.character)
num_comparisons <- length(comp)
signif_levels_bonferroni <- c(signif_levels[1:3]/num_comparisons, signif_levels[4])



gclassif <- ggplot(df_filt,
                   aes(x=group, y=Classified, col=group)) +
  geom_boxplot(width=0.5) +
  geom_jitter() +
  ylab("% Classified reads") +
  ggsignif::stat_signif(test="wilcox.test", na.rm=T, comparisons = comp,
                        step_increase=0.06,
                        tip_length = 0.01,
                        map_signif_level=signif_levels_bonferroni,
                        vjust=0.3,
                        color = "black"
  ) +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        strip.text = element_text(size=12),
        axis.title.x = element_text(size=12),
        axis.title.y = element_text(size=12),
        axis.text.y = element_text(size=12) )




library(patchwork)

load("/home/carlos/Escritorio/202311_DEPRESION/ReviewJune2025/Results_rstudio/results2/00_plot_rarecurve_all.RData")

rardf <- grar$df %>%
  merge(s_meta_read_data, by.x="sample", by.y="sampleID") %>%
  filter(!is.na(ob_o_sobrepeso)) %>%
  filter(Tanda == 1)

dfmax <- grar$df_max %>%
  merge(s_meta_read_data, by.x="sample", by.y="sampleID") %>%
  filter(!is.na(ob_o_sobrepeso))%>%
  filter(Tanda == 1)

min_reads <- min(dfmax$nreads) #<- grar$min_reads_all

(ggrare <- ggplot(rardf, aes(col=group, x=nreads, y=txcount, group=sample)) +
  geom_line() +
  geom_point(data=dfmax)+
  #geom_text_repel(data=dfmax, aes(label=as.character(sample)))+
  geom_vline(xintercept=min_reads, col="gray40", linetype=2) +
  xlab("Sample Size") +
  ylab("Species") +
  ggtitle("Rarefaction curve") +
  scale_x_continuous(labels = function(x) format(x, scientific = TRUE)) +
    theme_bw() +
    theme(axis.text.x = element_text(size=12),
          strip.text = element_text(size=12),
          axis.title.x = element_text(size=12),
          axis.title.y = element_text(size=12),
          axis.text.y = element_text(size=12),
          plot.title = element_text(hjust = 0.5))
  #theme(legend.position="none")
)

row1 <- greads + gclassif + plot_layout(widths = c(3, 1))
row1

final_plot <- row1 / ggrare +
  plot_layout(heights = c(1, 1.2)) +
  plot_annotation(tag_levels = 'A')
final_plot

ggsave(paste0(outdir, "num_reads_patchwork_noTanda2.pdf"), final_plot, width = 10, height = 10)

s_meta_print <- s_meta_read_data %>%
  select(sampleID, group, raw_reads, trimmed_reads, filt_reads, Classified, pct_trimmed_calc, pct_filt_calc, pct_filt_calc_rel2trim) %>%
  mutate(pct_trimmed_calc=100*pct_trimmed_calc,
         pct_filt_calc=100*pct_filt_calc,
         pct_filt_calc_rel2trim=100*pct_filt_calc_rel2trim) %>%
  dplyr::rename(pct_classified = Classified)

write_tsv(s_meta_print, paste0(outdir, "reads_per_sample_table_noTanda2.tsv"))

summ <- s_meta_print %>% group_by(group) %>%
  select(-sampleID) %>%
  summarise_all(.funs = list(mean=mean, median=median, min=min, max=max)) %>%
  dplyr::mutate(across(matches("reads"), ~ round(.x/1e6,1) )) %>%
  column_to_rownames("group") %>%
  as.matrix %>%
  t %>%
  as.data.frame %>%
  rownames_to_column("Var") %>%
  arrange(Var)

write_tsv(summ, paste0(outdir, "reads_per_sample_Summary_noTanda2.tsv"))
table(s_meta_read_data$group)


