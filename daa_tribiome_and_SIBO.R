
phseq_to_use <- names(all_phyloseq)[!grepl("tribiome", names(all_phyloseq))]

interestvar <- "Group"

opt <- restaurar(opt)
opt$reserva_0 <- opt$out
opt$out <- paste0(opt$out, "DeSEQ2/")
opt$minfreq <- 0
if(!dir.exists(opt$out)) dir.create(opt$out)

daa_all <- list()
for(phname in phseq_to_use){
  cat("Doing DESeq2 Analysys for: ", phname, "\n")
  phobj <- all_phyloseq[[phname]]
  daa_all[[phname]] <- list()
  
  name <- paste0(phname, '_', interestvar)
  res_tmp <-deseq_full_pipeline(phobj, name, interestvar, opt)
  daa_all[[phname]][[interestvar]] <- res_tmp
  
}

phseq_to_use <- names(all_phyloseq)[grepl("tribiome", names(all_phyloseq))]
interestvar <- "Weight"
for(phname in phseq_to_use){
  cat("Doing DESeq2 Analysys for: ", phname, "\n")
  phobj <- all_phyloseq[[phname]]
  daa_all[[phname]] <- list()
  
  name <- paste0(phname, '_', interestvar)
  res_tmp <-deseq_full_pipeline(phobj, name, interestvar, opt)
  daa_all[[phname]][[interestvar]] <- res_tmp
  
}

opt$out <- opt$reserva_0
save(daa_all, file=paste0(opt$out, "DeSEQ2/DESEQ2_tribiome_and_SIBO.RData"))
#load(paste0(opt$out, "DESeq2_ControlVarsAlone/DESEQ2_controlVarsAlone_all.RData"))


#### P faecium in thin vs overweight
## by ASV
phobj <- all_phyloseq$tribiome_raw
metad <- sample_data(phobj) %>% data.frame

pfaecium <- tax_table(phobj) %>% data.frame %>% 
  filter(grepl("faecium", Species)) %>% 
  filter(grepl("Phascolarcto", Genus))

asvdf <- daa_all$tribiome_raw$Weight$norm_counts_df %>% 
  dplyr::rename("ASV"="gene") %>% 
  filter(ASV %in% rownames(pfaecium))
nrow(asvdf)

rawdf <- daa_all$tribiome_raw$Weight$raw_df %>% 
  dplyr::rename("ASV"="gene") %>% 
  filter(ASV %in% rownames(pfaecium))

perdf <- daa_all$tribiome_raw$Weight$raw_df %>% 
  mutate(across(all_of(metad$sampleID), ~ 100*.x/sum(.x) )) %>% 
  dplyr::rename("ASV"="gene") %>% 
  filter(ASV %in% rownames(pfaecium))

asvsp <- daa_all$tribiome_raw_Species$Weight$norm_counts_df %>% 
  dplyr::rename("Species"="gene") %>% 
  filter(grepl("faecium", Species))

per_sp <- daa_all$tribiome_raw_Species$Weight$raw_df %>% 
  dplyr::rename("Species"="gene") %>% 
  mutate(across(all_of(metad$sampleID), ~ 100*.x/sum(.x) )) %>% 
  filter(grepl("faecium", Species)) %>% 
  filter(grepl("Phasco", Species))

nrow(rawdf)

prevalence <-rawdf %>% column_to_rownames("ASV") %>% 
  mutate_all(~ ifelse(.x > 0, 1, 0)) %>% 
  rowSums()

prev_any <- rawdf %>% column_to_rownames("ASV") %>% 
  mutate_all(~ ifelse(.x > 0, 1, 0)) %>% 
  colSums() %>% as.logical() %>% sum
  
  
prevalence_normal <-rawdf %>% 
  column_to_rownames("ASV") %>% 
  select(all_of(metad %>% filter(Weight=="Normal weight") %>% pull(sampleID))) %>% 
  mutate_all(~ ifelse(.x > 0, 1, 0)) %>% 
  rowSums()
prev_normal_any <- rawdf %>% 
  column_to_rownames("ASV") %>% 
  select(all_of(metad %>% filter(Weight=="Normal weight") %>% pull(sampleID))) %>% 
  mutate_all(~ ifelse(.x > 0, 1, 0)) %>% 
  colSums() %>% as.logical() %>% sum

prevalence_overweight <-rawdf %>% 
  column_to_rownames("ASV") %>% 
  select(all_of(metad %>% filter(Weight =="Overweight") %>% pull(sampleID))) %>% 
  mutate_all(~ ifelse(.x > 0, 1, 0)) %>% 
  rowSums()

prev_overweight_any <- rawdf %>% 
  column_to_rownames("ASV") %>% 
  select(all_of(metad %>% filter(Weight =="Overweight") %>% pull(sampleID))) %>% 
  mutate_all(~ ifelse(.x > 0, 1, 0)) %>% 
  colSums() %>% as.logical() %>% sum

prevalence <- prevalence[prevalence>0]

prev_df <- data.frame(ASV=names(prevalence), 
                      Prevalence=prevalence, 
                      Prevalence_normal = prevalence_normal[names(prevalence)],
                      Prevalence_overweight = prevalence_overweight[names(prevalence)],
                      pct_samples = 100*prevalence/(ncol(rawdf)-1) %>% round(2), 
                      pct_normal = 100*prevalence_normal[names(prevalence)]/length(which(metad$Weight=="Normal weight")) %>% round(2), 
                      pct_overweight = 100*prevalence_overweight[names(prevalence)]/length(which(metad$Weight=="Overweight")) %>% round(2) 
                      ) %>% 
  arrange(desc(prevalence)) %>% 
  rbind(data.frame(
    ASV = "All ASVs",
    Prevalence=prev_any, 
    Prevalence_normal = prev_normal_any,
    Prevalence_overweight = prev_overweight_any,
    pct_samples = 100*prev_any/(ncol(rawdf)-1) %>% round(2),
    pct_normal = prev_normal_any/length(which(metad$Weight=="Normal weight")) %>% round(2),
    pct_overweight = prev_overweight_any /length(which(metad$Weight=="Overweight"))%>% round(2)
    
  ))
  
write_tsv(prev_df, file = paste0(opt$out, "Prevalence_Pfaecium_tribiome.tsv"))

df2plot <- asvsp %>% gather("sampleID", "norm_abundance", all_of(metad$sampleID)) %>% 
  mutate(Percentage = per_sp[1, sampleID] %>% unlist) %>% 
  merge(metad, by="sampleID", all.x=T, all.y=F)

wilcox.test(df2plot$norm_abundance[df2plot$Weight=="Normal weight"], df2plot$norm_abundance[df2plot$Weight=="Overweight"])
t.test(df2plot$norm_abundance[df2plot$Weight=="Normal weight"], df2plot$norm_abundance[df2plot$Weight=="Overweight"])

wilcox.test(df2plot$Percentage[df2plot$Weight=="Normal weight"], df2plot$Percentage[df2plot$Weight=="Overweight"])
t.test(df2plot$Percentage[df2plot$Weight=="Normal weight"], df2plot$Percentage[df2plot$Weight=="Overweight"])


g0 <- ggplot(df2plot, aes(x=Weight, y=norm_abundance, col=Weight)) +
  geom_boxplot(aes_string(fill = v), alpha = 0.7, width=0.5, fill="white") +
  geom_point() +
  scale_color_uchicago() +
  labs(title = v, x = '') +
  theme_pubclean() +
  ylab("Normalized Abundance") +
  mytheme +
  ggsignif::stat_signif(test="wilcox.test", na.rm=T, comparisons = list(c("Normal weight", "Overweight")), 
                        step_increase=0.06,
                        tip_length = 0.01,
                        map_signif_level=signif_levels_bonferroni,
                        vjust=0.3,
                        color = "black"
  ) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(axis.text.x = element_blank()) 
ggsave(paste0(opt$out, "/Pfaecium_norm_abundance.pdf"), g0, width = 4, height = 4)

g1 <- ggplot(df2plot, aes(x=Weight, y=Percentage, col=Weight)) +
  geom_boxplot(aes_string(fill = v), alpha = 0.7, width=0.5, fill="white") +
  geom_point() +
  scale_color_uchicago() +
  labs(title = v, x = '') +
  theme_pubclean() +
  ylab("Relative Abundance") +
  mytheme +
  ggsignif::stat_signif(test="wilcox.test", na.rm=T, comparisons = list(c("Normal weight", "Overweight")), 
                        step_increase=0.06,
                        tip_length = 0.01,
                        map_signif_level=signif_levels_bonferroni,
                        vjust=0.3,
                        color = "black"
  ) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(axis.text.x = element_blank()) 
ggsave(paste0(opt$out, "/Pfaecium_percentage.pdf"), g1, width = 4, height = 4)
## by Species

### P. faecium prevalence
## ASVs

## Total species

