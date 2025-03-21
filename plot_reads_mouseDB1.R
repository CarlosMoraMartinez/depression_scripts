

multiqc_original <- read_tsv(opt$multiqcsummary) %>% arrange(Sample_Name)
multiqc_mousedb <- read_tsv(opt$multiqcsummary_mousedb) %>% arrange(Sample_Name)

assertthat::assert_that(all(multiqc_original$Sample_Name == multiqc_mousedb$Sample_Name))

multiqc <- multiqc_original %>% 
  mutate(Unclassified =multiqc_mousedb$Unclassified )

xmeta <- sample_data(all_phyloseq$filt) %>% data.frame
original_numreads <- read_tsv(opt$original_numreads) %>% 
  dplyr::mutate(original_reads = 2*(Unique_Reads + Duplicate_Reads))

assertthat::assert_that(all(xmeta$Sample.ID %in% multiqc$Sample_Name))
assertthat::assert_that(all(multiqc$Sample_Name %in% xmeta$Sample.ID))

assertthat::assert_that(all(xmeta$Sample.ID %in% original_numreads$Sample))
assertthat::assert_that(all(original_numreads$Sample %in% xmeta$Sample.ID))

met2plot <- merge(xmeta, multiqc, by.x= "Sample.ID", by.y="Sample_Name") %>% 
  merge(original_numreads, by.x= "Sample.ID", by.y="Sample") %>% 
  dplyr::mutate(pct_used = 100*num_reads/original_reads,
                pct_mapped_2total = Aligned*(1-0.01*Dropped),
                pct_unclassified_2total = Unclassified*(1-0.01*pct_mapped_2total-0.01*Dropped),
                total_discarded = pct_mapped_2total+Dropped + pct_unclassified_2total,
                total_calculated = total_discarded + pct_used
  ) %>% dplyr::mutate(
                pct_used = pct_used*100/total_calculated,
                pct_dropped = Dropped*100/total_calculated,
                pct_mapped_2total = pct_mapped_2total*100/total_calculated,
                pct_unclassified_2total = pct_unclassified_2total*100/total_calculated,
                total_discarded2 = pct_mapped_2total+pct_dropped + pct_unclassified_2total,
                total_calculated2 = total_discarded2 + pct_used
                
                )
  


outdir <- paste0(opt$out, "plot_reads")
if(!dir.exists(outdir)) dir.create(outdir)


(greads <- ggplot(met2plot, aes(x=Treatment,
                             y = log10(nreads), 
                             fill=Treatment,
                             col=Treatment))+
    facet_grid(Stress ~ Region_sequenced) +
    #geom_violin(alpha=0.6)+
    geom_boxplot(width=0.7, fill="white")+
    geom_point() +
    scale_color_npg() +
    theme_bw() +
    #scale_fill_npg() +
    ylab("log10(classified reads)")+
    theme(axis.text.x = element_text(vjust=1,hjust=1, angle = 45))
  
)
ggsave(filename = paste0(outdir, "/reads_filtPhylum_log10.pdf"), greads, width = 7, height = 4)


(greads <- ggplot(met2plot, aes(x=Treatment,
                                y = nreads, 
                                fill=Treatment,
                                col=Treatment))+
    facet_grid(Stress ~ Region_sequenced) +
    #geom_violin(alpha=0.6)+
    geom_boxplot(width=0.7, fill="white")+
    geom_point() +
    scale_color_npg() +
    theme_bw() +
    #scale_fill_npg() +
    theme(axis.text.x = element_text(vjust=1,hjust=1, angle = 45))
  
)
ggsave(filename = paste0(outdir, "/reads_filtPhylum_noLog.pdf"), greads, width = 7, height = 4)


(greads <- ggplot(met2plot, aes(x=Treatment,
                                y = Dropped, 
                                fill=Treatment,
                                col=Treatment))+
    facet_grid(Stress ~ Region_sequenced) +
    #geom_violin(alpha=0.6)+
    geom_boxplot(width=0.7, fill="white")+
    geom_point() +
    scale_color_npg() +
    theme_bw() +
    ylab("% of reads dropped during trimming") +
    #scale_fill_npg() +
    theme(axis.text.x = element_text(vjust=1,hjust=1, angle = 45))
  
)
ggsave(filename = paste0(outdir, "/reads_pct_trimmed.pdf"), greads, width = 7, height = 4)

###
(greads <- ggplot(met2plot, aes(x=Treatment,
                                y = Aligned, 
                                fill=Treatment,
                                col=Treatment))+
    facet_grid(Stress ~ Region_sequenced) +
    #geom_violin(alpha=0.6)+
    geom_boxplot(width=0.7, fill="white")+
    geom_point() +
    scale_color_npg() +
    theme_bw() +
    ylab("% of reads aligned to mouse genome") +
    #scale_fill_npg() +
    theme(axis.text.x = element_text(vjust=1,hjust=1, angle = 45))
  
)
ggsave(filename = paste0(outdir, "/reads_pct_aligned_to_mouse.pdf"), greads, width = 7, height = 4)

###
(greads <- ggplot(met2plot, aes(x=Treatment,
                                y = Unclassified, 
                                fill=Treatment,
                                col=Treatment))+
    facet_grid(Stress ~ Region_sequenced) +
    #geom_violin(alpha=0.6)+
    geom_boxplot(width=0.7, fill="white")+
    geom_point() +
    scale_color_npg() +
    theme_bw() +
    ylab("% of reads unclassified by Kraken2") +
    #scale_fill_npg() +
    theme(axis.text.x = element_text(vjust=1,hjust=1, angle = 45))
  
)
ggsave(filename = paste0(outdir, "/reads_pct_unclassified.pdf"), greads, width = 7, height = 4)

## total reads
(greads <- ggplot(met2plot, aes(x=Treatment,
                                y = log10(original_reads), 
                                fill=Treatment,
                                col=Treatment))+
    facet_grid(Stress ~ Region_sequenced) +
    #geom_violin(alpha=0.6)+
    geom_boxplot(width=0.7, fill="white")+
    geom_point() +
    scale_color_npg() +
    theme_bw() +
    #scale_fill_npg() +
    ylab("log10(raw reads)")+
    theme(axis.text.x = element_text(vjust=1,hjust=1, angle = 45))
  
)
ggsave(filename = paste0(outdir, "/original_reads_filtPhylum_log10.pdf"), greads, width = 7, height = 4)

###
met2plotb <- met2plot %>% mutate(used2=100*num_reads/original_reads)
(greads <- ggplot(met2plotb, aes(x=Treatment,
                                y = used2, 
                                fill=Treatment,
                                col=Treatment))+
    facet_grid(Stress ~ Region_sequenced) +
    #geom_violin(alpha=0.6)+
    geom_boxplot(width=0.7, fill="white")+
    geom_point() +
    scale_color_npg() +
    theme_bw() +
    ylab("% of reads used") +
    #scale_fill_npg() +
    theme(axis.text.x = element_text(vjust=1,hjust=1, angle = 45))
  
)
ggsave(filename = paste0(outdir, "/reads_pct_used.pdf"), greads, width = 7, height = 4)


met2long <- met2plot %>% 
  select(sampleID, Mouse, Treatment, Region_sequenced, Stress, pct_used, pct_mapped_2total, pct_unclassified_2total, pct_dropped) %>% 
  gather("PCT", "value", pct_used, pct_mapped_2total, pct_unclassified_2total, pct_dropped) %>% 
  dplyr::mutate(PCT=factor(PCT),
                PCT = fct_recode(PCT, Classified="pct_used", Unclassified = "pct_unclassified_2total", 
                                 Mouse="pct_mapped_2total", Dropped="pct_dropped" ))

met2plotlist <- met2long %>% 
  mutate(Treatment = ifelse(Treatment == "TRANSFER", "NO ABS", Treatment)) %>% 
  mutate(Mouse = ifelse(Mouse == "transfer", "T", Mouse)) %>% 
  arrange(Stress, Treatment) %>% 
  group_by(Treatment, Stress) %>% 
  group_split() %>% 
  map(\(xx){
    
    name <- paste(unique(xx$Stress), unique(xx$Treatment), sep=" - ", collapse=" - ")
    greads <- ggplot(xx, aes(x=Mouse,
                                y = value, 
                                fill=PCT,
                                col=PCT, 
                                group=PCT))+
    facet_grid( ~ Region_sequenced) +
    geom_col() +
    scale_color_npg() +
      scale_fill_npg() +
    theme_bw() +
      ylab("% of reads") + 
      xlab("Mouse") +
    #scale_fill_npg() +
      ggtitle(name) + 
      theme_classic()+
      theme(plot.title = element_text(hjust = 0.5, vjust=0.5))+
      theme(axis.text.x = element_text(size = 12))+
      theme(strip.text.x = element_text(size = 10))+
      theme(strip.text.y = element_text(size = 8))+
      theme(axis.title.y = element_text(size = 12))+
      theme(axis.title.x = element_text(size = 12)) +
      theme(legend.position = "none")
  

  })

pdf(paste0(outdir, "/reads_summary_all_nolegend.pdf"), width = 16, height = 7)
print(cowplot::plot_grid(plotlist = met2plotlist, ncol=5))
dev.off()

met2plotlist <- met2long %>% 
  mutate(Treatment = ifelse(Treatment == "TRANSFER", "NO ABS", Treatment)) %>% 
  mutate(Mouse = ifelse(Mouse == "transfer", "T", Mouse)) %>% 
  arrange(Stress, Treatment) %>% 
  group_by(Treatment, Stress) %>% 
  group_split() %>% 
  map(\(xx){
    
    name <- paste(unique(xx$Stress), unique(xx$Treatment), sep=" - ", collapse=" - ")
    greads <- ggplot(xx, aes(x=Mouse,
                             y = value, 
                             fill=PCT,
                             col=PCT, 
                             group=PCT))+
      facet_grid( ~ Region_sequenced) +
      geom_col() +
      scale_color_npg() +
      scale_fill_npg() +
      theme_bw() +
      ylab("% of reads") + 
      xlab("Mouse") +
      #scale_fill_npg() +
      ggtitle(name) + 
      theme_classic()+
      theme(plot.title = element_text(hjust = 0.5, vjust=0.5))+
      theme(axis.text.x = element_text(size = 12))+
      theme(strip.text.x = element_text(size = 10))+
      theme(strip.text.y = element_text(size = 8))+
      theme(axis.title.y = element_text(size = 12))+
      theme(axis.title.x = element_text(size = 12))
    
    
  })

pdf(paste0(outdir, "/reads_summary_all_legend.pdf"), width = 20, height = 7)
print(cowplot::plot_grid(plotlist = met2plotlist, ncol=5))
dev.off()

# ONLY CONTROLS

met2long_onlyNoAbs <- met2long %>% filter(Treatment == "NO ABS" & Stress == "Control") %>% 
  mutate(`Read assignment` = PCT)

(greads <- ggplot(met2long_onlyNoAbs, aes(x=Region_sequenced,
                                 y = value, 
                                 fill=`Read assignment`,
                                 col=`Read assignment`))+
    facet_grid(. ~ PCT) +
    #geom_violin(alpha=0.6)+
    geom_boxplot(width=0.7, fill="white")+
    geom_point() +
    scale_color_npg() +
    theme_bw() +
    ylab("% of reads") +
    #scale_fill_npg() +
    xlab("Region sequenced") +
    theme(axis.text.x = element_text(vjust=1,hjust=1, angle = 45))
  
)
ggsave(filename = paste0(outdir, "/reads_pct_only_Controls.pdf"), greads, width = 7, height = 3)


## Compare new and old Unclassified



met2plot_original <- merge(xmeta, multiqc_original, by.x= "Sample.ID", by.y="Sample_Name") %>% 
  merge(original_numreads, by.x= "Sample.ID", by.y="Sample") %>% 
  dplyr::mutate(pct_used = 100*num_reads/original_reads,
                pct_mapped_2total = Aligned*(1-0.01*Dropped),
                pct_unclassified_2total = Unclassified*(1-0.01*pct_mapped_2total-0.01*Dropped),
                total_discarded = pct_mapped_2total+Dropped + pct_unclassified_2total,
                total_calculated = total_discarded + pct_used
  ) %>% dplyr::mutate(
    pct_used = pct_used*100/total_calculated,
    pct_dropped = Dropped*100/total_calculated,
    pct_mapped_2total = pct_mapped_2total*100/total_calculated,
    pct_unclassified_2total = pct_unclassified_2total*100/total_calculated,
    total_discarded2 = pct_mapped_2total+pct_dropped + pct_unclassified_2total,
    total_calculated2 = total_discarded2 + pct_used
    
  )

assertthat::assert_that(all(met2plot_original$sampleID == met2plot$sampleID))

met2plot <- met2plot %>% mutate(
  Unclassified_StandardDB = met2plot_original$Unclassified,
  pct_used_StandarDB = met2plot_original$pct_used,
  pct_aligned = pct_used/(Unclassified+ pct_used),
  pct_aligned_StandarDB = pct_used_StandarDB/(Unclassified_StandardDB+ pct_used_StandarDB)
)


gm <- ggplot(met2plot, aes(x =pct_aligned_StandarDB, y =pct_aligned, 
                    col= Region_sequenced, fill=Region_sequenced,
                    shape=Treatment)) +
  geom_point(size=2)+
  geom_abline(intercept = 0, slope=1, col="gray", linetype=2) +
  theme_classic()+
  xlab("Classified / non-mouse - Standard DB") +
  ylab("Classified / non-mouse - Mouse DB") +
  scale_color_nejm()
ggsave(filename = paste0(outdir, "/PCT_classified_reads_ComparisonDatabases.pdf"), gm, width = 6, height = 4)
