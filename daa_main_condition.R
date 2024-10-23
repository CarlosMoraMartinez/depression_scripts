library(ggpmisc)
###########################
# Make combinations:
treatment_order <- c("REG3", "REG2", "REG1", "NO ABS", "ABS") %>% rev
region_order <- c("REG3", "REG2", "REG1") %>% rev
sd_order <- c("SD", "CONTROL") %>% rev

metad <- meta3 %>% filter(Treatment != "TRANSFER") #%>% 
  #dplyr::mutate(Treatment = factor(Treatment, levels = c("REG3", "REG2", "REG1", "NO ABS", "ABS")))

# Cominations treatment
combins_treatment <- combn(unique(metad$Treatment), 2, simplify = F) %>% 
  lapply(\(x)x[order(match(x, treatment_order))])
combins_treatment2 <- append(combins_treatment, lapply(combins_treatment, \(x) paste("SD", x, sep=" ")))

combins_treatment_all <- map(paste0("REG", 1:3), \(x) lapply(combins_treatment2, \(cc) paste(cc, x, sep=":"))) %>% flatten()

# Combinations Region
combins_region <- combn(unique(metad$Region_sequenced), 2, simplify = F) %>% 
  lapply(\(x)x[order(match(x, region_order))])
combins_region_all <- map(unique(unlist(combins_treatment2)), \(x) lapply(combins_region, \(cc) paste(x, cc, sep=":"))) %>% flatten()

# Combinations Stress
combins_stress_all <- metad %>% filter(Stress != "SD") %>% pull(Treatment_region) %>% unique %>% 
  lapply(\(x) c(paste0("SD ", x), x))

# Combinations to same level:

base_cond <- "ABS:REG1"
combins_tobaseline <- metad %>% filter(Treatment_region != base_cond) %>% 
  pull(Treatment_region) %>% unique %>% 
  lapply(\(x) c(base_cond, x))


ALL_COMBINS <- append(combins_treatment_all, combins_region_all) %>% 
  append(combins_stress_all) %>% 
  append(combins_tobaseline) %>% 
  lapply(\(x) c("Treatment_region", x))

ALL_COMBINS <- lapply(ALL_COMBINS, \(x) gsub(" |:", ".", x, perl=TRUE)) %>% unique
###########################

source(opt$r_functions)
opt <- restaurar(opt)

daa_all <- list()
vars2deseq <- c("Treatment_region")
vars2heatmap <- c("Treatment", "Region_sequenced", "Stress", "flowcell")
opt$mincount <- 10
opt$minsampleswithcount <- 3
phseq_to_use <- names(all_phyloseq)[1:2]
deseqname = "DeSEQ2_v8/"

for(phname in phseq_to_use){
  # cat("Doing DESeq2 Analysys for: ", phname, "\n")
  phobj <- all_phyloseq[[phname]] 
  sdata <- sample_data(phobj) %>% data.frame
  
  sample_data(phobj)$is_transfer <- ifelse(sdata$Treatment == "TRANSFER", TRUE, FALSE)
  sample_data(phobj)$Treatment <- ifelse(sdata$Treatment == "TRANSFER", "NO ABS", sdata$Treatment)
  sample_data(phobj)$Treatment_region <- gsub("TRANSFER", "NO ABS", sdata$Treatment_region)
  sample_data(phobj)$Treatment_region <- gsub(" |:", ".", sample_data(phobj)$Treatment_region, perl=TRUE)
  
  daa_all[[phname]] <-deseq_full_pipeline(phobj, 
                                           phname, vars2deseq, opt, 
                                           doPoscounts=FALSE, 
                                           all_combins=ALL_COMBINS,
                                          plot_all = TRUE,
                                          deseqname = deseqname,
                                          vars2heatmap = vars2heatmap
                                          )
  
}
opt <- restaurar(opt)
#save(daa_all, file = paste0(opt$out, deseqname, "DESEQ2_all.RData"))
#load(paste0(opt$out, deseqname, "DESEQ2_all.RData"))

combs_cp <- ALL_COMBINS
names(combs_cp) <- sapply(combs_cp, \(x) paste(x[1], x[3], "vs", x[2], sep="_"))
vars2heatmap2 <- vars2heatmap[1:2]

for(phname in phseq_to_use){
  phobj <- all_phyloseq[[phname]] 
  sdata <- sample_data(phobj) %>% data.frame
  sdata_simp <- sdata %>% select(Treatment_region, Treatment, Region_sequenced, Stress) %>% distinct
  by1 <- join_by(group1 == Treatment_region)
  by2 <- join_by(group2 == Treatment_region)
  
  sample_data(phobj)$is_transfer <- ifelse(sdata$Treatment == "TRANSFER", TRUE, FALSE)
  sample_data(phobj)$Treatment <- ifelse(sdata$Treatment == "TRANSFER", "NO ABS", sdata$Treatment)
  sample_data(phobj)$Treatment_region <- gsub("TRANSFER", "NO ABS", sdata$Treatment_region)
  sample_data(phobj)$Treatment_region <- gsub(" |:", ".", sample_data(phobj)$Treatment_region, perl=TRUE)
  sdata <- sample_data(phobj) %>% data.frame
  sdata_simp <- sdata %>% select(Treatment_region, Treatment, Region_sequenced, Stress) %>% distinct
  
  deseq_results <- daa_all[[phname]]$all_contrasts
  
  deseq_df <- tibble(
    phname = phname,
    cond_name = names(deseq_results),
    contrast_list = deseq_results
  ) %>% 
    dplyr::mutate(
      variable = sapply(combs_cp, \(x) x[1]),
      group1 = sapply(combs_cp, \(x) x[2]),
      group2 = sapply(combs_cp, \(x) x[3])
    ) %>% 
    left_join(sdata_simp, by1) %>% 
    left_join(sdata_simp, by2, , suffix = c("_group1", "_group2")) %>% 
    dplyr::mutate(
      diff_p05 = map(contrast_list, \(x)x$resdf %>% filter(padj < 0.05) %>% pull(taxon)),
      diff_p01 = map(contrast_list, \(x)x$resdf %>% filter(padj < 0.01) %>% pull(taxon)),
      samples_used = map2(group1, group2, \(g1, g2) sdata %>% filter(Treatment_region %in% c(g1, g2)) %>% pull(sampleID))
    )
  
  save(deseq_df, file = paste0(opt$out, deseqname, phname, "/DESEQ2_all_results_", phname, "_table.RData"))
  
  # Now make personalized heatmaps
  df2plot <- if(! nrow(daa_all[[phname]]$vst_counts_df)){daa_all[[phname]]$norm_counts_df}else{daa_all[[phname]]$vst_counts_df}
  
  # make_heatmap_subset(dearesult, df2plot, taxa, samples, vars2heatmap, dds, name, opt)
  opt$out <- paste0(opt$out, deseqname, phname, "/HeatmapsRegion/")
  if(!dir.exists(opt$out)) dir.create(opt$out)
  
  cond_df_reg <- deseq_df %>% 
    dplyr::filter(Treatment_group1 == Treatment_group2) %>% 
    dplyr::filter(Stress_group1 == Stress_group2) %>%
    dplyr::filter(Region_sequenced_group1 != Region_sequenced_group2) %>% 
    group_by(Treatment_group1, Stress_group1)
  cond_df_reg %>% 
    group_walk(~{
      hmname <- paste0(phname, "_Region_in_", unique(.x$Treatment_group2 %>% unlist), "_", unique(.x$Stress_group2 %>% unlist))
      print(hmname)
      write_tsv(.x, file=paste0(opt$out, hmname, ".tsv"))
      hmname <- paste0(phname, "_Region_in_", unique(.x$Treatment_group2 %>% unlist), "_", unique(.x$Stress_group2 %>% unlist), "_p05")
      make_heatmap_subset(.x$contrast_list, 
                          df2plot, 
                          .x$diff_p05%>% unlist %>% unique, 
                          .x$samples_used %>% unlist %>% unique, 
                          vars2heatmap, 
                          daa_all[[phname]]$dds, 
                          hmname, 
                          opt)
      hmname <- paste0(phname, "_Region_in_", unique(.x$Treatment_group2 %>% unlist), "_", unique(.x$Stress_group2 %>% unlist), "_p01")
      make_heatmap_subset(.x$contrast_list, 
                          df2plot, 
                          .x$diff_p01%>% unlist %>% unique, 
                          .x$samples_used %>% unlist %>% unique, 
                          vars2heatmap, 
                          daa_all[[phname]]$dds, 
                          hmname, 
                          opt)
    })
  
    opt <- restaurar(opt)
    
    ## By treatment, for each Region:Stress
    # make_heatmap_subset(dearesult, df2plot, taxa, samples, vars2heatmap, dds, name, opt)
    opt$out <- paste0(opt$out, deseqname, phname, "/HeatmapsTreatment/")
    if(!dir.exists(opt$out)) dir.create(opt$out)
    
    cond_df_treat <- deseq_df %>% 
      dplyr::filter(Region_sequenced_group1 == Region_sequenced_group2) %>% 
      dplyr::filter(Stress_group1 == Stress_group2) %>%
      dplyr::filter(Treatment_group1 != Treatment_group2) %>% 
      group_by(Region_sequenced_group1, Stress_group1)
    cond_df_treat %>% 
      group_walk(~{
        hmname <- paste0(phname, "_Treatment_in_", unique(.x$Region_sequenced_group2 %>% unlist), "_", unique(.x$Stress_group2 %>% unlist))
        print(hmname)
        write_tsv(.x, file=paste0(opt$out, hmname, ".tsv"))
        hmname <- paste0(phname, "_Treatment_in_", unique(.x$Region_sequenced_group2 %>% unlist), "_", unique(.x$Stress_group2 %>% unlist), "_p05")
        make_heatmap_subset(.x$contrast_list, 
                            df2plot, 
                            .x$diff_p05%>% unlist %>% unique, 
                            .x$samples_used %>% unlist %>% unique, 
                            vars2heatmap, 
                            daa_all[[phname]]$dds, 
                            hmname, 
                            opt)
        hmname <- paste0(phname, "_Treatment_in_", unique(.x$Region_sequenced_group2 %>% unlist), "_", unique(.x$Stress_group2 %>% unlist), "_p01")
        make_heatmap_subset(.x$contrast_list, 
                            df2plot, 
                            .x$diff_p01%>% unlist %>% unique, 
                            .x$samples_used %>% unlist %>% unique, 
                            vars2heatmap, 
                            daa_all[[phname]]$dds, 
                            hmname, 
                            opt)
      })
    
    
    ## Bigger heatmaps
    opt <- restaurar(opt)
    opt$out <- paste0(opt$out, deseqname, phname, "/HeatmapsRegion/")
    sdata2 <- sdata %>% dplyr::mutate(
      Treatment = factor(Treatment, levels = c("ABS", "NO ABS", "REG1", "REG2", "REG3")),
      Region_sequenced = factor(Region_sequenced, levels = c("REG1", "REG2", "REG3")), 
      Stress = factor(Stress, levels = c("Control", "SD"))
    ) %>% arrange(Treatment, Region_sequenced, Stress) %>% 
      mutate(sampleID = factor(sampleID, levels=sampleID)) %>% 
      group_by(Stress) %>% 
      group_split()
    
    df2plot_list <- sdata2 %>% 
      map( \(x) df2plot %>% dplyr::select(gene, all_of(x$sampleID)) %>% 
          column_to_rownames("gene") %>% 
          as.matrix
      )
    names(df2plot_list) <- sapply(sdata2, \(x)unique(x$Stress) %>% as.character) 
    
    df2plot_list %>% map2(sdata2, \(mat, sdata_parcial){
      
      taxalist <- deseq_df %>% 
        dplyr::filter(Stress_group1 == Stress_group2) %>% 
        dplyr::filter(Stress_group1 == unique(sdata_parcial$Stress)) %>% 
        pull(diff_p01) %>% unlist %>% unique
      name <- paste0("/heatmap_allRegions_", unique(sdata_parcial$Stress), "_p01.pdf")
      make_full_region_heatmap(mat, sdata_parcial, 
                               variables=vars2heatmap2, 
                               opt, name,
                               taxalist=taxalist,
                               logscale=FALSE, 
                               w=9, h=2,
                               trim_values=FALSE,
                               italics_rownames=TRUE)
      name <- paste0("/heatmap_allRegions_", unique(sdata_parcial$Stress), "_p01_trim02.pdf")
      make_full_region_heatmap(mat, sdata_parcial, 
                               variables=vars2heatmap2, 
                               opt, name,
                               taxalist=taxalist,
                               logscale=FALSE, 
                               w=9, h=2,
                               trim_values=TRUE,
                               italics_rownames=TRUE, trimquantile=0.02)
      name <- paste0("/heatmap_allRegions_", unique(sdata_parcial$Stress), "noABS_p01.pdf")
      sdata_parcial3 <- sdata_parcial %>% filter(Treatment != "ABS")
      mat <- mat[, sdata_parcial3$sampleID %>% as.character]
      make_full_region_heatmap(mat, 
                               sdata_parcial3, 
                               variables=vars2heatmap2, 
                               opt, name,
                               taxalist=taxalist,
                               logscale=FALSE, 
                               w=9, h=2,
                               trim_values=FALSE,
                               italics_rownames=TRUE)
      name <- paste0("/heatmap_allRegions_", unique(sdata_parcial$Stress), "noABS_p01_trim02.pdf")
      make_full_region_heatmap(mat, 
                               sdata_parcial3, 
                               variables=vars2heatmap2, 
                               opt, name,
                               taxalist=taxalist,
                               logscale=FALSE, 
                               w=9, h=2,
                               trim_values=TRUE,
                               italics_rownames=TRUE, trimquantile=0.02)
      
      })
      
    opt <- restaurar(opt)
    opt$out <- paste0(opt$out, deseqname, phname, "/CompareLFC/")
    if(! dir.exists(opt$out)) dir.create(opt$out)
    
    deseq_longdf <- deseq_df %>% dplyr::mutate(
      taxon = map(contrast_list, \(x) x$resdf$taxon),
      LFCshrink = map(contrast_list, \(x) x$resdf$log2FoldChangeShrink),
      padj = map(contrast_list, \(x) x$resdf$padj),
    ) %>% select(-contrast_list, -diff_p01, -diff_p05, -samples_used) %>% 
      unnest(cols = c(taxon, LFCshrink, padj))
    
    deseq2comp <- deseq_longdf %>% 
      filter(group1 == gsub(":", ".", base_cond)) %>% 
      filter(Treatment_group2 != "ABS")
    
    no_ab <- deseq2comp %>% dplyr::filter(Treatment_group2 == "NO ABS") %>% 
      dplyr::rename(NO_ABS_LFC = LFCshrink, NO_ABS_padj = padj) %>% 
      dplyr::select(-group2, -Treatment_group2, -cond_name)
    si_ab <-  deseq2comp %>% dplyr::filter(Treatment_group2 != "NO ABS")
    
    prep2plot <- merge(si_ab, no_ab)
    # test <- merge(deseq2comp, no_ab) %>% dplyr::filter(Treatment_group2 == "NO ABS")
    # plot(test$LFCshrink, test$NO_ABS_LFC)
    
    gx<-ggplot(prep2plot, aes(x=NO_ABS_LFC, y=LFCshrink, col=Treatment_group2)) +
      facet_grid(Stress_group2 ~ Region_sequenced_group2) +
      geom_smooth(method="lm", alpha=0) +
      geom_point(alpha=0.5, size=0.5)+
      ggpmisc::stat_poly_eq(use_label(c("R2", "p","eq")), #c("eq", "R2", "f", "p", "n")
                            method="lm", small.p=T, small.r=F, 
                            label.x = c("left","left", "left", "left"),
                            label.y=c(0.99, 0.93, 0.87, 0.81))+
      scale_color_npg()+
      xlab("LFC in the NO ABS condition") +
      ylab("LFC in the ABS + region condition") +
      theme(axis.text.x = element_text(size = 14))+
      theme(strip.text.x = element_text(size = 14))+
      theme(axis.title.y = element_text(size = 14))+
      theme(axis.title.x = element_text(size = 14))+
      theme(axis.text.y = element_text( size = 14)) +
      theme_classic()
    
    ggsave(paste0(opt$out, "/plot_vs_noABS_1.pdf"), gx, width = 8, height = 4)
    #This one works fine:
    gx<-ggplot(prep2plot, aes(x=NO_ABS_LFC, y=LFCshrink, col=Treatment_group2)) +
      facet_grid(Stress_group2+Treatment_group2 ~ Region_sequenced_group2) +
      geom_smooth(method="lm", alpha=0.5) +
      geom_point(alpha=0.5)+
      ggpmisc::stat_poly_eq(use_label(c("R2", "p","eq")), #c("eq", "R2", "f", "p", "n")
                            method="lm", small.p=T, small.r=F, 
                            label.x = c("left","left", "left", "left"),
                            label.y=c(0.99, 0.93, 0.87, 0.81))+
      scale_color_npg()+
      xlab("LFC in the NO ABS condition") +
      ylab("LFC in the ABS + region condition") +
      theme(axis.text.x = element_text(size = 14))+
      theme(strip.text.x = element_text(size = 14))+
      theme(axis.title.y = element_text(size = 14))+
      theme(axis.title.x = element_text(size = 14))+
      theme(axis.text.y = element_text( size = 14)) +
      theme_classic()
    ggsave(paste0(opt$out, "/plot_vs_noABS_2.pdf"), gx, width = 8, height = 7)
    
    gx<-ggplot(prep2plot, aes(x=NO_ABS_LFC, y=LFCshrink, col=Region_sequenced_group2, fill=Region_sequenced_group2)) +
      facet_grid(Stress_group2 ~ Treatment_group2 ) +
      geom_smooth(method="lm", alpha=0) +
      geom_point(alpha=0.5, size=0.5)+
      ggpmisc::stat_poly_eq(use_label(c("R2", "p","eq")), #c("eq", "R2", "f", "p", "n")
                            method="lm", small.p=T, small.r=F, 
                            label.x = c("left","left", "left", "left"),
                            label.y=c(0.99, 0.93, 0.87, 0.81))+
      scale_color_npg()+
      scale_fill_npg()+
      xlab("LFC in the NO ABS condition") +
      ylab("LFC in the ABS + region condition") +
      theme(axis.text.x = element_text(size = 14))+
      theme(strip.text.x = element_text(size = 14))+
      theme(axis.title.y = element_text(size = 14))+
      theme(axis.title.x = element_text(size = 14))+
      theme(axis.text.y = element_text( size = 14)) +
      theme_classic()
    ggsave(paste0(opt$out, "/plot_vs_noABS_3.pdf"), gx, width = 8, height = 4)
    
    abnlong_mean <- df2plot %>% gather("sampleID","abundance", 2:ncol(.)) %>% 
      inner_join(sdata) %>% 
      group_by(Treatment, Region_sequenced, Stress, gene) %>% 
      summarise(mean_abundance = mean(abundance))
    abn_noabs <- abnlong_mean %>% filter(Treatment == "NO ABS")%>% 
      dplyr::rename(Treatment_NoABS = Treatment, 
             mean_abundance_NoABS = mean_abundance)
    abn_other <- abnlong_mean %>% dplyr::filter(Treatment != "NO ABS") %>% 
      left_join(abn_noabs)
     
    gx<-ggplot(abn_other, aes(x=mean_abundance_NoABS, 
                              y=mean_abundance, 
                              col=Treatment, 
                              fill=Treatment, 
                              group = Treatment)) +
      facet_grid(Stress ~ Region_sequenced ) +
      geom_abline(intercept = 0, slope = 1, linetype=2, col="gray") +
      geom_smooth(method="lm", alpha=0) +
      geom_point(alpha=0.5, size=0.5)+
      ggpmisc::stat_poly_eq(use_label(c("R2", "p")), #c("eq", "R2", "f", "p", "n")
                   method="lm", small.p=T, small.r=F, 
                   label.x = c("left","left", "right", "right"),
                   label.y=c(0.99, 0.93, 0.99, 0.93))+
      scale_color_npg()+
      scale_fill_npg()+
      xlab("Abundance in the NO ABS condition") +
      ylab("Abundance in the ABS + region condition") +
      theme(axis.text.x = element_text(size = 14))+
      theme(strip.text.x = element_text(size = 16))+
      theme(strip.text.y = element_text(size = 16))+
      theme(axis.title.y = element_text(size = 16))+
      theme(axis.title.x = element_text(size = 16))+
      theme(axis.text.y = element_text( size = 14)) +
      theme_classic()
    ggsave(paste0(opt$out, "/plot_Abundance_vs_noABS_3.pdf"), gx, width = 11.5, height = 6)
    gx<-ggplot(abn_other, aes(x=mean_abundance_NoABS, 
                              y=mean_abundance, 
                              col=Treatment, 
                              fill=Treatment, 
                              group = Treatment)) +
      facet_grid(Stress ~ Region_sequenced ) +
      geom_abline(intercept = 0, slope = 1, linetype=2, col="gray") +
      geom_smooth(method="lm", alpha=0) +
      geom_point(alpha=0.5, size=0.5)+
      ggpmisc::stat_poly_eq(use_label(c("R2", "p","eq")), #c("eq", "R2", "f", "p", "n")
                            method="lm", small.p=T, small.r=F, 
                            label.x = c("left","left", "left", "left"),
                            label.y=c(0.99, 0.93, 0.87, 0.81))+
      scale_color_npg()+
      scale_fill_npg()+
      xlab("Abundance in the NO ABS condition") +
      ylab("Abundance in the ABS + region condition") +
      theme(axis.text.x = element_text(size = 14))+
      theme(strip.text.x = element_text(size = 16))+
      theme(strip.text.y = element_text(size = 16))+
      theme(axis.title.y = element_text(size = 16))+
      theme(axis.title.x = element_text(size = 16))+
      theme(axis.text.y = element_text( size = 14)) +
      theme_classic()
    ggsave(paste0(opt$out, "/plot_Abundance_vs_noABS_3b.pdf"), gx, width = 11.5, height = 6)
    #opt <- restaurar(opt)
    
    ## Models 
    allmodels <- prep2plot %>% group_by(Stress_group2, Treatment_group2, Region_sequenced_group2) %>% 
      group_modify(~ broom::tidy(lm(NO_ABS_LFC ~ LFCshrink, .x))) %>% 
      dplyr::mutate(term=gsub("[\\(\\)]", "", term, perl=TRUE)) %>% 
      select(-std.error, -statistic) %>% 
      gather("param", "value", estimate, p.value) %>% 
      unite("temp", param, term, sep="_") %>% 
      spread(temp, value)
    
    allmodels2 <- prep2plot %>% group_by(Stress_group2, Treatment_group2, Region_sequenced_group2) %>% 
      group_modify(~ broom::glance(lm(NO_ABS_LFC ~ LFCshrink, .x))) %>% 
      left_join(allmodels)
    write_tsv(allmodels2, file = paste0(opt$out, "models_LFC_against_baseline.tsv"))
    
    gbeta <- ggplot(allmodels2, aes(x = Region_sequenced_group2, y = adj.r.squared, 
                           col=Treatment_group2, 
                           fill=Treatment_group2, 
                           group=Treatment_group2))+
      geom_col(position="dodge")+
      facet_grid( ~ Stress_group2 ) +
      theme(axis.text.x = element_text(size = 14))+
      theme(strip.text.x = element_text(size = 16))+
      theme(strip.text.y = element_text(size = 16))+
      theme(axis.title.y = element_text(size = 16))+
      theme(axis.title.x = element_text(size = 16))+
      theme(axis.text.y = element_text( size = 14)) +
      ylab("Adjusted R2") + 
      xlab("Region Sequenced") +
      scale_fill_aaas()+
      scale_color_aaas() +
      theme_classic()
    ggsave(paste0(opt$out, "/plot_ModelsLFC_AdjRsquared.pdf"), gbeta, width = 8, height = 2.7) 
    
    gbeta <- ggplot(allmodels2, aes(x = Region_sequenced_group2, y = estimate_LFCshrink, 
                                    col=Treatment_group2, 
                                    fill=Treatment_group2, 
                                    group=Treatment_group2))+
      geom_col(position="dodge")+
      facet_grid( ~ Stress_group2 ) +
      theme(axis.text.x = element_text(size = 14))+
      theme(strip.text.x = element_text(size = 16))+
      theme(strip.text.y = element_text(size = 16))+
      theme(axis.title.y = element_text(size = 16))+
      theme(axis.title.x = element_text(size = 16))+
      theme(axis.text.y = element_text( size = 14)) +
      ylab("Slope of LFC shrink") + 
      xlab("Region Sequenced") +
      scale_fill_aaas()+
      scale_color_aaas() +
      theme_classic()
    ggsave(paste0(opt$out, "/plot_ModelsLFC_Slope.pdf"), gbeta, width = 8, height = 2.7) 
    
    ## BARPLOTS
    deseq2comp2 <- deseq_longdf %>% 
      dplyr::mutate(signif_05 = ifelse(is.na(padj) | padj>0.05, "NS", ifelse(LFCshrink<0, "Down", "UP")), 
                    signif_01 = ifelse(is.na(padj) | padj>0.01, "NS", ifelse(LFCshrink<0, "Down", "UP"))
      )
    
    aux <- deseq2comp2 %>% 
      dplyr::filter(Treatment_group1 == Treatment_group2) %>% 
      dplyr::filter( Treatment_group1 != "ABS") %>% 
      dplyr::filter(Stress_group1 == Stress_group2) %>% 
      dplyr::filter(Region_sequenced_group1 != Region_sequenced_group2) %>% 
      dplyr::mutate(taxon = gsub("_", " ", taxon) %>% gsub("[\\[\\]]", "", ., perl=TRUE))
    
    taxa2plot <- aux %>% 
      dplyr::filter(!is.na(padj)) %>% 
      group_by(taxon) %>% 
      summarise(minp = min(padj)) %>% 
      filter(minp < 0.01) %>% 
      pull(taxon)
    
    aux <- aux %>% dplyr::filter(taxon %in% taxa2plot) %>% 
      dplyr::mutate(comp_region = paste(Region_sequenced_group2, "vs", Region_sequenced_group1)) %>% 
      dplyr::mutate(comp_region = gsub("REG", "R", comp_region))
    
    tab2order <- aux %>% 
      dplyr::filter(taxon %in% taxa2plot) %>% 
      dplyr::filter(Stress_group1 == "Control") %>% 
      dplyr::filter(Treatment_group1 == "NO ABS") %>% 
      dplyr::filter(Region_sequenced_group2 == "REG3")%>% 
      dplyr::filter(Region_sequenced_group1 == "REG1") %>% 
      arrange(LFCshrink)
    
    aux <-  aux %>% 
      group_by(Stress_group1) %>% 
      group_split()
    
    auxplots <- aux %>% 
      map(\(aux2){
        
        stressmain = unique(aux2$Stress_group1)
        aux2 <- aux2 %>% dplyr::mutate(taxon = factor(taxon, levels=tab2order$taxon))
        
        gg <- ggplot(aux2, aes(x=taxon, y=LFCshrink, 
                               col=signif_05, fill=signif_05, group=signif_05))+
          facet_grid(~ Treatment_group1 + comp_region) +
          geom_col() +
          #thin_barplot_lines +
          #geom_hline(yintercept = 0, col="gray", linetype=2)+
          #geom_vline(xintercept = 0, col="gray", linetype=2)+
          scale_color_manual(values = c("steelblue2", "gray", "tomato")) +
          scale_fill_manual(values = c("steelblue2", "gray", "tomato")) +
          coord_flip() +
          theme_minimal() +
          theme(axis.text.x = element_text(size = 12))+
          theme(strip.text.x = element_text(size = 10))+
          theme(strip.text.y = element_text(size = 8))+
          theme(axis.title.y = element_text(size = 12))+
          theme(axis.title.x = element_text(size = 12))+
          theme(axis.text.y = element_text( size = 12, face = "italic")) +
          ggtitle(stressmain) +
          theme(plot.title = element_text(hjust = 0.5, vjust=0.5))
        h = 1.31 + 0.16*length(taxa2plot)
        ggsave(paste0(opt$out, paste0(stressmain, "_CompareRegionSeq_barplot_LFCShrink_p01_col05.pdf")), gg, width = 16.9, height = h)
          
        return(gg)
      })
    
    auxpoints <- aux %>% 
      map(\(aux2){
        
        stressmain = unique(aux2$Stress_group1)
        aux2 <- aux2 %>% dplyr::mutate(taxon = factor(taxon, levels=rev(tab2order$taxon)),
                                       comp_region = factor(comp_region, 
                                                            levels = c("R3 vs R1", "R3 vs R2", "R2 vs R1")))
        
        gg <- ggplot(aux2, aes(x=taxon, y=LFCshrink, 
                               col=Treatment_group1, 
                               fill=Treatment_group1, 
                               group=Treatment_group1))+
          facet_grid(comp_region~ .) +
          geom_hline(yintercept = 0, linetype=2, color="gray") +
          geom_point() +
          geom_line() +
          #thin_barplot_lines +
          #geom_hline(yintercept = 0, col="gray", linetype=2)+
          #geom_vline(xintercept = 0, col="gray", linetype=2)+
          #scale_color_manual(values = c("steelblue2", "gray", "tomato")) +
          #scale_fill_manual(values = c("steelblue2", "gray", "tomato")) +
          scale_color_npg() + 
          scale_fill_npg() +
          #coord_flip() +
          theme_classic() +
          theme(axis.text.x = element_text(size = 12, angle=45, vjust=1, hjust=1, face="italic"))+
          theme(strip.text.x = element_text(size = 14))+
          theme(strip.text.y = element_text(size = 14))+
          theme(axis.title.y = element_text(size = 12))+
          theme(axis.title.x = element_text(size = 12))+
          theme(axis.text.y = element_text( size = 12)) +
          ggtitle(stressmain) +
          theme(plot.title = element_text(hjust = 0.5, vjust=0.5)) +
          theme(plot.margin = unit(c(0,1,0,3), "cm"))
        w = 3.81 + 0.191*length(taxa2plot)
        ggsave(paste0(opt$out, paste0(stressmain, "_CompareRegionSeq_points_LFCShrink_p01.pdf")), gg, width = w, height = 6.1)
        
        return(gg)
      })
    pdf(paste0(opt$out, paste0("BothStress", "_CompareRegionSeq_points_LFCShrink_p01.pdf")), 
        height = 12.3, width = 3.81 + 0.191*length(taxa2plot))
    print(cowplot::plot_grid(plotlist = auxpoints, ncol=1))
    dev.off()
    
    ## Repeat line plot but now with ABS group
    aux <- deseq2comp2 %>% 
      dplyr::filter(Treatment_group1 == Treatment_group2) %>% 
      #dplyr::filter( Treatment_group1 != "ABS") %>% 
      dplyr::filter(Stress_group1 == Stress_group2) %>% 
      dplyr::filter(Region_sequenced_group1 != Region_sequenced_group2) %>% 
      dplyr::mutate(taxon = gsub("_", " ", taxon) %>% gsub("[\\[\\]]", "", ., perl=TRUE))
    
    taxa2plot <- aux %>% 
      dplyr::filter(!is.na(padj)) %>% 
      dplyr::filter(Treatment_group1 != "ABS") %>% 
      group_by(taxon) %>% 
      summarise(minp = min(padj)) %>% 
      filter(minp < 0.01) %>% 
      pull(taxon)
    
    aux <- aux %>% dplyr::filter(taxon %in% taxa2plot) %>% 
      dplyr::mutate(comp_region = paste(Region_sequenced_group2, "vs", Region_sequenced_group1)) %>% 
      dplyr::mutate(comp_region = gsub("REG", "R", comp_region))
    
    aux2cor <- aux %>% unite("full_cond", Stress_group1, Treatment_group1, comp_region, sep=":", remove = F) %>% 
      dplyr::select(full_cond, taxon, LFCshrink, Stress_group1, Treatment_group1, comp_region) 
    
    annrow <- aux2cor %>% dplyr::select(full_cond,Stress_group1, Treatment_group1, comp_region ) %>% 
      distinct() %>% 
      column_to_rownames("full_cond")
    
    aux2cor <- aux2cor %>% 
      dplyr::select(full_cond, taxon, LFCshrink) %>% 
      spread(full_cond, LFCshrink) %>% 
      column_to_rownames("taxon") %>% 
      as.matrix()
    
    numcols <- max(sapply(annrow, \(x)length(unique(x))))
    cc <- ggsci::pal_d3(palette = "category10")(numcols)
    user.colfn=colorRampPalette(cc)
    newcc <- user.colfn(numcols) # in case there are too many colors
    
    color_list_cols <- lapply(annrow, \(x) {y <-newcc[1:length(unique(x))]; names(y)<- unique(x); y})
    
    mat <- cor(aux2cor)  
    
    pheatmap(mat, annotation_col = annrow, annotation_row=annrow, annotation_colors = color_list_cols, 
             filename = paste0(opt$out, paste0("BothStress", "_CompareRegionSeq_heatmap_LFCShrink_p01_withABS.pdf")),
             width = 12, height = 8
    )
    annrowx <- annrow %>% rownames_to_column("cond") %>% group_by(Stress_group1) %>% group_split() %>% map(\(annrowx){
      stressmain <- annrowx$Stress_group1 %>% unlist %>% unique
      
      matx <- mat[annrowx$cond, annrowx$cond]
      rownames(matx) <- gsub(paste0(stressmain, ":"), "", rownames(matx))
      colnames(matx) <- gsub(paste0(stressmain, ":"), "", colnames(matx))
      annrowx <- annrowx %>% dplyr::select(-Stress_group1) %>% 
        dplyr::mutate(cond = gsub(paste0(stressmain, ":"), "", cond)) %>% 
        column_to_rownames("cond")
      pheatmap(matx, annotation_col = annrowx, annotation_row=annrowx, annotation_colors = color_list_cols[names(annrowx)], 
               filename = paste0(opt$out, paste0(stressmain, "_CompareRegionSeq_heatmap_LFCShrink_p01_withABS.pdf")),
               width = 9, height = 6
      )
      
    })
    
    tab2order <- aux %>% 
      dplyr::filter(taxon %in% taxa2plot) %>% 
      dplyr::filter(Stress_group1 == "Control") %>% 
      dplyr::filter(Treatment_group1 == "NO ABS") %>% 
      dplyr::filter(Region_sequenced_group2 == "REG3")%>% 
      dplyr::filter(Region_sequenced_group1 == "REG1") %>% 
      arrange(LFCshrink)
    
    aux <-  aux %>% 
      group_by(Stress_group1) %>% 
      group_split()
    auxpoints <- aux %>% 
      map(\(aux2){
        
        stressmain = unique(aux2$Stress_group1)
        aux2 <- aux2 %>% dplyr::mutate(taxon = factor(taxon, levels=rev(tab2order$taxon)),
                                       comp_region = factor(comp_region, 
                                                            levels = c("R3 vs R1", "R3 vs R2", "R2 vs R1")),
                                       Treatment_group1 = factor(Treatment_group1, 
                                                                 levels=c("NO ABS", "REG1", "REG2", "REG3", "ABS")))
        
        gg <- ggplot(aux2, aes(x=taxon, y=LFCshrink, 
                               col=Treatment_group1, 
                               fill=Treatment_group1, 
                               group=Treatment_group1))+
          facet_grid(comp_region~ .) +
          geom_hline(yintercept = 0, linetype=2, color="gray") +
          geom_point() +
          geom_line() +
          #thin_barplot_lines +
          #geom_hline(yintercept = 0, col="gray", linetype=2)+
          #geom_vline(xintercept = 0, col="gray", linetype=2)+
          #scale_color_manual(values = c("steelblue2", "gray", "tomato")) +
          #scale_fill_manual(values = c("steelblue2", "gray", "tomato")) +
          scale_color_npg() + 
          scale_fill_npg() +
          #coord_flip() +
          theme_classic() +
          theme(axis.text.x = element_text(size = 12, angle=45, vjust=1, hjust=1, face="italic"))+
          theme(strip.text.x = element_text(size = 14))+
          theme(strip.text.y = element_text(size = 14))+
          theme(axis.title.y = element_text(size = 12))+
          theme(axis.title.x = element_text(size = 12))+
          theme(axis.text.y = element_text( size = 12)) +
          ggtitle(stressmain) +
          theme(plot.title = element_text(hjust = 0.5, vjust=0.5)) +
          theme(plot.margin = unit(c(0,1,0,3), "cm"))
        w = 3.81 + 0.191*length(taxa2plot)
        ggsave(paste0(opt$out, paste0(stressmain, "_CompareRegionSeq_points_LFCShrink_p01_withABS.pdf")), gg, width = w, height = 6.1)
        
        return(gg)
      })
    pdf(paste0(opt$out, paste0("BothStress", "_CompareRegionSeq_points_LFCShrink_p01_withABS.pdf")), 
        height = 12.3, width = 3.81 + 0.191*length(taxa2plot))
    print(cowplot::plot_grid(plotlist = auxpoints, ncol=1))
    dev.off()
    
    
    ### By Stress
    aux <- deseq2comp2 %>% 
      dplyr::filter(Treatment_group1 == Treatment_group2) %>% 
      dplyr::filter( Treatment_group1 != "ABS") %>% 
      dplyr::filter(Stress_group1 != Stress_group2) %>% 
      dplyr::filter(Region_sequenced_group1 == Region_sequenced_group2) %>% 
      dplyr::mutate(taxon = gsub("_", " ", taxon) %>% gsub("[\\[\\]]", "", ., perl=TRUE))
    
    aux <-  aux %>% 
      group_by(Region_sequenced_group1) %>% 
      group_split()
    
    auxplots <- aux %>% 
      map(\(aux2){
        
        regmain = unique(aux2$Region_sequenced_group1)
        
        taxa2plot <- aux2 %>% 
          dplyr::filter(!is.na(padj)) %>% 
          group_by(taxon) %>% 
          summarise(minp = min(padj)) %>% 
          dplyr::filter(minp < 0.01) %>% 
          pull(taxon)
        
        aux2 <- aux2 %>% dplyr::filter(taxon %in% taxa2plot)
        
        tab2order <- aux2 %>% 
          dplyr::filter(taxon %in% taxa2plot) %>% 
          dplyr::filter(Treatment_group1 == "NO ABS") %>% 
          arrange(LFCshrink)
        
        
        aux2 <- aux2 %>% dplyr::mutate(taxon = factor(taxon, levels=tab2order$taxon))
        
        gg <- ggplot(aux2, aes(x=taxon, y=LFCshrink, 
                               col=signif_05, fill=signif_05, group=signif_05))+
          facet_grid(~ Treatment_group1) +
          geom_col() +
          #thin_barplot_lines +
          #geom_hline(yintercept = 0, col="gray", linetype=2)+
          #geom_vline(xintercept = 0, col="gray", linetype=2)+
          scale_color_manual(values = c("steelblue2", "gray", "tomato")) +
          scale_fill_manual(values = c("steelblue2", "gray", "tomato")) +
          coord_flip() +
          theme_minimal() +
          theme(axis.text.x = element_text(size = 12))+
          theme(strip.text.x = element_text(size = 10))+
          theme(strip.text.y = element_text(size = 8))+
          theme(axis.title.y = element_text(size = 12))+
          theme(axis.title.x = element_text(size = 12))+
          theme(axis.text.y = element_text( size = 12, face = "italic")) +
          ggtitle(regmain) +
          theme(plot.title = element_text(hjust = 0.5, vjust=0.5))
        h = 1.31 + 0.16*length(taxa2plot)
        ggsave(paste0(opt$out, paste0(regmain, "_CompareStress_barplot_LFCShrink_p01_col05.pdf")), gg, width = 10, height = h)
        
        return(gg)
      })
    
    auxpoints <- aux %>% 
      map(\(aux2){
        
        regmain = unique(aux2$Region_sequenced_group1)
        
        taxa2plot <- aux2 %>% 
          dplyr::filter(!is.na(padj)) %>% 
          group_by(taxon) %>% 
          summarise(minp = min(padj)) %>% 
          dplyr::filter(minp < 0.01) %>% 
          pull(taxon)
        
        aux2 <- aux2 %>% dplyr::filter(taxon %in% taxa2plot)
        
        tab2order <- aux2 %>% 
          dplyr::filter(taxon %in% taxa2plot) %>% 
          dplyr::filter(Treatment_group1 == "NO ABS") %>% 
          arrange(LFCshrink)
        
        
        aux2 <- aux2 %>% dplyr::mutate(taxon = factor(taxon, levels=rev(tab2order$taxon)))
        
        gg <- ggplot(aux2, aes(x=taxon, y=LFCshrink, 
                               col=Treatment_group1, fill=Treatment_group1, group=Treatment_group1))+
          #facet_grid(~ Treatment_group1) +
          #geom_col() +
          geom_hline(yintercept = 0, linetype=2, color="gray") +
          geom_point() +
          geom_line() +
          #thin_barplot_lines +
          #geom_hline(yintercept = 0, col="gray", linetype=2)+
          #geom_vline(xintercept = 0, col="gray", linetype=2)+
          scale_color_npg() + 
          scale_fill_npg() +
          #coord_flip() +
          theme_classic() +
          theme(axis.text.x = element_text(size = 12, angle=45, vjust=1, hjust=1, face="italic"))+
          theme(strip.text.x = element_text(size = 14))+
          theme(strip.text.y = element_text(size = 14))+
          theme(axis.title.y = element_text(size = 12))+
          theme(axis.title.x = element_text(size = 12))+
          theme(axis.text.y = element_text( size = 12)) +
          ggtitle(regmain) +
          theme(plot.title = element_text(hjust = 0.5, vjust=0.5))+
          theme(plot.margin = unit(c(0,1,0,3), "cm"))
        w = 3.81 + 0.191*length(taxa2plot)
        ggsave(paste0(opt$out, paste0(regmain, "_CompareStress_points_LFCShrink_p01.pdf")), gg, width = w, height = 5)
        
        return(gg)
      })
    
    pdf(paste0(opt$out, paste0("AllRegions", "_CompareStress_points_LFCShrink_p01.pdf")), 
        height = 5*length(auxpoints), width = 3.81 + 0.191*length(taxa2plot))
    print(cowplot::plot_grid(plotlist = auxpoints, ncol=1))
    dev.off()
    
    ## Again Stress, same thing but with ABS 
    
    aux <- deseq2comp2 %>% 
      dplyr::filter(Treatment_group1 == Treatment_group2) %>% 
      #dplyr::filter( Treatment_group1 != "ABS") %>% 
      dplyr::filter(Stress_group1 != Stress_group2) %>% 
      dplyr::filter(Region_sequenced_group1 == Region_sequenced_group2) %>% 
      dplyr::mutate(taxon = gsub("_", " ", taxon) %>% gsub("[\\[\\]]", "", ., perl=TRUE))
    
    aux <-  aux %>% 
      group_by(Region_sequenced_group1) %>% 
      group_split()
    
    auxpoints <- aux %>% 
      map(\(aux2){
        
        regmain = unique(aux2$Region_sequenced_group1)
        
        taxa2plot <- aux2 %>% 
          dplyr::filter(!is.na(padj)) %>% 
          dplyr::filter(Treatment_group1 != "ABS") %>% 
          group_by(taxon) %>% 
          summarise(minp = min(padj)) %>% 
          dplyr::filter(minp < 0.01) %>% 
          pull(taxon)
        
        aux2 <- aux2 %>% dplyr::filter(taxon %in% taxa2plot)
        
        tab2order <- aux2 %>% 
          dplyr::filter(taxon %in% taxa2plot) %>% 
          dplyr::filter(Treatment_group1 == "NO ABS") %>% 
          arrange(LFCshrink)
        
        
        aux2 <- aux2 %>% dplyr::mutate(taxon = factor(taxon, levels=rev(tab2order$taxon)),
                                       Treatment_group1 = factor(Treatment_group1, 
                                                                 levels=c("NO ABS", "REG1", "REG2", "REG3", "ABS")))
        
        gg <- ggplot(aux2, aes(x=taxon, y=LFCshrink, 
                               col=Treatment_group1, fill=Treatment_group1, group=Treatment_group1))+
          #facet_grid(~ Treatment_group1) +
          #geom_col() +
          geom_hline(yintercept = 0, linetype=2, color="gray") +
          geom_point() +
          geom_line() +
          #thin_barplot_lines +
          #geom_hline(yintercept = 0, col="gray", linetype=2)+
          #geom_vline(xintercept = 0, col="gray", linetype=2)+
          scale_color_npg() + 
          scale_fill_npg() +
          #coord_flip() +
          theme_classic() +
          theme(axis.text.x = element_text(size = 12, angle=45, vjust=1, hjust=1, face="italic"))+
          theme(strip.text.x = element_text(size = 14))+
          theme(strip.text.y = element_text(size = 14))+
          theme(axis.title.y = element_text(size = 12))+
          theme(axis.title.x = element_text(size = 12))+
          theme(axis.text.y = element_text( size = 12)) +
          ggtitle(regmain) +
          theme(plot.title = element_text(hjust = 0.5, vjust=0.5))+
          theme(plot.margin = unit(c(0,1,0,3), "cm"))
        w = 3.81 + 0.191*length(taxa2plot)
        ggsave(paste0(opt$out, paste0(regmain, "_CompareStress_points_LFCShrink_p01_withABS.pdf")), gg, width = w, height = 5)
        
        return(gg)
      })
    
    pdf(paste0(opt$out, paste0("AllRegions", "_CompareStress_points_LFCShrink_p01_withABS.pdf")), 
        height = 5*length(auxpoints), width = 3.81 + 0.191*length(taxa2plot))
    print(cowplot::plot_grid(plotlist = auxpoints, ncol=1))
    dev.off()
    
    ## Models from abundances

    allmodels <- abn_other %>% group_by(Stress, Treatment, Region_sequenced) %>% 
      group_modify(~ broom::tidy(lm(mean_abundance_NoABS ~ mean_abundance, .x))) %>% 
      dplyr::mutate(term=gsub("[\\(\\)]", "", term, perl=TRUE)) %>% 
      select(-std.error, -statistic) %>% 
      gather("param", "value", estimate, p.value) %>% 
      unite("temp", param, term, sep="_") %>% 
      spread(temp, value)
    
    allmodels2 <- abn_other %>% group_by(Stress, Treatment, Region_sequenced) %>% 
      group_modify(~ broom::glance(lm(mean_abundance_NoABS ~ mean_abundance, .x))) %>% 
      left_join(allmodels)
    write_tsv(allmodels2, file = paste0(opt$out, "models_VSTAbundances_sameRegionNoAbs.tsv"))
    
    gbeta <- ggplot(allmodels2, aes(x = Region_sequenced, y = adj.r.squared, 
                                    col=Treatment, 
                                    fill=Treatment, 
                                    group=Treatment))+
      geom_col(position="dodge")+
      facet_grid( ~ Stress ) +
      theme(axis.text.x = element_text(size = 14))+
      theme(strip.text.x = element_text(size = 16))+
      theme(strip.text.y = element_text(size = 16))+
      theme(axis.title.y = element_text(size = 16))+
      theme(axis.title.x = element_text(size = 16))+
      theme(axis.text.y = element_text( size = 14)) +
      ylab("Adjusted R2") + 
      xlab("Region Sequenced") +
      scale_fill_npg()+
      scale_color_npg() +
      theme_classic()
    ggsave(paste0(opt$out, "/plot_ModelsVSTAbundances_AdjRsquared.pdf"), gbeta, width = 8, height = 2.7) 
    
    gbeta <- ggplot(allmodels2, aes(x = Region_sequenced, y = estimate_mean_abundance, 
                                    col=Treatment, 
                                    fill=Treatment, 
                                    group=Treatment))+
      geom_col(position="dodge")+
      facet_grid( ~ Stress ) +
      theme(axis.text.x = element_text(size = 14))+
      theme(strip.text.x = element_text(size = 16))+
      theme(strip.text.y = element_text(size = 16))+
      theme(axis.title.y = element_text(size = 16))+
      theme(axis.title.x = element_text(size = 16))+
      theme(axis.text.y = element_text( size = 14)) +
      ylab("Slope of mean abundance predictor") + 
      xlab("Region Sequenced") +
      scale_fill_npg()+
      scale_color_npg() +
      theme_classic()
    ggsave(paste0(opt$out, "/plot_ModelsVSTAbundances_Slope.pdf"), gbeta, width = 8, height = 2.7) 
    opt <- restaurar(opt)
      
}
