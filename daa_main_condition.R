
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
    opt$out <- paste0(opt$out, deseqname, phname, "/HeatmapsLFC/")
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
      geom_smooth(method="lm") +
      geom_point(alpha=0.5)+
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
    gx<-ggplot(prep2plot, aes(x=NO_ABS_LFC, y=LFCshrink, col=Treatment_group2)) +
      facet_grid(Stress_group2+Treatment_group2 ~ Region_sequenced_group2) +
      geom_smooth(method="lm") +
      geom_point(alpha=0.5)+
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
      geom_smooth(method="lm") +
      geom_point(alpha=0.5)+
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
    ggsave(paste0(opt$out, "/plot_vs_noABS_3.pdf"), gx, width = 8, height = 7)
     
    opt <- restaurar(opt)
}
