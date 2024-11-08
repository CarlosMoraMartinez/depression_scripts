# Alpha 4 each
library(rstatix)

## Cualitativas

signif_levels=c("***"=0.001, "**"=0.01, "*"=0.05, "ns"=1.1)

alpha_indices <- c("Observed", "Chao1", "Shannon", "InvSimpson")
vars2test <- c("Treatment", "Region_sequenced", "Stress", "Treatment_full", "flowcell", "Treatment_region") #"Category_T0"

quant_vars <- c("nreads_filt", "num_reads")
vars2log <- c( "nreads_filt", "num_reads")

quant_vars_ext <- c(quant_vars, paste(vars2log, "_log", sep=""))
interestvar <- "Treatment"
extravars <- c(quant_vars, vars2test)
extravars <- extravars[extravars != interestvar]

outdir <- paste0(opt$out, "/AlphaDiversity/")
if(!dir.exists(outdir)) dir.create(outdir)

extravars2 <- extravars


phseq_to_use <- names(all_phyloseq)
#load(allphyloseqlist_fname)

for(phname in phseq_to_use[1:2]){
  cat("Alpha diversity in ", phname, "\n")
  phobj <- all_phyloseq[[phname]]
  phobj <- updatePsWithLogs(phobj, vars2log)
  

  divtab <- calculateAlphaDiversityTable(phobj, outdir, alpha_indices, paste0(phname, "_AlphaDiv") )
  divtab <- divtab %>% mutate(
    is_transfer = ifelse(Treatment == "TRANSFER", TRUE, FALSE),
    Treatment = ifelse(Treatment == "TRANSFER", "NO ABS", Treatment)
    
  )
  ## Statistics
  
  means_bygroups = divtab %>% 
    group_by(Treatment, Region_sequenced, Stress) %>% 
    summarise(across(all_of(alpha_indices), 
                     list(mean=mean, median=median, sd=sd, var=var, max=max, min=min)))
  
 
  
  means_bygroups = divtab %>% 
    group_by(Treatment, Region_sequenced, Stress) %>% 
    summarise(across(all_of(alpha_indices), 
                     list(mean=mean, median=median, sd=sd, var=var, max=max, min=min)))
  
  write_tsv(means_bygroups, paste0(outdir, "/", phname, "_AlphaDiv_meansByGroup.tsv"))
  
  ## Differences between stress and non-stress
  test_stress <- divtab %>% 
    filter(Treatment != "TRANSFER") %>% 
    group_by(Treatment, Region_sequenced) %>% 
    group_split() %>% 
    map(\(x){
      data.frame(
        alpha_ind = alpha_indices, 
      ttest = sapply(alpha_indices, \(ai) t.test( as.formula(paste0(ai,  "~ Stress")), x)$p.value),
      wilcox = sapply(alpha_indices, \(ai) wilcox.test( as.formula(paste0(ai,  "~ Stress")), x)$p.value),
      anova=summary(aov( as.formula(paste0(ai,  "~ Stress")), x))[[1]][1, "Pr(>F)"],
      shapiro = sapply(alpha_indices, \(ai) shapiro.test( unlist(x[, ai]))$p.value),
      bartlett = sapply(alpha_indices, \(ai) bartlett.test( unlist(x[, ai]), x$Stress)$p.value)
      ) %>% gather(key="test", value="pval", ttest:bartlett) %>% 
        unite("tmp", alpha_ind, test, sep="_") %>% 
        spread(tmp, pval) %>% 
        mutate(Treatment = unique(x$Treatment),
               Region_sequenced = unique(x$Region_sequenced),
               Comparison = "Stress:Control_vs_SD") %>% 
        select(Treatment, Region_sequenced, Comparison, everything())
    }) %>% 
    bind_rows() %>% 
    mutate_if(is.numeric, list("Adj" =\(x)p.adjust(x, method="BH"))) %>% 
    select(order(colnames(.))) %>% 
    select(Treatment, Region_sequenced, Comparison, everything())
  
  write_tsv(test_stress, paste0(outdir, "/", phname, "_AlphaDiv_testStress_long.tsv"))
  
  test_stress <- map(alpha_indices, \(ai){
    divtab %>% 
      filter(Treatment != "TRANSFER") %>% 
      group_by(Treatment, Region_sequenced) %>% 
      pairwise_t_test(as.formula(paste0(ai, " ~ Stress")), 
                      paired=FALSE,
                      p.adjust.method	= "BH"
      )}) %>% bind_rows()
  
  names(test_stress)[3] <- "alpha_index"
  write_tsv(test_stress, paste0(outdir, "/", phname, "_AlphaDiv_testStress.tsv"))
  
  ## Differences by region sequenced
  ## From https://www.datanovia.com/en/lessons/repeated-measures-anova-in-r/
  test_regionseq <- map(alpha_indices, \(ai){
    divtab %>% 
    filter(Treatment != "TRANSFER") %>% 
    group_by(Treatment, Stress) %>% 
    pairwise_t_test(as.formula(paste0(ai, " ~ Region_sequenced")), 
                    paired=FALSE,
                    p.adjust.method	= "BH"
    )}) %>% bind_rows()
  
  names(test_regionseq)[3] <- "alpha_index"
  write_tsv(test_regionseq, paste0(outdir, "/", phname, "_AlphaDiv_testRegionSequenced.tsv"))
  
  test_regionseq_paired <-  map(alpha_indices, \(ai){
    divtab %>% 
      filter(Treatment != "TRANSFER") %>% 
      group_by(Treatment, Stress) %>% 
      pairwise_t_test(as.formula(paste0(ai, " ~ Region_sequenced")), 
                      paired=TRUE,
                      p.adjust.method	= "BH"
      )}) %>% bind_rows()
  
  names(test_regionseq_paired)[3] <- "alpha_index"
  write_tsv(test_regionseq_paired, paste0(outdir, "/", phname, "_AlphaDiv_testRegionSequenced_paired.tsv"))
  
  ## test treatments
  test_treatments <- map(alpha_indices, \(ai){
    divtab %>% 
      filter(Region_sequenced != "TRANSFER") %>% 
      group_by(Region_sequenced, Stress) %>% 
      pairwise_t_test(as.formula(paste0(ai, " ~ Treatment")), 
                      paired=FALSE,
                      p.adjust.method	= "BH"
      )}) %>% bind_rows()
  
  names(test_treatments)[3] <- "alpha_index"
  write_tsv(test_treatments, paste0(outdir, "/", phname, "_AlphaDiv_testTreatment.tsv"))
  
  alphaplots_compRegions <- list()
  alphaplots_compTreatments <- list()
  alphaplots_compStress <- list()
  for(ind in alpha_indices){
    auxsig <- divtab %>% 
      filter(Treatment != "TRANSFER") %>% 
      group_by(Treatment, Stress) %>% 
      pairwise_t_test(as.formula(paste0(ind, " ~ Region_sequenced")), 
                      paired=FALSE,
                      p.adjust.method	= "BH"
      ) %>% 
      add_xy_position() 
    
    alphaplots_compRegions[[ind]] <- ggplot(divtab, aes(x=Region_sequenced ,
                                                        y = !!sym(ind), 
                                                        fill=Treatment,
                                                        col=Treatment))+
      facet_grid(Stress ~ Treatment) +
      #geom_violin(alpha=0.6)+
      geom_boxplot(width=0.7, fill="white")+
      geom_point() +
      scale_color_npg() +
      ggtitle(ind) +
      theme_bw() +
      xlab("Region Sequenced") +
      #scale_fill_npg() +
      theme(axis.text.x = element_text(vjust=1,hjust=1, angle = 45)) +
      theme(plot.title = element_text(hjust = 0.5)) +
      stat_pvalue_manual(auxsig, tip.length = 0.0, hide.ns = TRUE, vjust=0.6, bracket.nudge.y=0) +
      labs(
        caption = get_pwc_label(auxsig)
      )
    
    auxsig <- divtab %>% 
      filter(Treatment != "TRANSFER") %>% 
      group_by(Treatment, Region_sequenced) %>% 
      pairwise_t_test(as.formula(paste0(ind, " ~ Stress")), 
                      paired=FALSE,
                      p.adjust.method	= "BH"
      ) %>% 
      add_xy_position() 
    
    alphaplots_compStress[[ind]] <- ggplot(divtab %>% filter(Treatment != "TRANSFER"), aes(x=Stress,
                                                                                           y = !!sym(ind), 
                                                                                           fill=Treatment,
                                                                                           col=Treatment))+
      facet_grid(Region_sequenced ~  Treatment) +
      #geom_violin(alpha=0.6)+
      geom_boxplot(width=0.7, fill="white")+
      geom_point() +
      scale_color_npg() +
      ggtitle(ind) +
      theme_bw() +
      xlab("Stress") +
      theme(axis.text.x = element_text(vjust=1,hjust=1, angle = 45)) +
      theme(plot.title = element_text(hjust = 0.5))+
      stat_pvalue_manual(auxsig, tip.length = 0.0, hide.ns = TRUE, vjust=0.6, bracket.nudge.y=0) +
      labs(
        caption = get_pwc_label(auxsig)
      )
    
    auxsig <- divtab %>% 
      filter(Treatment != "TRANSFER") %>% 
      group_by(Stress, Region_sequenced) %>% 
      pairwise_t_test(as.formula(paste0(ind, " ~ Treatment")), 
                      paired=FALSE,
                      p.adjust.method	= "BH"
      ) %>% 
      add_xy_position() 
    
    alphaplots_compTreatments[[ind]] <- ggplot(divtab %>% filter(Treatment != "TRANSFER"), 
                                               aes(x=Treatment,
                                               y = !!sym(ind), 
                                               #fill=Treatment,
                                               col=Treatment))+
      facet_grid(Stress ~ Region_sequenced) +
      #geom_violin(alpha=0.6)+
      geom_boxplot(width=0.7, fill="white")+
      geom_point() +
      scale_color_npg() +
      ggtitle(ind) +
      theme_bw() +
      xlab("Treatment") +
      #scale_fill_npg() +
      theme(axis.text.x = element_text(vjust=1,hjust=1, angle = 45)) +
      theme(plot.title = element_text(hjust = 0.5))+
      stat_pvalue_manual(auxsig, tip.length = 0.00, hide.ns = TRUE, vjust=0.6) +
      labs(
        caption = get_pwc_label(auxsig)
      )
    
    
    ggsave(filename = paste0(outdir, "/", phname, "_", ind, "_AlphaDiv_byRegion.pdf"), alphaplots_compRegions[[ind]], width = 7, height = 4)
    ggsave(filename = paste0(outdir, "/", phname, "_", ind, "_AlphaDiv_byTreatment.pdf"), alphaplots_compTreatments[[ind]], width = 7, height = 4)
    ggsave(filename = paste0(outdir, "/", phname, "_", ind, "_AlphaDiv_byStress.pdf"), alphaplots_compStress[[ind]], width = 7, height = 6)
  }
  pdf( paste0(outdir, "/", phname, "_allIndices_AlphaDiv_byRegion.pdf"), width=14, height = 8)
  print(cowplot::plot_grid(plotlist = alphaplots_compRegions))
  dev.off()
  
  pdf( paste0(outdir, "/", phname, "_allIndices_AlphaDiv_byTreatment.pdf"), width=14, height = 8)
  print(cowplot::plot_grid(plotlist = alphaplots_compTreatments))
  dev.off()
  
  pdf( paste0(outdir, "/", phname, "_allIndices_AlphaDiv_byStress.pdf"), width=12, height = 9)
  print(cowplot::plot_grid(plotlist = alphaplots_compStress))
  dev.off()
  
  
  # Typical plots
  divplots <- getAlphaDiversity(phobj, vars2test, quant_vars_ext,
                                opt,
                                indices= alpha_indices,
                                correct_pvalues = T, correct_pvalues_indices = F,
                                name = paste0(phname, "_AlphaDiv"), w = 12, h = 4)
  divplots <- getAlphaDiversity(phobj, vars2test, quant_vars_ext,
                                opt,
                                indices= alpha_indices,
                                correct_pvalues = T, correct_pvalues_indices = T,
                                name = paste0(phname, "_AlphaDivAdjInd"), w = 12, h = 4)
}



# Beta 4 each
outdir <- paste0(opt$out, "/BetaDiversity/")
if(!dir.exists(outdir)) dir.create(outdir)

dists <- c("bray", "jaccard")
vars2pcoa <- c(vars2test, quant_vars_ext)
var2shape = "Stress"
ccaplots <- list()
for(phname in phseq_to_use){
  for(method in c("PCoA", "NMDS")){
    for(dist in dists){
      name <- paste0(phname, "_", dist, "_", method, "_5dims")
      cat("Beta diversity for ", name, "\n")
      
      ccaplots[[name]] <- makeAllPCoAs(all_phyloseq[[phname]], outdir,
                                       method = method,
                                       name = name, 
                                       dist_type = dist, 
                                       dist_name = dist,
                                       vars2plot = vars2pcoa, 
                                       var2shape = var2shape,
                                       extradims = 2:5, 
                                       create_pdfs = T, w=8, h = 8)
      
      
      name <- paste0(phname, "_", dist, "_", method)
      ccaplots[[name]] <- makeAllPCoAs(all_phyloseq[[phname]], outdir,
                                       method = method,
                                       name = name,
                                       dist_type = dist,
                                       dist_name = dist,
                                       vars2plot = vars2pcoa,
                                       var2shape = var2shape,
                                       extradims = 2:3,
                                       create_pdfs = T, w=8)

      ## Make better PCoA

      pcoa.bray <- ordinate(phobj, method = method, distance = dist)
      evals <- pcoa.bray$values$Eigenvalues

      df2plot <- pcoa.bray$points %>% data.frame %>%
        rownames_to_column()

      gg <- plot_ordination(phobj, pcoa.bray,
                            color = "Region_sequenced",
                            shape = "Treatment",
                            title = name, axes=c(1, 2)) +
        #coord_fixed(sqrt(evals[2] / evals[1])) +
        #scale_color_manual(values=palette2)+
        stat_ellipse(level=0.95, linetype=2, alpha = 0.8, na.rm = TRUE) +
        geom_point(size = 1.5) +
        #geom_point(size = 1, aes(col=Stress)) +
        #geom_text_repel(aes_string(label = labelsamples)) +
        theme_bw() +
        theme(axis.text.x = element_text(size = 14))+
        theme(strip.text.x = element_text(size = 14))+
        theme(axis.title.y = element_text(size = 14))+
        theme(axis.title.x = element_text(size = 14))+
        theme(axis.text.y = element_text( size = 14)) +
        scale_color_npg() +
        facet_grid( ~ Stress)
      ggsave(paste0(outdir, name, "_extra1.pdf"), gg, width = 8, height = 4)
      gg2 <- plot_ordination(phobj, pcoa.bray,
                            color = "Treatment",
                            #shape = "Stress",
                            title = name, axes=c(1, 2)) +
        #coord_fixed(sqrt(evals[2] / evals[1])) +
        #scale_color_manual(values=palette2)+
        stat_ellipse(level=0.95, linetype=2, alpha = 0.8, na.rm = TRUE) +
        geom_point(size = 1) +
        #geom_point(size = 1, aes(col=Stress)) +
        #geom_text_repel(aes_string(label = labelsamples)) +
        theme_bw() +
        theme(axis.text.x = element_text(size = 14))+
        theme(strip.text.x = element_text(size = 14))+
        theme(strip.text.y = element_text(size = 14))+
        theme(axis.title.y = element_text(size = 14))+
        theme(axis.title.x = element_text(size = 14))+
        theme(axis.text.y = element_text( size = 14)) +
        scale_color_npg() +
        facet_grid(Stress ~ Region_sequenced)
      ggsave(paste0(outdir, name, "_extra2.pdf"), gg2, width = 8, height = 5)
      
      
    }}}

# Composition 4 each

outdir <- paste0(opt$out, "/DescriptiveAbundances/")
if(!dir.exists(outdir)) dir.create(outdir)
tops <- c(5, 10)


library(plyr)
for(phname in phseq_to_use){
  phobj <- all_phyloseq[[phname]]
  # for(interestvar in vars2test){
  # cat("Doing Abundance Plots for: ", phname, ", ", interestvar, "\n")
  # #abund_plots <- plotAbundanceFullPipeline(all_phyloseq[[phname]], interestvar, outdir, phname, unique(meta3 %>% pull(!!sym(interestvar))), tops)
  #   
  #   
  #   oname <- paste0(outdir, "/relAbund_bySpecies_ColByGenus_", phname, "_", interestvar, "_", as.character(15), ".pdf")
  #   plotRelativeAbnBars_Fantaxtic(phobj, interestvar, topn = 15, tax_level="Genus", outname = oname)
  # }
  oname <- paste0(outdir, "/relAbund_bySpecies_ColByGenus_", phname, "_", "GRID", "_", as.character(15), ".pdf")
  plotRelativeAbnBars_Fantaxtic_grid(phobj, c("Treatment", "Region_sequenced","Stress" ), topn = 15, tax_level="Genus", outname = oname,
                                     height = 10, width = 10)
  }
