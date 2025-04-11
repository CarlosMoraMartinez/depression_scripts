# Alpha 4 each
library(rstatix)

## Cualitativas

signif_levels=c("***"=0.001, "**"=0.01, "*"=0.05, "ns"=1.1)

alpha_indices <- c("Observed", "Chao1", "Shannon", "InvSimpson")
vars2test <- c("Group") #"Category_T0"

quant_vars <- c()
vars2log <- c()
quant_vars_ext <- c()
interestvar <- "Group"


outdir <- paste0(opt$out, "/AlphaDiversity/")
if(!dir.exists(outdir)) dir.create(outdir)

phseq_to_use <- names(all_phyloseq)[!grepl("tribiome", names(all_phyloseq))]
#load(allphyloseqlist_fname)

for(phname in phseq_to_use){
  cat("Alpha diversity in ", phname, "\n")
  phobj <- all_phyloseq[[phname]]


  divtab <- calculateAlphaDiversityTable(phobj, outdir, alpha_indices, paste0(phname, "_AlphaDiv") )

  ## Statistics
  
  means_bygroups = divtab %>% 
    group_by(Group) %>% 
    summarise(across(all_of(alpha_indices), 
                     list(mean=mean, median=median, sd=sd, var=var, max=max, min=min)))
  
  write_tsv(means_bygroups, paste0(outdir, "/", phname, "_AlphaDiv_meansByGroup.tsv"))
  
  ## Differences between groups
  test_group <- data.frame(
        alpha_ind = alpha_indices, 
      ttest = sapply(alpha_indices, \(ai) t.test( as.formula(paste0(ai,  "~ Group")), divtab)$p.value),
      wilcox = sapply(alpha_indices, \(ai) wilcox.test( as.formula(paste0(ai,  "~ Group")), divtab)$p.value),
      anova=sapply(alpha_indices, \(ai) summary(aov( as.formula(paste0(ai,  "~ Group")), divtab))[[1]][1, "Pr(>F)"]),
      shapiro = sapply(alpha_indices, \(ai) shapiro.test( unlist(divtab[, ai]))$p.value),
      bartlett = sapply(alpha_indices, \(ai) bartlett.test( unlist(divtab[, ai]), divtab$Group)$p.value)
      ) %>% rownames_to_column("Index")
    
  
  write_tsv(test_group, paste0(outdir, "/", phname, "_AlphaDiv_testGroup_long.tsv"))
  
  # Typical plots
  divplots <- getAlphaDiversity(phobj, vars2test, quant_vars_ext,
                                opt,
                                indices= alpha_indices,
                                correct_pvalues = F, correct_pvalues_indices = F,
                                name = paste0(phname, "_AlphaDiv"), w = 8, h = 4)
}

## now between tribiome groups

vars2test <- c("Weight") #"Category_T0"
interestvar <- "Weight"

phseq_to_use <- names(all_phyloseq)[grepl("tribiome", names(all_phyloseq))]
#load(allphyloseqlist_fname)

for(phname in phseq_to_use){
  cat("Alpha diversity in ", phname, "\n")
  phobj <- all_phyloseq[[phname]]
  
  
  divtab <- calculateAlphaDiversityTable(phobj, outdir, alpha_indices, paste0(phname, "_AlphaDiv") )
  
  ## Statistics
  
  means_bygroups = divtab %>% 
    group_by(Group) %>% 
    summarise(across(all_of(alpha_indices), 
                     list(mean=mean, median=median, sd=sd, var=var, max=max, min=min)))
  
  write_tsv(means_bygroups, paste0(outdir, "/", phname, "_AlphaDiv_meansByGroup.tsv"))
  
  ## Differences between groups
  test_group <- data.frame(
    alpha_ind = alpha_indices, 
    ttest = sapply(alpha_indices, \(ai) t.test( as.formula(paste0(ai,  "~ Weight")), divtab)$p.value),
    wilcox = sapply(alpha_indices, \(ai) wilcox.test( as.formula(paste0(ai,  "~ Weight")), divtab)$p.value),
    anova=sapply(alpha_indices, \(ai) summary(aov( as.formula(paste0(ai,  "~ Weight")), divtab))[[1]][1, "Pr(>F)"]),
    shapiro = sapply(alpha_indices, \(ai) shapiro.test( unlist(divtab[, ai]))$p.value),
    bartlett = sapply(alpha_indices, \(ai) bartlett.test( unlist(divtab[, ai]), divtab$Weight)$p.value)
  ) %>% rownames_to_column("Index")
  
  
  write_tsv(test_group, paste0(outdir, "/", phname, "_AlphaDiv_testGroup_long.tsv"))
  
  # Typical plots
  divplots <- getAlphaDiversity(phobj, vars2test, quant_vars_ext,
                                opt,
                                indices= alpha_indices,
                                correct_pvalues = F, correct_pvalues_indices = F,
                                name = paste0(phname, "_AlphaDiv"), w = 8, h = 4)
}

# Beta 4 each
outdir <- paste0(opt$out, "/BetaDiversity/")
if(!dir.exists(outdir)) dir.create(outdir)

dists <- c("bray", "jaccard") # , "jaccard"
vars2pcoa <- c("Group")
var2shape = "Weight"
ccaplots <- list()
phseq_to_use <- names(all_phyloseq)[!grepl("tribiome", names(all_phyloseq))]
for(phname in phseq_to_use){
  for(method in c("PCoA", "NMDS")){
    for(dist in dists){
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
      phobj <- all_phyloseq[[phname]]
      pcoa.bray <- ordinate(phobj, method = method, distance = dist)
      evals <- pcoa.bray$values$Eigenvalues

      df2plot <- pcoa.bray$points %>% data.frame %>%
        rownames_to_column()

      gg <- plot_ordination(phobj, pcoa.bray,
                            color = "Group",
                            shape = "Weight",
                            title = name, axes=c(1, 2)) +
        #coord_fixed(sqrt(evals[2] / evals[1])) +
        #scale_color_manual(values=palette2)+
        #stat_ellipse(level=0.95, linetype=2, alpha = 0.8, na.rm = TRUE) +
        geom_point(size = 3) +
        #geom_point(size = 1, aes(col=Stress)) +
        #geom_text_repel(aes_string(label = labelsamples)) +
        theme_bw() +
        theme(axis.text.x = element_text(size = 14))+
        theme(strip.text.x = element_text(size = 14))+
        theme(axis.title.y = element_text(size = 14))+
        theme(axis.title.x = element_text(size = 14))+
        theme(axis.text.y = element_text( size = 14)) +
        scale_color_npg() 
      ggsave(paste0(outdir, name, "_extra1.pdf"), gg, width = 6, height = 4)
      gg2 <- plot_ordination(phobj, pcoa.bray,
                            color = "Group",
                            #shape = "Stress",
                            title = name, axes=c(1, 2)) +
        #coord_fixed(sqrt(evals[2] / evals[1])) +
        #scale_color_manual(values=palette2)+
        #stat_ellipse(level=0.95, linetype=2, alpha = 0.8, na.rm = TRUE) +
        geom_point(size = 3) +
        #geom_point(size = 1, aes(col=Stress)) +
        #geom_text_repel(aes_string(label = labelsamples)) +
        theme_bw() +
        theme(axis.text.x = element_text(size = 14))+
        theme(strip.text.x = element_text(size = 14))+
        theme(strip.text.y = element_text(size = 14))+
        theme(axis.title.y = element_text(size = 14))+
        theme(axis.title.x = element_text(size = 14))+
        theme(axis.text.y = element_text( size = 14)) +
        scale_color_npg() 
      ggsave(paste0(outdir, name, "_extra2.pdf"), gg2, width = 6, height = 5)
      
    }}}


# Beta 4 only tribiome

dists <- c("bray", "jaccard") # , "jaccard"
vars2pcoa <- c("Weight")
var2shape <- c()
ccaplots <- list()
phseq_to_use <- names(all_phyloseq)[grepl("tribiome", names(all_phyloseq))]

for(phname in phseq_to_use){
  for(method in c("PCoA", "NMDS")){
    for(dist in dists){
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
      phobj <- all_phyloseq[[phname]]
      pcoa.bray <- ordinate(phobj, method = method, distance = dist)
      evals <- pcoa.bray$values$Eigenvalues
      
      df2plot <- pcoa.bray$points %>% data.frame %>%
        rownames_to_column()
      
      gg2 <- plot_ordination(phobj, pcoa.bray,
                             color = "Weight",
                             #shape = "Stress",
                             title = name, axes=c(1, 2)) +
        #coord_fixed(sqrt(evals[2] / evals[1])) +
        #scale_color_manual(values=palette2)+
        #stat_ellipse(level=0.95, linetype=2, alpha = 0.8, na.rm = TRUE) +
        geom_point(size = 3) +
        #geom_point(size = 1, aes(col=Stress)) +
        #geom_text_repel(aes_string(label = labelsamples)) +
        theme_bw() +
        theme(axis.text.x = element_text(size = 14))+
        theme(strip.text.x = element_text(size = 14))+
        theme(strip.text.y = element_text(size = 14))+
        theme(axis.title.y = element_text(size = 14))+
        theme(axis.title.x = element_text(size = 14))+
        theme(axis.text.y = element_text( size = 14)) +
        scale_color_npg() 
      ggsave(paste0(outdir, name, "_extra2.pdf"), gg2, width = 6, height = 5)
      
    }}}

# Composition 4 each

outdir <- paste0(opt$out, "/DescriptiveAbundances/")
if(!dir.exists(outdir)) dir.create(outdir)
tops <- c(15)


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


# top taxa in controls only


phobj_controls <- subset_samples(all_phyloseq[["rarefied_min"]], Treatment == "NO ABS" & Stress == "Control")
otutab <- otu_table(phobj_controls) %>% data.frame %>% 
  rownames_to_column("Species") %>%
  mutate_if(is.numeric, \(x) 100*x/sum(x)) %>% 
  gather(key="Sample", value="Abundance", -Species) %>%
  mutate(Sample = gsub("^X", "", Sample, perl=T)) %>%
  merge(s_meta, by.x="Sample", by.y="sampleID", all.x=T, all.y=F)

toptable <- otutab %>% 
  group_by(Species, Region_sequenced) %>% 
  dplyr::summarise(mean = mean(Abundance), median = median(Abundance), sd = sd(Abundance), n = n()) %>% 
  ungroup() %>% 
  group_by(Region_sequenced) %>%
  group_split() %>% map(\(x){
    x %>% 
      arrange(-mean) %>% 
      head(15) %>%
      mutate(Species = factor(gsub("_", " ", Species), levels = gsub("_", " ", Species))) %>% 
      mutate(label = paste(as.character(round(mean, 2)), "%", sep=""))
  }) #%>% bind_rows()

max_xlim <- toptable %>% map(\(x) max(x$mean)) %>% unlist() %>% max()

colors <- ggsci::pal_npg()(3)
plots <- toptable %>% map2(colors, \(tab, cc){
  ggplot(tab, aes(x=Species, y=mean))+
    geom_segment(aes(x=Species, xend=Species, y=0, yend=mean), col=cc) +
    geom_point(size=3, col=cc) +
    geom_text(aes(label=label), col="gray30", nudge_y=4)+
    theme_classic()+
    coord_flip() +
    theme(axis.text.x = element_text(size=12))+
    theme(axis.text.y = element_text(size = 12))+
    theme(axis.text.x = element_text(size = 12))+
    theme(axis.title.y = element_text(size = 12))+
    theme(axis.title.x = element_text(size = 12))+
    theme(axis.text.y = element_text( size = 12, face = "italic")) +
    theme(legend.position = "none")+
    ylab("Mean Relative Abundance")+
    xlab("Species")+
    ylim(0, max_xlim + 8) +
    ggtitle(unique(tab$Region_sequenced))
})
library(cowplot)
pdf(paste0(outdir, "/top15Species_controls2.pdf"), width = 16, height = 4)
cowplot::plot_grid(plotlist = plots, ncol=3)
dev.off()


write_tsv(toptable %>% bind_rows(), paste0(outdir, "/relAbund_bySpecies_ColByGenus_controls.tsv"))
