library(tidyverse)
library(ggsci)
library(phyloseq)
library(ggvenn)
## GGPLOT THEMES

#options(ggplot2.discrete.fill = c("#1E90FF", "#00AA5A", "#F75A3F", "#8E7BFF","#00D1EE", "#00E6BB", "#F9F871", "#F45680", "#A5ABBD", "#B60E50"))
#options(ggplot2.discrete.colour = c("#1E90FF","#00AA5A", "#F75A3F",  "#8E7BFF","#00D1EE", "#00E6BB", "#F9F871", "#F45680", "#A5ABBD", "#B60E50"))

options(ggplot2.discrete.fill = c("#A1C6EA","#FD8B2F", "#00AA5A", "#8E7BFF","#00D1EE", "#00E6BB", "#F9F871", "#F45680", "#A5ABBD", "#B60E50"))
options(ggplot2.discrete.colour = c("#A1C6EA","#FD8B2F","#00AA5A",   "#8E7BFF","#00D1EE", "#00E6BB", "#F9F871", "#F45680", "#A5ABBD", "#B60E50"))
C_CASE = "#FD8B2F" #"rgba(200, 44, 44, 0.8)"
C_CASE2 = "tomato"
C_CASE_LINK = "#fBd895" #"#f9c784"
C_CTRL = "#A1C6EA" # "#006daa" #
C_CTRL2 = "steelblue2"
C_CTRL_LINK ="#DAE3E5" #"rgba(44, 44, 200, 0.8)"
C_CTRL_LINK2 ="#B8C1D4"
C_WHITE= "#DDDDDD"
C_NS =  "#A5ABBD" #"rgba(224, 224, 224, 0.8)"
C_OTHER = "gray30"

options(ggplot2.continuous.fill="viridis")
options(ggplot2.continuous.colour="viridis")

mytheme <-  theme_bw()+
  theme(plot.title = element_text(size = rel(1.5), hjust = 0.5,
                                  colour = "black", face = "bold")) +
  theme(legend.title = element_text(face = "bold")) +
  theme(axis.text.y = element_text(size = 12, 
                                   colour = "black", angle = 0, face = "bold")) +
  theme(strip.text.y = element_text(size = 12, 
                                    colour = "black", angle = 0, face = "bold")) +
  theme(strip.text.x = element_text(size = 12, 
                                    colour = "black", angle = 0, face = "bold")) +
  theme(axis.text.x = element_text(size = 11, 
                                   colour = "black", angle = 0, 
                                   face = "bold"))+
  theme(axis.title.x = element_text(vjust = 1, hjust = 0.5, 
                                    size = 12, colour = "black", 
                                    angle = 0, face = "bold")) +
  theme(axis.title.y= element_text(vjust = 1, hjust = 0.5, 
                                   size = 12, colour = "black", 
                                   angle = 90, face = "bold"))



mystyle <- theme_classic() +
  theme(axis.text = element_text(face="bold"), 
        axis.title = element_text(face="bold")) 

thin_barplot_lines <- theme(panel.grid.major.y = element_line(color = "lightgray",
                                                              size = 0.05,
                                                              linetype = 2))

restauraropt_mk <- function(opt){
  output <- opt$out
  restaurar <- function(optbad){
    optbad$out <- output
    return(opt)
  }
}
### Basic phlyloseq

MULTI_PAGE_PDFS = TRUE

updatePsWithLogs <- function(phobj, vars = c("nreads")){
  
  for(v in vars){
    newname <- paste0(v, "_log")
    sample_data(phobj)[[newname]] <- log(sample_data(phobj)[[v]] +1  )
  }
  return(phobj)
}

makeTotalCountPlot <- function(metadata, variable, outdir, name = ""){
  g1 <- ggplot(metadata, aes_string(x=variable, y="nreads", col=variable, fill=variable)) +
    geom_boxplot(alpha=0.6, width=0.5) +
    #scale_color_lancet() + 
    #scale_fill_lancet() +
    ylab("Total counts") +
    geom_jitter() + #size=3*samplesums$IMC/max(samplesums$IMC[!is.na(samplesums$IMC)])
    mytheme
  outname <- paste0(outdir,"/", name, '_', variable, ".pdf")
  ggsave(filename=outname, g1)
  return(g1)
}

makeRarefactionCurve <- function(phyloseq_rawdata, opt){
  otu2rare <- otu_table(phyloseq_rawdata)
  class(otu2rare) <- "matrix"
  png(paste0(opt$out, "/00_plot_rarecurve.png"), width = 15, height = 15, units = "cm", res = 300)
  
  rarplot <- rarecurve(t(otu2rare),
                       step = 100,
                       sample = 20000,
                       col = "blue",
                       cex = 0.75,
                       main = "Rarefaction curve")
  dev.off()
  
  #Table with all data
  rardf <- mapply(rarplot,  colnames(otu2rare), SIMPLIFY=FALSE, 
                  FUN=function(x, s) data.frame(sample=s, nreads = attr(x, "Subsample"), txcount=x)) %>%                   bind_rows() %>% 
    group_by(sample) %>% 
    dplyr::mutate(last_point = ifelse(nreads == max(nreads), 2, 0))
  
  #Table with max reads in each sample, to plot points
  rarmax <- rardf %>% filter(last_point == 2)
  
  #Table with the intersection between min reads of all samples and each of the curves, to plot hlines
  MIN_READS_ALLSAMPLES <- min(rarmax$nreads)
  
  intersection <- rardf %>% 
    dplyr::mutate(diff2min = abs( nreads - MIN_READS_ALLSAMPLES)) %>% 
    slice_min(n=2, order_by=diff2min) %>% 
    dplyr::mutate(prop_dist = 1- diff2min/sum(diff2min))
  intersection2 <- intersection %>% dplyr::summarise(nreads=MIN_READS_ALLSAMPLES,
                                                    txcount = sum(txcount*prop_dist))
  
  ggrare <- ggplot(rardf, aes(col=sample, x=nreads, y=txcount)) +
    geom_line() +
    geom_point(data=rarmax)+
    geom_text_repel(data=rarmax, aes(label=sample))+
    geom_hline(data=intersection2, aes(yintercept=txcount, col=sample), linetype=2, size=0.2) +
    geom_vline(xintercept=MIN_READS_ALLSAMPLES, col="gray40", linetype=2) +
    xlab("Species") +
    ylab("Sample Size") + 
    ggtitle("Rarefaction curve") +
    scale_x_continuous(labels = function(x) format(x, scientific = TRUE)) +
    mytheme +
    theme(legend.position="none")
  
  ggsave(paste0(opt$out, "/00_plot_rarecurve2.pdf"), ggrare)
  
  return(ggrare)
  
}

calculateAlphaDiversityTable <- function(phseq_obj, outdir, 
                                         indices = c("Observed", "Chao1", "Shannon", "InvSimpson", "Fisher"), 
                                         name="diversity"){
  library(ggpmisc) #Add regression formula
  div <- estimate_richness(phseq_obj, 
                       measures = indices)
  div$sampleID <- gsub("^X", "", rownames(div))
  div2 <- merge(data.frame(sample_data(phseq_obj)), div, by="sampleID")
  
  write_tsv(div2, file=paste0(outdir, "/", name, ".tsv") )
  return(div2)
}

getTestsForAllCombinations <- function(var, group, sep="-", paired=F){
  comp <- combn(unique(group[!is.na(group)]), 2, simplify = F)
  res <- lapply(comp, function(combo){
    a <- var[group==combo[1]]
    b <- var[group==combo[2]]
    res <- data.frame(groups_compared = paste(combo, sep=sep, collapse=sep),
                      t_pval = t.test(a, b, paired=paired)$p.value,
                      wilcox_pval = wilcox.test(a, b, paired=paired)$p.value,
                      shapiro_test =  shapiro.test(c(a, b))$p.value,
                      bartlett_test = bartlett.test(var[group %in% combo], 
                                                    group[group %in% combo])[[3]],
                      levene_test = leveneTest(var[group %in% combo], 
                                               factor(group[group %in% combo]))[1, 3]
    )
  }) %>% bind_rows  %>% dplyr::arrange(groups_compared)
  return(res)
}

getValidFactors <- function(df, max_nas=0.2, min_pct=0.1, max_classes=5, exclude=c()){
  is_factor <- function(x){
    (is.factor(x) | 
       is.character(x) | 
       (is.numeric(x) & length(unique(x)) <= max_classes))
  }
  meet_nas <- function(x){
    sum(is.na(x) | (x== "-"))/length(x) <= max_nas 
  }
  is_balanced <- function(x){
    tt <- table(x)
    tt <- tt/sum(tt)
    min(tt) >= min_pct & (length(tt) > 1)
  }
  
  
  res <- data.frame(var = names(df),
                    is_factor = sapply(df, is_factor),
                    meet_nas = sapply(df, meet_nas),
                    is_balanced = sapply(df, is_balanced),
                    not_excluded = (!names(df) %in% exclude)
  ) %>% dplyr::mutate(is_valid = is_factor & meet_nas & is_balanced & not_excluded)
  return(res)
}

getValidNumeric <- function(df, max_nas=0.0, max_classes=4){
  is_number <- function(x){
    (is.numeric(x) &
       (length(unique(x)) >= max_classes))
  }
  meet_nas <- function(x){
    sum(is.na(x) | (x== "-"))/length(x) <= max_nas 
  }
  
  res <- data.frame(var = names(df),
                    is_number = sapply(df, is_number),
                    meet_nas = sapply(df, meet_nas)
  ) %>% dplyr::mutate(is_valid = is_number & meet_nas)
  return(res)
}

# Test differences in several variables (vars) between groups in several grouping variables (groupvars)
testDiversityDifferences <- function(divtab, vars, groupvars, outdir, name="alpha_diversity"){
  library(car)
  res <- data.frame()
  for (v in vars){
    if(length(unique(divtab[, v])) < 2 ) next
    for(g in groupvars){
      if(length(unique(divtab[, g])) < 2 | min(table(divtab[!is.na(divtab[, v]), g]))<2) next
      meantab <- tapply(divtab[, v], divtab[, g], mean, na.rm=T)
      if(any(is.na(meantab))) next
      divtab[, g] <- as.factor(divtab[, g])
      divtab[, v] <- as.numeric(divtab[, v])
      num_groups <- length(unique(divtab[!is.na(divtab[, g]), g]))
      form <- as.formula(paste0(v, " ~ ", g))
      aovres <- aov(form, divtab)
      tres <- tryCatch(ifelse(num_groups==2, t.test(form, divtab)$p.value, NA), error=function(x)return(NA))
      wres <- tryCatch(ifelse(num_groups==2, wilcox.test(form, divtab)$p.value, NA), error=function(x)return(NA))
      swres <- shapiro.test(divtab[, v])$p.value
      bt <-  bartlett.test(form, divtab)[[3]]
      levt <-  leveneTest(form, divtab)[1, 3]
      aux <- data.frame(
        variable = v,
        groups = g, 
        comparison = "all",
        anova_F = summary(aovres)[[1]]["F value"][1, 1],
        anova_p = summary(aovres)[[1]]["Pr(>F)"][1, 1],
        t_test = tres,
        wilcox_test = wres,
        shapiro_normality_test = swres,
        bartlett_test = bt,
        levene_test = levt
      )
      res <- rbind(res, aux)
      if(num_groups > 2){
        aux_t <- getTestsForAllCombinations(divtab[, v], divtab[, g]) 
        aux2 <- data.frame(
          variable = v,
          groups = g, 
          comparison = aux_t$groups_compared,
          anova_F = NA,
          anova_p = NA,
          t_test = aux_t$t_pval,
          wilcox_test = aux_t$wilcox_pval,
          shapiro_normality_test = aux_t$shapiro_test,
          bartlett_test = aux_t$bartlett_test,
          levene_test = aux_t$levene_test          
        )
        res <- rbind(res, aux2)
      } ## if num_groups > 2
    } # For grouping variabbles
  } # For problem variables
  res$t_corrected <- p.adjust(res$t_test, method="BH")
  res$wilcox_corrected <- p.adjust(res$wilcox_test, method="BH")
  write_tsv(res, paste0(outdir, "/", name, "_QualitVarsTests.tsv"))
  return(res)
}

testDiversityDifferences_paired <- function(divtab, vars, groupvars, pairvar ,outdir, name="alpha_diversity"){
  library(car)
  divtab <- divtab[order(divtab[,pairvar]) , ]
  res <- data.frame()
  for (v in vars){
    for(g in groupvars){
      divtab[, g] <- as.factor(divtab[, g])
      divtab[, v] <- as.numeric(divtab[, v])
      gr_a <-divtab[divtab[, g] == levels(divtab[, g])[1] ,v]
      gr_b <-divtab[divtab[, g] == levels(divtab[, g])[2] ,v]
      num_groups <- length(unique(divtab[!is.na(divtab[, g]), g]))
      form2<- as.formula(paste0(v, " ~ ", g, " + Error(", pairvar, "/", g, ")"))
      form <- as.formula(paste0(v, " ~ ", g))
      aovres <- aov(form2, divtab)
      sum_aovres <- summary(aovres)[[2]] %>% unlist
      tres <- tryCatch(ifelse(num_groups==2, t.test(gr_a, gr_b, paired = T)$p.value, NA), error=function(x)return(NA))
      wres <- tryCatch(ifelse(num_groups==2, wilcox.test(gr_a, gr_b, paired = T)$p.value, NA), error=function(x)return(NA))
      swres <- shapiro.test(divtab[, v])$p.value
      bt <-  bartlett.test(form, divtab)[[3]]
      levt <-  leveneTest(form, divtab)[1, 3]
      aux <- data.frame(
        variable = v,
        groups = g, 
        comparison = "all",
        anova_F =  sum_aovres["F value1"],
        anova_p = sum_aovres["Pr(>F)1"],
        t_test = tres,
        wilcox_test = wres,
        shapiro_normality_test = swres,
        bartlett_test = bt,
        levene_test = levt
      )
      res <- rbind(res, aux)
      if(num_groups > 2){
        aux_t <- getTestsForAllCombinations(divtab[, v], divtab[, g]) 
        aux2 <- data.frame(
          variable = v,
          groups = g, 
          comparison = aux_t$groups_compared,
          anova_F = NA,
          anova_p = NA,
          t_test = aux_t$t_pval,
          wilcox_test = aux_t$wilcox_pval,
          shapiro_normality_test = aux_t$shapiro_test,
          bartlett_test = aux_t$bartlett_test,
          levene_test = aux_t$levene_test          
        )
        res <- rbind(res, aux2)
      } ## if num_groups > 2
    } # For grouping variabbles
  } # For problem variables
  res$t_corrected <- p.adjust(res$t_test, method="BH")
  res$wilcox_corrected <- p.adjust(res$wilcox_test, method="BH")
  write_tsv(res, paste0(outdir, "/", name, "_QualitVarsTestsPaired.tsv"))
  return(res)
}

#Makes a linear regression with vars as dependent variables and qvars as predictors
testDiversityWithQuantVars <- function(divtab, vars, qvars,
                                       outdir,
                                       name = "AlphaDiv_quant_vars_regressions"){
  library(broom)
  res <- data.frame()
  
  getLm <- function(v, g){
    aux <- divtab %>% dplyr::mutate_at(c(v, g), scale )
    
    form <- as.formula(paste0(v, " ~ ", g))
    tryCatch(mod <- lm(form, data=aux), error=function(x)cat("Error fitting ", as.character(form)))
    return(mod)
  }
  
  df <- expand.grid(variable=vars, predictor=qvars) %>% 
    dplyr::mutate_all(as.character) %>% 
    dplyr::mutate(model = map2(variable, predictor,getLm)) %>% 
    dplyr::mutate(model2 = map(model, broom::glance)) %>% 
    dplyr::mutate(intercept = sapply(model, function(mod) summary(mod)$coefficients[1] ),
           slope = sapply(model, function(mod) summary(mod)$coefficients[2] ),
           std_error_intercept = sapply(model, function(mod) summary(mod)$coefficients[3] ),
           std_error_slope = sapply(model, function(mod) summary(mod)$coefficients[4] )
           ) %>% 
    unnest(model2)
  fname <- paste0(outdir, "/", name, ".tsv")
  write_tsv(df, file=fname)
  return(df)
}

#Transform a vector of p-values into a vector of characters, "***", "**", etc
getSignif <- function(vector){
  x <- ifelse(vector < 0.001, "***",
              ifelse(vector < 0.01, "**",
                     ifelse(vector< 0.05, "*", "NS")
                     )
              )
  return(x)
}

WriteManyPlots <- function(plotlist, name, outdir, w=12, h=7, separate=F, opt=list()){
  if(! dir.exists(outdir)){dir.create(outdir)}
  if(separate){
    for(p in names(plotlist)){
      fname = paste0(outdir, "/", name, "_", p, ".pdf")
      pdf(fname, width=w, height=g)
      tryCatch({print(plotlist[[p]])}, error = function(x){print(getNullPlot(opt, p))})
      dev.off()
    }
  }else{
    fname = paste0(outdir, "/", name, "_all.pdf")
    pdf(fname, width=w, height=h)
    for(p in names(plotlist)){
      tryCatch({print(plotlist[[p]])}, error = function(x){print(getNullPlot(opt, p))})
    }
    dev.off()
  }
}

getAlphaDiversity <- function(phseq_obj, vars, qvars= c(), 
                              opt = list(), 
                              indices=c("Observed", "Chao1", "Shannon", "InvSimpson", "Fisher"),
                              name="AlphaDiversity", 
                              signif_levels=c("***"=0.001, "**"=0.01, "*"=0.05, "ns"=1.1),
                              correct_pvalues = TRUE, correct_pvalues_indices = FALSE,
                              test2show="wilcox.test", w=8, h=4){
  outdir <- paste0(opt$out, "AlphadivPlots")
  if(! dir.exists(outdir)){dir.create(outdir)}
  
  #Get Alpha Diversity values
  divtab <- calculateAlphaDiversityTable(phseq_obj, outdir, indices, name)
  vars <-  map_vec(divtab[, vars], \(x)length(unique(x[!is.na(x)]))) %>% 
    base::subset(. > 1) %>% names
  # Get statistical tests (not used later)
  alphadif <- testDiversityDifferences(divtab, indices, vars, outdir, name)
  plots <- list()
  sign_df <- tibble()
  for(v in vars){
    #Si no lo convertimos a caracter, geom_signif/stat_signif fallan
    standa1 <- sample_data(phseq_obj) %>% data.frame
    standa1 <- standa1$sampleID[!is.na(standa1[, v])]
    phseq_obj_filt <- phyloseq::prune_samples(standa1, phseq_obj) 
    sample_data(phseq_obj_filt)[, v] <- sample_data(phseq_obj_filt)[, v] %>% unlist %>% as.character

    divtab2 <- divtab[!is.na(divtab[, v]), ]
    divtab2[, v] <- as.character(divtab2[, v])
    
    comp <- combn(unique(divtab2[, v]), 2, simplify = F)
    num_comparisons <- alphadif %>% 
      filter(comparison != "all" & 
               groups == v & 
               variable == "Observed") %>% 
      nrow
    if(correct_pvalues & num_comparisons>1){
      signif_levels_bonferroni <- c(signif_levels[1:3]/(num_comparisons*length(indices)), signif_levels[4])
    }else{
      signif_levels_bonferroni <- signif_levels
    }
    if(correct_pvalues_indices ){
      signif_levels_bonferroni <- c(signif_levels_bonferroni[1:3]/length(indices), 
                                    signif_levels_bonferroni[4])
    }
    
    plots[[v]] <- plot_richness(phseq_obj_filt, x = v,
                                color = v, 
                                measures = indices) +
      #ggplot(divtab, aes(x=Psoriasis, y=Shannon, col=Psoriasis, fill=Psoriasis)) +
      #facet_wrap(. ~ Sex, scales = "free") +
      geom_boxplot(aes_string(fill = v), alpha = 0.7, width=0.5) +
      #scale_color_manual(values = c("#ffafcc", "#90DBF4")) + 
      #scale_fill_manual(values = c("#ffafcc", "#90DBF4")) +
      #scale_color_lancet() + 
      #scale_fill_lancet() +
      labs(title = v, x = '') +
      theme_pubclean() +
      mytheme +
      ggsignif::stat_signif(test=test2show, na.rm=T, comparisons = comp, 
                  step_increase=0.06,
                  tip_length = 0.01,
                  map_signif_level=signif_levels_bonferroni,
                  vjust=0.3,
                  color = "black"
      ) +
      theme(plot.title = element_text(hjust = 0.5)) +
      theme(axis.text.x = element_blank()) 
    
    
    #theme(axis.text.x = element_text(angle = 360, hjust = 0.5, size = 10))
  } #Plots qualitative variables
  
  regressions <- testDiversityWithQuantVars(divtab, indices, qvars, outdir, paste0(name, "QuantVarsRegression"))
  
  for(v in qvars){
    auxtext <- regressions %>% filter(predictor == v) %>% 
      dplyr::mutate(text = paste0("R^2=", as.character(round(r.squared, 2)), ", p=", as.character(round(p.value, 3)) ))
    
    plots[[v]] <- plot_richness(phseq_obj, x = v,
                                color = v, 
                                measures = c("Observed", "Chao1", "Shannon", "InvSimpson")) +
      geom_point(aes_string(fill = v), alpha = 0.7) +
      stat_poly_eq(use_label(c("R2", "p")), #c("eq", "R2", "f", "p", "n")
                   method="lm", small.p=T, small.r=F, label.y=0.99)+
      #stat_poly_line(method = "lm") +
      geom_smooth(method="lm", fullrange = TRUE, linetype=1) +
      #scale_color_manual(values = c("#ffafcc", "#90DBF4")) + 
      #scale_fill_manual(values = c("#ffafcc", "#90DBF4")) +
      labs(title = v, x = '') +
      #geom_text(data=auxtext, aes(y=3, x = 20, labels=text)) +
      
      theme_pubclean() +
      mytheme 
    theme(plot.title = element_text(hjust = 0.5)) +
      #theme(axis.text.x = element_blank())
      theme(axis.text.x = element_text(angle = 360, hjust = 0.5, size = 10))
  }
  names(plots) <- gsub("_", "-", names(plots))
  WriteManyPlots(plots, name, outdir, w=w, h=h, separate=F)
  return(plots)
}


makePCoA <- function(phobj, pcoa.bray, evals, var2color="Condition", name = "Bray-Curtis", extradims= 2:5, labelsamples="sampleID"){
  gs <- lapply(extradims, FUN=function(axis){
    gg <- plot_ordination(phobj, pcoa.bray, color = var2color, 
                          title = name, axes=c(1, axis)) + 
      #coord_fixed(sqrt(evals[2] / evals[1])) +
      #scale_color_manual(values=palette2)+ 
      stat_ellipse(level=0.95, linetype=2, alpha = 0.8, na.rm = TRUE) +
      geom_point(size = 2) +
      geom_text_repel(aes_string(label = labelsamples)) +
      mytheme 
    # if(!is.numeric(gg$data[, var2color])){
    #   gg <- gg + scale_color_lancet() + scale_fill_lancet()
    # }
    return(gg)
  })
  cw <- cowplot::plot_grid(plotlist=gs)
  return(cw)
}
makeAllPCoAs <- function(phobj, outdir, 
                         method = "PCoA",
                         name="PCoAs", 
                         dist_type = "bray", 
                         dist_name = "Bray-Curtis", 
                         vars2plot = c(),
                         extradims = 2:5,
                         create_pdfs = MULTI_PAGE_PDFS,
                         labelsamples="sampleID", w=12, h=5){
  #palette2 <- RColorBrewer::brewer.pal(n = 9, name = 'Set1')[8:9]
  pcoa.bray <- ordinate(phobj, method = method, distance = dist_type)
  
  evals <- pcoa.bray$values$Eigenvalues
  
  if(create_pdfs){
    # Plot all the variables and save to PDF
    if(length(vars2plot) == 0){
      vars2plot <- phobj %>% sample_data() %>% names()
      vars2plot <- vars2plot[! vars2plot %in% c("Codigo", "Num_paciente")]
    }
    all_pcoas_plots <- lapply(vars2plot, FUN=function(vv, phobj, pcoa.bray, evals){
      makePCoA(phobj, pcoa.bray, evals, vv, dist_name, extradims, labelsamples = labelsamples)
    },phobj, pcoa.bray, evals)
    names(all_pcoas_plots) <- vars2plot
    
    WriteManyPlots(all_pcoas_plots, name, outdir, w=w, h=h, separate=F, opt)
    
  }else{
    #Plot only one variable and return the plot without saving it
    all_pcoas_plots <- list("Condition"=makePCoA(phobj, pcoa.bray, evals, vars2plot[1], dist_name, extradims))
  }
  
  return(all_pcoas_plots)
}


makeConstrainedOrdination<-function(phobj, variables = c("Condition"),  dist_type="bray", color_by=""){
  phobj_not_na <- phobj
  metad <- sample_data(phobj_not_na) %>% data.frame
  for(v in variables){
    na_samples_2Remove <- metad$sampleID[! is.na(metad[, v]) ]
    phobj_not_na <- prune_samples(na_samples_2Remove, phobj_not_na)
    #phobj_not_na <- subset_samples(phobj_not_na, sampleID %in% na_samples_2Remove)
    metad <- sample_data(phobj_not_na) %>% data.frame
  }
  form <- paste0("~ ", paste(variables, sep=" + ", collapse = " + ")) %>% as.formula
  bray_not_na <- phyloseq::distance(physeq = phobj_not_na, method = dist_type)
  
  # CAP ordinate
  cap_ord <<- ordinate(
    physeq = phobj_not_na, 
    method = "CAP",
    distance = bray_not_na,
    formula = form
  )
  
  # CAP plot
  if(color_by == ""){
    color_by <- variables[1]
  }
  cap_plot <- plot_ordination(
    physeq = phobj_not_na, 
    ordination = cap_ord, 
    color = color_by, 
    axes = c(1,2)
  ) + 
    aes_string(col = color_by) + 
    geom_point( size = 4)
    #geom_point(colour = "grey90", size = 1.5, type=1) 
  #scale_color_manual(values = c("#a65628", "red", "#ffae19", "#4daf4a", 
  #                              "#1919ff", "darkorchid3", "magenta"))
  
  # if(! is.numeric(metad[, color_by])){
  #   cap_plot <- cap_plot + scale_color_lancet()
  # } 
  # Now add the environmental variables as arrows
  arrowmat <- vegan::scores(cap_ord, display = "bp")
  
  # Add labels, make a data.frame
  arrowdf <- data.frame(labels = rownames(arrowmat), arrowmat)
  text_colors <- rep("gray19", length(variables))
  if(color_by %in% variables){
    text_colors[color_by == variables] <- "green4"
  }
  
  # Define the arrow aesthetic mapping
  arrow_map <- aes(xend = CAP1, 
                   yend = CAP2, 
                   x = 0, 
                   y = 0, 
                   shape = NULL, 
                   color = NULL, 
                   label = labels)
  
  label_map <- aes(x = 1.3 * CAP1, 
                   y = 1.3 * CAP2, 
                   shape = NULL, 
                   color = NULL, 
                   label = labels)
  
  arrowhead = arrow(length = unit(0.02, "npc"))
  
  # Make a new graphic
  g1 <- cap_plot + 
    geom_text_repel( aes(label=sampleID))+
    geom_segment(
      mapping = arrow_map, 
      size = 1, 
      data = arrowdf, 
      color = text_colors, #"gray", 
      arrow = arrowhead
    ) + 
    geom_text(
      mapping = label_map, 
      size = 6,  
      color = text_colors,
      data = arrowdf, 
      show.legend = FALSE
    ) + mytheme
  
  return(g1)
}

makeConstrainedOrdinationSingleVar<-function(phobj, variables = c("Condition"),  dist_type="bray"){
  phobj_not_na <- phobj
  metad <- data.frame(sample_data(phobj_not_na))
  na_samples_2Remove <- metad$sampleID[! is.na(metad[, variables[1]])]
  phobj_not_na <- prune_samples(na_samples_2Remove, phobj_not_na)
  #phobj_not_na <- subset_samples(phobj_not_na, ! sampleID %in% na_samples_2Remove)
  
  form <- paste0("~ ", paste(variables, sep=" + ", collapse = " + ")) %>% as.formula
  bray_not_na <- phyloseq::distance(physeq = phobj_not_na, method = dist_type)
  
  # CAP ordinate
  cap_ord <<- ordinate(
    physeq = phobj_not_na, 
    method = "CAP",
    distance = bray_not_na,
    formula = form
  )
  
  # CAP plot
  cap_plot <- plot_ordination(
    physeq = phobj_not_na, 
    ordination = cap_ord, 
    color = variables[1], 
    axes = c(1,2)
  ) + 
    aes_string(col = variables[1]) + 
    geom_point( size = 4) 
    #geom_point(colour = "grey90", size = 1.5, type=1) 
  
    # if(! is.numeric(metad[, variables[1]])){
    #   cap_plot <- cap_plot + scale_color_lancet()
    # } 
    
  #scale_color_manual(values = c("#a65628", "red", "#ffae19", "#4daf4a", 
  #                              "#1919ff", "darkorchid3", "magenta"))
  
  
  # Now add the environmental variables as arrows
  arrowmat <- vegan::scores(cap_ord, display = "bp")
  
  # Add labels, make a data.frame
  arrowdf <- data.frame(labels = rownames(arrowmat), arrowmat)
  
  # Define the arrow aesthetic mapping
  arrow_map <- aes(xend = CAP1, 
                   yend = 0, 
                   x = 0, 
                   y = 0, 
                   shape = NULL, 
                   color = NULL, 
                   label = labels)
  
  label_map <- aes(x = 1.3 * CAP1, 
                   y = 1.3 * 0, 
                   shape = NULL, 
                   color = NULL, 
                   label = labels)
  
  arrowhead = arrow(length = unit(0.02, "npc"))
  
  # Make a new graphic
  g1 <- cap_plot + 
    geom_text_repel( aes(label=sampleID))+
    geom_segment(
      mapping = arrow_map, 
      size = 1, 
      data = arrowdf, 
      color = "gray", 
      arrow = arrowhead
    ) + 
    geom_text(
      mapping = label_map, 
      size = 6,  
      data = arrowdf, 
      show.legend = FALSE
    ) + mytheme
  
  return(g1)
}

makeConstrainedOrdinationAll <- function(phobj, opt, name="CCA", outdir= "CCA", variableList = list(c("Condition")),  dist_type="bray"){
  plots <- list()
  i = 1
  for(vars in variableList){
    if(length(vars) ==1 ){
      plots[[i]] <- makeConstrainedOrdinationSingleVar(phobj, variables =vars,  dist_type=dist_type)
    }else{
      plots[[i]] <- makeConstrainedOrdination(phobj, variables =vars,  dist_type=dist_type)
    }
    i <- i+1
  }
  names(plots)<- sapply(variableList, paste, sep="-", collapse="-")
  WriteManyPlots(plots, name, outdir, opt=opt)
  return(plots)
}

### Abundance barplots etc

getRelAbundanceTab <- function(phobj){
  ottmp <- phyloseq::otu_table(phobj)
  
  pre_prevalence <- apply(X = ottmp,
                          MARGIN = ifelse(taxa_are_rows(phobj), yes = 1, no = 2),
                          FUN = function(x){sum(x > 0)})
  pre_prevalence = data.frame(Prevalence = pre_prevalence,
                              TotalAbundance = phyloseq::taxa_sums(phobj),
                              tax_table(phobj), 
                              relative_prevalence = pre_prevalence/ nsamples(phobj)
  )
  return(pre_prevalence)
}

getMeanRelAbundancesByPhylum <- function(phobj){
  #Mean prevalence (absolute, not relative) of ASVs by Phylum
   pre_prevalence <- getRelAbundanceTab(phobj)
   df_prevalence <- ddply(pre_prevalence, "Phylum", function(df1){cbind(mean(df1$Prevalence),
                                                                        sum(df1$Prevalence))})
  
  names(df_prevalence) <- c("Phylum", "prevalence", "um of prevalences")
  return(df_prevalence)
}

getMeanRelAbundancesByGenus <- function(phobj){
  #Mean prevalence of ASVs by Genus (absolute prevalence, not relative)
  pre_prevalence <- getRelAbundanceTab(phobj)
  
  ps_rel_abund = phyloseq::transform_sample_counts(phobj, function(x){x / sum(x)})
  
  df_prevalence <- ddply(pre_prevalence, "Genus", function(df1){cbind(mean(df1$Prevalence),
                                                                      sum(df1$Prevalence))})
  
  names(df_prevalence) <- c("Genus", "prevalence", "Sum of prevalences")
  return(df_prevalence)
}

getRelAbundancesByPhylumAndVariable <- function(phobj, variable="Condition", 
                                                outname="prevalence_by_phylum.tsv",
                                                oldlevs=c("Control", "Depression")){
  #Mean prevalence of Phyla
  levs=c("no", "yes")
  aa <- psmelt(phobj) %>% data.frame 
  aa$Condition <- levs[match(aa[, variable], oldlevs)]

  asvmean <- aa %>% group_by(OTU) %>% 
    dplyr::summarise(Phylum = unique(Phylum),
                     mean_abund = mean(Abundance),
                     prevalence = sum(Abundance>0),
                     num_samples = n(),
                     rel_prevalence = prevalence/n()) %>% 
    group_by(Phylum) %>% 
    dplyr::summarise(mean_rel_prev = 100*mean(rel_prevalence)
    )
  
  aa <- aa %>% 
    group_by(sampleID, Condition, Phylum) %>% 
    dplyr::summarise(
      Total_Abundance = sum(Abundance),
      ASV_prevalence = sum(Abundance > 0),
      num_ASVs = n(),
      rel_ASV_prevalence = ASV_prevalence/num_ASVs
    ) %>% 
    ungroup() %>% 
    group_by(Condition, Phylum) %>% 
    dplyr::summarise(prevalence = sum(Total_Abundance>0),
                     mean_abundance = mean(Total_Abundance),
                     n_samples = length(unique(sampleID)),
                     rel_prevalence = prevalence/n_samples,
                     mean_asv_prevalence = mean(rel_ASV_prevalence)
    )
  
  
  bb <- aa %>% gather("variable", "value", prevalence, 
                      mean_abundance, rel_prevalence, 
                      n_samples, mean_asv_prevalence) %>% 
    unite("tmp" , Condition, variable) %>% 
    tidyr::spread(tmp, value) %>% 
    dplyr::mutate(total_n_samples = no_n_samples + yes_n_samples, 
                  total_prevalence = yes_prevalence + no_prevalence,
                  total_rel_prevalence = total_prevalence/total_n_samples, 
                  mean_asv_prevalence = (yes_mean_asv_prevalence*yes_n_samples + no_mean_asv_prevalence*no_n_samples)/total_n_samples)
  if(outname != "") write_tsv(bb, file = outname)
  
  return(bb)
}


getRelAbundancesByPhylumAndVariable_beforeAfter <- function(phobj, variable="Condition", 
                                                            outname="prevalence_by_phylum.tsv"){
  #Mean prevalence of Phyla
  
  aa <- psmelt(phobj) %>% data.frame 
  aa$Condition <- aa[, variable]
  
  asvmean <- aa %>% group_by(OTU) %>% 
    dplyr::summarise(Phylum = unique(Phylum),
                     mean_abund = mean(Abundance),
                     prevalence = sum(Abundance>0),
                     num_samples = n(),
                     rel_prevalence = prevalence/n()) %>% 
    group_by(Phylum) %>% 
    dplyr::summarise(mean_rel_prev = 100*mean(rel_prevalence)
    )
  
  aa <- aa %>% 
    group_by(sampleID, Condition, Phylum) %>% 
    dplyr::summarise(
      Total_Abundance = sum(Abundance),
      ASV_prevalence = sum(Abundance > 0),
      num_ASVs = n(),
      rel_ASV_prevalence = ASV_prevalence/num_ASVs
    ) %>% 
    ungroup() %>% 
    group_by(Condition, Phylum) %>% 
    dplyr::summarise(prevalence = sum(Total_Abundance>0),
                     mean_abundance = mean(Total_Abundance),
                     n_samples = length(unique(sampleID)),
                     rel_prevalence = prevalence/n_samples,
                     mean_asv_prevalence = mean(rel_ASV_prevalence)
    )
  
  
  bb <- aa %>% gather("variable", "value", prevalence, 
                      mean_abundance, rel_prevalence, 
                      n_samples, mean_asv_prevalence) %>% 
    unite("tmp" , Condition, variable) %>% 
    tidyr::spread(tmp, value) %>% 
    dplyr::mutate(total_n_samples = Before_n_samples + After_n_samples, 
                  total_prevalence = Before_prevalence + After_prevalence,
                  total_rel_prevalence = total_prevalence/total_n_samples, 
                  mean_asv_prevalence = (Before_mean_asv_prevalence*Before_n_samples + Before_mean_asv_prevalence*Before_n_samples)/total_n_samples)
  if(outname != "") write_tsv(bb, file = outname)
  
  return(bb)
}

getRelAbundancesByGenusAndVariable <- function(phobj, variable="Condition", 
                                               outname="prevalence_by_genus.tsv",
                                               oldlevs=c("Control", "Depression")){
  #Mean prevalence of Phyla
  
  
  aa <- psmelt(phobj) %>% data.frame 
  aa$Condition <- aa[, variable]
  
  levs=c("no", "yes")
  aa <- psmelt(phobj) %>% data.frame 
  aa$Condition <- levs[match(aa[, variable], oldlevs)]
  
  asvmean <- aa %>% group_by(OTU) %>% 
    dplyr::summarise(Genus = unique(Genus),
              mean_abund = mean(Abundance),
              prevalence = sum(Abundance>0),
              num_samples = length(unique(sampleID)),
              rel_prevalence = prevalence/num_samples) %>% 
    group_by(Genus ) %>% 
    dplyr::summarise(mean_rel_prev = mean(rel_prevalence)
    )
  
  aa <- aa %>% 
    group_by(sampleID, Condition, Genus) %>% 
    dplyr::summarise(
      Total_Abundance = sum(Abundance),
      ASV_prevalence = sum(Abundance > 0),
      num_ASVs = n(),
      rel_ASV_prevalence = ASV_prevalence/num_ASVs
    ) %>% 
    ungroup() %>% 
    group_by(Condition, Genus) %>% 
    dplyr::summarise(prevalence = sum(Total_Abundance>0),
              mean_abundance = mean(Total_Abundance),
              n_samples = length(unique(sampleID)),
              rel_prevalence = prevalence/n_samples,
              mean_asv_prevalence = mean(rel_ASV_prevalence)
    )
  
  
  bb <- aa %>% tidyr::gather("variable", "value", prevalence, 
                      mean_abundance, rel_prevalence, 
                      n_samples, mean_asv_prevalence) %>% 
    unite("tmp" , Condition, variable) %>% 
    tidyr::spread(tmp, value) %>% 
    dplyr::mutate(total_n_samples = no_n_samples + yes_n_samples, 
           total_prevalence = yes_prevalence + no_prevalence,
           total_rel_prevalence = total_prevalence/total_n_samples, 
           mean_asv_prevalence = (yes_mean_asv_prevalence*yes_n_samples + no_mean_asv_prevalence*no_n_samples)/total_n_samples)
  if(outname != "") write_tsv(bb, file = outname)
  
  return(bb)
}

plotRelativeAbnBarsPhylum <- function(phobj, variable="Condition", 
                                      outname="phylumBarplot.pdf", 
                                      height=8, width=12, ocluster=F){
  library(RColorBrewer)
  ps_rel_abund = phyloseq::transform_sample_counts(phobj, function(x){x / sum(x)})
  facet_form <- as.formula(paste0(". ~ ", variable))
  
  
  g1 <- phyloseq::plot_bar(ps_rel_abund, fill = "Phylum") +
    geom_bar(aes(color = Phylum, 
                 fill = Phylum), 
             stat = "identity", position = "stack") +
    labs(x = "", y = "Relative Abundance\n") +
    facet_wrap(facet_form, scales = "free") +

    #scale_fill_uchicago() +
    #scale_color_uchicago() +
    #theme_pubclean() +
    mytheme + 
    scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
    theme(panel.background = element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank())
  
  df_prevalence <- getMeanRelAbundancesByPhylum(phobj)  %>% 
    dplyr::arrange(prevalence)
  g1$data$Phylum <- factor(g1$data$Phylum, levels = df_prevalence$Phylum) #Sort phylums from most to least prevalent
  
  #Sort samples according to prevalence of the most prevalent taxon or by clustering
  if(! ocluster){
    most_prevalent <- df_prevalence$Phylum[nrow(df_prevalence)]
    sorted_samples <- g1$data %>% filter(Phylum == most_prevalent) %>% 
      group_by(sampleID) %>% dplyr::summarise(Abundance = sum(Abundance)) %>% 
      dplyr::arrange(Abundance) %>% pull(sampleID)
    g1$data$Sample <- factor(g1$data$Sample, levels = sorted_samples)
  }else{
    dmat <- ps_rel_abund %>% otu_table() %>% t %>% dist
    hcl <- hclust(dmat)
    sorted_samples <- gsub("^X", "", hcl$labels[hcl$order])
    g1$data$Sample <- factor(g1$data$Sample, levels = sorted_samples)
  }
  
  mycolors <- colorRampPalette(brewer.pal(8, "Spectral"))(nrow(df_prevalence)) # %>% rev
  #mycolors[is.na(df_prevalence$Phylum) | df_prevalence$Phylum == "NA"] <- "#000000"
  
  g1 <- g1 + 
    scale_fill_manual(values = mycolors) +
    scale_color_manual(values = mycolors) 
  ggsave(outname, g1, width=width, height = height)
  return(g1)
}

#Shows the most prevalent genera (aggregating ASVs/OTUs)
plotRelativeAbnBarsGenus <- function(phobj, variable="Condition", topn = 15, 
                                     outname="phylumBarplot.pdf", height=8, 
                                     width=12, ocluster=T, 
                                     oldlevs=c("Control", "Depression")){
  library(RColorBrewer)

  df_prevalence <- getRelAbundancesByGenusAndVariable(phobj, variable, outname="", oldlevs=oldlevs)
  df_top <- df_prevalence %>% top_n(topn, total_rel_prevalence) %>% 
    filter(!is.na(Genus)) %>% 
    dplyr::arrange(desc(total_rel_prevalence))
  taxdata <- tax_table(phobj) %>% data.frame() %>% filter(Genus %in% df_top$Genus)
  
  gendata <- taxdata %>% group_by(Genus) %>% dplyr::summarise(Phylum=unique(Phylum),
                                                       Family=unique(Family))
  
  otus_gen <- data.frame(otu_table(phobj) ) %>% 
    rownames_to_column("ASV") %>% 
    dplyr::mutate(is_top = ASV %in% rownames(taxdata),
           Genus = taxdata$Genus[match(ASV, rownames(taxdata))]
           ) %>% 
    dplyr::mutate(Genus = ifelse(is.na(Genus), "Other", Genus)) %>%
    group_by(Genus) %>% 
    dplyr::summarise_if(is.numeric, sum) %>% 
    filter(Genus != "Other") %>% 
    dplyr::mutate(
      Phylum = gendata$Phylum[match(Genus, gendata$Genus)],
      Family = gendata$Family[match(Genus, gendata$Genus)],
    )
  
  otus_gen_gen_rel <- otus_gen %>% 
    dplyr::mutate_if(is.numeric, function(x)x/sum(x))
  
  otus_df <- otus_gen_gen_rel %>% 
    tidyr::gather(key="sampleID", value="Abundance", -Genus, -Phylum, -Family) %>% 
    dplyr::mutate(sampleID = gsub("^X", "", sampleID))
  
  ranked_genera <- otus_gen_gen_rel %>% 
    dplyr::mutate(mean_rel = rowMeans(across(where(is.numeric)))) %>% 
    #top_n(topn, mean_rel) %>% 
    dplyr::arrange(mean_rel) #desc(mean_rel)
  most_prev_genus <- ranked_genera$Genus[nrow(ranked_genera)]
  
  df_merged <- merge(otus_df, data.frame(sample_data(phobj)), by = "sampleID") %>% 
    dplyr::mutate(Genus = factor(Genus, levels = ranked_genera$Genus))
  
  if(! ocluster){
    ranked_samples <- df_merged %>% 
      filter(Genus == most_prev_genus) %>% 
      dplyr::arrange(Abundance) %>% 
      pull(sampleID)
  }else{
    dmat <- otus_gen_gen_rel %>% select_if(is.numeric) %>% t %>% dist
    hcl <- hclust(dmat)
    ranked_samples <- gsub("^X", "", hcl$labels[hcl$order])
  }
  
  df_merged <- df_merged %>% dplyr::mutate(sampleID = factor(sampleID, levels=ranked_samples))
  
  #cols <- brewer.pal(n=nrow(ranked_genera), name="Set2") # No more colors than those in palette
  mycolors <- colorRampPalette(brewer.pal(8, "Spectral"))(nrow(ranked_genera))
  facet_form <- paste0( ". ~ ", variable) %>% as.formula
  
  g1 <-ggplot(df_merged, aes(x=sampleID, y=Abundance, color = Genus, 
                             fill = Genus)) +
    geom_bar(stat = "identity", position = "stack") +
    labs(x = "", y = "Relative Abundance\n") +
    facet_wrap(facet_form, scales = "free") +
    #scale_fill_uchicago() +
    #scale_color_uchicago() +
    #theme_pubclean() +
    scale_fill_manual(values = mycolors) +
    scale_color_manual(values = mycolors) +
    mytheme +
    theme(panel.background = element_blank(),
          # axis.text.x=element_text(angle=40, vjust = -0.5, hjust=0.5)
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          legend.text=element_text(face="italic")
          ) +
    scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) 
    
  ggsave(outname, g1, width=width, height = height)
  return(g1)
}

##This one should be equivalent to the fantastix one, 
# Shows the most prevalent species/OTUs/ASVs but coloreda by genus
plotRelativeAbnBarsSpecies_ColByGenus <- function(phobj, variable="Condition", topn = 15, 
                                     outname="SpeciescolByGenusBarplot.pdf", height=8, 
                                     width=12, docluster=T){
  library(RColorBrewer)
  ps_rel_abund = phyloseq::transform_sample_counts(phobj, function(x){x / sum(x)})
  aa <- otu_table(ps_rel_abund)
  bb <- rowSums(aa) %>% sort(decreasing = T)
  top_asv <- names(bb)[1:topn]
  taxdata <- tax_table(phobj) %>% data.frame() %>% 
    rownames_to_column("ASV") %>% 
    dplyr::filter(ASV %in% top_asv) %>% 
    dplyr::arrange(Genus) 
  taxdata$Genus2 <- sapply(1:nrow(taxdata), FUN=function(i, lista){
    aux <- table(lista[1:i])
    aux2 <- table(lista)[names(aux)]
    aux <- ifelse(aux2 > 1, as.character(aux), "")
    paste0(lista[i], as.character(aux[lista[i]] ))
  }, taxdata$Genus)
  
  
  otus_gen <- data.frame(otu_table(phobj) ) %>% 
    rownames_to_column("ASV") %>% 
    dplyr::mutate(is_top = ASV %in% taxdata$ASV,
           Genus = taxdata$Genus2[match(ASV, taxdata$ASV)]
    ) %>% 
    dplyr::mutate(Genus = ifelse(is.na(Genus), "Other", Genus)) %>%
    # group_by(Genus) %>% 
    # summarise_if(is.numeric, sum) %>% 
    filter(Genus != "Other") %>% 
    dplyr::mutate(
      Phylum = taxdata$Phylum[match(ASV, taxdata$ASV)],
      Family = taxdata$Family[match(ASV, taxdata$ASV)],
    ) 
  
  otus_gen_gen_rel <- otus_gen %>% 
    dplyr::mutate_if(is.numeric, function(x)x/sum(x))
  
  otus_df <- otus_gen_gen_rel %>% 
    tidyr::gather(key="sampleID", value="Abundance", -ASV, -Genus, -Phylum, -Family, -is_top) %>% 
    dplyr::mutate(sampleID = gsub("^X", "", sampleID))
  
  
  df_merged <- merge(otus_df, data.frame(sample_data(phobj)), by = "sampleID") 
  
  if(! ocluster){
    ranked_genera <- otus_gen_gen_rel %>% 
      dplyr::mutate(mean_rel = rowMeans(across(where(is.numeric)))) %>% 
      top_n(topn, mean_rel) %>% dplyr::arrange(mean_rel) #desc(mean_rel)
    most_prev_genus <- ranked_genera$Genus[nrow(ranked_genera)]
    
    ranked_samples <- df_merged %>% 
      filter(Genus == most_prev_genus) %>% 
      dplyr::arrange(Abundance) %>% 
      pull(sampleID)
  }else{
    # dmat <- otus_gen_gen_rel 
    # rownames(dmat) <- dmat$Genus
    # dmat <- dmat %>% 
    #   select_if(is.numeric) %>% dist
    # hcl <- hclust(dmat)
    # ranked_genera <-hcl$labels[hcl$order]
    ranked_genera <- taxdata$Genus2
  
    dmat <- otus_gen_gen_rel %>% select_if(is.numeric) %>% t %>% dist
    hcl <- hclust(dmat)
    ranked_samples <- gsub("^X", "", hcl$labels[hcl$order])
  }
  
  df_merged <- df_merged %>% dplyr::mutate(sampleID = factor(sampleID, levels=ranked_samples),
                                    Genus = factor(Genus, levels = ranked_genera)
                                    )
  
  #cols <- brewer.pal(n=nrow(ranked_genera), name="Set2") # No more colors than those in palette
  mycolors <- colorRampPalette(brewer.pal(8, "Spectral"))(length(ranked_genera))
  g1 <-ggplot(df_merged, aes(x=sampleID, y=Abundance, color = Genus, 
                             fill = Genus)) +
    geom_bar(stat = "identity", position = "stack") +
    labs(x = "", y = "Relative Abundance\n") +
    facet_wrap(facet_form, scales = "free") +
    #scale_fill_uchicago() +
    #scale_color_uchicago() +
    #theme_pubclean() +
    scale_fill_manual(values = mycolors) +
    scale_color_manual(values = mycolors) +
    mytheme +
    theme(panel.background = element_blank(),
          # axis.text.x=element_text(angle=40, vjust = -0.5, hjust=0.5)
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          legend.text=element_text(face="italic")
          ) +
    scale_y_continuous(limits = c(0, 1), expand = c(0, 0))
  ggsave(outname, g1, width=width, height = height)
  return(g1)
}

plotRelativeAbnBars_Fantaxtic <- function(phobj, variable="Condition", topn = 15, 
                                               tax_level = "Genus",
                                     outname="GenusBarplotFx.pdf", height=7, width=12){
  library(fantaxtic)
  topntx <- get_top_taxa(physeq_obj = phobj, n = topn, relative = T,
                        discard_other = T, other_label = "Other")
  
  topntx <- name_taxa(topntx, label = "", species = F, other_label = "Other")
  topntx <- fantaxtic_bar(topntx, color_by = tax_level, label_by = tax_level, 
                        facet_by = variable, grid_by = NULL, 
                        other_color = "Grey") +
    mytheme +
    theme(axis.text.x = element_text(size = 10, 
                                     colour = "black", angle = 90, 
                                     face = "plain", hjust=1, vjust=1))
    #theme(strip =element_rect(fill="white"))+
  ggsave(filename = outname, topntx, height = height, width = width) 
  return(topntx)
}


plotPrevalenceVsAbundance <- function(phobj, outname="phylumBarplot.pdf", height=10, width=12){
  pre_prevalence <- getRelAbundanceTab(phobj)
  g1 <- ggplot(pre_prevalence, aes(TotalAbundance, relative_prevalence, color = Phylum)) +
    # Agregamos una línea para nuestro umbral
    geom_hline(yintercept = 0.05, alpha = 0.5, linetype = 2) +  geom_point(size = 2, alpha = 0.7) +
    scale_x_log10() +  xlab("Total Abundance") +
    ylab("Prevalence/num. samples") +
    facet_wrap(~Phylum) +
    theme_pubclean() +
    guides(color = FALSE)
  ggsave(outname, g1, width=width, height = height)
  return(g1)
}

plotPhylumBoxplots <- function(phobj, var="Condition", outname="phylumBarplot.pdf", height=10, width=8, paired=F){
  ## Total abundance by phylum, apparently by summing over all ASVs in phylum
  ps_phylum <- phyloseq::tax_glom(phobj, "Phylum")
  phyloseq::taxa_names(ps_phylum) <- phyloseq::tax_table(ps_phylum)[, "Phylum"]
  
  sample_data(ps_phylum)[, var] <- unlist(sample_data(ps_phylum)[, var]) %>% as.character
  comps <- combn(unique(unlist(sample_data(ps_phylum)[, var])), 2, simplify = F)
  signif_codes <- list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, Inf), symbols = c("****", "***", "**", "*", "ns"))
  
  g1 <- phyloseq::psmelt(ps_phylum) %>%
    ggplot(data = ., aes_string(x = var, y = "Abundance")) +
    facet_wrap(~ OTU, scales = "free")+ 
    geom_boxplot(outlier.shape  = NA) +
    geom_jitter(aes(color = OTU), height = 0, width = .2) +
    stat_compare_means(method="wilcox.test", comparisons = comps, 
                       symnum.args = signif_codes, paired=paired) +
    labs(x = var, y = "Abundance\n") +
    theme_pubclean() + 
    guides(color = FALSE)
  if(paired){
    g1 <- g1 + geom_line(aes(group = pacienteID), color = "gray20", linetype=1, size=0.1, alpha=0.25)
  }
  ggsave(outname, g1, width=width, height = height)
  return(g1)
}

getPhylumTests <- function(phobj, var="Condition", outname="phylumTests", paired=F){
  ## Total abundance by phylum, apparently by summing over all ASVs in phylum
  ps_phylum <- phyloseq::tax_glom(phobj, "Phylum")
  phyloseq::taxa_names(ps_phylum) <- phyloseq::tax_table(ps_phylum)[, "Phylum"]
  
  comps <- combn(unique(unlist(sample_data(ps_phylum)[, var])), 2, simplify = F)
  signif_codes <- list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, Inf), symbols = c("****", "***", "**", "*", "ns"))
  
  dfmelt <- phyloseq::psmelt(ps_phylum) %>% 
    group_by(Phylum)
  
  pvals <- dfmelt %>% 
    dplyr::group_modify(~ broom::tidy(wilcox.test(data=., as.formula(paste0( "Abundance ~", var)), paired = paired)))
  pvals$padj = p.adjust(pvals$p.value, method="BH")
  write_tsv(pvals, file=outname)
  return(pvals)
}

clusterBrayCurtis <- function(phobj, variable="Condition", renamestr="", 
                              dist_type="bray", clust_method="ward.D2"){
  ps_rel_abund = phyloseq::transform_sample_counts(phobj, function(x){x / sum(x)})
  ps_rel_otu <- data.frame(phyloseq::otu_table(ps_rel_abund))
  ps_rel_otu <- t(ps_rel_otu)
  rownames(ps_rel_otu) <- gsub("^X", renamestr, rownames(ps_rel_otu))
  #bc_dist <- vegan::vegdist(ps_rel_otu, method = dist_type)
  bc_dist <- phyloseq::distance(physeq = phobj, method = dist_type)
  
  #Save as dendrogram
  ward <- as.dendrogram(hclust(bc_dist, method = clust_method))
  #Provide color codes
  meta <- data.frame(phyloseq::sample_data(ps_rel_abund))
  levv <- unique(meta[, variable])
  colorCode <- c( "firebrick3", "dodgerblue3")
  names(colorCode) <- levv
  labels_colors(ward) <- colorCode[meta[, variable]][order.dendrogram(ward)]
  
  hmplot <- plot_heatmap(ps_rel_abund, sample.label=variable)
  return(list("ward"=ward, "hmplot"=hmplot))
}

testHMP <- function(phobj, seed=123){
  ps_phylum <- phyloseq::tax_glom(pre_phyloseq_filt, "Phylum")
  controls <- phyloseq::subset_samples(ps_phylum, Condition == "Control")
  cf <- phyloseq::subset_samples(ps_phylum, Condition == "Depression")
  #Output OTU tables
  control_otu <- data.frame(phyloseq::otu_table(controls))
  cf_otu <- data.frame(phyloseq::otu_table(cf))
  
  group_data <- list(control_otu, cf_otu)
  (xdc <- HMP::Xdc.sevsample(group_data))  
}

######################## PERMANOVA
adonis2table <- function(mod1){
  aux <- data.frame(
    variable = rownames(mod1)[1],
    DF_var=mod1$Df[1],
    DF_Residual = mod1$Df[2],
    DF_Total = mod1$Df[3],
    SumOfSQs_var = mod1$SumOfSqs[1],
    SumOfSQs_Residual =mod1$SumOfSqs[2],
    SumOfSQs_Total = mod1$SumOfSqs[3],
    R2_var = mod1$R2[1],
    R2_Residual = mod1$R2[2],
    R2_Total = mod1$R2[3],
    F_statistic = mod1$F[1],
    P = mod1$`Pr(>F)`[1]
  )
  return(aux)
}

makePermanova <- function(phobj, dist_method = "bray", seed = 123, 
                          exclude_vars = c("sampleID"), outname = "permanovas.tsv"){
  ## From https://deneflab.github.io/MicrobeMiseq/demos/mothur_2_phyloseq.html#permanova
  library(stringi)
  
  meta <- sample_data(phobj)
  braydist <- phyloseq::distance(phobj, method = dist_method)
  sampledf <- data.frame(sample_data(phobj))
  names(sampledf) <- stri_trans_general(str = names(sampledf), id = "Latin-ASCII") %>% 
    gsub(" ", "_", .)
  
  vars2test <- names(sampledf)[! names(sampledf) %in% exclude_vars]
  
  res <- data.frame()
  for(var in vars2test){
    set.seed(seed)
    form <- paste0("braydist ~ ", var) %>% as.formula
    if(length(unique(sampledf[!is.na(sampledf[, var]) , var])) > 1){
      # Adonis test
      #cat(var, "\n")
      mod1 <- adonis2(form, data = sampledf, na.action=na.exclude)
      res <-rbind(res, adonis2table(mod1))
    } 
  }
  res <- res %>% dplyr::arrange(P)
  res$padj <- p.adjust(res$P, method="BH")
  write_tsv(res, file=outname)
  return(res)
  
}

makeBetaDispTests<- function(phobj, dist_method = "bray", 
                  exclude_vars = c("sampleID"), 
                  outname){
  
  braydist <- phyloseq::distance(phobj, method = dist_method)
  sampledf <- data.frame(sample_data(phobj))
  vars2test <- names(sampledf)[! names(sampledf) %in% exclude_vars]
  res <- data.frame()
  for(var in vars2test){
    beta <- betadisper(braydist, sampledf$Psoriasis)
    betaper <- permutest(beta)
    aux <- data.frame(variable = var, 
                      F = betaper$tab$F[1],
                      p = betaper$tab$`Pr(>F)`[1])
    res <- rbind(res, aux)
  }
  write_tsv(res, file=outname)
  return(res)
}

makePermanovaSeveralFactors <- function(phobj, 
                                        dist_method = "bray",
                                        seed = 123, 
                                        modlist = c(), 
                                        outname = "permanovas_mult.RData"){
  ## From https://deneflab.github.io/MicrobeMiseq/demos/mothur_2_phyloseq.html#permanova
  library(stringi)
  
  meta <- sample_data(phobj)
  braydist <- phyloseq::distance(phobj, method = dist_method)
  sampledf <- data.frame(sample_data(phobj))
  
  res <- tibble()
  for(vars in modlist){
    set.seed(seed)
    form <- paste0("braydist ~ ", paste(vars, sep= " + ", collapse = " + "))
    mod1 <- adonis2(as.formula(form), data = sampledf, na.action=na.exclude)
    aux <- tibble(formula = form, model = list(mod1), explained=1 - mod1["Residual", "R2"])
    res <-rbind(res, aux)
  }
  save(res, file=outname)
  return(res)
  
}


##################################################################
##################################################################
### DESEQ

getDeseqContrastFromCategorical <- function(dds, lev_combin, opt, name){
  contrastvec <- c(lev_combin[1], lev_combin[3], lev_combin[2])
  contrast_name <- paste0(lev_combin[1],'_', lev_combin[3], '_vs_', lev_combin[2]) %>% gsub(" ", ".", .)
  res <- results(dds, contrast = contrastvec)
  resLFC <- lfcShrink(dds, contrast = contrastvec, type="normal", lfcThreshold = log2(opt$fc)) #apeglm gives weird results
  resLFC_ape <- tryCatch(lfcShrink(dds, coef = contrast_name, type="apeglm", lfcThreshold = log2(opt$fc)),error=\(x)data.frame() ) #apeglm gives weird results
  resLFC_ashr <- lfcShrink(dds, contrast = contrastvec, type="ashr", lfcThreshold = log2(opt$fc)) #apeglm gives weird results
  resdf <- defWriteDEAResults(res, resLFC, opt, paste0(name, "_",contrast_name, "_DAAshrinkNormal.tsv"))
  resdf_ape <- defWriteDEAResults(res, resLFC_ape, opt, paste0(name, "_", contrast_name, "_DAAshrinkApe.tsv"))
  resdf_ashr <- defWriteDEAResults(res, resLFC_ashr, opt, paste0(name, "_", contrast_name, "_DAAshrinkAshr.tsv"))
  
  return(list("res"=res,
              "resLFC"=resLFC,
              "resLFC_ape"=resLFC_ape, 
              "resLFC_ashr"=resLFC_ashr,
              "resdf"=resdf,
              "resdf_ape"=resdf_ape,
              "resdf_shr"=resdf_ashr, 
              contrast_vec = contrastvec,
              contrast_name = contrast_name))
}

getDeseqContrastFromNumerical <- function(dds, nvarname, opt, name){
  res <- results(dds, name = nvarname)
  resLFC <- lfcShrink(dds, coef = nvarname, type="normal", lfcThreshold = log2(opt$fc)) #apeglm gives weird results
  resLFC_ape <- tryCatch(lfcShrink(dds, coef = nvarname, type="apeglm", lfcThreshold = log2(opt$fc)),error=\(x)data.frame() ) #apeglm gives weird results
  resLFC_ashr <- lfcShrink(dds, coef = nvarname, type="ashr", lfcThreshold = log2(opt$fc)) #apeglm gives weird results
  resdf <- defWriteDEAResults(res, resLFC, opt, paste0(name, "_",nvarname, "_DAAshrinkNormal.tsv"))
  resdf_ape <- defWriteDEAResults(res, resLFC_ape, opt, paste0(name, "_", nvarname, "_DAAshrinkApe.tsv"))
  resdf_ashr <- defWriteDEAResults(res, resLFC_ashr, opt, paste0(name, "_", nvarname, "_DAAshrinkAshr.tsv"))
  
  return(list("res"=res,
              "resLFC"=resLFC,
              "resLFC_ape"=resLFC_ape, 
              "resLFC_ashr"=resLFC_ashr,
              "resdf"=resdf,
              "resdf_ape"=resdf_ape,
              "resdf_shr"=resdf_ashr,
              "nvarname"=nvarname))
}

getDeseqResults <- function(phobj, opt, name="", variables = c("Condition")){
  formula <- paste0("~ ", paste(variables, sep=" + ", collapse=" + ")) %>% 
    as.formula
  dds <- phyloseq_to_deseq2(phobj, design= formula)
  raw_counts <- counts(dds)
  if(opt$minsampleswithcount == 0){ 
    dds <- dds[rowSums(counts(dds)) >= opt$mincount,] 
  }else{
    dds <- dds[rowSums(counts(dds) >= opt$mincount) >= opt$minsampleswithcount,] 
  }  
  filt_counts <- counts(dds)
  
  ##Add pseudocount if necessary
  anyNonZero <- raw_counts %>% apply(MAR=1, all) %>% any
  if(!anyNonZero){
    do_poscounts = TRUE
    dds <- DESeq(dds, betaPrior = F, sfType = "poscounts")
  }else{
    do_poscounts = FALSE
    dds <- DESeq(dds, betaPrior = F)
  }
  
  write_file(paste(resultsNames(dds), collapse="\t" ), file=paste0(opt$out, name, "_", "DEA_resultsNames.tsv"))

  design <- dds@colData %>% as.data.frame()
  all_combos_done <- TRUE
  
  all_combins <- map(variables, \(x){
    if(is.numeric(design[, x])){
      return(list(c(x, "NUMERIC")))
    }
    levs <- levels(design[, x] %>% unlist) 
    combins <- lapply(combn(1:length(levs), 2, simplify = F), \(y)c(x, levs[y]))
  }) %>% flatten

  all_contrasts <- map(all_combins, \(lev_combin){
    if(lev_combin[2] == "NUMERIC"){
      getDeseqContrastFromNumerical(dds, lev_combin[1], opt, name)
    }else{
      getDeseqContrastFromCategorical(dds, lev_combin, opt, name)
    }
  })  
  names(all_contrasts) <- lapply(all_combins, \(x)ifelse(x[2]=="NUMERIC", x[1], paste0(x[1],'_', x[3], '_vs_', x[2]) %>% gsub(" ", ".", .)))
  all_combos_done <- TRUE
  # Write raw counts  
  rawc_df <- defWriteMatAsDF(raw_counts, opt, paste0(name, "_", "raw_counts.tsv") )
  #filtx_df <- defWriteMatAsDF(filt_counts, opt, "raw_counts_filtered.tsv") 

  # Normalized counts  
  norm_counts <- counts(dds, normalized = T)
  norm_counts_df <- defWriteMatAsDF(norm_counts, opt, paste0(name, "_", "norm_counts.tsv") )
  
  tryCatch({
    vstds <- varianceStabilizingTransformation(raw_counts, blind=F)
    vst_counts_df <- defWriteMatAsDF(vstds, opt, paste0(name, "_", "vst_counts.tsv") )
  }, error= function(x){
    vstds <<- NULL
    vst_counts_df <<- data.frame()
  })
  
  all_results_list <- list(
    "dds"=dds, 
    "raw_counts"=raw_counts,
    "all_contrasts"=all_contrasts, 
    "res"=all_contrasts[[1]]$res,
    "resLFC"=all_contrasts[[1]]$resLFC,
    "resLFC_ape"=all_contrasts[[1]]$resLFC_ape, 
    "resLFC_ashr"=all_contrasts[[1]]$resLFC_ashr,
    "resdf"=all_contrasts[[1]]$resdf,
    "resdf_ape"=all_contrasts[[1]]$resdf_ape,
    "resdf_shr"=all_contrasts[[1]]$resdf_ashr,
    "raw_df" = rawc_df, 
    "norm_counts"=norm_counts,
    "norm_counts_df"=norm_counts_df,
    "vstds"=vstds, 
    "vst_counts_df"=vst_counts_df,
    "all_combos_done"=all_combos_done,
    "options"= list(mincount=opt$mincount, minsampleswithcount=opt$minsampleswithcount, 
                    minfreq=opt$minfreq, poscount=do_poscounts)
  )
  save(file = paste0(opt$out, "DESEQ2_all_results_", name, ".R"), all_results_list)
  return(all_results_list)
}

getDeseqResults_paired <- function(phobj, opt, name="", variables = c("Condition"), individual = "pacienteID"){
  formula <- paste0("~ ", paste(variables, sep=" + ", collapse=" + "), " + ",
                          paste(individual, sep=" + ", collapse=" + ")) %>% 
                          as.formula
  reduced_formula = paste0("~ ", paste(individual, sep=" + ", collapse=" + ")) %>% 
    as.formula
  
  counts_table <- otu_table(phobj) %>% data.frame %>% as.matrix
  colnames(counts_table) <- gsub("^X", "", colnames(counts_table))
  design <- sample_data(phobj) %>% data.frame
  design <- design[colnames(counts_table) ,]
  
  dds <- DESeqDataSetFromMatrix(countData = counts_table, design = formula, colData = design)
  
  #dds <- phyloseq_to_deseq2(phobj, design= formula)
  raw_counts <- counts(dds)
  if(opt$minsampleswithcount == 0){ 
    dds <- dds[rowSums(counts(dds)) >= opt$mincount,] 
  }else{
    dds <- dds[rowSums(counts(dds) >= opt$mincount) >= opt$minsampleswithcount,] 
  }  
  filt_counts <- counts(dds)
  
  ##Add pseudocount if necessary
  anyNonZero <- raw_counts %>% apply(MAR=1, all) %>% any
  if(!anyNonZero){
    do_poscounts = TRUE
    dds_null <- DESeq(dds, betaPrior = F, sfType = "poscounts")
    dds <- DESeq(dds, betaPrior = F, sfType = "poscounts", reduced=reduced_formula, test="LRT")
  }else{
    do_poscounts = FALSE
    dds_null <- DESeq(dds, betaPrior = F)
    dds <- DESeq(dds, betaPrior = F, reduced=reduced_formula, test="LRT")
  }
  
  levelscond <- levels(design[, variables[1]])
  res <- tryCatch(results(dds, contrast=c(variables[1], levelscond[2], levelscond[1])),
                  error=results(dds, contrast=list(paste0(variables[1], "_" ,levelscond[2], "_vs_",levelscond[1]) ))
                  )
  res_null <- tryCatch(results(dds_null, contrast=c(variables[1], levelscond[2], levelscond[1])),
                       error=results(dds_null, contrast=list(paste0(variables[1], "_" ,levelscond[2], "_vs_",levelscond[1]) ))
  )
  
  resLFC <- lfcShrink(dds, coef=resultsNames(dds)[2], type="normal", lfcThreshold = log2(opt$fc)) #apeglm gives weird results
  resLFC_null <- lfcShrink(dds_null, coef=resultsNames(dds)[2], type="normal", lfcThreshold = log2(opt$fc)) #apeglm gives weird results
  resLFC_ape <- lfcShrink(dds, coef=resultsNames(dds)[2], type="apeglm", lfcThreshold = log2(opt$fc)) #apeglm gives weird results
  resLFC_ashr <- lfcShrink(dds, coef=resultsNames(dds)[2], type="ashr", lfcThreshold = log2(opt$fc)) #apeglm gives weird results
  resdf <- defWriteDEAResults(res, resLFC, opt, paste0(name, "_", "DEA_results_shrinkNormal.tsv"))
  resLFC_null <- defWriteDEAResults(res, resLFC, opt, paste0(name, "_", "DEA_results_shrinkNormal_againstNullModel.tsv"))
  resdf_ape <- defWriteDEAResults(res, resLFC_ape, opt, paste0(name, "_", "DEA_results_shrinkApe.tsv"))
  resdf_ashr <- defWriteDEAResults(res, resLFC_ashr, opt, paste0(name, "_", "DEA_results_shrinkAshr.tsv"))
  
  # Write raw counts  
  rawc_df <- defWriteMatAsDF(raw_counts, opt, paste0(name, "_", "raw_counts.tsv") )
  #filtx_df <- defWriteMatAsDF(filt_counts, opt, "raw_counts_filtered.tsv") 
  
  # Normalized counts  
  #  dds <- estimateSizeFactors(dds)##  already done
  norm_counts <- counts(dds, normalized = T)
  norm_counts_df <- defWriteMatAsDF(norm_counts, opt, paste0(name, "_", "norm_counts.tsv") )
  
  tryCatch({
    vstds <- varianceStabilizingTransformation(raw_counts, blind=F)
    vst_counts_df <- defWriteMatAsDF(vstds, opt, paste0(name, "_", "vst_counts.tsv") )
  }, error= function(x){
    vstds <<- NULL
    vst_counts_df <<- data.frame()
  })
  
  all_results_list <- list(
    "dds"=dds, 
    "dds_null"=dds_null,
    "raw_counts"=raw_counts,
    "res"=res,
    "res_null"=res_null,
    "resLFC"=resLFC,
    "resLFC_ape"=resLFC_ape, 
    "resLFC_ashr"=resLFC_ashr,
    "resdf"=resdf,
    "resdf_ape"=resdf_ape,
    "resdf_shr"=resdf_ashr,
    "raw_df" = rawc_df, 
    "norm_counts"=norm_counts,
    "norm_counts_df"=norm_counts_df,
    "vstds"=vstds, 
    "vst_counts_df"=vst_counts_df,
    "options"= list(mincount=opt$mincount, minsampleswithcount=opt$minsampleswithcount, 
                    minfreq=opt$minfreq, poscount=do_poscounts)
  )
  save(file = paste0(opt$out, "DESEQ2_all_results_", name, ".R"), all_results_list)
  return(all_results_list)
}

defWriteMatAsDF <- function(mat, opt, name){
  fname <-paste(opt$out, name, sep="/", collapse="/") 
  df <- mat %>% as.data.frame(row.names = rownames(.)) %>% 
    rownames_to_column("taxon")
  write_tsv(df, fname)
  return(df)
}

defWriteDEAResults <- function(res, resLFC, opt, name){
  fname <-paste(opt$out, name, sep="/", collapse="/") 
  resdf <- res %>% as.data.frame(row.names = rownames(.)) %>% 
    rownames_to_column("taxon") %>% 
    dplyr::mutate(log2FoldChangeShrink = resLFC$log2FoldChange,
           lfcSE_Shrink = resLFC$lfcSE, svalue = resLFC$svalue)
  write_tsv(resdf, fname)
  return(resdf)
}
plotCountsDeseq <- function(dds, raw_counts, norm_counts){ 
  d <- plotCounts(dds, gene=which.min(res$padj), intgroup="Condition", normalized = T,
                  returnData=TRUE)
  d <- raw_counts %>% colSums() %>% as.data.frame()
  d$sample <- rownames(d)
  names(d)[1] <- "count"
  d$Condition <- colData(dds)$Condition[match( d$sample, colData(dds)$sampleID)]  
  
  g1 <- ggplot(d, aes(x=Condition, y=count, col=Condition, fill=Condition)) + 
    geom_point(position=position_jitter(w=0.1,h=0), size=2) + 
    #  scale_y_log10(breaks=c(25,100,400)) + 
    geom_text_repel(aes(label=sample)) +
    ggtitle("Total raw counts by sample") +
    theme_bw()
  
  d2 <- norm_counts %>% colSums() %>% as.data.frame()
  d2$sample <- rownames(d)
  names(d2)[1] <- "count"
  d2$Condition <- colData(dds)$Condition[match( d2$sample, colData(dds)$sampleID)]  
  
  g2 <- ggplot(d2, aes(x=Condition, y=count, col=Condition, fill=Condition)) + 
    geom_point(position=position_jitter(w=0.1,h=0), size=2) + 
    #  scale_y_log10(breaks=c(25,100,400)) + 
    geom_text_repel(aes(label=sample)) +
    ggtitle("Total normalized counts by sample") +
    theme_bw()
  
  nc_exp <- raw_counts %>% as.data.frame() %>%
    dplyr::mutate(gene=rownames(raw_counts)) %>%
    gather(key=sample, value='norm_counts', colnames(norm_counts)) %>% 
    filter(norm_counts > 0)
  
  g3 <- ggplot(nc_exp, aes(x=norm_counts, col=sample))+
    geom_density() +
    xlim(x=c(0,  quantile(nc_exp$norm_counts, c(0.95))[1] ))+
    ggtitle("Distribution of normalized counts by sample") +
    theme_bw() +
    scale_x_continuous(trans='log2') +
    guides(color=FALSE)
  
  grid.arrange(g1, g2, g3, nrow = 1)
  
} 
make_maplot<- function(res, opt, name="maplot.pdf"){
  outname <- paste(opt$out, name, sep="/", collapse="/")
  pdf(outname)
  ma <- DESeq2::plotMA(res, ylim=c(-3, 3), alpha=opt$pval)
  abline(h = c(-1*log2(opt$fc), log2(opt$fc)), col = "dodgerblue3", lwd =2, lty =2)
  tmp <- dev.off()  
  DESeq2::plotMA(res, ylim=c(-3, 3), alpha=opt$pval)
  abline(h = c(-1*log2(opt$fc), log2(opt$fc)), col = "dodgerblue3", lwd =2, lty =2)
  
} 
make_volcano <- function(res, opt, name, pcol="pvalue"){
  outname <- paste(opt$out, name, sep="/", collapse="/")
  pdf(outname, width = 12, height = 8)
  volc <- NULL
  try({ 
    volc <- EnhancedVolcano(res,
                            lab = rownames(res),
                            x = 'log2FoldChange', title = "", subtitle="",
                            y = pcol,
                            pCutoff = opt$pval,
                            pCutoffCol = pcol,
                            FCcutoff = log2(opt$fc),
                            legendPosition = 'right',
                            legendLabSize = 12,
                            legendIconSize = 3.0,
                            legendLabels=c('Not sig.',
                                           paste("FC > ", as.character(opt$fc), sep="", collapse=""),
                                           paste(pcol, " < ", as.character(opt$pval), sep="", collapse=""),
                                           paste("FC > ", as.character(opt$fc), " & ", pcol, " < ", as.character(opt$pval), sep="", collapse=""))
    )
  })  
  print(volc)
  tmp <- dev.off()
  volc
} 

defWriteMatAsDF <- function(mat, opt, name){
  fname <-paste(opt$out, name, sep="/", collapse="/") 
  df <- mat %>% as.data.frame(row.names = rownames(.)) %>% 
    rownames_to_column("gene")
  write_tsv(df, fname)
  return(df)
}

makeHeatmap <- function(resdf, dds, df2plot, 
                        variable = "condition",
                        opt, 
                        name = "heatmap.pdf", 
                        logscale=FALSE, 
                        ptype = "padj", w=5, h=4, 
                        trim_values=FALSE,
                        italics_rownames=TRUE, taxalist=c()){
  outname <- paste(opt$out, name, sep="/", collapse="/")
  annot <- as.data.frame(colData(dds)[variable])
  names(annot) <- c(variable)
  rownames(annot) <- colData(dds)$sampleID
  
  if(length(taxalist)==0){
    if(ptype == "padj"){
      taxa <- resdf %>% filter(padj <= opt$pval & 
                             !is.na(padj) &
                             abs(log2FoldChangeShrink) >= log2(opt$fc) ) %>% 
      pull(taxon)
    }else{
      taxa <- resdf %>% filter(pvalue <= opt$pval & 
                               abs(log2FoldChangeShrink) >= log2(opt$fc) ) %>% 
                        pull(taxon)
    }
    if(length(taxa) < .GlobalEnv$opt$num_genes_default){
      #get only first n genes  
      taxa <- resdf[order(resdf$pvalue), ]  %>% 
        head(.GlobalEnv$opt$num_genes_default) %>% 
        pull(taxon)
    } 
  }else{
    taxa <- resdf$taxon[resdf$taxon %in% taxalist]
  }
  
  mat <- df2plot %>%  
    dplyr::filter(gene %in% taxa) %>% 
    column_to_rownames("gene") %>% 
    as.matrix %>% 
    t %>% scale %>% t
  if(logscale){
    mat <- log(mat + 1)
  }
  
  fontsize_row = 10 - nrow(mat) / 15
  
  mat_trim <- mat
  if(trim_values){
    mat_trim[mat> quantile(mat, 0.99)] <- quantile(mat, 0.99)
    mat_trim[mat < quantile(mat, 0.01)] <- quantile(mat, 0.01)
  }
  
  labels_row <- gsub("_", " ", rownames(mat_trim))
  if(italics_rownames){
    labels_row <- lapply(labels_row, \(x){
      bquote(italic(.(x)))
    }) %>% as.expression()
  }
  hm <- pheatmap(mat_trim, cluster_rows=T, 
                 show_rownames=nrow(mat) < 120,
                 annotation_col = annot,
                 fontsize_row = fontsize_row,
                 border_color = NA,
                 labels_row = labels_row,
                 cluster_cols=T, 
                 annotcluster_rowsation_col=annot
                 )
  w <- w+0.05*ncol(mat)
  h <- h+0.05*nrow(mat)
  pdf(outname, width = w, height = h)
  print(hm)
  tmp <- dev.off()
}


getPropVar <- function(pca, PC){
  xx <- summary(pca)
  percc<- round(100*xx$importance["Proportion of Variance" , PC],2) %>% as.character
  labb <- paste(PC, " (", percc, "%)", sep="", collapse="")
  return(labb)
}

getSignASVBoxplot_paired<-function(resdf_annot, raw_counts, outdir="./", name="sign_asvs_boxplot.pdf", mode="PROP", max_n=10){
  ASVs <- resdf_annot %>% dplyr::arrange(pvalue) %>% filter(pvalue < 0.05)  %>% unite("taxon2", taxon, Genus, sep=":", remove=F)
  if(nrow(ASVs) > max_n){ASVs <- ASVs[1:10,]}
  
  if(mode=="PROP"){
    prop_counts <- apply(raw_counts, MAR=2,function(x)x/sum(x))
  }else{
    prop_counts <- raw_counts
  }
  df2plot <- prop_counts[ASVs$taxon, ] %>% 
    as.data.frame %>% rownames_to_column("ASV") %>%
    dplyr::mutate(Genus = ASVs$Genus,taxon=ASVs$taxon2) %>% 
    tidyr::gather("sampleID", "raw_counts", -ASV, -Genus, -taxon ) %>% 
    dplyr::mutate(patient = metadata$pacienteID[match(sampleID, metadata$sampleID)],
                  Condition = metadata$Condition[match(sampleID, metadata$sampleID)])
  df2plot$taxon=gsub(":$", "", df2plot$taxon)
  
  g1 <- ggplot(df2plot, aes(x = Condition, y=raw_counts, col=Condition, fill=Condition)) + 
    facet_wrap(. ~ taxon , scales = "free")+
    geom_violin(alpha=0.1)+
    geom_boxplot(width=0.1, fill="darkgray", notchwidth = 0.5, notch=F)+
    geom_point(size=1)+
    geom_line(aes(group = patient), color = "gray20", linetype=1, size=0.1, alpha=0.25) +
    # scale_color_lancet()+
    # scale_fill_lancet() +
    theme_pubclean() +
    ylab("Proportion")
  
  oname <- paste0(outdir, "/", name)
  ggsave(filename=oname, g1, width=6, height = 8)
  return(g1)
  
}

getPropVar_princomp <- function(pca){
  pcts <- round(100*pca$sdev^2/sum(pca$sdev^2), 2) %>% as.character
  nnames <- gsub("Comp.", "PC", names(pca$sdev))
  labb <- paste(nnames, " (", pcts, "%)", sep="")
  names(labb)<- nnames
  return(labb)
}

plotPCA <- function(counts, design, genes, plotvar="Condition"){
  
  mat <- counts %>%  
    dplyr::filter(gene %in% genes) %>% 
    column_to_rownames("gene") %>% 
    dplyr::select(design$sampleID) %>% 
    as.matrix %>% 
    t %>% scale
  
  if(nrow(mat) < 3){
    cw <- getNullPlot()
    return(cw)
  }
  pca <- prcomp(mat)
  df <- pca$x %>% as.data.frame %>% 
    rownames_to_column("sample")
  df$group <- design[match(df$sample, design$sampleID), plotvar] %>% unlist
  
  legname <- gsub("_", " ", plotvar)
  legname <- gsub("\\.", "-", legname)
  
  g1 <- ggplot(df, aes(x=PC1, y=PC2, col=group)) +
    geom_point() + 
    geom_text_repel(aes(label = sample)) + 
    xlab(getPropVar(pca, "PC1")) +
    ylab(getPropVar(pca, "PC2")) +
    mystyle +
    guides(col=guide_legend(title=legname))
  g2 <- ggplot(df, aes(x=PC1, y=PC3, col=group)) +
    geom_point() + 
    geom_text_repel(aes(label = sample)) + 
    xlab(getPropVar(pca, "PC1")) +
    ylab(getPropVar(pca, "PC3")) +
    mystyle +
    guides(col=guide_legend(title=legname))
  g3 <- ggplot(df, aes(x=PC1, y=PC4, col=group)) +
    geom_point() + 
    geom_text_repel(aes(label = sample)) + 
    xlab(getPropVar(pca, "PC1")) +
    ylab(getPropVar(pca, "PC4")) +
    mystyle +
    guides(col=guide_legend(title=legname))
  g4 <- ggplot(df, aes(x=PC1, y=PC5, col=group)) +
    geom_point() + 
    geom_text_repel(aes(label = sample)) + 
    xlab(getPropVar(pca, "PC1")) +
    ylab(getPropVar(pca, "PC5")) +
    mystyle +
    guides(col=guide_legend(title=legname))
  cw <- cowplot::plot_grid(g1, g2, g3, g4)
  
  return(list(pca=pca, plots=cw))
}

getNullPlot <- function(opt, name="var", error=FALSE){
  if(!error){
    emptylab <- paste("Error plotting ", name, sep="", collapse="")
  }else{
    emptylab = "Error while plotting"
  }
  
  gempty <- ggplot() + 
    geom_text(aes(x=1, y=1),label=emptylab, size=14)+ 
    theme_light() +
    theme(axis.text = element_blank(), axis.title = element_blank())
  return(gempty)
}

makeAllPCAs <- function(phobj, counts_df, genes, vars2pca, opt, name = "PCAs"){
  design <- sample_data(phobj)
  design<- data.frame(design)
  nreads <- otu_table(phobj) %>% colSums()

  design <- design %>% 
    dplyr::mutate(reads_log10_current = log10(nreads[as.character(sampleID)]))
  
  vars2pca <- c(vars2pca,"reads_log10_current")
  
  pca_plots <- lapply(vars2pca,  FUN=function(vv, counts, design, genes){
     plotPCA(counts,design, genes, vv)
  },counts_df, design, genes)
  names(pca_plots) <- vars2pca
  
  save(pca_plots, file = paste0(opt$out, name, "PCA.RData"))
  write_tsv(as.data.frame(pca_plots[[1]]$pca$x), file = paste0(opt$out, name, "PCAnewvars.tsv"))
  write_tsv(as.data.frame(pca_plots[[1]]$pca$rotation), file = paste0(opt$out, name, "PCArotation.tsv"))
  
  pdf(paste0(opt$out, name, ".pdf"), width=16, height = 12)
  for(pl in vars2pca){
    tryCatch({print(pca_plots[[pl]][["plots"]])}, error = function(x){print(getNullPlot(opt, pl))})
  }
  tmp <- dev.off()

  return(pca_plots)
}


plotPCA2 <- function(dfraw, vars2pc, labvar = "PACIENTE", plotvars=c("Condition"), 
                    transform = "scale", dims=2:5, name = "PCAs.pdf", outdir="", w=12, h=8){
  
  dfraw <- dfraw %>% column_to_rownames(labvar) 
  mat <- dfraw[, vars2pc] %>%  
    as.matrix
  
  if(transform=="scale"){
    mat <- mat %>% scale %>% t
  }else if(transform == "log"){
    mat <- log(mat + 1) %>% t
  }else if(transform == "sqrt"){
    mat <- mat %>% sqrt %>% t
  }else{
    mat <- mat %>% t
  }
  
  dfraw$sample <- rownames(dfraw)
  pca <- prcomp(mat)
  df <- pca$rotation %>% as.data.frame %>% 
    rownames_to_column("sample")
  df2 <- merge(df, dfraw, by="sample")
  
  plotlist <- list()
  for(plotvar in plotvars){
    
    legname <- gsub("_", " ", plotvar)
    legname <- gsub("\\.", "-", legname)
    gs <- lapply(paste0("PC", dims), FUN=function(PC){
      g1 <- ggplot(df2, aes_string(x="PC1", y=PC, col=plotvar)) +
        geom_point() + 
        geom_text_repel(aes(label = sample)) + 
        xlab(getPropVar(pca, "PC1")) +
        ylab(getPropVar(pca, PC)) +
        mystyle +
        guides(col=guide_legend(title=legname))
    # if(! is.numeric(df2[, plotvar])){
    #   g1 <- g1 + scale_color_lancet() 
    # }
    return(g1)
    })
  plotlist[[plotvar]] <- cowplot::plot_grid(plotlist=gs)
  }
  
  WriteManyPlots(plotlist,name = name, outdir = outdir, w = w, h=h)
  dev.off()
  return(list(pca=pca, plots=plotlist))
}

getGTTableFromRes <- function(res, genes, name){
  gttable <- res %>% as.data.frame() %>% 
    rownames_to_column("Taxa") %>% 
    filter(Taxa %in% genes) %>% 
    gt() %>% 
    tab_caption(name) %>% 
    data_color(
      method = "numeric",
      palette = c("firebrick3", "dodgerblue2"), 
      columns = `pvalue` 
    ) %>% 
    data_color(
      method = "numeric",
      palette = c("firebrick3", "dodgerblue2"), 
      columns = `padj` 
    ) %>% 
    data_color(
      method = "numeric",
      palette = c("dodgerblue2", "firebrick3"), 
      columns = `log2FoldChange` 
    ) %>% 
    tab_options(
      table.background.color = "white",
      column_labels.background.color = "lightgray",
      column_labels.font.size = px(16),
      table.font.size = px(12),
      data_row.padding = px(16),
      table.width = px(600)
    )
  return(gttable)
}



getSummaryTablesDeseq <- function(res, opt){
  namestab <- c(paste("p < ", as.character(opt$pval), sep="", collapse="") ,
                paste("LFC > ", as.character(log2(opt$fc)), sep="", collapse="")
  )
  restab <- res %>% as.data.frame() %>% 
    dplyr::mutate(a = pvalue <= opt$pval, 
           b = ifelse(log2FoldChange <= -log2(opt$fc),"less frequent",
                      ifelse( log2FoldChange >= log2(opt$fc), "more frequent", "equal"))
    ) %>% 
    dplyr::select(a:b) %>% 
    set_names(namestab) %>% 
    table() 
  rownames(restab) <- ifelse(rownames(restab) == "TRUE", 
                             paste("p <= ", as.character(opt$pval), sep="", collapse=""), 
                             paste("p > ", as.character(opt$pval), sep="", collapse=""))
  
  #restab %>% kable(caption="Number of taxons(species) per category using raw p-values")
  
  restab_adj <- res %>% as.data.frame() %>% 
    dplyr::mutate(a = padj <= opt$pval, 
           b = ifelse(log2FoldChange <= -log2(opt$fc),"less frequent",
                      ifelse( log2FoldChange >= log2(opt$fc), "more frequent", "equal"))
    ) %>% 
    dplyr::select(a:b) %>% 
    set_names(namestab) %>% 
    table() 
  rownames(restab_adj) <- ifelse(rownames(restab_adj) == "TRUE", 
                                 paste("p <= ", as.character(opt$pval), sep="", collapse=""), 
                                 paste("p > ", as.character(opt$pval), sep="", collapse=""))
  
  #restab_adj %>% kable(caption="Number of taxons (species) per category using adjusted p-values")  
  write_tsv(as.data.frame(restab), paste0(opt$out, "/num_diff_rawpval.tsv"))
  write_tsv(as.data.frame(restab_adj), paste0(opt$out, "/num_diff_adjpval.tsv"))
  return(list("restab"=restab, "restab_adj"=restab_adj))
}

getCorrelations <- function(mat1, mat2){
  res <- data.frame()
  sp_pmat <-matrix(nrow=ncol(mat1), ncol=ncol(mat2))
  colnames(sp_pmat) <- colnames(mat2)
  rownames(sp_pmat) <- colnames(mat1)
  pr_pmat <-sp_pmat
  
  for(na in colnames(mat1)){
    for(nb in colnames(mat2)){
      ct <- cor.test(mat1[, na], mat2[, nb], method="pearson")
      ct2 <- cor.test(mat1[, na], mat2[, nb], method="spearman")
      sp1 <- shapiro.test(mat1[, na])
      sp2 <- shapiro.test(mat2[, nb])
      mod <- lm(mat1[, na] ~ mat2[, nb]) %>% summary() 
      
      aux <- data.frame(
        var1 = na,
        var2 = nb,
        pearson_cor = ct$estimate,
        pearson_pval = ct$p.value,
        spearman_cor = ct2$estimate,
        spearman_pval = ct2$p.value,
        shapiro_v1 = sp1$p.value,
        shapiro_v2 = sp2$p.value,
        rsquared = mod$r.squared,
        slope_pval = mod$coefficients[2, 4]
      )
      res <- rbind(res, aux)
      #Store results also in matrices
      sp_pmat[na, nb] <- ct2$p.value
      pr_pmat[na, nb] <- ct$p.value
    }
  }
  res$pearson_adj_pval <- p.adjust(res$pearson_pval, method="BH")
  res$spearman_adj_pval <- p.adjust(res$spearman_pval, method="BH")
  
  spcors <- cor(mat1, mat2, method="spearman")
  pcors <-  cor(mat1, mat2, method="pearson")
  return(list(res=res, spcors=spcors, pcors=pcors, sp_pvals = sp_pmat, pr_pmat = pr_pmat))
}

getCorsByGroup <- function(mat1, mat2, group1, group2, outdir, name, levnames=c("PSO", "CTRL")){
  
  cors <- getCorrelations(mat1, mat2)
  cors_pso <- getCorrelations(mat1[group1,], mat2[group1,])
  cors_ct <- getCorrelations(mat1[group2,], mat2[group2,])
  
  cors_pso_df <- cors_pso[[1]] %>% dplyr::select(-var1, -var2)
  names(cors_pso_df) <- paste(levnames[1], names(cors_pso_df), sep="_")
  
  cors_ct_df <- cors_ct[[1]] %>% dplyr::select(-var1, -var2)
  names(cors_ct_df) <- paste(levnames[2], names(cors_ct_df), sep="_")
  
  
  cordf <- cbind(cors[[1]], cors_pso_df, cors_ct_df) 
  
  oname <- paste0(outdir, name)
  write_tsv(cordf, file = oname)
  return(list(alldf=cordf, all=cors, gr1=cors_pso, gr2=cors_ct))
}

getCLRdata <- function(phobj, prevalence_lim = 0.05){
  library(CoDaSeq)
  #1) Imputate zeros
  # cmultRepl expects samples in rows and taxa in columns
  tax_matrix <- getZCompositionImputedData(phobj, prevalence_lim)
  
  #2) CLR transform
  tax_matrix_clr <- codaSeq.clr(tax_matrix, samples.by.row=F, IQLR=TRUE)
  
  return(tax_matrix_clr)
}

getZCompositionImputedData <- function(phobj, prevalence_lim = 0.05){
  library(zCompositions)
  otutab <- otu_table(phobj)
  #Umbral de prevalencia
  prevalence <- apply(otutab, MAR=1, FUN=function(x)sum(x>0)/length(x))
  otutab_filt <- otutab[prevalence > prevalence_lim ,]
  
  # CLR transform as in Xia, Sun & Xen book, pag 350
  
  #1) Imputate zeros
  # cmultRepl expects samples in rows and taxa in columns
  tax_matrix <- zCompositions::cmultRepl(X = t(otutab_filt), output="p-counts") %>% t 
  
  return(tax_matrix)
}

plotRegressionsASV_vs_vars<- function(df3, asvs, vars2cor, groupvar="Psoriasis", 
                                      group_levels = c("no", "yes"), 
                                      groupnames = c("Control", "Psoriasis"),
                                      xlabel = "log(pg/mL)",
                                      outdir = "", name = "regr.pdf", w=7, h=14, opt=list()){
  library(ggpmisc)
  df4 <- df3 %>% gather("IL", "log_val", vars2cor) 
  df4[, groupvar] <- factor(df4[, groupvar], levels=group_levels)
  df4[, "IL"] <- factor(df4$IL, levels=vars2cor)
  
  df5 <- df4 %>% filter(df4[, metaname]==group_levels[1])
  df6 <-  df4 %>% filter(df4[, metaname]==group_levels[2])
  
  plotlist <- list()
  for(ASV in asvs){
    gno<-ggplot(df5, aes_string(x="log_val", y=ASV)) + 
      facet_wrap( IL ~ . , scales="free", ncol=1) + 
      geom_point(col="dodgerblue3", alpha=1, size=0.5) + 
      geom_smooth(method="lm", col="dodgerblue3") +
      stat_poly_eq(use_label(c("R2", "p", "n")), #c("eq", "R2", "f", "p", "n")
                   method="lm", small.p=T, small.r=F, label.y=0.99, 
                   col="dodgerblue3") +
      theme_pubr() + 
      #scale_color_lancet() +
      ylab(ASV) +
      xlab(xlabel)+
      ggtitle(groupnames[1])
    
    
    gyes<-ggplot(df6, aes_string(x="log_val", y=ASV))+ 
      facet_wrap( IL ~ . , scales="free", ncol=1) + 
      geom_smooth(method="lm", col="firebrick3") +
      geom_point( col="firebrick3", alpha=1, size=0.5) + 
      stat_poly_eq(use_label(c("R2", "p", "n")), #c("eq", "R2", "f", "p", "n")
                   method="lm", small.p=T, small.r=F, label.y=0.99, 
                   col="firebrick3") +
      theme_pubr() + 
      #scale_color_lancet() +
      ylab(ASV)+
      xlab(xlabel) +
      ggtitle(groupnames[2])
    
    plotlist[[ASV]] <- cowplot::plot_grid(gno, gyes)
    
  }
  WriteManyPlots(plotlist, name, outdir, w=6, h=12, separate = F, opt=opt)
  return(plotlist)
}

plotRegressionsASV_vs_vars2<- function(df3, corrdf_annot, vars2cor, groupvar="Psoriasis", 
                                      group_levels = c("no", "yes"), 
                                      top_n=10,
                                      groupnames = c("Control", "Psoriasis"),
                                      ylabel = "log(pg/mL)",
                                      outdir = "", name = "regr.pdf", w=7, h=14, opt=list()){
  library(ggpmisc)
  df4 <- df3 %>% gather("ASV", "Abundance", matches("^ASV")) 
  df4[, groupvar] <- factor(df4[, groupvar], levels=group_levels)
  
  corrdf_annot <- corrdf_annot %>% 
    dplyr::mutate(var1_g = paste0(sapply(Genus, FUN=function(x){
      gsub("\\[|\\]", "", strsplit(x, " ")[[1]][1], perl=T)
      }), ":", var1))
  df4$ASV_g <- corrdf_annot$var1_g[match(df4$ASV, corrdf_annot$var1)]
  plotlist <- list()
  for(var2cor in vars2cor){
    best_asvs <- corrdf_annot %>% filter(var2==var2cor) %>% 
      arrange(dplyr::desc(PSO_spearman_cor + CTRL_spearman_cor)) %>% 
      head(top_n) %>% pull(var1_g)
    df4b <- df4 %>%  filter(ASV_g %in% best_asvs) %>% dplyr::mutate(ASV_g = factor(ASV_g, levels=best_asvs))
    df5 <- df4b %>% filter(df4b[, metaname]==group_levels[1])
    df6 <-  df4b %>% filter(df4b[, metaname]==group_levels[2])
    
    gno<-ggplot(df5, aes_string(x="Abundance", y=paste0("`", var2cor, "`"))) + 
      facet_wrap( ASV_g ~ . , scales="free", ncol=1) + 
      geom_point(col="dodgerblue3", alpha=1, size=0.5) + 
      geom_smooth(method="lm", col="dodgerblue3") +
      stat_poly_eq(use_label(c("R2", "p", "n")), #c("eq", "R2", "f", "p", "n")
                   method="lm", small.p=T, small.r=F, label.y=0.99, 
                   col="dodgerblue3") +
      theme_pubr() + 
      #scale_color_lancet() +
      ylab(ylabel) +
      xlab("CLR-transformed Abundance")+
      ggtitle(groupnames[1])
    
    
    gyes<-ggplot(df6, aes_string(x="Abundance", y=paste0("`", var2cor, "`")))+ 
      facet_wrap( ASV_g ~ . , scales="free", ncol=1) + 
      geom_smooth(method="lm", col="firebrick3") +
      geom_point( col="firebrick3", alpha=1, size=0.5) + 
      stat_poly_eq(use_label(c("R2", "p", "n")), #c("eq", "R2", "f", "p", "n")
                   method="lm", small.p=T, small.r=F, label.y=0.99, 
                   col="firebrick3") +
      theme_pubr() + 
      #scale_color_lancet() +
      ylab(ylabel)+
      xlab("CLR-transformed Abundance") +
      ggtitle(groupnames[2])
    
    plotlist[[var2cor]] <- cowplot::plot_grid(gno, gyes)
    
  }
  WriteManyPlots(plotlist, name, outdir, w=w, h=h, separate = F, opt=opt)
  return(plotlist)
}


calc_scImpute <- function(otutab){
  xx <- otutab %>% as.data.frame()
  write.table(xx, file="test_input.csv", sep=",", quote=F, row.names = T)
  library(scImpute)
  outlier_samples <- scimpute(# full path to raw count matrix
    count_path ="test_input.csv", 
    infile = "csv",           # format of input file
    outfile = "csv",          # format of output file
    out_dir = "./",           # full path to output directory
    labeled = FALSE,          # cell type labels not available
    drop_thre = 0.5,          # threshold set on dropout probability
    Kcluster = 2,             # 2 cell subpopulations
    ncores = 4) 
  scimputed <- read.csv("scimpute_count.csv")
  return(list(outliers = outlier_samples, imputed_counts=scimputed))
}


plotPCsCorwithASVs <- function(tab, contrres, num_model=1, component="Comp.2"){
  asvs <- strsplit(contrres$asvs_p05[num_model], "_")[[1]]
  if(length(asvs)>10) asvs <- asvs[1:10]
                   
  tab2 <- tab %>% filter(OTU %in% asvs) %>% 
    group_by(OTU) %>% 
    dplyr::mutate(Abundance2= log( (Abundance + 1)/exp(mean(log(Abundance+1))) ),
           Genus2 = sapply(Genus, FUN=function(x){y<-strsplit(x, " ")[[1]]; return(y[length(y)])} ),
           OTU2 = paste(OTU, Genus2, Psoriasis, sep=":"))
  
  g0 <- ggplot(tab2 , aes_string(x=component, y="Abundance2", col="Psoriasis", fill="Psoriasis"))+
    facet_wrap(OTU2 ~ ., scales = "free", ncol=4) + 
    geom_point() + 
    geom_smooth(method="lm") +
    stat_poly_eq(use_label(c("R2", "p")), #c("eq", "R2", "f", "p", "n")
                 method="lm", small.p=T, small.r=F, label.y=0.99) +
    theme_pubr() + 
    #scale_color_lancet() +
    #scale_fill_lancet() +
    ylab("CLR Abundance")
  

  return(g0)
}

plotPCsCorwithASVs_PSO <- function(tab, contrres, num_model=1){
  asvs <- strsplit(contrres$asvs_p05[num_model], "_")[[1]]
  #if(length(asvs)>10) asvs <- asvs[1:10]
  
  tab2 <- tab %>% filter(OTU %in% asvs) %>% 
    group_by(OTU) %>% 
    dplyr::mutate(Abundance2= log( (Abundance + 1)/exp(mean(log(Abundance+1))) ),
           Genus2 = sapply(Genus, FUN=function(x){y<-strsplit(x, " ")[[1]]; return(y[length(y)])} ),
           OTU2 = paste(OTU, Genus2, Psoriasis, sep=":"))
  
  g0 <- ggplot(tab2 , aes_string(x="Psoriasis", y="Abundance2", col="Psoriasis"))+
    facet_wrap(OTU ~ ., scales = "free", ncol=4) + 
    geom_boxplot() +
    geom_jitter(alpha=0.6) +
    #geom_signif(test="wilcox.test", na.rm=T, comparisons = c("yes", "no"), 
                # step_increase=0.03,
                # tip_length = 0.01,
                # vjust=0.4,
                # color = "black"
    #) +
    theme_pubr() + 
    #scale_color_lancet() +
    #scale_fill_lancet() +
    ylab("CLR Abundance")
  
  
  return(g0)
}

basicLimma <- function(phobj, opt, variables =  c("Condition"), individual="pacienteID", outdir="", name="basic_limma"){
  formula <- paste0("~ 0 +", paste(variables, sep=" + ", collapse=" + ")) %>% 
    as.formula
  #formula_ind <- paste0("~ ", paste(variables, sep=" + ", collapse=" + "), " + (1|", individual, ")") %>% 
  as.formula
  count_mat <- otu_table(phobj) %>% data.frame %>% as.matrix()
  count_mat <- count_mat[rowSums(count_mat  > opt$mincount) > opt$minsampleswithcount, ] 
  
  # Standard usage of limma/voom
  genes = DGEList( count_mat )
  genes = calcNormFactors( genes )
  
  snames <- gsub("^X", "", colnames(count_mat))
  metadata <- metadata[snames, ]
  design = model.matrix(formula, metadata)
  vobj_tmp = voom( genes, design, plot=FALSE)
  
  # apply duplicateCorrelation 
  dupcor <- duplicateCorrelation(vobj_tmp,design,block=metadata[, individual])
  
  # run voom considering the duplicateCorrelation results
  # in order to compute more accurate precision weights
  # Otherwise, use the results from the first voom run
  vobj = voom( genes, design, plot=FALSE, block=metadata[, individual], correlation=dupcor$consensus)
  
  # But this step uses only the genome-wide average for the random effect
  fitDupCor <- lmFit(vobj, design, block=metadata$Individual, correlation=dupcor$consensus)
  
  # Fit Empirical Bayes for moderated t-statistics
  #fitDupCor_bay <- eBayes( fitDupCor )
  contr <- makeContrasts(ConditionAfter -ConditionBefore, levels =  colnames(coef(fitDupCor)))
  tmp <- contrasts.fit(fitDupCor, contr)
  tmp <- eBayes(tmp)
  result <- topTable(tmp, sort.by = "P", n = Inf)
  
  write.table(result, file=paste0(outdir, name, ".tsv"))
  reslist <- list(result=result, vobj=vobj, vobj_tmp = vobj_tmp, fitmm=tmp, contrast=contr)
  oname <- paste0(outdir, "/", name, ".RData")
  save(reslist, file=oname)
  return(reslist)
}

dreamLimma_old2 <- function(phobj, opt, variables =  c("Condition"), individual="pacienteID"){
  formula <- paste0("~ 0 +", paste(variables, sep=" + ", collapse=" + ")) %>% 
    as.formula
  formula_ind <- paste0("~ ", paste(variables, sep=" + ", collapse=" + "), " + (1|", individual, ")") %>% 
    as.formula
  count_mat <- otu_table(phobj) %>% data.frame %>% as.matrix()
  count_mat <- count_mat[rowSums(count_mat  > opt$mincount) > opt$minsampleswithcount, ] 
  
  # Standard usage of limma/voom
  genes = DGEList( count_mat )
  genes = calcNormFactors( genes )
  
  snames <- gsub("^X", "", colnames(count_mat))
  metadata <- metadata[snames, ]
  design = model.matrix(formula, metadata)
  vobj_tmp = voom( genes, design, plot=FALSE)
  dupcor <- duplicateCorrelation(vobj_tmp,design,block=metadata[, individual])
  vobj = voom( genes, design, plot=FALSE, block=metadata[, individual], correlation=dupcor$consensus)
  ##DREAM
  cl <- makeCluster(4)
  registerDoParallel(cl)
  
  # Get the contrast matrix for the hypothesis test
  L = getContrast( dupcor, formula_ind, metadata,  
                   colnames(vobj_tmp$design)[2])
  fitmm = dream( vobj, form, metadata, L)
  fitmm_ebayes = eBayes( fitmm )
  
  # get results
  result <- topTable(fitmm_ebayes, sort.by = "P", n = Inf)
  return(list(result=result, vobj_tmp=vobj_tmp, fitmm=fitmm, L=L))
}
dreamLimma <- function(phobj, opt, variables =  c("Condition"), individual="pacienteID", outdir = "", name = "dream"){
  formula <- paste0("~ 0 +", paste(variables, sep=" + ", collapse=" + ")) %>% 
    as.formula
  formula_ind <- paste0("~ ", paste(variables, sep=" + ", collapse=" + "), " + (1|", individual, ")") %>% 
    as.formula
  count_mat <- otu_table(phobj) %>% data.frame %>% as.matrix()
  count_mat <- count_mat[rowSums(count_mat  > opt$mincount) > opt$minsampleswithcount, ] 
  colnames(count_mat)<- gsub("^X", "", colnames(count_mat))
  # Standard usage of limma/voom
  genes = DGEList( count_mat )
  genes = calcNormFactors( genes )
  
  metadata <- metadata[colnames(count_mat), ]
  #design = model.matrix(formula_ind, metadata)
  vobj_tmp = voomWithDreamWeights(genes, formula_ind, metadata)
  
  contrvec <- paste(variables[1] , levels(metadata[, variables[1]]), sep="")
  L = getContrast( vobj_tmp, formula_ind, metadata,  
                   contrvec[2])
  
  #fitmm_1 = dream(vobj_tmp, formula_ind, metadata)
  fitmm_dream = dream(vobj_tmp, formula_ind, metadata, L)
  #result_dream1 <- variancePartition::topTable(fitmm_1, sort.by = "P", n = Inf, coef=contrvec[2])
  result_dream <- variancePartition::topTable(fitmm_dream, sort.by = "P", n = Inf, coef=contrvec[2])
  
  write.table(result_dream, file=paste0(outdir, "/", name, "_contrast.tsv"))
  #write.table(result_dream1, file=paste0(outdir, "/", name, ".tsv"))
  
  reslist <- list(result=result_dream, vobj_tmp=vobj_tmp, fitmm=fitmm_dream, L=L)
  oname <- paste0(outdir, "/" ,name, ".RData")
  save(reslist, file=oname)
  return(reslist)
}

dreamLimma_old <- function(phobj, opt, variables =  c("Condition"), individual="pacienteID", outdir = "", name = "dream"){
  formula <- paste0("~ ", paste(variables, sep=" + ", collapse=" + ")) %>% 
    as.formula
  formula_ind <- paste0("~ ", paste(variables, sep=" + ", collapse=" + "), " + (1|", individual, ")") %>% 
    as.formula
  
  count_mat <- otu_table(phobj) %>% data.frame %>% as.matrix()
  count_mat <- count_mat[rowSums(count_mat  > opt$mincount) > opt$minsampleswithcount, ] 
  colnames(count_mat)<- gsub("^X", "", colnames(count_mat))
  # Standard usage of limma/voom
  genes = DGEList( count_mat )
  genes = calcNormFactors( genes )
  
  #snames <- gsub("^X", "", colnames(count_mat))
  metadata <- metadata[colnames(count_mat), ]
  design = model.matrix(formula, metadata)
  vobj_tmp = voom( genes, design, plot=FALSE)
  #colnames(vobj_tmp) <- gsub("^X", "", colnames(vobj_tmp))
  dupcor <- duplicateCorrelation(vobj_tmp,design,block=metadata[, individual])
  vobj = voom( genes, design, plot=FALSE, block=metadata[, individual], correlation=dupcor$consensus)
  ##DREAM
  cl <- makeCluster(4)
  registerDoParallel(cl)
  
  # Get the contrast matrix for the hypothesis test
  L = getContrast( dupcor, formula_ind, metadata,  
                   colnames(vobj_tmp$design)[2])
  fitmm = dream( vobj, formula_ind, metadata, L)
  fitmm_ebayes = eBayes( fitmm )
  
  
  # get results
  result <- variancePartition::topTable(fitmm_ebayes, sort.by = "P", n = Inf)
  
  write.table(result, file=paste0(outdir, "/", name, ".tsv"))
  reslist <- list(result=result, vobj_tmp=vobj_tmp, fitmm=fitmm, L=L)
  oname <- paste0(outdir, name, ".RData")
  save(reslist, file=oname)
  return(reslist)
}



makeCCAfromMatrixVegan <- function(datamat, tax_matrix_clr, pcvar2retain=0.9, outdir = "", name="CCA_from_CLR"){
  #Assumes matrices have the same number of rows (subjects) and that are ordered in the same way
  txpca <- prcomp(tax_matrix_clr %>% t)
  txsum <- txpca %>% summary()
  pc2use <- colnames(txsum$importance)[txsum$importance["Cumulative Proportion",] > pcvar2retain][1]
  txpcdata <- txpca$rotation[, 1:(which(colnames(txpca$rotation) == pc2use))]
  ccares <- CCorA(datamat, txpcdata)
  
  save(ccares, file = paste0(outdir, name, ".RData"))
  pdf(paste0(outdir, name, "_BiplotVegan.pdf"), width = 12, height = 12)
  biplot(ccares)
  dev.off()
  return(ccares)
}

makeCCA <- function(datamat, tax_matrix_clr, metadata,metadatavar="Psoriasis",
                    pcvar2retain=0.9, outdir = "", name="CCA_from_CLR"){
  library(C)
  #Assumes matrices have the same number of rows (subjects) and that are ordered in the same way
  txpca <- prcomp(tax_matrix_clr %>% t)
  txsum <- txpca %>% summary()
  pc2use <- colnames(txsum$importance)[txsum$importance["Cumulative Proportion",] > pcvar2retain][1]
  txpcdata <- txpca$rotation[, 1:(which(colnames(txpca$rotation) == pc2use))]
  
  datamat <- cbind(datamat, ifelse(metadata[rownames(datamat), "Psoriasis"] == "yes" , 10, 0))
  colnames(datamat)[ncol(datamat)] <- "Psoriasis"
  #ccares <- CCorA(datamat, txpcdata)
  #ccares <- cca(txpcdata, datamat)
  ccares <- cancor(txpcdata, datamat)
  
  CC1_X <- txpcdata %*% ccares$xcoef[, 1]
  CC1_Y <- datamat %*% ccares$ycoef[, 1]
  
  CC2_X <- txpcdata %*% ccares$xcoef[, 2]
  CC2_Y <- datamat %*% ccares$ycoef[, 2]
  
  cca_df <- data.frame(CC1_X=scale(CC1_X), 
                       CC1_Y=scale(CC1_Y), 
                       CC2_X=scale(CC2_X), 
                       CC2_Y=scale(CC2_Y)) %>% 
    rownames_to_column("sampleID") 
  cca_df[, metadatavar] <- metadata[match(cca_df$sampleID, metadata$sampleID), metadatavar]
  arrows_y = data.frame(x=0, 
                        xend=ccares$ycoef[,1]*2, 
                        y=0, 
                        yend=ccares$ycoef[,2]*2) %>% 
    rownames_to_column("var2name")
  arrows_x = data.frame(x=0, 
                        xend=ccares$xcoef[,1]*2, 
                        y=0, 
                        yend=ccares$xcoef[,2]*2) %>% 
    rownames_to_column("var2name")
  arrows_x <- arrows_x[1:5,]
  
  g1 <- ggplot(cca_df, aes_string(x=metadatavar, y="CC1_X", col=metadatavar)) + 
    geom_boxplot()+
     mytheme 
    # scale_color_lancet()
  g2 <- ggplot(cca_df, aes_string(x=metadatavar, y="CC1_Y", col=metadatavar)) + 
    geom_boxplot()+
    mytheme 
    #scale_color_lancet()
  # Y points and Y variables
  g3 <- ggplot(cca_df, aes_string(x="CC1_Y", y="CC2_Y", col=metadatavar, fill=metadatavar))+
    geom_point()+
    geom_segment(aes(x=x, xend=xend, y=y, yend=yend),
                 data= arrows_y, inherit.aes = F, arrow = arrow(length=unit(0.3, "cm")))+
    geom_text_repel(aes(x=xend, y=yend, label=var2name), data=arrows_y,
                    inherit.aes = F)+
    stat_ellipse()+
    mytheme 
    #scale_color_lancet()
  # X points and Y variables
  g4 <- ggplot(cca_df, aes_string(x="CC1_X", y="CC2_X", col=metadatavar, fill=metadatavar))+
    geom_point()+
    geom_segment(aes(x=x, xend=xend, y=y, yend=yend),
                 data= arrows_y, inherit.aes = F, arrow = arrow(length=unit(0.3, "cm")))+
    geom_text_repel(aes(x=xend, y=yend, label=var2name), data=arrows_y,
                    inherit.aes = F)+
    stat_ellipse()+
    mytheme 
    #scale_color_lancet()
  # Y points and X variables
  g5 <- ggplot(cca_df, aes_string(x="CC1_Y", y="CC2_Y", col=metadatavar, fill=metadatavar))+
    geom_point()+
    geom_segment(aes(x=x, xend=xend, y=y, yend=yend),
                 data= arrows_x, inherit.aes = F, arrow = arrow(length=unit(0.3, "cm")))+
    geom_text_repel(aes(x=xend, y=yend, label=var2name), data=arrows_x,
                    inherit.aes = F)+
    stat_ellipse()+
    mytheme 
    #scale_color_lancet()
  # X points and X variables
  g6 <- ggplot(cca_df, aes_string(x="CC1_X", y="CC2_X", col=metadatavar, fill=metadatavar))+
    geom_point()+
    geom_segment(aes(x=x, xend=xend, y=y, yend=yend),
                 data= arrows_x, inherit.aes = F, arrow = arrow(length=unit(0.3, "cm")))+
    geom_text_repel(aes(x=xend, y=yend, label=var2name), data=arrows_x,
                    inherit.aes = F)+
    stat_ellipse()+
    mytheme 
    #scale_color_lancet()
  
  save(ccares, file = paste0(outdir, name, ".RData"))
  pdf(paste0(outdir, "/", name, "_plots.pdf"), width = 12, height = 12)
  g1
  g2
  g3
  g4
  g5
  g6
  dev.off()
  return(list(ccares=ccares, g1=g1, g2=g2, g3=g3, g4=g4, g5=g5, g6=g6))
}

getSignif2 <- function(vector, limits=c(0.05, 0.01, 0.001), 
                       legends = c("", "*", "**", "***")){
  x <- ifelse(vector < limits[3], legends[4],
              ifelse(vector < limits[2], legends[3],
                     ifelse(vector< limits[1], legends[2], legends[1])
              )
  )
  return(x)
}

makeCorrelationHeatmap <- function(mat1, mat2, cormethod="pearson",
                                   pval=0.05, select_asvs=c(), outdir="", name="corrHeatmap", 
                                   clust=T, order_rows = c(), order_cols = c()){
  # select_asvs: plot these ASVs. Overrides pval
  # pval: plot ASVs with any correlation >  than this threshold
  cormat <- cor(mat1, mat2, method=cormethod)
  cor_pval <- expand.grid(colnames(mat1), colnames(mat2)) %>% 
    rowwise() %>% 
    mutate(pval = cor.test(mat1[, Var1], mat2[, Var2], method=cormethod)$p.value) %>% 
    ungroup() %>% 
    tidyr::spread(Var2, pval) %>% 
    column_to_rownames("Var1") %>% 
    as.matrix()
  cor_pval <- cor_pval[rownames(cormat), colnames(cormat)]
  

  if(length(select_asvs) > 0){
    cormat <- cormat[select_asvs, ]
    cor_pval <- cor_pval[select_asvs, ]
  }else if(pval < 1.0){
    select_asvs <- cor_pval %>% apply(MAR=1, function(x) any(x < pval)) %>% which() %>% names
    cormat <- cormat[select_asvs, ]
    cor_pval <- cor_pval[select_asvs, ]
  }
  cor_ptext <- getSignif2(cor_pval)
  if(! clust){
    cormat <- cormat[order_rows, order_cols]
    cor_ptext <- cor_ptext[order_rows, order_cols]
    cor_pval <- cor_pval[order_rows, order_cols]
  }
  
  cormat <- cormat %>% t
  cor_ptext <- cor_ptext %>% t
  cor_pval <- cor_pval %>% t
  
  oname <- paste0(outdir, "/", name, ".pdf")
  fontsize_row = 16 - nrow(cormat) / 15
  fontsize_col = 16 - ncol(cormat) / 15
  
  pheatmap(cormat, display_numbers = cor_ptext, filename=oname,
           width = 20, height = 6,
           fontsize_row = fontsize_row,
           fontsize_number = 16,
           fontsize_col = fontsize_col,
           cluster_rows = clust, cluster_cols = clust
           )
  tmp <- tryCatch(dev.off(), error=function(x){})
  hm <- pheatmap(cormat, display_numbers = cor_ptext,
           filename=oname,
           width = 20, height = 6,
           fontsize_row = fontsize_row,
           fontsize_number = 16,
           fontsize_col = fontsize_col,
           cluster_rows = clust, cluster_cols = clust
           )
  return(list(hm=hm, asvs=select_asvs, cormat=cormat, 
              pvals=cor_pval, cor_ptext=cor_ptext))
}


makeCorrelationHeatmapByGroup <- function(mat1, mat2, metadata,var2plot="Psoriasis",
                                          varwithnames = "sampleID",
                                          cormethod="pearson",
                                          pval=0.05, select_asvs=c(),
                                          outdir="", name="corrHeatmap", 
                                          clust=T, order_rows = c(), 
                                          order_cols = c(),
                                          annot_asvs = NULL){
  
  hmall <- makeCorrelationHeatmap(mat1, mat2, cormethod=cormethod,
                                  pval=0.05, select_asvs=select_asvs, outdir=outdir, 
                                  name=paste0(name, "_allsamples"), clust=T)
  default_annot_colors <- c("gray40", "dodgerblue3", "firebrick3", "green4", "pink4", "skyblue")
  levs <- levels(metadata[, var2plot])
  hmlevs <- lapply(levs, FUN=function(lev){
    samples2select <- metadata[ metadata[, var2plot]==lev, varwithnames]
    mat1_tmp <- mat1[rownames(mat1) %in% samples2select, ]
    mat2_tmp <- mat2[rownames(mat2) %in% samples2select, ]
    hmlist <- makeCorrelationHeatmap(mat1_tmp, mat2_tmp, cormethod=cormethod,
                                  pval=0.05, select_asvs=hmall$asvs, outdir=outdir, 
                                  name=paste0(name,"_", var2plot,"_", lev ), clust=T)
    return(hmlist)
  })
  names(hmlevs)<-levs
  
  hm1 <- hmall$hm
  asvorder <- hm1$tree_col$labels[hm1$tree_col$order]
  varorder <- hm1$tree_row$labels[hm1$tree_row$order]
  
  newcormat <- hmall$cormat[varorder, asvorder]
  newptext <- hmall$cor_ptext[varorder, asvorder]
  annot_rows <- rep("All samples", nrow(newcormat))
  vspaces <- c()
  for(ll in levs){
    vspaces <- c(vspaces, nrow(newcormat))
    newcormat <- rbind(newcormat, hmlevs[[ll]]$cormat[varorder, asvorder])
    newptext <- rbind(newptext, hmlevs[[ll]]$cor_ptext[varorder, asvorder])
    annot_rows <- c(annot_rows,rep(paste0(var2plot, "-", ll), nrow(hmlevs[[ll]]$cor_ptext)))
  }
  
  # select_asvs: plot these ASVs. Overrides pval
  # pval: plot ASVs with any correlation >  than this threshold
  
  oname <- paste0(outdir, "/", name, ".pdf")
  fontsize_row = 16 - nrow(newcormat) / 15
  fontsize_col = 16 - ncol(newcormat) / 15
  
  annot_rows <- data.frame( Var=rownames(newcormat), Condition = annot_rows)
  rownames(annot_rows) <- paste(annot_rows$Condition, annot_rows$Var, sep="-")
  rownames(newcormat) <- rownames(annot_rows)
  rownames(newptext) <- rownames(annot_rows)
  
  cond_colors <- default_annot_colors[1:length(unique(annot_rows$Condition))]
  names(cond_colors) <- unique(annot_rows$Condition)
  annot_colors <- list(Condition=cond_colors)
  if(! is.null(annot_asvs)){annot_colors[["LFC"]] <- c("Up"="red", "-"="white", "Down"="blue")}
  
  #annot_asvs <- annot_asvs[colnames(newcormat), ]

  tmp <-1
  while(!is.null(tmp)){tmp <- tryCatch(dev.off(), error=function(x){})}
  hm <- pheatmap(newcormat, display_numbers = newptext, 
           #filename=oname,
           #width = 20, height = 12,
           fontsize_row = fontsize_row,
           fontsize_number = 16,
           fontsize_col = fontsize_col,
           cluster_rows = F, cluster_cols = F,
           gaps_row = vspaces,
           annotation_row = annot_rows %>% dplyr::select(Condition),
           labels_row = annot_rows$Var,
           annotation_colors = annot_colors,
           annotation_col = annot_asvs
  )
  #tmp <- tryCatch(dev.off(), error=function(x){})
  pdf(oname, width = 20, height = 10)
  print(hm)
  tmp <- dev.off()
  return(hm)
  #return(list(hm=hm, hm_allonly=hmall, hmbygroup=hmlevs))
}

#### depression extra
makeLinearModelsSingleVariable <- function(divtab, 
                                           interestvar, 
                                           extravars, 
                                           alphaindices =c("Observed", "Chao1", "Shannon", "InvSimpson"), 
                                           combos=1:3,
                                           outdir = "", name = "linearmodels" ){
  if(length(unique(divtab[, interestvar])) < 2) return(list()) #Error: only 1 level, not possible to fit model
  extravars <-  map_vec(divtab[, extravars], \(x)length(unique(x[!is.na(x)]))) %>% 
    base::subset(. > 1) %>% names # remove variables without 2 or more levels
  
  models <- list()
  anovas_singlevar <- data.frame()
  anovas <- data.frame()
  for(aind in alpha_indices){
    formulanull_ch <- paste0(aind, " ~ 1")
    formulanull <- as.formula(formulanull_ch)
    models[[formulanull_ch]] <- lm(formula = formulanull, data = divtab )
    # Variable por variable
    combolist <- c(interestvar, extravars)
    for(varg in combolist){
      formula_ch <- paste0(aind, " ~ ", varg)
      formula <- as.formula(formula_ch)
      models[[formula_ch]] <-  lm(formula = formula, data = divtab )
      auxanova <- anova(models[[formula_ch]]) %>% as.data.frame
      aux <- cbind(data.frame(nvars = 0,
                              Index=aind, 
                              model=formula_ch, 
                              reduced_model = formulanull_ch
                              #mod1=list(models[[formula_ch]]),
                              #mod2=list(models[[formulared_ch]])
      ),
      auxanova[1, ])
      anovas_singlevar <- rbind(anovas_singlevar, aux)
    }
    for(n_comb in combos){
      if(n_comb == 0)next
      combolist <- combn(extravars, n_comb, simplify = F)
      for(varg in combolist){
        formulared_ch <- paste0(aind, " ~ ", paste(varg, sep=" + ", collapse=" + "))
        formulared <- as.formula(formulared_ch)
        models[[formulared_ch]] <-  lm(formula = formulared, data = divtab )
        formula_ch <- paste0(aind, " ~ ", paste(varg, sep=" + ", collapse=" + "), " + ", interestvar)
        formula <- as.formula(formula_ch)
        models[[formula_ch]] <-  lm(formula = formula, data = divtab )
        
        auxanova <- anova(models[[formula_ch]], 
                          models[[formulared_ch]]) %>% as.data.frame
        aux <- cbind(data.frame(nvars = n_comb,
                                Index=aind, 
                                model=formula_ch, 
                                reduced_model = formulared_ch
                                #mod1=list(models[[formula_ch]]),
                                #mod2=list(models[[formulared_ch]])
        ),
        auxanova[2, ])
        anovas <- rbind(anovas, aux)
      }#combos
    }#combo length
  }#alpha dif measures
  results <- list(single_anovas =anovas_singlevar, anovas=anovas,models=models)
  save(results, file=paste0(outdir, "/", name, "linearmodels.RData"))
  write_tsv(anovas_singlevar, file=paste0(outdir, "/", name, "linearmodels_anovas_singlevar.tsv"))
  write_tsv(anovas, file=paste0(outdir, "/", name, "linearmodels_anovas_allvars.tsv"))
  return(results)
}


deseq_full_pipeline <- function(phobj, name, vars2deseq, opt){
  if(!dir.exists(paste0(opt$out, "DeSEQ2"))) dir.create(paste0(opt$out, "DeSEQ2"))
  outdir <- paste0(opt$out, "DeSEQ2/", name, "/")
  opt$reserva <- opt$out
  opt$out <- outdir
  if(!dir.exists(opt$out)) dir.create(opt$out)
  if(opt$minfreq > 0){
    opt$minsampleswithcount <- opt$minfreq*nsamples(phobj)
    cat("Minfreq: ", opt$minfreq, ", setting minsampleswithcount to ", opt$minsampleswithcount)
  }
  dearesults <- getDeseqResults(phobj, opt, name, variables = vars2deseq)
  
  list2env(dearesults, envir = environment())
  tax2annot <- tax_table(phobj)
  resdf_annot <- resdf %>% 
    dplyr::mutate(Genus = data.frame(tax2annot[taxon, "Genus"])$Genus) %>%
    dplyr::select(Genus, everything()) %>% 
    dplyr::arrange(pvalue) 
  write_tsv(resdf_annot, file=paste0(opt$out, "/DEA_annot.tsv"))
  # resdf_annot %>% filter(pvalue < 0.05) %>% dplyr::select(Genus, taxon) %>% 
  #   kable(caption="Differentially abundant ASVs at adjusted p-value < 0.05")
  tryCatch(make_maplot(res, opt, paste0(name, "_MAPlot-rawFC.pdf")),  error=\(x)cat("Error make_maplot"))
  tryCatch(make_maplot(resLFC, opt,  paste0(name, "_MAPlot-rawFC.pdf")),  error=\(x)cat("Error make_maplot"))
  tryCatch(make_maplot(resLFC_ape, opt,  paste0(name, "_MAPlot-rawFC-ape.pdf")),  error=\(x)cat("Error make_maplot"))
  tryCatch(make_maplot(resLFC_ashr, opt,  paste0(name, "_MAPlot-rawFC-ashr.pdf")),  error=\(x)cat("Error make_maplot"))
  plotDispEsts(dds, CV=T , ylim = c(1e-6, 1e1))
  rtabs <- getSummaryTablesDeseq(res, opt)
  tryCatch(make_volcano(resLFC, opt, paste0(name, "volcano_rawfc_rawpval.pdf"), "pvalue"), error=\(x)cat("Error make_volcano"))
  tryCatch(make_volcano(res, opt, paste0(name, "volcano_rawfc_adjpval.pdf"), "padj"),  error=\(x)cat("Error make_volcano"))
  df2plot <- if(! nrow(vst_counts_df)){norm_counts_df}else{vst_counts_df}
  
  if(all_combos_done & length(all_contrasts) > 1){
    cat("All contrasts TRUE, intersecting Taxon list")
    taxalist_praw <- map(all_contrasts, \(x){
      x$resdf %>% dplyr::filter(pvalue < opt$pval) %>% pull(taxon)
    }) %>% unlist %>% unique
    taxalist_padj <- map(all_contrasts, \(x){
      x$resdf %>% dplyr::filter(padj < opt$pval) %>% pull(taxon)
    }) %>% unlist %>% unique
  }else{
    taxalist_praw = taxalist_padj = c()
  }
  tryCatch(makeHeatmap(resdf, dds, df2plot, vars2deseq,
              opt, name = paste0(name, "diff_ab_heatmap_rawpval.pdf"), 
              logscale = F, ptype="pvalue", trim_values = TRUE, taxalist=taxalist_praw), 
           error=\(x) cat("Error makeHeatmap praw"))
  tryCatch(makeHeatmap(resdf, dds, df2plot, vars2deseq,
              opt, name = paste0(name, "diff_ab_heatmap_adjpval.pdf"), 
              logscale = F, ptype="padj", trim_values = TRUE, taxalist=taxalist_padj), 
              error=\(x) cat("Error makeHeatmap padj"))
  opt$out <- opt$reserva
  return(dearesults)
}


plotAbundanceFullPipeline <- function(phobj, interestvar, outdir, phname, levs, tops=c(5,10,15,20)){
  oname <- paste0(outdir, phname, "_phylumBarplot.pdf")
  g1 <- plotRelativeAbnBarsPhylum(phobj, interestvar, oname)
  
  oname <- paste0(outdir, phname, "_phylumBoxplots.pdf")
  g2 <-plotPhylumBoxplots(phobj, interestvar, oname, paired=F)
  
  oname <- paste0(outdir, phname, "_phylumBivariateTests.tsv")
  tests_tab <- getPhylumTests(phobj, interestvar, oname, paired=F)
  
  outname <- paste0(outdir, phname, "_PrevalenceOfPhyla.tsv")
  df_prevalence <- getRelAbundancesByPhylumAndVariable(phobj, interestvar, outname, 
                                                       oldlevs=levs)
  togenplots <- list()
  for(topn in tops){
    oname <- paste0(outdir, phname, "_relAbund_byGenus_top", as.character(topn), ".pdf")
    togenplots[[as.character(topn)]] <- plotRelativeAbnBarsGenus(phobj, interestvar, topn=topn, outname=oname, 
                                                                 oldlevs=levs, 
                                                                 width = 20, height = 8)
  }
  
  oname <- paste0(outdir, phname, "_prevalence_byGenus.tsv")
  pre_prevalence <- getRelAbundancesByGenusAndVariable(phobj, interestvar, oname,
                                                       oldlevs=levs)
  
  oname <- paste0(outdir, phname, "_relAbund_byASV_ColByGenus_top", as.character(n_species), ".pdf")
  g3 <- plotRelativeAbnBars_Fantaxtic(phobj, interestvar, topn = n_species, tax_level="Genus", outname = oname)
  
  oname <- paste0(outdir, phname, "_PrevVsAbund.pdf")
  g4 <- plotPrevalenceVsAbundance(phobj, oname)
  
  resall <- list(
    prev_vs_abn=g4,
    sp_col_by_gen=g3,
    prevalence_tab_genus=pre_prevalence,
    prevalence_tab_phylum=df_prevalence,
    top_genus_plots=togenplots,
    phylum_abundance_tests = tests_tab,
    phylum_boxplots = g2,
    phylum_rel_bars = g1
  )
  return(resall)
}

makeHeatmapsFromMultipleDeseqResults <- function(deseq_results_list, 
                                                 daa_main, 
                                                 main_name,
                                                 vars2plot,
                                                 italics_rownames=T, 
                                                 pfilt=0.05, pplot=0.05,
                                                 name, outdir, w=5, h=4){
  library(RColorBrewer)
  deseq_results_list <- deseq_results_list[vars2plot]
  deseq_results_list[[main_name]] <- daa_main
  all_daa_tab <- names(deseq_results_list) %>% 
    map(\(vname) {
      deseq_results_list[[vname]] %>% 
        dplyr::mutate(variable_controlled=vname)
    }) %>% 
    bind_rows() %>% 
    dplyr::mutate(variable_controlled = gsub("tratamiento_coagul_betabloq_etc", "tratamiento_coag", .$variable_controlled),
                  pvalue = ifelse(pvalue>pplot, NA, -10*log10(pvalue)), #non-significant p-values will be shown in gray
                  padj=ifelse(padj > pplot, NA, -10*log10(padj))
                  ) %>% 
    dplyr::select(taxon, variable_controlled, pvalue, padj, log2FoldChange, log2FoldChangeShrink) %>% 
    pivot_wider( names_from = variable_controlled, 
                 values_from = c(pvalue, padj, log2FoldChange, log2FoldChangeShrink))
  write_tsv(all_daa_tab, paste0(outdir, name, "_MergedCorrectedDAA.tsv"))
  lapply(c("padj", "pvalue"), \(vname){
    mat <- all_daa_tab %>% 
      column_to_rownames("taxon") %>% 
      dplyr::select(starts_with(vname)) %>% 
      as.matrix 
    colnames(mat) <- gsub(paste0(vname, "_"), "", colnames(mat)) 
    mat <- mat[mat[, main_name] >= -10*log10(pfilt) & !is.na(mat[, main_name]) , ]
    mat <- mat[order(mat[, main_name], decreasing = T), ]

    #mat <- mat[ !apply(mat, MAR=1, \(x) all(is.na(x))) ,]
    annot <- data.frame(main_name = mat[, main_name])
    names(annot) <- main_name
    mat <- mat[, colnames(mat) != main_name]
    colnames(mat) <- colnames(mat) %>% gsub("_", " ", .) %>% gsub("tratamiento", "treat.", .) 
    
    labels_row <- gsub("_", " ", rownames(mat))
    if(italics_rownames){
      labels_row <- lapply(labels_row, \(x){
        bquote(italic(.(x)))
      }) %>% as.expression()
    }
    annoCol<-list(Group= colorRampPalette(rev(brewer.pal(n = 7, name =
                                                           "RdYlBu")))(100))#colorRamp(c("dodgerblue3", "firebrick4")))
    names(annoCol) <- main_name
    fontsize_row = 10 - nrow(mat) / 15
    fname <- paste0(outdir, name, '_',main_name, '_', vname, "_.pdf")
    pheatmap(mat, cluster_rows=F,
             show_rownames=nrow(mat) < 120,
             annotation_row = annot,
             fontsize_row = fontsize_row,
             filename = fname,
             width =  w+0.05*ncol(mat),
             height =  h+0.05*nrow(mat),
             color =  colorRampPalette(rev(brewer.pal(n = 7, name =
                                                        "RdYlBu")))(100),
             #legend_breaks = c(-10*log10(c(0.05, 0.001, 0.000001, 0.00000001)), max(mat)),
             #legend_names = as.character(c(0.05, 0.001, 0.000001, 0.00000001, max(mat[!is.na(mat)]))),
             #border_color = NA,
             labels_row = labels_row,
             na_col = "#DDDDDD",
             cluster_cols=T,
             annotation_colors = annoCol
             #annotcluster_rowsation_col=annot
    )
    
  })
  
}

## Se queda con los taxones con p < plim_select en cualquiera de los contrastes
# Colorea según plim_plot (gris si p > plim_plot)
compareLFCContrats <- function(contrastlist, firstContrast, 
                               contrastNamesOrdered, mainContrastName, 
                               plim_select= 0.000001, plim_plot=0.05,
                               name2remove = "",
                               resdfname="resdf", outdir = "./", name="LFC_compare", w=12, h=8,
                               scale_mode="fixed"){
  alldeatables <- map(names(contrastlist), 
                      \(x) contrastlist[[x]][[resdfname]] %>% mutate(comparison=x)) %>% 
    bind_rows()
  alldeatables <- rbind(alldeatables, firstContrast[[resdfname]] %>% 
                          mutate(comparison=mainContrastName)) %>% 
    mutate(taxon = gsub("_", " ", taxon))
  
  tax2plot <- alldeatables %>% filter(padj<=plim_select) %>% pull(taxon) %>% unique
  tax2plot %>% length
  taxorder <- firstContrast[[resdfname]] %>% arrange(desc(log2FoldChangeShrink)) %>% 
    mutate(taxon = gsub("_", " ", taxon)) %>% 
    filter(taxon %in% tax2plot) %>% pull(taxon)
  
  alldeatables_filt <- alldeatables %>% 
    dplyr::filter(taxon %in% taxorder) %>% 
    dplyr::mutate(taxon=factor(taxon, levels=taxorder),
                  comparison=gsub(name2remove, "", comparison),
                  comparison=gsub("_", " ", comparison),
                  comparison=factor(comparison, 
                                    levels=contrastNamesOrdered),
                  UpOrDown = ifelse(log2FoldChangeShrink<0, "Down", "Up"),
                  UpOrDownSig = ifelse(padj<=plim_plot & !is.na(padj), UpOrDown, "NS"),
                  UpOrDownSig = factor(UpOrDownSig, levels=c("Up", "NS", "Down"))
    )
  
  g1<-ggplot(alldeatables_filt, aes(x=taxon, y=log2FoldChangeShrink, fill=UpOrDownSig)) +
    facet_grid(comparison ~ ., scales=scale_mode) +
    geom_hline(yintercept = 0, linetype=3, col="gray80") +
    geom_bar(stat="identity") +
    scale_fill_manual(values=c(C_CASE, C_NS, C_CTRL))+ #c("firebrick1", "lightgray","steelblue2")
    mytheme +
    theme_classic() +
    theme(strip.text.y = element_text(size = 10, 
                                      colour = "black", angle = 0, face = "bold"))+
    theme(axis.text.x = element_text(size = 6, 
                                     colour = "black", angle = 45, 
                                     face = "italic", hjust=1, vjust=1))
    
  
  ggsave(filename = paste0(outdir, "/", name, ".pdf"), g1, width = w, height = h)
  return(g1)
}


compareLFCContrats2 <- function(contrastlist, firstContrast, 
                               contrastNamesOrdered, mainContrastName, 
                               plim_select= 0.000001, plim_plot=0.05,
                               name2remove = "",
                               resdfname="resdf", outdir = "./", name="LFC_compare", w=12, h=8,
                               scale_mode="fixed"){
  alldeatables <- map(names(contrastlist), 
                      \(x) contrastlist[[x]][[resdfname]] %>% mutate(comparison=x)) %>% 
    bind_rows()
  alldeatables <- rbind(alldeatables, firstContrast[[resdfname]] %>% 
                          mutate(comparison=mainContrastName)) %>% 
    mutate(taxon = gsub("_", " ", taxon))
  
  tax2plot <- alldeatables %>% filter(padj<=plim_select) %>% pull(taxon) %>% unique
  tax2plot %>% length
  taxorder <- firstContrast[[resdfname]] %>% arrange(desc(log2FoldChangeShrink)) %>% 
    mutate(taxon = gsub("_", " ", taxon)) %>% 
    filter(taxon %in% tax2plot) %>% pull(taxon)
  
  alldeatables_filt <- alldeatables %>% 
    dplyr::filter(taxon %in% taxorder) %>% 
    dplyr::mutate(taxon=factor(taxon, levels=taxorder),
                  comparison=gsub(name2remove, "", comparison),
                  comparison=gsub("_", " ", comparison),
                  comparison=factor(comparison, 
                                    levels=contrastNamesOrdered),
                  UpOrDown = ifelse(log2FoldChangeShrink<0, "Down", "Up"),
                  UpOrDownSig = ifelse(padj<=plim_plot & !is.na(padj), UpOrDown, "NS"),
                  UpOrDownSig = factor(UpOrDownSig, levels=c("Up", "NS", "Down"))
    )
  
  g1<-ggplot(alldeatables_filt, aes(x=taxon, y=log2FoldChangeShrink, fill=UpOrDownSig)) +
    facet_grid(. ~ comparison, scales=scale_mode) +
    geom_hline(yintercept = 0, linetype=3, col="gray80") +
    geom_bar(stat="identity") +
    scale_fill_manual(values=c(C_CASE, C_NS,C_CTRL))+
    mytheme +
    theme_classic()+
    theme(strip.text.y = element_text(size = 8, 
                                      colour = "black", angle = 0, face = "bold"))+
    theme(axis.text.y = element_text(size = 8, 
                                     colour = "black", angle = 0, 
                                     face = "italic", hjust=1, vjust=0.3))+
    coord_flip() +
    thin_barplot_lines
  
  ggsave(filename = paste0(outdir, "/", name, "2.pdf"), g1, width = w, height = h)
  return(g1)
}


compareLFCContrats3 <- function(contrastlist, firstContrast, 
                                contrastNamesOrdered, mainContrastName, 
                                plim_select= 0.000001, plim_plot=0.05,
                                name2remove = "",
                                resdfname="resdf", outdir = "./", name="LFC_compare3", w=8, h=4,
                                scale_mode="fixed"){
  alldeatables <- map(names(contrastlist), 
                      \(x) contrastlist[[x]][[resdfname]] %>% mutate(comparison=x)) %>% 
    bind_rows()
  alldeatables <- rbind(alldeatables, firstContrast[[resdfname]] %>% 
                          mutate(comparison=mainContrastName)) %>% 
    mutate(taxon = gsub("_", " ", taxon))
  
  tax2plot <- alldeatables %>% dplyr::select(taxon, padj, comparison) %>% 
    spread(key=comparison, value = padj) %>% 
    filter_if(is.numeric, all_vars(. < plim_select)) %>% 
    pull(taxon)

  tax2plot %>% length
  taxorder <- firstContrast[[resdfname]] %>% arrange(desc(log2FoldChangeShrink)) %>% 
    mutate(taxon = gsub("_", " ", taxon)) %>% 
    filter(taxon %in% tax2plot) %>% pull(taxon)
  
  alldeatables_filt <- alldeatables %>% 
    dplyr::filter(taxon %in% taxorder) %>% 
    dplyr::mutate(taxon=factor(taxon, levels=taxorder),
                  comparison=gsub(name2remove, "", comparison),
                  comparison=gsub("_", " ", comparison),
                  comparison=factor(comparison, 
                                    levels=contrastNamesOrdered),
                  UpOrDown = ifelse(log2FoldChangeShrink<0, "Down", "Up"),
                  UpOrDownSig = ifelse(padj<=plim_plot & !is.na(padj), UpOrDown, "NS"),
                  UpOrDownSig = factor(UpOrDownSig, levels=c("Up", "NS", "Down"))
    )
  
  g1<-ggplot(alldeatables_filt, aes(x=taxon, y=log2FoldChangeShrink, fill=UpOrDownSig)) +
    facet_grid(. ~ comparison, scales=scale_mode) +
    geom_hline(yintercept = 0, linetype=3, col="gray80") +
    geom_bar(stat="identity") +
    scale_fill_manual(values=c(C_CASE,C_CTRL))+
    mytheme +
    theme_classic()+
    theme(strip.text.x = element_text(size = 12, 
                                      colour = "black", angle = 0, face = "bold"))+
    theme(axis.text.y = element_text(size = 12, 
                                     colour = "black", angle = 0, 
                                     face = "italic", hjust=1, vjust=0.3))+
    theme(axis.title.x = element_text(size = 12))+
    coord_flip() +
    thin_barplot_lines
  
  ggsave(filename = paste0(outdir, "/", name, "3.pdf"), g1, width = w, height = h)
  return(g1)
}


makeVennLocal <- function(vars2venn, name="VennDiagram", outdir, w=5, h=5){
  
  gv <- ggvenn(
    vars2venn, columns = names(vars2venn),
    stroke_size = 0.5,
    stroke_color = C_NS,
    fill_color = c(C_CASE, C_CASE2, C_CTRL2, C_CTRL),show_elements = F
  )
  ggsave(filename = paste0(outdir, name, ".pdf"), gv, width = w, height = h)
  return(gv)
}

compareLFCContratsNumeric <- function(contrastlist, firstContrast, 
                                contrastNamesOrdered, mainContrastName, 
                                plim_select= 0.000001, plim_plot=0.05,
                                name2remove = "",
                                resdfname="resdf", outdir = "./", name="LFC_compare", w=12, h=8,
                                scale_mode="fixed"){
  alldeatables <- map(names(contrastlist), 
                      \(x) contrastlist[[x]][[resdfname]] %>% mutate(comparison=x)) %>% 
    bind_rows()
  alldeatables <- rbind(alldeatables, firstContrast[[resdfname]] %>% 
                          mutate(comparison=mainContrastName)) %>% 
    mutate(taxon = gsub("_", " ", taxon))
  
  tax2plot <- alldeatables %>% dplyr::filter(comparison==i & padj<=plim_select) %>% pull(taxon) %>% unique
  
  # 1) Calculate correlations
  mat2cor <- alldeatables %>% dplyr::filter(taxon %in% tax2plot) %>% 
    dplyr::select(taxon, log2FoldChange, comparison) %>% 
    spread(key=comparison, value=log2FoldChange) %>% 
    column_to_rownames("taxon") 
  cordf <- data.frame()
  i <- gsub(" ", "_", contrastNamesOrdered[1])
  for(j in names(mat2cor[, names(mat2cor)!=i])){
    mat <- mat2cor[, c(i, j)]
    mat <- mat[!apply(mat, MAR=1, \(x)any(is.na(x))), ]
    aux <- data.frame(
      Contrast1 = i,
      Contrast2 = j,
      pearson_cor = cor(mat[, i], mat[, j], method= "pearson"),
      pearson_pval = cor.test(mat[, i], mat[, j], method= "pearson")$p.value,
      spearman_cor = cor(mat[, i], mat[, j], method= "spearman"),
      spearman_pval = cor.test(mat[, i], mat[, j], method= "spearman")$p.value,
      kendall_cor = cor(mat[, i], mat[, j], method= "kendall"),
      kendall_pval = cor.test(mat[, i], mat[, j], method= "kendall")$p.value
    )
    cordf <- rbind(cordf, aux)
  }
  cordf <- cordf %>% dplyr::mutate(across(matches("_pval$"), list(ajd=\(x)p.adjust(x, method="BH"))))
  write_tsv(cordf, file = paste0(outdir, '/', name, "_Correlations.tsv"))
  
  # 2) Make Venn Diagram
  tax2plot <- alldeatables %>% dplyr::filter(comparison==i & padj<=plim_select) %>% pull(taxon) %>% unique
  tax2plot_up <- alldeatables %>% dplyr::filter(comparison==i & padj<=plim_select & log2FoldChangeShrink >0) %>% pull(taxon) %>% unique
  tax2plot_down <- alldeatables %>% dplyr::filter(comparison==i & padj<=plim_select & log2FoldChangeShrink <0) %>% pull(taxon) %>% unique
  venplots <- map(names(mat2cor[, names(mat2cor)!=i]), \(x){
    tmp <-  alldeatables %>% dplyr::filter(comparison == x & padj < plim_select & ! is.na(padj))
    vars2venn <- list(
      "More in Depr" = tax2plot_up,
      "More in Ctrl" = tax2plot_down,
      " Up" = tmp %>% dplyr::filter(log2FoldChangeShrink > 0) %>% pull(taxon),
      " Down" = tmp %>% dplyr::filter(log2FoldChangeShrink < 0) %>% pull(taxon)
    )
    names(vars2venn)[3:4] <- paste0(x, names(vars2venn)[3:4])
    gvenn <- makeVennLocal(vars2venn, name=paste0(name, "_", x, "_VennDiagram"), outdir, w=8, h=8)
    return(list(groups=vars2venn, plot = gvenn))
  })
  names(venplots)<- names(mat2cor[, names(mat2cor)!=i])
   ## 3) Mosaic plot
  library(ggmosaic)
  
  df2mosaic <- alldeatables %>% 
    dplyr::mutate(direction = ifelse(padj > plim_select | is.na(padj), "NS",
                                     ifelse(log2FoldChangeShrink < 0, "Down", "Up"))) %>%
   
    dplyr::select(taxon, comparison,direction) %>% 
    tidyr::spread(key = comparison, value = direction) %>% 
    dplyr::filter(taxon %in% tax2plot) %>% 
    dplyr::mutate(across(-taxon, \(x)ifelse(is.na(x), "NS", x)),
     across(-taxon, \(x)factor(x, levels=c("Down", "NS", "Up" ))),
     across(!!sym(mainContrastName), \(x)factor(x, levels=c("Down", "NS", "Up")))
    )
  
  gmosplots <- map(names(df2mosaic)[!names(df2mosaic) %in% c("taxon", mainContrastName)], \(contr){
    gmos <- ggplot(data = df2mosaic) +
      geom_mosaic(aes(x = product(!!sym(mainContrastName), !!sym(contr)), fill=!!sym(mainContrastName)), na.rm=T) +
      geom_mosaic_text(aes(x = product(!!sym(mainContrastName), !!sym(contr)), 
                         fill=!!sym(mainContrastName), label=after_stat(.wt)),
      na.rm=T, as.label=T, size=6) +
      scale_fill_manual(values = c(C_CASE, C_NS, C_CTRL))+
      theme_classic() +
      mytheme +
      theme(axis.text=element_text(size=10),
          axis.title=element_text(size=10,face="bold"),
          legend.position = 'none')+
      xlab(gsub("_", " ", contr))+
      ylab(gsub("_", " ", mainContrastName))
    ggsave(filename = paste0(outdir, '/', name,'_',contr, "_MosaicPlots.pdf"), width = 4, height = 4) 
    return(gmos)
  })
  return(list(correlations=cordf, vennplots=venplots, mosaicplots=gmosplots))
}


filter_taxa_padj <- function(resdf, plim=0.05, fc=1){
  taxa <- resdf %>% 
    dplyr::filter(padj <= plim & abs(log2FoldChangeShrink) >= log2(fc) ) %>% 
    pull(taxon)
  return(taxa)
}

filter_taxa_praw <- function(resdf, plim=0.05, fc=1){
  taxa <- resdf %>% 
    dplyr::filter(pvalue <= plim & abs(log2FoldChangeShrink) >= log2(fc) ) %>% 
    pull(taxon)
  return(taxa)
}

makeContingencyPlot <- function(df, var1, var2, outdir, name, w=8, h=6){
  gmos <- ggplot(data = df) +
    geom_mosaic(aes(x = product(!!sym(var1), !!sym(var2)), fill=!!sym(var1)), na.rm=T) +
    geom_mosaic_text(aes(x = product(!!sym(var1), !!sym(var2)), 
                         fill=!!sym(var1), label=after_stat(.wt)),
                     na.rm=T, as.label=T, size=6) +
    theme_classic() +
    mytheme +
    theme(axis.text.x=element_text(size=14),
          axis.title=element_text(size=14,face="bold"),
          axis.text.y = element_text(size = 14, 
                                     colour = "black", angle = 90, face = "bold"),
          legend.position = 'none')+
    xlab(gsub("_", " ", var2))+
    ylab(gsub("_", " ", var1))
  fname <- paste0(outdir, "/", name, "_mosaic.pdf")
  ggsave(filename = fname, gmos, width = w, height = h)
  return(gmos)
}
