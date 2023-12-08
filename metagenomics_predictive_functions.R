library(tidyverse)
library(ggsci)
library(phyloseq)


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


plotPCA <- function(dfraw, vars2pc, labvar = "PACIENTE", plotvars=c("Condition"), 
                    transform = "scale", dims=2:5, name = "PCAs.pdf", outdir="", w=12, h=8){
  
  dfraw <- dfraw %>% column_to_rownames(labvar) 
  mat <- dfraw[, vars2pc] %>%  
    as.matrix
  
  if(transform=="scale"){
    mat <- mat %>% scale 
  }else if(transform == "log"){
    mat <- log(mat + 1) 
  }else if(transform == "sqrt"){
    mat <- mat %>% sqrt 
   }#else{
  #   mat <- mat 
  # }
  
  dfraw$sample <- rownames(dfraw)
  pca <- prcomp(mat, center = T, scale. = T)
  df <- pca$x %>% as.data.frame %>% 
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
    if(! is.numeric(df2[, plotvar])){
      g1 <- g1 + scale_color_lancet() 
    }
    return(g1)
    })
  plotlist[[plotvar]] <- cowplot::plot_grid(plotlist=gs)
  }
  
  WriteManyPlots(plotlist,name = name, outdir = outdir, w = w, h=h)
  dev.off()
  return(list(pca=pca, plots=plotlist))
}



plotPCA_princomp <- function(dfraw, vars2pc, labvar = "PACIENTE", plotvars=c("Condition"), 
                    transform = "scale", dims=2:5, name = "PCAs.pdf", outdir="", w=12, h=8){
  
  dfraw <- dfraw %>% column_to_rownames(labvar) 
  mat <- dfraw[, vars2pc] %>%  
    as.matrix
  
  if(transform=="scale"){
    mat <- mat %>% scale 
  }else if(transform == "log"){
    mat <- log(mat + 1) 
  }else if(transform == "sqrt"){
    mat <- mat %>% sqrt
  }else{
    mat <- mat 
  }
  
  dfraw$sample <- rownames(dfraw)
  pca <- princomp(mat, cor=TRUE)
  df <- pca$scores %>% as.data.frame 
  names(df) <- gsub("Comp.", "PC", names(df))
  df2 <- df %>%  rownames_to_column("sample") %>% 
    merge( dfraw, by="sample")
  
  pcts <- getPropVar_princomp(pca)
  plotlist <- list()
  for(plotvar in plotvars){
    
    legname <- gsub("_", " ", plotvar)
    legname <- gsub("\\.", "-", legname)
    gs <- lapply(paste0("PC", dims), FUN=function(PC){
      g1 <- ggplot(df2, aes_string(x="PC1", y=PC, col=plotvar)) +
        geom_point() + 
        geom_text_repel(aes(label = sample)) + 
        xlab(pcts["PC1"]) +
        ylab(pcts[PC]) +
        mystyle +
        guides(col=guide_legend(title=legname))
      if(! is.numeric(df2[, plotvar])){
        g1 <- g1 + scale_color_lancet() 
      }
      return(g1)
    })
    plotlist[[plotvar]] <- cowplot::plot_grid(plotlist=gs)
  }
  
  WriteManyPlots(plotlist,name = name, outdir = outdir, w = w, h=h) 
  return(list(pca=pca, plots=plotlist))
}



makeLogisticRegressions <- function(depvar, posval, df, vars2test){
  df[, depvar] <- ifelse(df[, depvar]==posval, 1, 0)
  names(df) <- gsub("-", ".", names(df))
  vars2test <- gsub("-", ".", vars2test)
  mods_alone <- lapply(vars2test, FUN=function(var){
    form <- paste0(depvar, "~ ", var) %>% as.formula
    fm1 <- glm(form, data = df, family = binomial(link = "logit"))
    return(fm1)
  }) 
  mods <- lapply(vars2test, FUN=function(var){
    form <- paste0(depvar, "~ ", var) %>% as.formula
    fm1 <- stats::glm(form, data = df, family = binomial(link = "logit"))
    ss <- summary(fm1)
    pred <- ifelse(predict(fm1, type="response") > 0.5, 1, 0)
    confmat <- confusionMatrix(as.factor(pred), as.factor(df$Psoriasis))
    tab <- cbind(
      variable = var,
      broom::glance(fm1),
      coef_intercept=ss$coefficients[1, 1],
      coef_slope=ss$coefficients[2, 1],
      stderr_intercept=ss$coefficients[1, 2],
      stderr_slope=ss$coefficients[2, 2],      
      zval_intercept=ss$coefficients[1, 3],
      zval_slope=ss$coefficients[2, 3],      
      pval_intercept=ss$coefficients[1, 4],
      pval_slope=ss$coefficients[2, 4],
      Sensitivity=confmat$byClass["Sensitivity"],
      Specificity=confmat$byClass["Specificity"],
      PosPredValue=confmat$byClass["Pos Pred Value"],
      NegPredValue=confmat$byClass["Neg Pred Value"],
      Accuracy=confmat$overall["Accuracy"]
    )
    return(tab)
  }) %>% bind_rows()
  rownames(mods) <- names_cyt
  mods <- mods %>% arrange(dplyr::desc(Accuracy))
  return(list(mods_alone=mods_alone, mods=mods))
}



plotBoxplots <- function(df, variable, grouping_var, wrap_var, 
                              fname, 
                              signif_levels=c("***"=0.001, "**"=0.01, "*"=0.05, "ns"=1.1),
                              num_comparisons = 1,
                              ylabel = "pg/mL",
                              correct_pvalues = TRUE, write=TRUE){
   if(correct_pvalues){
     signif_levels_bonferroni <- c(signif_levels[1:3]/num_comparisons, signif_levels[4])
   }else{
     signif_levels_bonferroni <- signif_levels
   }
    wrap_formula <- paste0(". ~ ",  wrap_var) %>% as.formula
    df[, grouping_var] <- as.factor(df[, grouping_var])
    comp <- combn(unique(as.character(df[, grouping_var])), 2, simplify = F)
    g1 <-  ggplot(df, aes_string(x=grouping_var, y=variable, col=grouping_var)) +
      facet_wrap(wrap_formula, scales = "free") +
      geom_boxplot(alpha = 0.7, width=0.5, outlier.alpha=0) +
      geom_jitter(alpha=0.5)+
      #scale_color_manual(values = c("#ffafcc", "#90DBF4")) + 
      #scale_fill_manual(values = c("#ffafcc", "#90DBF4")) +
      scale_color_lancet() + 
      scale_fill_lancet() +
      labs(title = variable, x = '') +
      #mytheme +
      theme_pubclean() +
      ggtitle("") +
      geom_signif(test="wilcox.test", na.rm=T, comparisons = comp, 
                  step_increase=0.03,
                  tip_length = 0.01,
                  map_signif_level=signif_levels_bonferroni,
                  vjust=0.4,
                  color = "black"
      ) +
      theme(plot.title = element_text(hjust = 0.5)) +
      theme(axis.text.x = element_blank()) +
      ylab(ylabel)
    
    if(write){
      ggsave(filename=fname, g1, width = 6, height = 8)
    }
    return(g1)
    #theme(axis.text.x = element_text(angle = 360, hjust = 0.5, size = 10))
  } #Plots qualitative variables


plotBarplot <- function(df, variable, grouping_var, wrap_var, 
                         fname, 
                         signif_levels=c("***"=0.001, "**"=0.01, "*"=0.05, "ns"=1.1),
                         num_comparisons = 1,
                         ylabel = "pg/mL",
                         correct_pvalues = TRUE, write=TRUE){
  if(correct_pvalues){
    signif_levels_bonferroni <- c(signif_levels[1:3]/num_comparisons, signif_levels[4])
  }else{
    signif_levels_bonferroni <- signif_levels
  }
  wrap_formula <- paste0(". ~ ",  wrap_var) %>% as.formula
  df[, grouping_var] <- as.factor(df[, grouping_var])
  comp <- combn(unique(as.character(df[, grouping_var])), 2, simplify = F)
  
  g1 <-  ggplot(df, aes_string(x=grouping_var, y=variable, fill=grouping_var, col=grouping_var)) +
    facet_wrap(wrap_formula, scales = "free", nrow=1) +
    #geom_bar(stat="identity") +
    stat_summary(fun.data=mean_sdl, geom="bar", width=0.8, alpha=0.8) +
    stat_summary(fun.data=mean_cl_boot, geom="errorbar", width=0.3, col="black") +
    geom_point()+
    #scale_color_manual(values = c("#ffafcc", "#90DBF4")) + 
    #scale_fill_manual(values = c("#ffafcc", "#90DBF4")) +

    labs(title = variable, x = '') +
    #mytheme +
    theme_pubclean() +
    geom_signif(test="wilcox.test", na.rm=T, comparisons = comp, 
                step_increase=0.03,
                tip_length = 0.01,
                map_signif_level=signif_levels_bonferroni,
                vjust=0.4,
                color = "black"
    ) +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(axis.text.x = element_blank()) +
    mytheme +
    ggtitle("Cytokine levels per group") +
    ylab(ylabel)
  if(length(unique(df[, grouping_var] )) == 2){
    g1 <- g1  +
      scale_color_manual(values =c("firebrick4", "dodgerblue4")) + 
      scale_fill_manual(values =c("firebrick2", "dodgerblue2")) 
  }else{
    g1 <- g1  +
      scale_color_lancet() + 
      scale_fill_lancet()
  }
  if(write){
    ggsave(filename=fname, g1, width = 12, height = 6)
  }
  return(g1)
  #theme(axis.text.x = element_text(angle = 360, hjust = 0.5, size = 10))
} #Plots qualitative variables


plotAgainstAllFactors<- function(alldf, groupvars, name, outdir, opt, correct_bonferroni=F, names_cyt=c()){
  
  plotlist <- list()
  for(grouping_var in groupvars){
    if(length(unique(alldf[, grouping_var])) < 2 | min(table(alldf[, grouping_var])) <2) next
    df <- alldf[, c(grouping_var, names_cyt)] %>% gather(IL, pg_mL, names_cyt)
    oname <- paste0(outdir, name)
    plotlist[[grouping_var]] <- plotBoxplots(df, "pg_mL", grouping_var, "IL", 
                                             fname=oname, 
                                             signif_levels=c("***"=0.001, "**"=0.01, "*"=0.05, "ns"=1.1),
                                             num_comparisons = nrow(test_res),
                                             ylabel = "pg/mL",
                                             correct_pvalues = correct_bonferroni, 
                                             write = F)
  }
  WriteManyPlots(plotlist, name = name, outdir =outdir, w = 12, h = 12, separate = F, opt)
  return(plotlist)
}


plotAgainstAllFactors_Barplot<- function(alldf, groupvars, name, outdir, opt, correct_bonferroni=F, names_cyt=c()){
  
  plotlist <- list()
  for(grouping_var in groupvars){
    if(length(unique(alldf[, grouping_var])) < 2 | min(table(alldf[, grouping_var])) <2) next
    df <- alldf[, c(grouping_var, names_cyt)] %>% gather(IL, pg_mL, names_cyt)
    oname <- paste0(outdir, name)
    plotlist[[grouping_var]] <- plotBarplot(df, "pg_mL", grouping_var, "IL", 
                                             fname=oname, 
                                             signif_levels=c("***"=0.001, "**"=0.01, "*"=0.05, "ns"=1.1),
                                             num_comparisons = nrow(test_res),
                                             ylabel = "pg/mL",
                                             correct_pvalues = correct_bonferroni, 
                                             write = F)
  }
  WriteManyPlots(plotlist, name = name, outdir =outdir, w = 12, h = 12, separate = F, opt)
  return(plotlist)
}


plotAgainstAllQuant<- function(alldf, groupvars, name, outdir, col_var, opt, names_cyt=c()){
  
  plotlist <- list()
  for(grouping_var in groupvars){
    df <- alldf[, c(grouping_var, names_cyt, col_var)] %>% tidyr::gather(IL, pg_mL, names_cyt)
    oname <- paste0(outdir, name)
    plotlist[[grouping_var]] <- plotRegr(df, variable="pg_mL", 
                                             x_var = grouping_var, 
                                             wrap_var = "IL",
                                             col_var ="Psoriasis",
                                             fname=oname, 
                                             write = F)
  }
  WriteManyPlots(plotlist, name = name, outdir =outdir, w = 12, h = 12, separate = F, opt)
  return(plotlist)
}


plotRegr <- function(df, variable, x_var, wrap_var, col_var,
                         fname, 
                         write=TRUE){

  wrap_formula <- paste0(". ~ ",  wrap_var) %>% as.formula
  df[, col_var] <- as.factor(df[, col_var])
  points_col <- ifelse(df[, col_var] == "no", "dodgerblue3", "firebrick3") #levels(df[, col_var])[1]
  g1 <-  ggplot(df, aes_string(x=x_var, y=variable)) +
    facet_wrap(wrap_formula, scales = "free") +
    geom_smooth(method="lm") +
    geom_point(alpha=0.6, col=points_col)+
    stat_poly_eq(use_label(c("R2", "p", "n")), #c("eq", "R2", "f", "p", "n")
                 method="lm", small.p=T, small.r=F, label.y=0.99)+
    #scale_color_manual(values = c("#ffafcc", "#90DBF4")) + 
    #scale_fill_manual(values = c("#ffafcc", "#90DBF4")) +
    scale_color_lancet() + 
    scale_fill_lancet() +
    labs(title = variable, x = '') +
    #mytheme +
    theme_pubclean() +
    ggtitle("") +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(axis.text.x = element_blank()) +
    ggtitle(x_var)
  
  if(write){
    ggsave(filename=fname, g1, width = 6, height = 8)
  }
  return(cowplot::plot_grid(g1))
  #theme(axis.text.x = element_text(angle = 360, hjust = 0.5, size = 10))
} #Plots qualitative variables
