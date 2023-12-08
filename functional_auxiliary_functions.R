
readFunctionalMatrix <- function(opt, fname){
  #ftab <- read.table(paste0(opt$input_funcional, fname), comment.char = "", sep="\t", head=T)
  ftab <- read_delim(paste0(opt$input_funcional, fname), delim="\t")
  names(ftab) <- sapply(names(ftab), function(x)strsplit(x, "_")[[1]][1]) %>% gsub("G4M0", "G4M", .)
  names(ftab)[1] <- "Pathway"
  return(ftab)
}

filterWholeProcessesAndFreq <- function(df, opt){
  df2 <- df %>% dplyr::filter(!grepl("\\|", Pathway)) %>% 
    dplyr::filter(!grepl("UNINTEGRATED|UNMAPPED|UNGROUPED", Pathway, perl=T))
  keep <- rowSums(df2 > opt$mincount) > opt$minsampleswithcount
  df2 <- df2[keep, ]
  return(df2)
}

filterBySpeciesAndFreq <- function(df, opt){
  df2 <- df %>% dplyr::filter(grepl("\\|", Pathway))%>% 
    dplyr::filter(!grepl("UNINTEGRATED|UNMAPPED|UNGROUPED", Pathway, perl=T))
  keep <- rowSums(df2 > opt$mincount) > opt$minsampleswithcount
  return(df2[keep, ])
}

getCazyClass <- function(cazy_tt){
  cazytypes <- c( "GT"="GlycosylTransferases",
                  "GH"="Glycoside Hydrolases",
                  "PL"="Polysaccharide Lyases",
                  "CE"="Carbohydrate Esterases",
                  "AA"="Auxiliary Activities",
                  "CB"="Carbohydrate-Binding Modules")
  split_list <- str_split(cazy_tt, pattern = "\\|") %>% unlist %>% substr(1,2)
  cazy_class <- cazytypes[split_list]
  return(cazy_class)
}

limma4functional <- function(df2, metad2, interestvar = "Condition", covars=c()){
  library(limma)
  for(covar in covars){
    metad2 <- metad2 %>% dplyr::filter(!is.na(metad2[, covar]))
  }
  rownames(df2) <- NULL
  expr <- df2 %>% column_to_rownames("Pathway") %>% as.matrix
  expr <- expr[, metad2$sampleID]
  expr <- log(expr+1)
  metad2[, interestvar] <- as.factor(metad2[, interestvar])
  if(length(covars) > 0){
    form <- paste("~0 ", interestvar, paste(covars, collapse = ' + '), 
                  sep = ' + ', collapse=" + ") %>% 
      as.formula
  }else{
    form <- paste0("~0 + ", interestvar) %>% as.formula()
  }
  print(form)
  design <- model.matrix(form, metad2)
  colnames(design) <- gsub(interestvar, "", colnames(design), perl=F)
  #colnames(design) <- gsub("metad2\\$Condition", "", colnames(design), perl=F)
  fit <- lmFit(expr, design)
  cont.matrix <- makeContrasts(case_vs_control = Depression - Control,
                               levels = design)
  fit2 <- contrasts.fit(fit, cont.matrix)
  fit2 <- eBayes(fit2)
  tt <- topTable(fit2, n=Inf, sort.by = "P", adjust.method = "BH")
  names(tt) <- c("log2FoldChange", "AveExpr", "t", "pvalue", "padj", "B") 
  dim(tt)
  return(tt)
}


makeHeatmapFunctional <- function(resdf, met2use, df2plot, 
                                  variable = "condition",
                                  opt, 
                                  name = "heatmap.pdf", 
                                  logscale=FALSE, 
                                  ptype = "padj", w=20, h=14){
  default_annot_colors <- c("gray80", "dodgerblue3", "firebrick3", "green4", "pink4", "skyblue")
  outname <- paste(opt$out, name, sep="/", collapse="/")
  annot <- as.data.frame(met2use[, variable])
  names(annot) <- c(variable)
  rownames(annot) <- met2use$sampleID
  
  if(ptype == "padj"){
    taxa <- resdf %>% dplyr::filter(padj <= opt$pval & 
                                      abs(log2FoldChange) >= log2(opt$fc) ) %>% 
      rownames
  }else{
    taxa <- resdf %>% dplyr::filter(pvalue <= opt$pval & 
                                      abs(log2FoldChange) >= log2(opt$fc) ) %>% 
      rownames
  }
  
  if(length(taxa) < .GlobalEnv$opt$num_genes_default){
    #get only first n genes  
    taxa <- resdf[order(resdf$pvalue), ]  %>% 
      head(.GlobalEnv$opt$num_genes_default) %>% 
      rownames
  } 
  mat <- df2plot %>%  
    filter(Pathway %in% taxa) %>% 
    column_to_rownames("Pathway") %>% 
    as.matrix 
  if(logscale){
    mat <- log(mat + 1)
  }
  mat <- mat %>% 
    t %>% scale %>% t
  
  cond_colors <- default_annot_colors[1:length(unique(annot$Condition))]
  names(cond_colors) <- unique(annot$Condition)
  annot_colors <- list(Condition=cond_colors)
  
  fontsize_row = 14 - nrow(mat) / 15
  fontsize_col = 14 - ncol(mat) / 15
  hm <- pheatmap(mat, cluster_rows=T, 
                 show_rownames=nrow(mat) < 100,
                 annotation_col = annot,
                 fontsize_row = fontsize_row,
                 fontsize_col = fontsize_col,
                 cluster_cols=T, annotcluster_rowsation_col=annot,
                 annotation_colors = annot_colors
  )
  pdf(outname, width = w, height = h)
  print(hm)
  tmp <- dev.off()
  return(hm)
}

nameProcesses<- function(tab, namedtab){
  if(is.null(namedtab)){return(tab)}
  tab$process_name <- namedtab$long[match(rownames(tab), namedtab$short)]
  return(tab)
}

makeBarplotFunctional <- function(tab, plim=0.001, name, outdir, 
                                  include_longnames = FALSE, 
                                  w=12, h=12){
  tab <- tab %>% rownames_to_column("Pathway")
  if(include_longnames & "process_name" %in% names(tab)){
    pnames <- gsub(" \\[[a-zA-Z0-9\\.\\]:\\- ]+", "", tab$process_name, perl=T) %>% gsub(" / ", "/", .)
    tab$Pathway <- paste(tab$Pathway, pnames, sep=": ")
  }
  tab2 <- tab %>% dplyr::filter(padj < plim) %>% 
    dplyr::arrange(log2FoldChange) %>% 
    dplyr::mutate(Pathway= factor(Pathway, levels = Pathway),
                  direction = ifelse(log2FoldChange<0, "less abundant", "more abundant"),
                  #direction = factor(direction, levels=legendnames),
                  `-10log(padj)` = -10*log10(padj)
    )
  
  
  g1 <- ggplot(tab2, aes(x=Pathway, y = log2FoldChange, fill=direction, col=direction ))+ #`-10log(padj)`
    geom_hline(yintercept = 0, linetype=3, col="gray80") +
    geom_bar(stat="identity") +
    coord_flip()+
    mytheme +
    scale_fill_manual(values = c("#A1C6EA", "#FD8B2F")) +
    scale_color_manual(values = c("#A1C6EA", "#FD8B2F")) +
    #scale_fill_gradient2(low="steelblue1", "#DDDDDD", high="tomato")+
    #scale_color_gradient2(low="steelblue1", "#DDDDDD", high="tomato")+
    theme(strip.text.y = element_text(size = 10, 
                                      colour = "black", angle = 0, face = "bold"))+
    guides(fill=guide_legend(title="LFC"), colour=guide_legend(title="LFC"))+
    theme(panel.grid = element_blank())+
    ggtitle(paste0(name, " - proc. with padj<", as.character(plim)))+
    theme(axis.text.y = element_text(size = 12, 
                                     colour = "black", angle = 0, face = "plain"))
  
  ggsave(filename = paste0(outdir, "/", name, ".pdf"), g1, width = w, height = h)
  return(g1)
}
