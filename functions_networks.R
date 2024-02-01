library(tidyverse)
library(ggrepel)
library(geomtextpath)
library(cowplot)

C_CASE = "#FD8B2F" #"rgba(200, 44, 44, 0.8)"
C_CASE2 = "tomato"
C_CASE_LINK = "#fBd895" #"#f9c784"
C_CTRL = "#A1C6EA" # "#006daa" #
C_CTRL2 = "steelblue2"
C_CTRL_LINK ="#DAE3E5" #"rgba(44, 44, 200, 0.8)"
C_WHITE= "#DDDDDD"
C_NS =  "#A5ABBD" #"rgba(224, 224, 224, 0.8)"
C_OTHER = "gray30"

C_CASE3 = "#00AA5A"
C_CTRL3 = "#8E7BFF"

process_names <- function(namelist){
  gsub("-", ".", namelist) %>% 
    gsub("\\[|\\]|\\(|\\)", "", ., perl=T) 
}

process_names2 <- function(namelist){
  
  res <- gsub("_", " ", namelist) %>% 
    sapply(\(x){
      if(x=="") return("")
      x <- strsplit(x, " ")[[1]]
      y <- toupper(strsplit(x[1], "")[[1]][1])
      y <- paste0(y, ". ", x[2])
      return(y)
    })
}

getGraphFromSamples <- function(phobj, samples, daataxa, vstdf, 
                                net_estimator, net_method, 
                                daatab, daatab2=NULL,
                                outdir="./",
                                filter_empty=FALSE, 
                                filt_quantile=0.95, 
                                name="testnet", w=12, h=12){
  # https://github.com/ryanjw/co-occurrence
  # https://rdrr.io/bioc/minet/man/minet.html  
  auxf <- function(x){sd(x)>0}
  df2net <- vstdf %>% 
    dplyr::filter(gene %in% daataxa) %>% 
    dplyr::mutate(gene = process_names(gene)) %>% 
    column_to_rownames("gene") %>% 
    as.matrix %>% t %>% data.frame %>% 
    rownames_to_column("sample") %>% 
    dplyr::filter(sample %in% samples) %>% 
    column_to_rownames("sample") %>% 
    dplyr::select_if(auxf)

  netres <- minet(df2net, 
                  method=net_method, 
                  estimator=net_estimator, 
                  disc="none", 
                  nbins=sqrt(NROW(df2net)))
  
  vert_cols <- ifelse(is.na(daatab$padj), C_NS,
                      ifelse(daatab$padj < opt$pval, 
                             ifelse(daatab$log2FoldChangeShrink > 0, C_CASE, C_CTRL), C_NS))
  if(! is.null(daatab2)){
    daatab2 <- daatab2 %>% dplyr::filter(taxon %in% daatab$taxon)
    daatab2 <- daatab2[match(daatab2$taxon, daatab$taxon), ]
    vert_cols <- ifelse(is.na(daatab2$padj) | vert_cols != C_NS, vert_cols,
                        ifelse(daatab2$padj < opt$pval, 
                               ifelse(daatab2$log2FoldChangeShrink > 0, C_CASE3, C_CTRL3), 
                               vert_cols))
  }
 
  names(vert_cols) <- process_names(daatab$taxon)
  assertthat::assert_that(all(names(df2net) %in% names(vert_cols) ), msg = "Error! names in dfnet and vert colors are different")

  print(table(vert_cols))
  
  if(filt_quantile > 0){
    filtthres <- quantile(netres, filt_quantile)
    netres <- ifelse(netres < filtthres, 0, 1)
  }
  if(filter_empty){
    keep <- apply(netres, MAR=1, \(x)any(x>0)) %>% which %>% names
    netres <- netres[, keep]
  }
  
  
  net <- graph.adjacency(netres, mode="undirected")
  vert_cols <- vert_cols[names(V(net))] 
  print(length(V(net)))
  
  fname <- paste0(outdir, "/",name, "net.pdf")
  pdf(fname, width=12, height=12)
  p1 <- plot(net, vertex.label.size=0.1, vertex.label=NA,
             vertex.size=4, vertex.color=vert_cols)
  dev.off()
  return(list(net=net, cols=vert_cols))
}

getGraphProps <- function(net_obj){
  net <- net_obj$net
  gprop <- data.frame(
    vertex = names(V(net)),
    closeness = closeness(net), 
    betweenness = betweenness(net),
    degree = degree(net),
    color = net_obj$cols
  ) %>%
    dplyr::mutate(closeness = ifelse(is.na(closeness),0,closeness)) %>% 
    dplyr::mutate(across(all_of(c("closeness", "betweenness", "degree")), \(x){
        scaled <- (x-min(x))/(max(x)-min(x))
        scaled[is.na(scaled)] <- 0
        return(scaled)
        },.names = "scaled_{.col}")) %>% 
    dplyr::mutate(sum_measures = scaled_closeness+scaled_betweenness+scaled_degree) %>% 
    dplyr::arrange(desc(sum_measures))
  return(gprop)
}

make_boxplot_nodeprops <- function(nodeprops, outdir, name, w=12, h=8,
                                   signif_levels=c("***"=0.001, "**"=0.01, "*"=0.05, "ns"=1.1),
                                   correct_pvals = TRUE){
  propslong <- nodeprops %>% 
    dplyr::mutate(closeness = ifelse(closeness==1, NA, closeness)) %>% 
    gather("property", "value", closeness, betweenness, degree)
  
  tests <- map(c("closeness", "betweenness", "degree"),
               \(x)getTestsForAllCombinations(nodeprops[, x], nodeprops[, "class"]) %>% 
                 dplyr::mutate(measure=x) %>% 
                 dplyr::select(measure, everything())
  ) %>% bind_rows
  fname <- paste0(outdir, "/", name, "_stats.tsv")
  write_tsv(tests, file = fname)
  
  num_comparisons <- length(unique(tests$groups_compared)>1)
  if(correct_pvals & num_comparisons>1){
    signif_levels <- c(signif_levels[1:3]/num_comparisons, signif_levels[4])
  }
  
  comp <- combn(unique(as.character(propslong$class)), 2, simplify = F)
  
  (g1 <- ggplot(propslong, aes(x=class, y=value, fill=class))+
      facet_wrap(~ property, scales="free")+
      geom_violin(alpha=0.5, outlier.shape = NA)+
      geom_boxplot(width=0.2, fill="darkgray", 
                   outlier.shape = NA,
                   notchwidth = 0.5, notch=F)+
      ggsignif::stat_signif(test="wilcox.test", na.rm=T, comparisons = comp, 
                            step_increase=0.03,
                            tip_length = 0.01,
                            map_signif_level=signif_levels,
                            vjust=0.4,
                            color = "black"
      ) +
      #coord_flip() +
      theme_pubclean()+
      theme(axis.text.x = element_text(size = 12, 
                                       colour = "black", angle = 45, 
                                       face = "plain", hjust=1, vjust=1))
  )
  
  fname <- paste0(outdir, name, "_boxplot.pdf")
  ggsave(filename = fname, g1, width = w, height = h)
  return(list(tests=tests, plot=g1))
}


plotNetWithTopNames <- function(net_obj, nodeprops, net_class, ntop=10, outdir, name, add_labels=TRUE,
                                w=8, h=8){
  
  top_taxa <- nodeprops %>% dplyr::filter(class==net_class) %>% top_n(sum_measures, n=10) %>% 
    pull(vertex)
  
  net_vnames <- names(V(net_obj$net))
  if(add_labels){
    ver_names <- ifelse(net_vnames %in% top_taxa, net_vnames, "") %>% 
      process_names %>% process_names2
  }else{
    ver_names <- ""
  }
  fname <- paste0(outdir, name, ".pdf")
  V(net_obj$net)$label.cex <- 1.2
  V(net_obj$net)$label.face <- "italic"
  pdf(fname, width=w, height=h)
  p1 <- plot(net_obj$net, vertex.label.size=6, vertex.label=ver_names,
             vertex.label.color="black",
             vertex.label.dist=0.5,
             vertex.size=ifelse(net_vnames %in% top_taxa, 8, 4), 
             vertex.color=net_obj$cols)
  dev.off()
  
}


makeTopTaxaPlot <- function(nodeprops, ntop, outdir, name, 
                            w=8, h=8, manual_scale=c("#A1C6EA", "#FD8B2F", "#A5ABBD")){
  top_taxa <- nodeprops %>% dplyr::group_by(class) %>% top_n(sum_measures, n=ntop) %>% 
    gather("measure", "value", closeness, betweenness, degree) %>% 
    dplyr::mutate(vertex = gsub("_", " ", vertex))
  fname <- paste0(outdir, name, "_top", as.character(ntop), "Table.tsv")
  write_tsv(x = top_taxa, file = fname)
  gprops <- ggplot(top_taxa, aes(x=vertex, y=value, fill=color))+
    geom_col()+
    coord_flip()+
    facet_grid(class~measure, scales = "free")+
    theme_pubclean()+
    scale_fill_manual(values=manual_scale)+
    theme(axis.text.y=element_text(face="italic"))
  fname <- paste0(outdir, name, "_top", as.character(ntop), "Barplot.pdf")
  ggsave(filename = fname, gprops, width = w, height = h)
  return(gprops)
}
