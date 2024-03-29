library(plotly)
library(rjson)
library(reticulate)

# Install venv and virtualenv from terminal before reticulate:
# sudo apt-get install python3-virtualenv
# sudo apt-get install python3-venv

# Create virtual environment with plotly:
# virtualenv_create("r-reticulate")
# virtualenv_install("r-reticulate", "plotly")
# use_virtualenv("r-reticulate")

use_python("/usr/bin/python")

C_CASE = "#FD8B2F" #"rgba(200, 44, 44, 0.8)"
C_CASE_LINK = "#fBd895" #"#f9c784"
C_CTRL = "#A1C6EA" # "#006daa" #
C_CTRL_LINK ="#DAE3E5" #"rgba(44, 44, 200, 0.8)"
C_WHITE= "#DDDDDD"
C_NS =  "#A5ABBD" #"rgba(224, 224, 224, 0.8)"




prepareSankeyPlot <- function(daa_proc_all, 
                              daa_by_sp_all, 
                              quant_by_proc, 
                              metadf, 
                              plim = 0.01, 
                              plim_sp = 0.01, 
                              include_others = FALSE, 
                              include_longnames = FALSE,
                              others_name = "Aggr. taxa",
                              case_name = "Depression", 
                              control_name = "Control", 
                              outdir = "~/", 
                              name = "test"){
  ## Returns a list with 
  ##  1) a data frame of nodes and 
  ##  2) a data frame of links between nodes.
  # Get processes to plot
  daa_proc1 <- daa_proc_all %>% 
    dplyr::filter(padj< plim) %>% 
    rownames_to_column("Pathway")
  
  # Get process-species pairs that are present in selected processes, 
  # with or without aggregating non-significant species
  if(include_others){
    daa_sp <- daa_by_sp_all %>% 
      dplyr::filter(Pathway %in% daa_proc1$Pathway) %>% 
      #dplyr::filter(padj < plim_sp) #%>% 
      dplyr::mutate(Species=ifelse(padj < plim_sp, Species, others_name)) %>% 
      group_by(Pathway, Genus, Species) %>% 
      dplyr::summarise_if(is.numeric, mean) %>% 
      ungroup()
  }else{
    daa_sp <- daa_by_sp_all %>% 
      dplyr::filter(Pathway %in% daa_proc1$Pathway) %>% 
      dplyr::filter(padj < plim_sp) #%>% 
  }
  #all(daa_proc1$Pathway %in% daa_sp$Pathway)
  
  # Get color scale for species nodes (limits according to process:species LFC)
  scale_species = scales::gradient_n_pal(colours=c(C_CTRL, C_WHITE, C_CASE), 
                                         values=c(min(c(-1,min(daa_sp$log2FoldChange))), 
                                                  0,
                                                  max(c(1,max(daa_sp$log2FoldChange)))))
  # Get color scale for process nodes (limits according to process LFC)
  scale_proc = scales::gradient_n_pal(colours=c(C_CTRL, C_WHITE, C_CASE), 
                                      values=c(min(c(-1,min(daa_proc1$log2FoldChange))), 
                                               0,
                                               max(c(1,max(daa_proc1$log2FoldChange)))))
  # Get color scale for process-species link (limits according to process:species LFC)
  scale_species_link = scales::gradient_n_pal(colours=c(C_CTRL_LINK, C_WHITE, C_CASE_LINK), 
                                              values=c(min(c(-1,min(daa_sp$log2FoldChange))), 
                                                       0,
                                                       max(c(1,max(daa_sp$log2FoldChange)))))
  
  part1 <- daa_sp %>% dplyr::select(Pathway, Species, AveExpr, log2FoldChange, padj) %>% 
    #mutate(col = ifelse(padj > plim, C_NS, ifelse(log2FoldChange < 0, C_CTRL, C_CASE))) %>% 
    dplyr::mutate(col=scale_species_link(log2FoldChange)) %>% 
    dplyr::select(-padj) %>% 
    group_by(Pathway) %>% 
    dplyr::mutate(AveExpr = AveExpr/(sum(AveExpr)))
  names(part1) <- c("node1", "node2", "Size", "LFC", "color")
  
  tab_quant_all <- quant_by_proc %>% dplyr::filter(Pathway %in% daa_proc1$Pathway)
  
  part2 <- tab_quant_all %>% 
    pivot_longer(2:ncol(tab_quant_all)) %>% 
    dplyr::mutate(Condition = metadf$Condition[match(name, metadf$sampleID)]) %>% 
    group_by(Condition, Pathway) %>% 
    dplyr::summarise(value=mean(value)) %>% 
    dplyr::mutate(LFC= 0,
           color = ifelse(Condition == case_name, C_CASE_LINK, C_CTRL_LINK)) %>% 
    ungroup() %>% 
    group_by(Condition) %>% 
    dplyr::mutate(value = value/sum(value)) %>% 
    ungroup() %>% group_by(Pathway) %>% dplyr::mutate(value = value/sum(value))
  
  
  part2 <- part2 %>% select(Condition, Pathway, value, LFC, color) %>% 
    dplyr::arrange(Condition)
  names(part2) <-  c("node1", "node2", "Size", "LFC", "color")
  
  nodes2fill <- daa_proc1$Pathway[!daa_proc1$Pathway %in% part1$node1]
  if(length(nodes2fill)>0){
    part3 <- data.frame(node1 = nodes2fill, 
                      node2 = others_name, 
                      Size=OTHERS_SIZE, 
                      LFC=daa_proc1$log2FoldChange[match(nodes2fill, daa_proc1$Pathway)]) %>% 
    dplyr::mutate(color = scale_species_link(LFC))
  
    full_v1 <- rbind(part2, part1, part3)
  }else{
    full_v1 <- rbind(part2, part1)  
  }
  
  
  nodelist_1 <- data.frame(value = c(control_name, case_name), 
                           value2match = c(control_name, case_name), 
                           xpos = c(0,0),
                           color = c(C_CTRL, C_CASE),
                           nodesize = c(prop.table(table(metad$Condition))[control_name],
                                        prop.table(table(metad$Condition))[case_name]))
  
  tb_norm <- tab_quant_all %>% mutate_if(is.numeric, \(x)x/sum(x)) %>% 
    dplyr::select(-Pathway) %>% rowSums() 
  tb_norm <- tb_norm/sum(tb_norm)
  
  include_longnames = (include_longnames & "process_name" %in% names(daa_proc1) & "process_name" %in% names(daa_by_sp_all))
  nodelist_2 <- rbind(nodelist_1, 
                      data.frame(
                        value = {if(include_longnames)paste(daa_proc1$Pathway, daa_proc1$process_name, sep=":") else daa_proc1$Pathway},
                        value2match = daa_proc1$Pathway,
                        xpos = 0.25,
                        color = scale_proc(daa_proc1$log2FoldChange),
                        #color = ifelse(daa_proc1$log2FoldChange < 0, C_CASE, C_CTRL), 
                        nodesize = tb_norm[match(daa_proc1$Pathway, tab_quant_all$Pathway)]
                      ))
  species <- part1 %>% group_by(node2) %>% 
    dplyr::summarise(LFC= mean(LFC),
              SumExpr = sum(Size)) #%>% 
  # dplyr::mutate(color = ifelse(LFC < 0, C_CASE, C_CTRL))
  scale_species_node = scales::gradient_n_pal(colours=c(C_CTRL, C_WHITE, C_CASE), 
                                              values=c(min(species$LFC), 
                                                       0,
                                                       max(species$LFC)))
  nodelist_3 <- data.frame(
    value = species$node2,
    value2match = species$node2,
    xpos = 1,
    #color = species$color, 
    color = scale_species_node(species$LFC),
    nodesize = species$SumExpr
  )
  
  nodelist <- rbind(nodelist_2, nodelist_3) %>% dplyr::mutate(nodenum=0:(n()-1))
  if(others_name %in% full_v1$node2 & ! others_name %in% nodelist$value){
    nodelist <- rbind(nodelist, 
                      data.frame(value = others_name, 
                                 value2match = others_name,
                                 xpos = 1,
                                 color=C_NS, 
                                 nodesize=1, nodenum=nrow(nodelist))
    )
  }
  nodelist$color[nodelist$value==others_name] <- C_NS
  
  # nodelist <- full_v1 %>% select(node1, node2, color) %>% pivot_longer(node1:node2) %>% 
  #   select(-name) %>% dplyr::distinct() %>% 
  #   mutate(nodenum = 0:(n()-1))
  full_v1 <- full_v1 %>% dplyr::mutate(source=nodelist$nodenum[match(node1, nodelist$value2match)],
                                       target = nodelist$nodenum[match(node2, nodelist$value2match)])
  write_tsv(nodelist, file=paste0(outdir, "/", name, "_nodelist4sankey.tsv"))
  write_tsv(full_v1, file=paste0(outdir, "/", name, "_linklist4sankey.tsv"))
  return(list(nodelist=nodelist, linklist=full_v1))
}

makeSankeyPlot <- function(nodelist, linklist, name, fname, outdir, fontsize=18){
  # From https://plotly.com/r/sankey-diagram/
  print(head(nodelist))
  fig <- plot_ly(
    type = "sankey",
    arrangement="fixed",
    domain = list(
      x =  c(0,1),
      y =  c(0,1)
    ),
    orientation = "h",
    valueformat = ".0f",
    valuesuffix = "TWh",
    node = list(
      label = nodelist$value,
      color = nodelist$color,
      value = nodelist$nodesize,
      x = nodelist$xpos,
      pad = 15,
      thickness = 15,
      line = list(
        color = "black",
        width = 0.5
      )
    ),
    link = list(
      source = linklist$source,
      target = linklist$target,
      value =  linklist$Size,
      label =  as.character(round(linklist$LFC, 2)),
      color =  linklist$color
    )
  )
  
  fig <- fig %>% layout(
    title = name,
    font = list(
      size = fontsize
    ),
    xaxis = list(showgrid = F, zeroline = F),
    yaxis = list(showgrid = F, zeroline = F)
  )
  
  #fig
  #tmp <- paste0(outdir, "/", fname, ".png")
  #save_image(fig, tmp)
  #use_virtualenv("r-reticulate")
  htmlwidgets::saveWidget(fig, file=paste0(outdir, fname, ".html"))
  return(fig)
}
