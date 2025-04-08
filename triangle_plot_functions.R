library(Ternary)
library(ggsci)
library(Cairo)


make_ternary_plot_byTreatment <- function(df2plot, phobj, outdir="./", name=""){
  phmeta <- sample_data(phobj) %>% data.frame
  
  sim <- df2plot %>% column_to_rownames("rowname") %>% 
    as.matrix %>% 
    dist(method = "euclidean") %>% as.matrix %>% as.data.frame %>% 
    rownames_to_column("rowname") %>%
    gather(key = "colname", value = "dist", -rowname) %>%
    mutate(norm_dist = dist/max(dist),
           sim = 100*(1 - norm_dist)) %>%
    dplyr::rename(sampleA = rowname, sampleB=colname) %>% 
    merge(phmeta, 
          by.x = "sampleA", by.y = "sampleID") %>%
    merge(phmeta, 
          by.x = "sampleB", by.y = "sampleID", suffixes = c("_sA","_sB")) %>% 
    filter( Treatment_sB == "NO ABS") %>%
    filter(Stress_sA == Stress_sB) %>%
    group_by(sampleA, Region_sequenced_sA, Treatment_sA, Stress_sA, Region_sequenced_sB, Treatment_sB, Stress_sB) %>% 
    summarise(mean_sim = mean(sim)) %>% ungroup() %>% 
    #mutate(mean_sim = 100* (mean_sim - min(mean_sim)) / (max(mean_sim) - min(mean_sim)) ) %>%   
    spread(key = "Region_sequenced_sB", value = "mean_sim") 
    
  
  CairoPDF(paste0(outdir, name, "_ternary_plot.pdf"), width = 14, height = 7, family = "Arial")
  par(mfrow = c(2, 4), mar = c(0.3, 0.3, 1.3, 0.3))
  
  curr_treat <- "REG3"
  curr_sd <- "Control"
  
  for(curr_sd in c("Control", "SD")){
    
    limits <- sim %>% filter( Stress_sA == curr_sd) %>% ungroup %>% 
      select(sampleA, starts_with("REG", ignore.case=F)) %>% 
      group_by(sampleA) %>%
      group_split() %>%
      map( \(x) x %>% select(-1) %>% unlist )
    
    for(curr_treat in c("ABS", "REG1", "REG2", "REG3")){
    
    
    this_df <- sim %>% 
      filter(Treatment_sA == curr_treat & Stress_sA == curr_sd) 
    #limits =  sim %>% 
    #  filter(Treatment_sA == "NO ABS") %>% 
    #  group_by(Region_sequenced_sA) %>%
    #  summarise(REG1 = mean(REG1), REG2 = mean(REG2), REG3 = mean(REG3)) %>%
    #  group_by(Regisimon_sequenced_sA) %>% 
    #  group_split()%>%
    #  map( \(x) x %>% select(-1) %>% unlist )
    
    data_points <- this_df %>% ungroup %>% 
      select(sampleA, starts_with("REG", ignore.case=F)) %>% 
      group_by(sampleA) %>%
      group_split() %>%
      map( \(x) x %>% select(-1) %>% unlist )
    names(data_points) <- this_df$sampleA 
    
    
    cols <- ggsci::pal_npg()(3)
    
    names(cols) <- c("REG1", "REG2", "REG3")
    # Initial plot
    TernaryPlot(alab = "REG1 \u2192", blab = "\u2190 REG2", clab = "REG3 \u2192",
                region = limits,
                #xlim = c(30, 70), ylim = c(30, 70),
                lab.col =cols,
                main = paste0(curr_sd, '-', curr_treat), # Title
                point = "right", 
                lab.cex = 1.8, 
                grid.minor.lines = 0,
                grid.lty = "solid", 
                col = rgb(1, 1, 1), 
                grid.col = "lightgray", 
                axis.col =cols, 
                ticks.col = cols,
                axis.rotate = FALSE,
                padding = 0.1)
    # Colour the background:
    #cols <- TernaryPointValues(cols)
    #ColourTernary(cols, spectrum = NULL)
    
    
    AddToTernary(graphics::points, data_points, pch = 16, cex = 2, 
                 bg = cols[this_df$Region_sequenced_sA],
                 col = cols[this_df$Region_sequenced_sA],
    )
  
  }}
 dev.off()
}


make_ternary_plot_byRegion <- function(df2plot, phobj, outdir="./", name=""){
  phmeta <- sample_data(phobj) %>% data.frame
  
  sim <- df2plot %>% column_to_rownames("rowname") %>% 
    as.matrix %>% 
    dist(method = "euclidean") %>% as.matrix %>% as.data.frame %>% 
    rownames_to_column("rowname") %>%
    gather(key = "colname", value = "dist", -rowname) %>%
    mutate(norm_dist = dist/max(dist),
           sim = 100*(1 - norm_dist)) %>%
    dplyr::rename(sampleA = rowname, sampleB=colname) %>% 
    merge(phmeta, 
          by.x = "sampleA", by.y = "sampleID") %>%
    merge(phmeta, 
          by.x = "sampleB", by.y = "sampleID", suffixes = c("_sA","_sB")) %>% 
    filter( Treatment_sB == "NO ABS") %>%
    filter(Stress_sA == Stress_sB) %>%
    group_by(sampleA, Region_sequenced_sA, Treatment_sA, Stress_sA, Region_sequenced_sB, Treatment_sB, Stress_sB) %>% 
    summarise(mean_sim = mean(sim)) %>% ungroup() %>% 
    #mutate(mean_sim = 100* (mean_sim - min(mean_sim)) / (max(mean_sim) - min(mean_sim)) ) %>%   
    spread(key = "Region_sequenced_sB", value = "mean_sim") 
  
  
  CairoPDF(paste0(outdir, name, "_ternary_plot2.pdf"), width = 12, height = 7, family = "Arial")
  par(mfrow = c(2, 3), mar = c(0.3, 0.3, 1.3, 0.3))
  
  curr_treat <- "REG3"
  curr_sd <- "Control"
  
  for(curr_sd in c("Control", "SD")){
    
    limits <- sim %>% filter( Stress_sA == curr_sd) %>% ungroup %>% 
      select(sampleA, starts_with("REG", ignore.case=F)) %>% 
      group_by(sampleA) %>%
      group_split() %>%
      map( \(x) x %>% select(-1) %>% unlist )
    
    for(curr_treat in c("REG1", "REG2", "REG3")){
      
      
      this_df <- sim %>% 
        filter(Region_sequenced_sA == curr_treat & Stress_sA == curr_sd) 
      #limits =  sim %>% 
      #  filter(Treatment_sA == "NO ABS") %>% 
      #  group_by(Region_sequenced_sA) %>%
      #  summarise(REG1 = mean(REG1), REG2 = mean(REG2), REG3 = mean(REG3)) %>%
      #  group_by(Regisimon_sequenced_sA) %>% 
      #  group_split()%>%
      #  map( \(x) x %>% select(-1) %>% unlist )
      
      data_points <- this_df %>% ungroup %>% 
        select(sampleA, starts_with("REG", ignore.case=F)) %>% 
        group_by(sampleA) %>%
        group_split() %>%
        map( \(x) x %>% select(-1) %>% unlist )
      names(data_points) <- this_df$sampleA 
      
      
      cols <- ggsci::pal_npg()(5)[2:5]
      
      names(cols) <- c("ABS", "REG1", "REG2", "REG3")
      # Initial plot
      TernaryPlot(alab = "REG1 \u2192", blab = "\u2190 REG2", clab = "REG3 \u2192",
                  region = limits,
                  #xlim = c(30, 70), ylim = c(30, 70),
                  lab.col =cols[2:5],
                  main = paste0(curr_sd, '-', curr_treat), # Title
                  point = "right", 
                  lab.cex = 1.8, 
                  grid.minor.lines = 0,
                  grid.lty = "solid", 
                  col = rgb(1, 1, 1), 
                  grid.col = "lightgray", 
                  axis.col =cols[2:5], 
                  ticks.col = cols[2:5],
                  axis.rotate = FALSE,
                  padding = 0.1)
      # Colour the background:
      #cols <- TernaryPointValues(cols)
      #ColourTernary(cols, spectrum = NULL)
      
      
      AddToTernary(graphics::points, data_points, pch = 16, cex = 2, 
                   bg = cols[this_df$Treatment_sA],
                   col = cols[this_df$Treatment_sA],
      )
      
    }}
  dev.off()
}


#make_ternary_plot_byTreatment(df2plot, phobj, outdir, name )
#make_ternary_plot_byRegion(df2plot, phobj, outdir, name )
