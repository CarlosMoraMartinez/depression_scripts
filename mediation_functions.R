library(tidyverse)
library(bmem)
library(sem)
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


# #So amazing
# library(gtExtras)
# metadata %>% gt_plt_summary()
# #Also cool
# library(gtsummary)
# metadata %>% tbl_summary()

makeMediationSimple <- function(df, xname, yname, medname){
  library(bmem)
  library(sem)
  
  eq1 = paste0(yname, " = b*", medname, " + cp*", xname)
  eq2 = paste0(medname, " = a*", xname)
  eqs <- paste(eq1, eq2, sep="\n", collapse="\n")
  
  iris.model<-specifyEquations(exog.variances=T, text=eqs)
  
  effects<-c('a*b', 'cp+a*b') 
  nlsy.res<-bmem.sobel(df, iris.model, effects) 
  return(nlsy.res)
  
}

makeMediationSimple_mergedX <- function(df, xnames, yname, medname){
  
  a_params <- paste("a", 1:length(xnames), sep="")
  cp_params <- paste("cp", 1:length(xnames), sep="")
  
  eq1 = paste(yname, " = b*", medname, " + ", cp_params, "*", xnames, sep="") %>% 
    paste(collapse = "\n")
  eq2 = paste0(medname, " = ", a_params, "*", xnames) %>% 
    paste(collapse = "\n")
  eqs <- paste(eq1, eq2, sep="\n", collapse="\n")
  
  iris.model<-specifyEquations(exog.variances=T, text=eqs)
  
  effects<-c(paste(a_params, "*b", sep="" ), paste(cp_params, '+', a_params,'*b', sep="")) 
  nlsy.res<-bmem.sobel(df, iris.model, effects) 
  oldnames2 <- rownames(nlsy.res$estimates) %>% 
    sapply(\(x)strsplit(x, "\\*|\\+",perl = TRUE)[[1]][1]) 
  # newnames <- ifelse(oldnames %in% a_params, 
  #        paste(oldnames, xnames[match(oldnames, a_params)], sep="_"),
  #        ifelse(oldnames %in% cp_params, 
  #               paste(oldnames, xnames[match(oldnames, cp_params)], sep="_"), 
  #               oldnames))
  resdf <- nlsy.res$estimate %>% 
    rownames_to_column("param")
  resdf$x_labels <- ifelse(oldnames2 %in% a_params, 
                           xnames[match(oldnames2, a_params)],
                           ifelse(oldnames2 %in% cp_params, 
                                  xnames[match(oldnames2, cp_params)], 
                                  gsub("V\\[|\\]", "", oldnames2)))
  #rownames(nlsy.res$estimates) <- newnames
  return(resdf)
}

makeMediationSimple_mergedMediat <- function(df, x_name, y_name, mednames){
  
  a_params <- paste("a", 1:length(mednames), sep="")
  b_params <- paste("b", 1:length(mednames), sep="")
  
  eq1 = paste(y_name, " = ", b_params, "*", mednames, " + ", "cp*", x_name, sep="") %>% 
    paste(collapse = "\n")
  eq2 = paste0(mednames, " = ", a_params, "*", x_name) %>% 
    paste(collapse = "\n")
  eqs <- paste(eq1, eq2, sep="\n", collapse="\n")
  
  iris.model<-specifyEquations(exog.variances=T, text=eqs)
  
  effects<-c(paste(a_params, "*", b_params, sep="" ), paste('cp+', a_params,'*', b_params, sep="")) 
  nlsy.res<-bmem.sobel(df, iris.model, effects) 
  
  oldnames2 <- rownames(nlsy.res$estimates) %>% 
    sapply(\(x){y<-strsplit(x, "\\*|\\+",perl = TRUE)[[1]];y[length(y)]}) 
  # newnames <- ifelse(oldnames %in% a_params, 
  #        paste(oldnames, xnames[match(oldnames, a_params)], sep="_"),
  #        ifelse(oldnames %in% cp_params, 
  #               paste(oldnames, xnames[match(oldnames, cp_params)], sep="_"), 
  #               oldnames))
  resdf <- nlsy.res$estimate %>% 
    rownames_to_column("param")
  resdf$x_labels <- ifelse(oldnames2 %in% a_params, 
                           mednames[match(oldnames2, a_params)],
                           ifelse(oldnames2 %in% b_params, 
                                  mednames[match(oldnames2, b_params)], 
                                  gsub("V\\[|\\]", "", oldnames2)))
  
  #rownames(nlsy.res$estimates) <- newnames
  return(resdf)
}

makeMediationComplex <- function(df, xname="Edad", yname="Condition_bin", 
                                 medname="IMC", medname2="Faecalibacterium_prausnitzii"){
  
  eq1 = paste0(yname, " = b*", medname, " + cp*", medname2, " + fp*", xname)
  eq2 = paste0(medname, " = a*", medname2, " + d*", xname)
  eq3 = paste0(medname2, " = e*", xname)
  eqs <- paste(eq1, eq2, eq3, sep="\n", collapse="\n")
  
  iris.model<-specifyEquations(exog.variances=T, text=eqs)
  
  effects<-c('a*b', 'cp+a*b', 'd*b', 'fp+d*b', 'e*cp', 'fp+e*cp', 'd*b+fp+e*cp+e*a*b') 
  nlsy.res<-bmem.sobel(df, iris.model, effects) 
  return(nlsy.res)
  
}

makeMediationComplex_mergedMEdiat2 <- function(df, xname="Edad", 
                                               yname="Condition_bin", 
                                               medname="IMC", 
                                               mednames2=c("Faecalibacterium_prausnitzii")){
  
  a_params <- paste("a", 1:length(mednames2), sep="")
  cp_params <- paste("cp", 1:length(mednames2), sep="")
  e_params <- paste("e", 1:length(mednames2), sep="")
  
  eq1 = paste(yname, " = b*", medname, " + ", cp_params, "*", mednames2, " + fp*", xname, sep="") %>% 
    paste(collapse = "\n")
  eq2 = paste0(medname, " = ", a_params, "*", mednames2, " + d*", xname) %>% 
    paste(collapse = "\n")
  eq3 = paste0(mednames2, " = ", e_params, "*", xname)%>% 
    paste(collapse = "\n")
  eqs <- paste(eq1, eq2, eq3, sep="\n", collapse="\n")
  
  iris.model<-specifyEquations(exog.variances=T, text=eqs)
  
  effects<-c(paste(a_params, "*b", sep="" ), 
             paste(cp_params, '+', a_params,'*b', sep=""),
             'fp+d*b',
             paste(e_params, '*', cp_params, sep=""),
             paste('fp + ', e_params, '*', cp_params, sep=""),
             paste('d*b + fp + ', e_params, '*', cp_params, ' + ', e_params, '*',a_params, "*b", sep="")
  ) 
  nlsy.res<-bmem.sobel(df, iris.model, effects) 
  
  resdf <- nlsy.res$estimate %>% 
    rownames_to_column("param")
  resdf$x_labels <- resdf$param %>% 
    sapply(\(x){
      y <-strsplit(x, "\\*|\\+",perl = TRUE)[[1]] %>% gsub(" ", "", .) %>% sort 
      
      if(any(stringr::str_starts(y, "a"))){
        pthis <- y[stringr::str_starts(y, "a")]
        return(paste(mednames2[a_params==pthis], collapse=":"))
      }else if(any(stringr::str_starts(y, "cp"))){
        pthis <- y[stringr::str_starts(y, "cp")]
        return(paste(mednames2[cp_params==pthis], collapse=":"))
      }else if(any(stringr::str_starts(y, "e"))){
        pthis <- y[stringr::str_starts(y, "e")]
        return(paste(mednames2[e_params==pthis], collapse=":"))
      }else{
        return(gsub("V\\[|\\]", "", x))
      }
      
    }) 
  return(resdf)
}

plotIMCMediationSimple <- function(res, vars2test, outdir, outname, 
                                   plim_plot = 0.05, use_color_scale=FALSE, w=14, h=10){
  
  edge_width_factor <- 10 ##Constants to adjust edge width. However, not used anymore
  minq <- 0.2
  rangefactor <- 10
  
  bacorder <- res %>% filter(grepl("^cp[0-9]+\\+a[0-9]+\\*b$", param)) %>% arrange(Estimate) #^cp[0-9]+$
  assertthat::assert_that(all(vars2test %in% bacorder$x_labels) & nrow(bacorder)==length(vars2test))
  
  vnames <- c(bacorder$x_labels, "D", "BMI")
  vertices <- data.frame(names = gsub("_", " ", vnames),
                         xpos1 = c( rep(1, length(vnames)-2), #1:(length(vnames)-2),
                                    15, #length(vnames)+5,
                                    10),#length(vnames)+5),
                         ypos1 = c(length(vnames):3, 
                                   as.integer(length(vnames)*0.9),
                                   as.integer(length(vnames)*0.2) 
                         ),
                         size = c(rep(4, length(vars2test)), 10, 10),
                         totaleffect = c(bacorder$Estimate, NA, NA),
                         total_pvalue = c(bacorder$p.value, NA, NA),
                         colors = c(ifelse(bacorder$p.value < plim_plot, 
                                           ifelse(bacorder$Estimate< 0, C_CTRL, C_CASE), 
                                           C_NS),
                                    C_CASE, C_OTHER),
                         xoffset = ifelse(vnames == "D", 0.5, ifelse(vnames== "BMI", 1, 0))
  )
  
  
  edgenames <- expand.grid(c("a", "cp"), as.character(1:length(vars2test)) ) %>% 
    as.matrix %>% apply(MAR=1,paste,sep="", collapse="")
  edgenames <- c("b", edgenames)
  assertthat::assert_that(all(edgenames %in% res$param))  
  
  minval <- quantile(abs(res$Estimate)*edge_width_factor, minq) #to trim edge widths
  
  
  edgetab <- res %>% dplyr::filter(param %in% edgenames) %>% 
    dplyr::mutate(
      from=ifelse(param=="b", "BMI", gsub("_", " ", x_labels)),
      to = ifelse(grepl("^a", param), "BMI", "D"),
      x0 = vertices$xpos1[match(from, vertices$names)],
      y0 = vertices$ypos1[match(from, vertices$names)],
      x1 = vertices$xpos1[match(to, vertices$names)],
      y1 = vertices$ypos1[match(to, vertices$names)],
      labelnames = ifelse(to=="D" & from == "BMI", paste0("b=", as.character(round(Estimate*10,2))),
                          ifelse(to=="D", paste0("c'=", as.character(round(Estimate*10,2))), 
                                 paste0("a=", as.character(round(Estimate*10,2))))
      ),
      color = ifelse(Estimate < 0, C_CTRL, C_CASE),
      color = ifelse(p.value > plim_plot, C_NS, color),
      linetype = ifelse(p.value > plim_plot, 2, 1),
      width = abs(Estimate)*edge_width_factor,
      width = ifelse(width <  minval, minval, 
                     ifelse(width > rangefactor*minval, 
                            rangefactor*minval, width))
    ) %>% dplyr::select(from, to, everything())
  
  edgetab$y1[edgetab$to=="D" & edgetab$from=="BMI"] = edgetab$y1[edgetab$to=="D" & edgetab$from=="BMI"]-0.75
  edgetab$x1[edgetab$to=="D" & edgetab$from=="BMI"] = edgetab$x1[edgetab$to=="D" & edgetab$from=="BMI"]+0.1
  
  if(use_color_scale){
    scale_link = scales::gradient_n_pal(colours=c(C_CTRL2, C_CASE_LINK), 
                                        values=c(min(c(-1,min(edgetab$Estimate))), 
                                                 0,
                                                 max(c(1,max(edgetab$Estimate)))))
    edgetab <- edgetab %>% dplyr::mutate(color = ifelse(p.value > plim_plot, C_NS, scale_link(Estimate)))
    vertices <- vertices %>% dplyr::mutate(colors = c(ifelse(bacorder$p.value < plim_plot, 
                                                             scale_link(totaleffect), 
                                                             C_NS),
                                                      C_CASE, C_OTHER))
    
  }
  
  vertices$labelnames <- mapply(vertices$names,vertices$totaleffect, FUN=function(x, val) {
    if( x %in% c("D", "BMI"))return(x)else 
      return(paste0("italic('", gsub("sp ", "sp. ", x), 
                    "')~(", as.character(round(10*val, 2)), ")"))}, SIMPLIFY = T)
  (g1 <- ggplot()+
      # geom_point(data=vertices, aes(x=xpos1, y=ypos1),
      #            size=5, col=vertices$colors)+
      # geom_segment(data = edgetab, aes(x=x0, y=y0, xend=x1, yend=y1),
      #              size = edgetab$width,
      #              col = edgetab$color,
      #              linetype = edgetab$linetype,
      #              arrow=arrow(length=unit(0.01, "npc"), type="closed")) +
      geom_textsegment(data = edgetab, 
                       aes(label=labelnames,x=x0, y=y0, xend=x1, yend=y1),
                       #size = edgetab$width,
                       col = edgetab$color,
                       linetype = edgetab$linetype,
                       arrow=arrow(length=unit(0.01, "npc"), type="closed")) +
      geom_label(data=vertices, aes(label=labelnames, x=xpos1, y=ypos1), 
                 parse=TRUE,
                 size=vertices$size,
                 hjust=1,
                 nudge_x = vertices$xoffset,
                 col=vertices$colors, inherit.aes = T)+
      xlim(c(min(vertices$xpos1-7), max(vertices$xpos1+7)))+
      theme_void()
  )
  ggsave(filename = paste0(opt$out, "/", name,'_p', as.character(plim_plot), ".pdf"), g1, 
         width = w, height = h)
  return(list(plot=g1, bacorder=bacorder, vertices=vertices, edges=edgetab))
}

plotAgeMediationSimple <- function(res, vars2test, outdir, outname, 
                                   plim_plot = 0.05, use_color_scale=FALSE, w=14, h=6){
  
  edge_width_factor <- 10 ##Constants to adjust edge width. However, not used anymore
  minq <- 0.2
  rangefactor <- 10
  
  bacorder <- res %>% filter(grepl("^+b[0-9]+$", param)) %>% arrange(Estimate) #^cp[0-9]+$
  assertthat::assert_that(all(vars2test %in% bacorder$x_labels) & nrow(bacorder)==length(vars2test))
  
  vnames <- c(bacorder$x_labels, "D", "Age")
  total_age_effect <- res %>% dplyr::filter(grepl("^cp\\+", param, perl=T)) %>% 
    pull(Estimate) %>% 
    mean
  vertices <- data.frame(names = gsub("_", " ", vnames),
                         xpos1 = c( rep(1, length(vnames)-2), #1:(length(vnames)-2),
                                    20, #length(vnames)+5,
                                    14),#length(vnames)+5),
                         ypos1 = c(length(vnames):3, 
                                   7,
                                   3 ),
                         size = c(rep(4, length(vars2test)), 10, 10),
                         totaleffect = c(bacorder$Estimate, NA, total_age_effect),
                         total_pvalue = c(bacorder$p.value, NA, NA),
                         colors = c(ifelse(bacorder$p.value < plim_plot, 
                                           ifelse(bacorder$Estimate< 0, C_CTRL, C_CASE), 
                                           C_NS),
                                    C_CASE, C_OTHER),
                         xoffset = ifelse(vnames == "D", 0.5, ifelse(vnames== "BMI", 1, 0))
  )
  
  
  edgenames <- expand.grid(c("a", "b"), as.character(1:length(vars2test)) ) %>% 
    as.matrix %>% apply(MAR=1,paste,sep="", collapse="")
  edgenames <- c("cp", edgenames)
  assertthat::assert_that(all(edgenames %in% res$param))  
  
  minval <- quantile(abs(res$Estimate)*edge_width_factor, minq) #to trim edge widths
  
  
  edgetab <- res %>% dplyr::filter(param %in% edgenames) %>% 
    dplyr::mutate(
      from=ifelse(grepl("^b", param), gsub("_", " ", x_labels), "Age"),
      to = ifelse(grepl("^a", param), gsub("_", " ", x_labels), "D"),
      x0 = vertices$xpos1[match(from, vertices$names)],
      y0 = vertices$ypos1[match(from, vertices$names)],
      x1 = vertices$xpos1[match(to, vertices$names)],
      y1 = vertices$ypos1[match(to, vertices$names)],
      labelnames = ifelse(to=="D" & from == "Age", paste0("c'=", as.character(round(Estimate*10,2))),
                          ifelse(to=="D", paste0("b=", as.character(round(Estimate*10,2))), 
                                 paste0("a=", as.character(round(Estimate*10,2))))
      ),
      color = ifelse(Estimate < 0, C_CTRL, C_CASE),
      color = ifelse(p.value > plim_plot, C_NS, color),
      linetype = ifelse(p.value > plim_plot, 2, 1),
      width = abs(Estimate)*edge_width_factor,
      width = ifelse(width <  minval, minval, 
                     ifelse(width > rangefactor*minval, 
                            rangefactor*minval, width))
    ) %>% dplyr::select(from, to, everything())
  
  edgetab$y1[edgetab$to=="D" & edgetab$from=="Age"] = edgetab$y1[edgetab$to=="D" & edgetab$from=="Age"]-0.3
  edgetab$x1[edgetab$to=="D" & edgetab$from=="Age"] = edgetab$x1[edgetab$to=="D" & edgetab$from=="Age"]-0.3
  
  if(use_color_scale){
    scale_link = scales::gradient_n_pal(colours=c(C_CTRL2, C_CASE_LINK), 
                                        values=c(min(c(-1,min(edgetab$Estimate))), 
                                                 0,
                                                 max(c(1,max(edgetab$Estimate)))))
    edgetab <- edgetab %>% dplyr::mutate(color = ifelse(p.value > plim_plot, C_NS, scale_link(Estimate)))
    vertices <- vertices %>% dplyr::mutate(colors = c(ifelse(bacorder$p.value < plim_plot, 
                                                             scale_link(totaleffect), 
                                                             C_NS),
                                                      C_CASE, C_OTHER))
    
  }
  
  vertices$labelnames <- mapply(vertices$names,vertices$totaleffect, FUN=function(x, val) {
    if( x == "D"){
      return(x)
    }else if(x=="Age"){
      return("Age") #paste0("Age (", as.character(round(10*val, 2)), ")")
    }else{
      return(paste0("italic('", gsub("sp ", "sp. ", x), "')"))
                                  }}, SIMPLIFY = T)
  (g1 <- ggplot()+
      # geom_point(data=vertices, aes(x=xpos1, y=ypos1),
      #            size=5, col=vertices$colors)+
      # geom_segment(data = edgetab, aes(x=x0, y=y0, xend=x1, yend=y1),
      #              size = edgetab$width,
      #              col = edgetab$color,
      #              linetype = edgetab$linetype,
      #              arrow=arrow(length=unit(0.01, "npc"), type="closed")) +
      geom_textsegment(data = edgetab, 
                       aes(label=labelnames,x=x0, y=y0, xend=x1, yend=y1),
                       #size = edgetab$width,
                       col = edgetab$color,
                       linetype = edgetab$linetype,
                       arrow=arrow(length=unit(0.01, "npc"), type="closed")) +
      geom_label(data=vertices, aes(label=labelnames, x=xpos1, y=ypos1), 
                 parse=TRUE,
                 size=vertices$size,
                 hjust=1,
                 nudge_x = vertices$xoffset,
                 col=vertices$colors, inherit.aes = T)+
      xlim(c(min(vertices$xpos1-7), max(vertices$xpos1+7)))+
      ylim(c(-1, max(vertices$ypos1)+2))+
      theme_void()
  )
  ggsave(filename = paste0(opt$out, "/", name,'_p', as.character(plim_plot), ".pdf"), g1, 
         width = w, height = h)
  return(list(plot=g1, bacorder=bacorder, vertices=vertices, edges=edgetab))
}


plotAgeAndIMCMediationComplex <- function(res, vars2test, outdir, outname, 
                                   plim_plot = 0.05, use_color_scale=FALSE, w=14, h=10){
  
  edge_width_factor <- 10 ##Constants to adjust edge width. However, not used anymore
  minq <- 0.2
  rangefactor <- 10
  
  bacorder <- res %>% filter(grepl("^d\\*b\\+fp\\+e", gsub(" ", "", param))) %>% arrange(Estimate)#^cp[0-9]+$
  assertthat::assert_that(all(vars2test %in% bacorder$x_labels) & nrow(bacorder)==length(vars2test))
  
  vnames <- c(bacorder$x_labels, "D", "Age", "BMI")
  total_age_effect <- res %>% dplyr::filter(grepl("^d\\*b\\+fp\\+e", param, perl=T)) %>% 
    pull(Estimate) %>% 
    mean
  total_age_pval <- res %>% dplyr::filter(grepl("^d\\*b\\+fp\\+e", param, perl=T)) %>% 
    pull(p.value) %>% 
    mean
  
  vertices <- data.frame(names = gsub("_", " ", vnames),
                         xpos1 = c( rep(1, length(bacorder$x_labels)), #1:(length(vnames)-2),
                                    24, #length(vnames)+5,
                                    12, 
                                    18),#length(vnames)+5),
                         ypos1 = c(1.1*(length(vnames):4-2), 
                                   as.integer(length(vnames)*0.7),
                                   as.integer(length(vnames)*1.1),
                                   1),#as.integer(length(vnames)*0.2)),
                         size = c(rep(4, length(vars2test)), 10, 10, 10),
                         totaleffect = c(bacorder$Estimate, NA, total_age_effect, NA),
                         total_pvalue = c(bacorder$p.value, NA, total_age_pval, NA),
                         colors = c(ifelse(bacorder$p.value < plim_plot, 
                                           ifelse(bacorder$Estimate< 0, C_CTRL, C_CASE), 
                                           C_NS),
                                    C_CASE, C_OTHER, C_OTHER),
                         xoffset = ifelse(vnames == "D", 0.5, ifelse(vnames== "BMI", 1, 0))
  )
  
  
  edgenames <- expand.grid(c("a", "cp", "e"), as.character(1:length(vars2test)) ) %>% 
    as.matrix %>% apply(MAR=1,paste,sep="", collapse="")
  edgenames <- c("d", "b", "fp", edgenames)
  assertthat::assert_that(all(edgenames %in% res$param))  
  
  minval <- quantile(abs(res$Estimate)*edge_width_factor, minq) #to trim edge widths
  
  fromdict <- list(e="Age", d="Age", b="BMI", fp="Age")
  todict <- list(a="BMI", d="BMI", b="D", cp="D", fp="D")
  selectfun <- function(x, pname, plist){
    return(if(gsub("[0-9]+", "", x, perl=T) %in% names(plist)) 
      plist[[gsub("[0-9]+", "", x, perl=T)]] else 
        gsub("_", " ", pname))
  }
  edgetab <- res %>% dplyr::filter(param %in% edgenames) %>% 
    dplyr::mutate(
      from=map2_vec(param, x_labels, selectfun, fromdict),
      to = map2_vec(param, x_labels, selectfun, todict),
      x0 = vertices$xpos1[match(from, vertices$names)],
      y0 = vertices$ypos1[match(from, vertices$names)],
      x1 = vertices$xpos1[match(to, vertices$names)],
      y1 = vertices$ypos1[match(to, vertices$names)],
      labelnames = paste0(gsub("p", "'", gsub("[0-9]+", "", param)), '=',as.character(round(Estimate*10,2))),
      color = ifelse(Estimate < 0, C_CTRL, C_CASE),
      color = ifelse(p.value > plim_plot, C_NS, color),
      linetype = ifelse(p.value > plim_plot, 2, 1),
      width = abs(Estimate)*edge_width_factor,
      width = ifelse(width <  minval, minval, 
                     ifelse(width > rangefactor*minval, 
                            rangefactor*minval, width))
    ) %>% dplyr::select(from, to, everything())
  
  edgetab$y1[edgetab$to=="D" & edgetab$from=="Age"] = edgetab$y1[edgetab$to=="D" & edgetab$from=="Age"]+0.7
  edgetab$x1[edgetab$to=="D" & edgetab$from=="Age"] = edgetab$x1[edgetab$to=="D" & edgetab$from=="Age"]-0.3
  
  edgetab$y1[edgetab$to=="D" & edgetab$from=="BMI"] = edgetab$y1[edgetab$to=="D" & edgetab$from=="BMI"]-0.7
  edgetab$x1[edgetab$to=="D" & edgetab$from=="BMI"] = edgetab$x1[edgetab$to=="D" & edgetab$from=="BMI"]-0.7
  
  edgetab$x1[edgetab$to=="D" & !edgetab$from %in% c("Age", "BMI")] = edgetab$x1[edgetab$to=="D" & !edgetab$from %in% c("Age", "BMI")] -0.7
  
  edgetab$y0[edgetab$from=="Age" ] = edgetab$y0[edgetab$from=="Age" ]-0.7
  edgetab$x0[edgetab$from=="Age" ] = edgetab$x0[edgetab$from=="Age" ]-0.7
  
  edgetab$y1[edgetab$to=="BMI"] = edgetab$y1[edgetab$to=="BMI"]+0.7
  edgetab$x1[edgetab$to=="BMI" & edgetab$from !="Age"] = edgetab$x1[edgetab$to=="BMI" & edgetab$from !="Age"]-0.7
  
  if(use_color_scale){
    scale_link = scales::gradient_n_pal(colours=c(C_CTRL2, C_CASE_LINK), 
                                        values=c(min(c(-1,min(edgetab$Estimate))), 
                                                 0,
                                                 max(c(1,max(edgetab$Estimate)))))
    edgetab <- edgetab %>% dplyr::mutate(color = ifelse(p.value > plim_plot, C_NS, scale_link(Estimate)))
    vertices <- vertices %>% dplyr::mutate(colors = c(ifelse(bacorder$p.value < plim_plot, 
                                                             scale_link(totaleffect), 
                                                             C_NS),
                                                      C_CASE, C_OTHER, C_OTHER))
    
  }
  
  vertices$labelnames <- mapply(vertices$names,vertices$totaleffect, FUN=function(x, val) {
    if( x %in% c("D", "BMI", "Age")){
      return(x)
    }else{
      return(paste0("italic('", gsub("sp ", "sp. ", x), "')"))
    }}, SIMPLIFY = T)
  (g1 <- ggplot()+
      # geom_point(data=vertices, aes(x=xpos1, y=ypos1),
      #            size=5, col=vertices$colors)+
       # geom_segment(data = edgetab, aes(x=x0, y=y0, xend=x1, yend=y1),
       #              size = edgetab$width,
       #              col = edgetab$color,
       #              linetype = edgetab$linetype,
       #              arrow=arrow(length=unit(0.01, "npc"), type="closed")) +
      geom_textsegment(data = edgetab, 
                       aes(label=labelnames,x=x0, y=y0, xend=x1, yend=y1),
                       #size = edgetab$width,
                       col = edgetab$color,
                       linetype = edgetab$linetype,
                       arrow=arrow(length=unit(0.01, "npc"), type="closed")) +
      geom_label(data=vertices, aes(label=labelnames, x=xpos1, y=ypos1), 
                 parse=TRUE,
                 size=vertices$size,
                 hjust=1,
                 nudge_x = vertices$xoffset,
                 col=vertices$colors, inherit.aes = T)+
      xlim(c(min(vertices$xpos1-7), max(vertices$xpos1+7)))+
      ylim(c(-1, max(vertices$ypos1)+2))+
      theme_void()
  )
  ggsave(filename = paste0(opt$out, "/", name,'_p', as.character(plim_plot), ".pdf"), g1, 
         width = w, height = h)
  
  return(list(plot=g1, bacorder=bacorder, vertices=vertices, edges=edgetab))
}
makePlotBySpecies <- function(bacnames, df_all, outdir, name, quantvar="IMC_log", 
                              quantvar_name = "log(IMC)",
                              corrmethod="pearson", w=8, h=10){
  assertthat::assert_that(all(bacnames %in% names(df_all)))
  vars2factor <- bacnames %>% gsub("_", " ", .) %>% 
    gsub("sp ", "sp. ", .) #%>% 
  # paste0("italic('", ., "')")
  
  vars2factor <- factor(vars2factor, levels = vars2factor)
  
  df2box <- df_all %>% dplyr::select(all_of(c(bacnames, "Condition", mediator_name))) %>% 
    gather(key="taxon", value="abundance", bacnames) %>% 
    mutate(taxon = gsub("_", " ", taxon) %>% 
             gsub("sp ", "sp. ", .) %>% 
             #paste0("italic('", ., "')") %>% 
             factor(levels=vars2factor)) %>% 
    mutate(Condition = ifelse(Condition == "Control", "C", "D"),
           Condition = factor(Condition, levels=c("C", "D")))
  dfmeans <- df2box %>% group_by(Condition, taxon) %>% summarise(abundance=mean(abundance))
  g1 <- ggplot(df2box, aes(x=Condition, y=abundance, col=Condition))+
    facet_grid( taxon ~ .)+
    geom_boxplot(outlier.alpha = 0, notch = F, width=0.1, size=1.5, varwidth = T)+
    geom_point(data=dfmeans, aes(x=Condition, y=abundance, size=1.8)) +
    coord_flip() +
    theme_minimal() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
    theme(strip.text.y = element_text(size = 10, 
                                      colour = "black", angle = 0, face = "italic")) +
    theme(axis.text.x = element_text(size = 8, 
                                     colour = "black", angle = 0, 
                                     face = "plain"))+
    theme(axis.text.y = element_text(size = 8, 
                                     colour = "black", angle = 0, 
                                     face = "plain")) +
    #theme(axis.text.y = element_blank())+
    ylab("Taxon abundance")
  ggsave(filename = paste0(opt$out, "/", name, ".pdf"), g1, width = w, height = h)
  
  mods <- do.call('rbind', map(bacnames, \(x){
    form <- as.formula(paste0(quantvar, " ~ ", x))
    mm <- lm(data=df_all, formula = form) %>% summary() 
    mm2 <- cbind(broom::glance(mm), 
                 data.frame(intercept=mm$coefficients[1], 
                            slope=mm$coefficients[2], 
                            taxon=x,
                            correlation = cor(df_all[, quantvar], df_all[, x], method = corrmethod, use="complete.obs"),
                            correlation_pvalue = cor.test(df_all[, quantvar], df_all[, x], method = corrmethod, use="complete.obs")$p.value
                            ))
    return(mm2)
  } 
  )) %>% mutate(taxon = gsub("_", " ", taxon) %>% 
                  gsub("sp ", "sp. ", .) %>% 
                  #paste0("italic('", ., "')") %>% 
                  factor(levels=vars2factor),
                color = ifelse(correlation_pvalue > 0.05, C_NS, 
                               ifelse(slope <0 , C_CTRL, C_CASE))
                )
  
  g2 <- ggplot(mods, aes(y=correlation, x=taxon))+
    facet_grid( taxon ~ ., scales="free")+
    geom_col(fill=mods$color)+
    coord_flip()+ 
    theme_minimal() +
    geom_hline(yintercept = 0, col=C_NS, linetype=2)+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
    #theme(strip.text.y = element_text(size = 10, 
    #                                  colour = "black", angle = 0, face = "italic")) +
    theme(strip.text.y = element_blank())+
    theme(axis.text.x = element_text(size = 8, 
                                     colour = "black", angle = 0, 
                                     face = "plain")) +
    theme(axis.text.y = element_text(size = 12, 
                                     colour = "black", angle = 0, 
                                     face = "italic"))+
    ylab(paste0(make_clean_names(corrmethod, case="big_camel"), " corr. ", quantvar_name))
    
  ggsave(filename = paste0(opt$out, "/", name,'_',quantvar ,"_barplot.pdf"), g2, width = w, height = h)
  
  g3 <- g1 + theme(strip.text.y = element_blank()) + theme(legend.position = "none")
  cw <- cowplot::plot_grid(plotlist=list(g2, g3), nrow = 1, rel_widths = c(2,1))
  pdf(paste0(opt$out, name, "_merged.pdf"), width = h, height = h)
  print(cw)
  dev.off()
  write_tsv(mods, file = paste0(opt$out, name,'_',quantvar ,"_correlations.tsv"))
  return(list(box=g1, bars=g2, cw=cw, correlations=mods))
}
