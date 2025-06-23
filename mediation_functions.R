library(tidyverse)
library(bmem)
library(sem)
library(ggrepel)
library(geomtextpath)
library(cowplot)
library(janitor)
library("wesanderson")


C_CASE = "#FD8B2F" #"rgba(200, 44, 44, 0.8)"
C_CASE2 = "tomato"
C_CASE_LINK = "#fBd895" #"#f9c784"
C_CTRL = "#A1C6EA" # "#006daa" #
C_CTRL2 = "steelblue2"
C_CTRL_LINK ="#DAE3E5" #"rgba(44, 44, 200, 0.8)"
C_WHITE= "#DDDDDD"
C_NS =  "#A5ABBD" #"rgba(224, 224, 224, 0.8)"
C_OTHER = "gray30"


pal1 <- wes_palette("Darjeeling2", 2)
pal <- wes_palette("AsteroidCity2")
C_POS_EFF = pal1[2]
C_NEG_EFF = pal[4]

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
  #nlsy.res<-bmem.bs(df, iris.model, effects)
  return(nlsy.res)

}

makeMediationSimpleLavaan <- function(df, xname, yname, medname){
  library(lavaan)

  eq1 = paste0(medname, " ~ a*", xname)
  eq2 = paste0(yname, " ~ b*", medname, " + cp*", xname)
  eq3 = "ab := a*b"
  eq4 = "total := cp + (a*b)"

  if(length(unique(df[, yname])) < 2){
    ss <- data.frame(
      Estimate = numeric(5),
      S.E. = numeric(5),
      `z-score` = numeric(5),
      p.value = rep(0.05001, 5)
    )
    names(ss)[3] <- "z-score"
    rownames(ss) <- c("a", "b", "cp", "a*b", "cp+a*b")
    return(list(sem=list(), estimates=ss))
  }
    
  df <- df %>% dplyr::mutate(!!yname := ordered(!!sym(yname)))
    

  eqs <- paste(eq1, eq2, eq3, eq4, sep="\n", collapse="\n")

  fit <- lavaan::sem(eqs, data = df)
  ss <- summary(fit)$pe %>%
    dplyr::select(label, est, se, z, pvalue) %>%
    dplyr::filter(label %in% c("a", "b", "cp", "ab", "total")) %>%
    dplyr::mutate(label = ifelse(label == "ab", "a*b", label),
                  label = ifelse(label == "total", "cp+a*b", label)) %>%
    column_to_rownames("label")
  names(ss) <- c("Estimate", "S.E.", "z-score", "p.value")

  return(list(sem=fit, estimates=ss))

}

makeMediationSimplePowerCurve <- function(df, xname, yname, medname, param_estimates = c(a=0.1, b=0.1, cp=0.1),
                                          nboot=1000, nrep=1000, nboot_curve=100, nrep_curve=100,
                                          min_n=100, max_n=500, inc_n=20){
  library(bmemLavaan)

  eq1 = paste0(yname, " ~ b*start(", param_estimates["b"], ")*", medname, " + cp*start(", param_estimates["a"], ")*", xname)
  eq2 = paste0(medname, " ~ a*start(", param_estimates["a"], ")*", xname)
  eq3 = "ab := a*b"
  #eq4 = "total := cp + (a*b)"
  eqs <- paste(eq1, eq2, eq3, sep="\n", collapse="\n")

  eq1 = paste0(medname, " ~ a*", xname)
  eq2 = paste0(yname, " ~ b*", medname, " + cp*", xname)
  eq3 = "ab:=a*b"
  eq4 = "total := cp + (a*b)"
  eqs <- paste(eq1, eq2, eq3, eq4, sep="\n", collapse="\n")

  pw <- power.basic(eqs, indirect = eq3, nobs = seq(100, 1000, 100), nrep=100)
  #eqs<-'
  #Condition_bin ~ b*start(5)*IMC_log + cp*start(0.1)*Actinomyces_naeslundii
  #IMC_log ~ a*start(0.2)*Actinomyces_naeslundii
  #ab := a*b
  #'

  power_res <- bmemLavaan::power.bmem(model = eqs,
                         nobs = nrow(df),
                         nrep = nrep,
                         nboot = nboot,
                         alpha=0.95, verbose = F)
  powsum <- summary(power_res)

  check_ssize <- seq(min_n, max_n, inc_n)
  pwc <- power.curve(model = eqs,
              nobs=check_ssize,
              method='normal', nrep=nrep_curve,
              nboot=nboot_curve, alpha=.95,interactive=FALSE)
  pwcnames <- colnames(pwc)
  pwcmeanings <- c("b", "cp", "a", "Ycov", "Mcov", "Xcov", "Indir")
  names(pwcmeanings) <- pwcnames
  pwcdf <- pwc %>% as.data.frame %>%
    dplyr::mutate(N = check_ssize) %>%
    gather("param_name", "power", -N) %>%
    dplyr::mutate(param = pwcmeanings[param_name],
                  Taxon = xname)

  return(list(pwres = power_res, pwsum=powsum, pwcurve = pwcdf, estimates=param_estimates))

}


makeMediationSimplePowerCurveLavaan <- function(df, xname, yname, medname, med_res,
                                          nrep=1000,
                                          min_n=100,
                                          max_n=500,
                                          inc_n=20,
                                          error=0.1,
                                          ncores=12){
  library(parallel)
  library(doParallel)
  library(foreach)

  cat("Doing Power Test for ", xname, "\n")
  cl <- makeCluster(ncores)
  registerDoParallel(cl)

  valsMed <- df %>% pull(!!sym(medname))
  valsX <- df %>% pull(!!sym(xname))

  ssizes <- seq(min_n, max_n, inc_n)

  eq1_lm <- paste0(medname, " ~ ", xname)
  mod1 = lm(as.formula(eq1_lm), data=df)
  mod1res <- resid(mod1)

  ressum_all <-  foreach(ssize = ssizes, .combine = bind_rows,
                         .export = c("makeMediationSimpleLavaan"),
                         .packages = c("dplyr", "tidyr", "magrittr", "tibble", "bmem", "sem", "lavaan")) %dopar%{

               cat("simulating ", xname, ", sample size=", ssize, "\n")

    #res <- foreach(numrep = 1:nrep, .combine = bind_rows,
    #               .packages = c("dplyr", "tidyr", "magrittr", "tibble", "bmem", "sem")) %dopar% {

      res <- data.frame()
      for(numrep in 1:nrep){
        X <- sample(valsX, ssize, replace=TRUE)
        M = med_res["a", "Estimate"]*X + rnorm(n = ssize , mean = 0, sd = sd(mod1res))
        M = med_res["a", "Estimate"]*X + sample(mod1res, ssize, replace=TRUE)

        logit_Y <- med_res["cp", "Estimate"] * X + med_res["b", "Estimate"] * M
        p_Y <- 1 / (1 + exp(-logit_Y))
        Y <- rbinom(ssize, 1, p_Y)

        #simdf <- data.frame(
        #  !!xname := X,
        #  !!medname := M,
        #  !!yname := Y
        #)

        simdf <- data.frame(
          X = X,
          M = M,
          Y = Y
        )
        colnames(simdf) <- c(xname, medname, yname)

        sim_sobelres <- tryCatch({

          makeMediationSimpleLavaan(simdf, xname = xname, yname = yname, medname = medname)$estimates %>%
              rownames_to_column("param") %>%
              filter(param %in% c("a", "b", "cp", "a*b", "cp+a*b ")) %>%
              dplyr::mutate(numrep = numrep)

        }, error = function(e) {
          warning(e)
          xx <- data.frame(param=character(0),
                     Estimate=numeric(0),
                     S.E.=numeric(0),
                     `z-score`=numeric(0),
                     p.value=numeric(0),
                     numrep=numeric(0))
          names(xx)[4] <- "z-score"
        }) 

        res <- res %>% bind_rows(sim_sobelres)
      }

      ressum <- res %>% group_by(param) %>%
        dplyr::mutate(
          p.value = ifelse(is.na(p.value), 1, p.value)
        ) %>%
        dplyr::summarise(
          ssize = ssize,
          mean_est = mean(Estimate, na.rm=T),
          mean_se = mean(`S.E.`, na.rm=T),
          mean_zscore = mean(`z-score`, na.rm=T),
          n = n(),
          power = sum(p.value < 0.05)/n()
        )
      ressum
  }
  stopCluster(cl)
  return(ressum_all)

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
  C_CTRL = C_CTRL2

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

getTotalColors <- function(bacorder, plim_plot, custom_colors,
                           col_ctrl = C_CTRL, col_case = C_CASE, col_ns = C_NS, col_other = C_OTHER){
  cols <- ifelse(bacorder$Estimate < 0, col_ctrl, col_case)
  if(is.null(custom_colors)){
    cols <- ifelse(bacorder$p.value <= plim_plot, cols, col_ns)
  }else{
    if(is.na(custom_colors$total)){
      cols <- ifelse(bacorder$p.value <= plim_plot, cols, col_ns)
    }else if(!custom_colors$total){
      cols <- rep(col_ns, nrow(bacorder))
    }
  }
  cols <- c(cols, col_case, col_other)
  return(cols)
}


getEdgeColors <- function(param, p.value, Estimate, plim_plot, custom_colors,
                          col_ctrl = C_CTRL, col_case = C_CASE, col_ns = C_NS){
  cols <- ifelse(Estimate < 0, col_ctrl, col_case)
  if(is.null(custom_colors)){
    cols <- ifelse(p.value > plim_plot, col_ns, cols)
  }else{
    cols <- ifelse(param=="b",
                   cols,
                   ifelse(grepl("^a", param),
                          ifelse(is.na(custom_colors$indirect),
                                 ifelse(p.value > plim_plot, col_ns, cols),
                                 ifelse(custom_colors$indirect, cols, col_ns)
                          ),
                          ifelse(is.na(custom_colors$direct),
                                 ifelse(p.value > plim_plot, col_ns, cols),
                                 ifelse(custom_colors$direct, cols, col_ns)
                          )
                   )
    )
  }
  return(cols)
}


plotIMCMediationSimple <- function(res, vars2test, outdir, outname,
                                   plim_plot = 0.05, use_color_scale=FALSE, w=14, h=10,custom_colors=NULL,
                                   name = "plotIMCMediationSimple"){
  C_CTRL = C_CTRL2
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
                         colors = getTotalColors(bacorder, plim_plot, custom_colors, col_ctrl = C_POS_EFF, col_case = C_NEG_EFF),
                         xoffset = ifelse(vnames == "D", 0.5, ifelse(vnames== "BMI", 1, 0))
  )

  #print(vertices$colors)
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
      color = getEdgeColors(param, p.value, Estimate, plim_plot, custom_colors, col_ctrl = C_POS_EFF, col_case = C_NEG_EFF),
      #linetype = ifelse(color == C_NS, 2, 1),
      width = abs(Estimate)*edge_width_factor,
      width = ifelse(width <  minval, minval,
                     ifelse(width > rangefactor*minval,
                            rangefactor*minval, width))
    ) %>%
    dplyr::mutate(linetype = ifelse(color == C_NS, 2, 1)) %>%
    dplyr::select(from, to, everything())
  #print(edgetab$color)
  edgetab$y1[edgetab$to=="D" & edgetab$from=="BMI"] = edgetab$y1[edgetab$to=="D" & edgetab$from=="BMI"]-0.75
  edgetab$x1[edgetab$to=="D" & edgetab$from=="BMI"] = edgetab$x1[edgetab$to=="D" & edgetab$from=="BMI"]+0.1

  edgetab$x1[edgetab$to=="D" & edgetab$from!="BMI"] = edgetab$x1[edgetab$to=="D" & edgetab$from=="BMI"]-0.5

  edgetab$y1[edgetab$to=="BMI"] =  edgetab$y1[edgetab$to=="BMI"] + 0.7
  edgetab$x1[edgetab$to=="BMI"] = edgetab$x1[edgetab$to=="BMI"] - 0.7

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
  ofname <- paste0(outdir, "/", name,'_p', as.character(plim_plot), ifelse(is.null(custom_colors),"", "CS"), ".pdf")
  ggsave(filename = ofname, g1,
         width = w, height = h)
  return(list(plot=g1, bacorder=bacorder, vertices=vertices, edges=edgetab))
}

plotAgeMediationSimple <- function(res, vars2test, outdir, outname,
                                   plim_plot = 0.05, use_color_scale=FALSE, w=14, h=6){

  C_CTRL = C_CTRL2
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
  ggsave(filename = paste0(outdir, "/", name,'_p', as.character(plim_plot), ".pdf"), g1,
         width = w, height = h)
  return(list(plot=g1, bacorder=bacorder, vertices=vertices, edges=edgetab))
}


plotAgeAndIMCMediationComplex <- function(res, vars2test, outdir, outname,
                                   plim_plot = 0.05, use_color_scale=FALSE, w=14, h=10){

  C_CTRL = C_CTRL2
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
  ggsave(filename = paste0(outdir, "/", name,'_p', as.character(plim_plot), ".pdf"), g1,
         width = w, height = h)

  return(list(plot=g1, bacorder=bacorder, vertices=vertices, edges=edgetab))
}
makePlotBySpecies <- function(bacnames, df_all, outdir, name, quantvar="IMC_log",
                              quantvar_name = "log(IMC)",
                              corrmethod="pearson", plim_plot=0.05, w=8, h=10){
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
  dfmeans <- df2box %>% group_by(Condition, taxon) %>% dplyr::summarise(abundance=mean(abundance))
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
  ggsave(filename = paste0(outdir, "/", name, ".pdf"), g1, width = w, height = h)

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
                color = ifelse(correlation_pvalue > plim_plot, C_NS,
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
                                     vjust = 0.5,
                                     face = "italic"))+
    ylab(paste0(make_clean_names(corrmethod, case="big_camel"), " corr. ", quantvar_name))+
    thin_barplot_lines

  ggsave(filename = paste0(outdir, "/", name,'_',quantvar ,"_barplot.pdf"), g2, width = w, height = h)

  g3 <- g1 + theme(strip.text.y = element_blank()) + theme(legend.position = "none")
  cw <- cowplot::plot_grid(plotlist=list(g2, g3), nrow = 1, rel_widths = c(2,1))
  pdf(paste0(outdir, name, "_merged.pdf"), width = w, height = h)
  print(cw)
  dev.off()
  write_tsv(mods, file = paste0(outdir, name,'_',quantvar ,"_correlations.tsv"))
  return(list(box=g1, bars=g2, cw=cw, correlations=mods))
}


makePlotBySpeciesEffects_BMI <- function(bacnames, df_all, res_final, outdir, name,
                               w=8, h=10, plim_plot=0.05, fix_limits=FALSE){
  assertthat::assert_that(all(bacnames %in% names(df_all)))
  vars2factor <- bacnames %>% gsub("_", " ", .) %>%
    gsub("sp ", "sp. ", .) #%>%
  # paste0("italic('", ., "')")

  vars2factor <- factor(vars2factor, levels = vars2factor)

  df2box <- df_all %>% dplyr::select(all_of(c(bacnames, "Condition", mediator_name))) %>%
    gather(key="taxon", value="abundance", bacnames) %>%
    dplyr::mutate(taxon = gsub("_", " ", taxon) %>%
             gsub("[\\[\\]]", "", ., perl=T) %>%
             gsub("sp ", "sp. ", .) %>%
             #paste0("italic('", ., "')") %>%
             factor(levels=vars2factor)) %>%
    dplyr::mutate(Condition = ifelse(Condition == "Control", "C", "D"),
           Condition = factor(Condition, levels=c("C", "D")))
  dfmeans <- df2box %>% group_by(Condition, taxon) %>% dplyr::summarise(abundance=mean(abundance))
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
  ggsave(filename = paste0(outdir, "/", name, ".pdf"), g1, width = w, height = h)


  res_mod <- res_final %>% dplyr::mutate(
    param = gsub("[0-9]+", "", param),
    Effect = ifelse(param=="a*b", "Indirect (a*b)", ifelse(param=="cp", "Direct (c')", ifelse(param=="cp+a*b", "Total (c' + a*b)", "Other"))),
  ) %>% dplyr::filter(Effect!="Other") %>%
    dplyr::mutate(taxon = gsub("_", " ", x_labels) %>%
                    gsub("[\\[\\]]", "", ., perl=T) %>%
                    gsub("sp ", "sp. ", .) %>%
                    #paste0("italic('", ., "')") %>%
                    factor(levels=vars2factor),
                  color = ifelse(p.value <= plim_plot, ifelse(Estimate< 0, "Neg", "Pos") , "NS"),
                  color = factor(color, levels = c("Neg", "Pos", "NS")),
                  color2 = ifelse(p.value <= plim_plot, ifelse(Estimate< 0, C_CTRL, C_CASE) , C_NS),
                  Estimate = 10*Estimate
                  )
  val_lims = c(min(res_mod$Estimate), max(res_mod$Estimate))
  plots <- map(c("Total (c' + a*b)", "Direct (c')", "Indirect (a*b)"), \(xx){
    ptab <- res_mod %>% dplyr::filter(Effect == xx) %>%
      dplyr::mutate(taxon = factor(taxon, levels = rev(vars2factor)))
    g2 <- ggplot(ptab, aes(y=Estimate, x=taxon))+ #, fill=color
      geom_col(fill=ptab$color2)+
      coord_flip()+
      scale_fill_manual(values = c(C_CTRL, C_CASE, C_NS))+
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
                                     vjust = 0.5,
                                      face = "italic"))+
      ylab("")+
      thin_barplot_lines
    if(fix_limits){
      g2 <- g2 + ylim(val_lims)
    }
    return(g2)
  })

  ga = plots[[1]] + theme(legend.position = "none") + ylab("Total effects (c' + a*b)")
  gb = plots[[2]] + theme(legend.position = "none")  + theme(axis.text.y = element_blank()) + xlab("") + ylab("Direct Effects (c')")
  gc = plots[[3]] + theme(legend.position = "none")  + theme(axis.text.y = element_blank())+ xlab("") + ylab("Indirect Effects (a*b)")
  gd <- g1 + theme(legend.position = "none") + theme(strip.text.y = element_blank())
  cw <- cowplot::plot_grid(plotlist=list(ga, gb, gc, gd), nrow = 1, rel_widths = c(2,1, 1,1))
  pdf(paste0(outdir, name, "_BarplotEffectsmerged.pdf"), width = w, height = h)
  print(cw)
  dev.off()
  write_tsv(res_mod, file = paste0(outdir, name ,"_BarplotEffectsmerged.tsv"))
  return(list(box=g1, bars=plots, cw=cw, tab=res_mod))
}


makePlotBySpeciesEffects_BMI2 <- function(bacnames, df_all, res_final, outdir, name,
                                         w=8, h=10, plim_plot=0.05, fix_limits=FALSE){

  titlesize = ifelse(h <=10, 12, h+2)
  assertthat::assert_that(all(bacnames %in% names(df_all)))
  vars2factor <- bacnames %>% gsub("_", " ", .) %>%
    gsub("sp ", "sp. ", .) #%>%
  # paste0("italic('", ., "')")

  vars2factor <- factor(vars2factor, levels = vars2factor)

  df2box <- df_all %>% dplyr::select(all_of(c(bacnames, "Condition", mediator_name))) %>%
    gather(key="taxon", value="abundance", bacnames) %>%
    dplyr::mutate(taxon = gsub("_", " ", taxon) %>%
                    gsub("[\\[\\]]", "", ., perl=T) %>%
                    gsub("sp ", "sp. ", .) %>%
                    #paste0("italic('", ., "')") %>%
                    factor(levels=vars2factor)) %>%
    dplyr::mutate(Condition = ifelse(Condition == "Control", "C", "D"),
                  Condition = factor(Condition, levels=c("C", "D")))
  dfmeans <- df2box %>% group_by(Condition, taxon) %>% dplyr::summarise(abundance=mean(abundance))
  g1 <- ggplot(df2box, aes(x=Condition, y=abundance, col=Condition))+
    facet_grid( taxon ~ .)+
    geom_boxplot(outlier.alpha = 0, notch = F, width=0.1, size=1.5, varwidth = T)+
    geom_point(data=dfmeans, aes(x=Condition, y=abundance, size=1.8)) +
    coord_flip() +
    theme_minimal() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
    theme(strip.text.y = element_text(size = 10,
                                      colour = "black", angle = 0, face = "italic")) +
    theme(axis.text.x = element_text(size = ifelse(h<=8, 8, 10),
                                     colour = "black", angle = 0,
                                     face = "plain"))+
    theme(axis.text.y = element_text(size = ifelse(h<=8, 8, 10),
                                     colour = "black", angle = 0,
                                     face = "plain")) +
    #theme(axis.text.y = element_blank())+
    theme(axis.title = element_text(size = titlesize)) +
    ylab("Taxon abundance")
  ggsave(filename = paste0(outdir, "/", name, ".pdf"), g1, width = w, height = h)


  res_mod <- res_final %>% dplyr::mutate(
    param = gsub("[0-9]+", "", param),
    Effect = ifelse(param=="a*b", "Indirect (a*b)",
                    ifelse(param=="cp", "Direct (c')",
                           ifelse(param=="cp+a*b", "Total (c' + a*b)",
                                  ifelse(param=="a", "Effects on BMI (a)", "Other")))),
  ) %>% dplyr::filter(Effect!="Other") %>%
    dplyr::mutate(taxon = gsub("_", " ", x_labels) %>%
                    gsub("[\\[\\]]", "", ., perl=T) %>%
                    gsub("sp ", "sp. ", .) %>%
                    #paste0("italic('", ., "')") %>%
                    factor(levels=vars2factor),
                  color = ifelse(p.value <= plim_plot, ifelse(Estimate< 0, "Neg", "Pos") , "NS"),
                  color = factor(color, levels = c("Neg", "Pos", "NS")),
                  color2 = ifelse(p.value <= plim_plot, ifelse(Estimate< 0, C_POS_EFF, C_NEG_EFF) , C_NS),
                  Estimate = 10*Estimate
    )
  val_lims = c(min(res_mod$Estimate), max(res_mod$Estimate))
  plots <- map(c("Total (c' + a*b)", "Direct (c')", "Effects on BMI (a)", "Indirect (a*b)"), \(xx){
    ptab <- res_mod %>% dplyr::filter(Effect == xx) %>%
      dplyr::mutate(taxon = factor(taxon, levels = rev(vars2factor)))
    g2 <- ggplot(ptab, aes(y=Estimate, x=taxon))+ #, fill=color
      geom_col(fill=ptab$color2)+
      coord_flip()+
      scale_fill_manual(values = c(C_CTRL, C_CASE, C_NS))+
      theme_minimal() +
      geom_hline(yintercept = 0, col=C_NS, linetype=2)+
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
      #theme(strip.text.y = element_text(size = 10,
      #                                  colour = "black", angle = 0, face = "italic")) +
      theme(strip.text.y = element_blank())+
      theme(axis.text.x = element_text(size =  ifelse(h<=10, 8, 10), # 8,
                                       colour = "black", angle = 0,
                                       face = "plain")) +
      theme(axis.text.y = element_text(size = titlesize, # 12,
                                       colour = "black", angle = 0,
                                       vjust = 0.5,
                                       face = "italic"))+
      theme(axis.title = element_text(size = titlesize))+
      ylab("")+
      thin_barplot_lines
    if(fix_limits){
      g2 <- g2 + ylim(val_lims)
    }
    return(g2)
  })

  ga = plots[[1]] + theme(legend.position = "none") + ylab("Total effects (c' + a*b)")
  gb = plots[[2]] + theme(legend.position = "none")  + theme(axis.text.y = element_blank()) + xlab("") + ylab("Direct Effects (c')")
  gc = plots[[3]] + theme(legend.position = "none")  + theme(axis.text.y = element_blank())+ xlab("") + ylab("Effects on BMI (a)")
  gd = plots[[4]] + theme(legend.position = "none")  + theme(axis.text.y = element_blank())+ xlab("") + ylab("Indirect Effects (a*b)")
  ge <- g1 + theme(legend.position = "none") + theme(strip.text.y = element_blank())
  cw <- cowplot::plot_grid(plotlist=list(ga, gb, gc, gd, ge), nrow = 1, rel_widths = c(2.5,1, 1, 1,1))
  pdf(paste0(outdir, name, "_BarplotEffectsmerged.pdf"), width = w, height = h)
  print(cw)
  dev.off()
  write_tsv(res_mod, file = paste0(outdir, name ,"_BarplotEffectsmerged.tsv"))
  return(list(box=g1, bars=plots, cw=cw, tab=res_mod))
}


makeBarplotDAA <- function(daalist, outdir, plim=0.05, name=""){
  daatab <- lapply(names(daalist), \(x){daalist[[x]]$contrast <- x; return(daalist[[x]])}) %>% bind_rows() %>%
    dplyr::select(taxon, log2FoldChangeShrink, padj, contrast) %>%
    dplyr::mutate(padj = ifelse(is.na(padj), 1, padj)) %>%
    gather(key="param", value ="value", padj, log2FoldChangeShrink) %>%
    tidyr::unite(col="tmp", contrast, param, sep="__") %>%
    spread(tmp, value) %>%
    dplyr::mutate(
      CondAdjusted = D_vs_C_adj_BMIplusAge__padj <= plim & D_vs_C__padj <= plim,
      CondOnly = D_vs_C_adj_BMIplusAge__padj > plim & D_vs_C__padj <= plim,
      CondAdjustedOnly = D_vs_C_adj_BMIplusAge__padj <= plim & D_vs_C__padj > plim,
      BMIAdjusted = BMI_adj_DeprplusAge__padj <= plim,
      AgeAdjusted = Age_adj_DeprplusBMI__padj <= plim
    ) %>%
    gather(key="param", value="value", -taxon, -CondAdjusted, -CondOnly, -CondAdjustedOnly, -BMIAdjusted, -AgeAdjusted ) %>%
    separate(param, into=c("Contrast", "param"), sep="__") %>%
    spread(key=param, value = value) %>%
    dplyr::mutate(Contrast= gsub("plus", "+", Contrast),
                  Contrast = gsub("_", " ", Contrast),
                  Contrast = gsub("adj", "adj.", Contrast),
                  taxon = gsub("_", " ", taxon),
                  taxon = gsub("[\\[\\]]", "", taxon),
                  color = ifelse(padj < plim, ifelse(log2FoldChangeShrink < 0, "Down", "Up"), "NS"),
                  color = factor(color, levels=c("Down", "Up", "NS")))

  orderedlevels <- daatab %>% filter(Contrast == "D vs C") %>% arrange(log2FoldChangeShrink) %>% pull(taxon)
  contrastlevels <- c("D vs C", "D vs C adj. BMI+Age", "BMI adj. Depr+Age", "Age adj. Depr+BMI")
  daatab <- daatab %>% mutate(taxon=factor(taxon, levels=orderedlevels),
                              Contrast = factor(Contrast, levels = contrastlevels))

  vars2separate <- c("CondAdjusted", "CondOnly", "CondAdjustedOnly", "BMIAdjusted", "AgeAdjusted")
  vars2separatenames <- c("Different in D vs C after adjusting", "Not significant after BMI+Age correction",
                          "Significant only after correction", "Significant BMI after Depr+Age correction",
                          "Significant Age after Depr+BMI correction")
  plots <- map2(vars2separate, vars2separatenames,\(x, y){
    tmptab <- daatab %>% dplyr::filter(!!sym(x))
    g1 <- ggplot(tmptab, aes(x=taxon, y=log2FoldChangeShrink, fill=color))+
      facet_grid( . ~ Contrast) +
      geom_hline(yintercept = 0, linetype=2, col=C_NS) +
      scale_fill_manual(values=c(C_CTRL, C_CASE, C_NS)) +
      geom_col() +
      coord_flip() +
      theme_minimal() +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
      theme(strip.text.y = element_text(size = 10,
                                      colour = "black", angle = 0, face = "italic")) +
      theme(axis.text.x = element_text(size = 8,
                                     colour = "black", angle = 0,
                                     face = "plain"))+
      theme(axis.text.y = element_text(size = 8,
                                     colour = "black", angle = 0,
                                     vjust = 0.5,
                                     face = "italic")) +
      theme(legend.position = 'none')+
    #theme(axis.text.y = element_blank())+
    xlab(y)+
    ylab("") +
    ylim(c(-9, 9))+
      thin_barplot_lines
  return(list(plot=g1, size=length(unique(tmptab$taxon))))
})
  names(plots) <- vars2separate
  cw <- cowplot::plot_grid(plotlist = list(plots$CondAdjusted$plot,
                                     plots$CondOnly$plot + theme(strip.text.x = element_blank()),
                                     plots$CondAdjustedOnly$plot+ theme(strip.text.x = element_blank()) + ylab("LFC")
                                     ),
                       rel_heights = c(plots$CondAdjusted$size, plots$CondOnly$size, plots$CondAdjustedOnly$size),
                     ncol=1)

  filename = paste0(outdir, "/", name, "_Barplot_LFC_all1.pdf")
  pdf(filename, w=8, h=16)
  print(cw)
  dev.off()
  ggsave(filename = paste0(outdir, "/", name, "_Barplot_LFC_all_BMI.pdf"), plots$BMIAdjusted$plot, w=8, h=12)
  ggsave(filename = paste0(outdir, "/", name, "_Barplot_LFC_all_Age.pdf"), plots$AgeAdjusted$plot, w=8, h=6)
  write_tsv(daatab, file =  paste0(outdir, "/", name, "_Barplot_LFC_all1.tsv"))
  return(list(plots=plots, tab=daatab))
}

makeBarplotDAA2 <- function(daalist, outdir, plim=0.05, name=""){
  daatab <- lapply(names(daalist), \(x){daalist[[x]]$contrast <- x; return(daalist[[x]])}) %>% bind_rows() %>%
    dplyr::select(taxon, log2FoldChangeShrink, padj, contrast) %>%
    dplyr::mutate(padj = ifelse(is.na(padj), 1, padj)) %>%
    gather(key="param", value ="value", padj, log2FoldChangeShrink) %>%
    tidyr::unite(col="tmp", contrast, param, sep="__") %>%
    spread(tmp, value) %>%
    dplyr::mutate(
      CondAdjusted = D_vs_C_adj_BMI__padj <= plim & D_vs_C__padj <= plim,
      CondOnly = D_vs_C_adj_BMI__padj > plim & D_vs_C__padj <= plim,
      CondAdjustedOnly = D_vs_C_adj_BMI__padj <= plim & D_vs_C__padj > plim,
      BMIAdjusted = BMI_adj_Depr__padj <= plim & D_vs_C_adj_BMI__padj > plim &  D_vs_C__padj > plim
    ) %>%
    gather(key="param", value="value", -taxon, -CondAdjusted, -CondOnly, -CondAdjustedOnly, -BMIAdjusted) %>%
    separate(param, into=c("Contrast", "param"), sep="__") %>%
    spread(key=param, value = value) %>%
    dplyr::mutate(Contrast= gsub("plus", "+", Contrast),
                  Contrast = gsub("_", " ", Contrast),
                  Contrast = gsub("adj", "adj.", Contrast),
                  Contrast = gsub("Depr", "Depr.", Contrast),
                  taxon = gsub("_", " ", taxon),
                  taxon = gsub("[\\[\\]]", "", taxon),
                  color = ifelse(padj < plim, ifelse(log2FoldChangeShrink < 0, "Down", "Up"), "NS"),
                  color = factor(color, levels=c("Down", "Up", "NS")))

  orderedlevels <- daatab %>% filter(Contrast == "D vs C") %>% arrange(log2FoldChangeShrink) %>% pull(taxon)
  contrastlevels <- c("D vs C", "D vs C adj. BMI", "BMI", "BMI adj. Depr.")
  daatab <- daatab %>% mutate(taxon=factor(taxon, levels=orderedlevels),
                              Contrast = factor(Contrast, levels = contrastlevels))

  vars2separate <- c("CondAdjusted", "CondOnly", "CondAdjustedOnly", "BMIAdjusted")
  vars2separatenames <- c("Different in D vs C after adjusting", "Not significant after BMI correction",
                          "Significant only after correction", "Significant BMI after Depr correction")
  plots <- map2(vars2separate, vars2separatenames,\(x, y){
    tmptab <- daatab %>% dplyr::filter(!!sym(x))
    g1 <- ggplot(tmptab, aes(x=taxon, y=log2FoldChangeShrink, fill=color))+
      facet_grid( . ~ Contrast) +
      geom_hline(yintercept = 0, linetype=2, col=C_NS) +
      scale_fill_manual(values=c(C_CTRL, C_CASE, C_NS)) +
      geom_col() +
      coord_flip() +
      theme_minimal() +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
      theme(strip.text.y = element_text(size = 10,
                                        colour = "black", angle = 0, face = "italic")) +
      theme(axis.text.x = element_text(size = 8,
                                       colour = "black", angle = 0,
                                       face = "plain"))+
      theme(axis.text.y = element_text(size = 8,
                                       colour = "black", angle = 0,
                                       vjust = 0.5,
                                       face = "italic")) +
      theme(legend.position = 'none')+
      #theme(axis.text.y = element_blank())+
      xlab(y)+
      ylab("") +
      ylim(c(-9, 9))+
      thin_barplot_lines
    return(list(plot=g1, size=length(unique(tmptab$taxon))))
  })
  names(plots) <- vars2separate
  cw <- cowplot::plot_grid(plotlist = list(plots$CondAdjusted$plot,
                                           plots$CondOnly$plot + theme(strip.text.x = element_blank()),
                                           plots$CondAdjustedOnly$plot+ theme(strip.text.x = element_blank()) + ylab("LFC")
  ),
  rel_heights = c(plots$CondAdjusted$size, plots$CondOnly$size, plots$CondAdjustedOnly$size),
  ncol=1)

  filename = paste0(outdir, "/", name, "_Barplot_LFC_all1.pdf")
  pdf(filename, w=8, h=16)
  print(cw)
  dev.off()
  ggsave(filename = paste0(outdir, "/", name, "_Barplot_LFC_all_BMI.pdf"), plots$BMIAdjusted$plot, w=8, h=12)
  write_tsv(daatab, file =  paste0(outdir, "/", name, "_Barplot_LFC_all1.tsv"))

  ## Cowplot including BMI
  cw <- cowplot::plot_grid(plotlist = list(plots$CondOnly$plot + theme(strip.text.x = element_blank()),
                                           plots$CondAdjustedOnly$plot+ theme(strip.text.x = element_blank()),
                                           plots$BMIAdjusted$plot + theme(strip.text.x = element_blank()) + ylab("LFC")
  ),
  rel_heights = c(plots$CondOnly$size, plots$CondAdjustedOnly$size, plots$BMIAdjusted$size),
  ncol=1)

  filename = paste0(outdir, "/", name, "_Barplot_LFC_all2.pdf")
  pdf(filename, w=8, h=20)
  print(cw)
  dev.off()
  return(list(plots=plots, tab=daatab))
}

###################################################
makeBarplotDAA3_Int <- function(daalist, outdir, plim=0.05, name="",
                                cond_names=list()){

  if(length(cond_names) == 0) cond_names <- names(daalist)
  cond_names <- paste0(cond_names, "__padj")

  daatab <- lapply(names(daalist), \(x){daalist[[x]]$contrast <- x; return(daalist[[x]])}) %>% bind_rows() %>%
    dplyr::select(taxon, log2FoldChangeShrink, padj, contrast) %>%
    dplyr::mutate(padj = ifelse(is.na(padj), 1, padj)) %>%
    gather(key="param", value ="value", padj, log2FoldChangeShrink) %>%
    tidyr::unite(col="tmp", contrast, param, sep="__") %>%
    spread(tmp, value) %>%
    dplyr::mutate(
      AllSig = !!sym(cond_names[1]) <= plim & !!sym(cond_names[2]) <= plim & !!sym(cond_names[3]) <= plim & !!sym(cond_names[4]) <= plim,
      CondAdjusted = !!sym(cond_names[1]) <= plim & !!sym(cond_names[2]) <= plim,
      CondOnly = !!sym(cond_names[2])  > plim & !!sym(cond_names[1]) <= plim,
      CondAdjustedOnly = !!sym(cond_names[2]) <= plim & !!sym(cond_names[1]) > plim,
      BMIAdjusted = !!sym(cond_names[4]) <= plim & !!sym(cond_names[2]) > plim &  !!sym(cond_names[1]) > plim,

      #CondAdjusted = CondAdjusted & !AllSig,
      #CondOnly = CondAdjusted & !AllSig,
      #CondAdjustedOnly = CondAdjusted & !AllSig,
      #BMIAdjusted = CondAdjusted & !AllSig

    ) %>%
    gather(key="param", value="value", -taxon, -AllSig, -CondAdjusted, -CondOnly, -CondAdjustedOnly, -BMIAdjusted) %>%
    separate(param, into=c("Contrast", "param"), sep="__") %>%
    spread(key=param, value = value) %>%
    dplyr::mutate(Contrast= gsub("plus", "+", Contrast),
                  Contrast= gsub("Inter", "*", Contrast),
                  Contrast = gsub("_", " ", Contrast),
                  Contrast = gsub("adj", "adj.", Contrast),
                  Contrast = gsub("Depr", "Depr.", Contrast),
                  taxon = gsub("_", " ", taxon),
                  taxon = gsub("[\\[\\]]", "", taxon),
                  color = ifelse(padj < plim, ifelse(log2FoldChangeShrink < 0, "Down", "Up"), "NS"),
                  color = factor(color, levels=c("Down", "Up", "NS")))

  orderedlevels <- daatab %>% filter(Contrast == "D vs C") %>% arrange(log2FoldChangeShrink) %>% pull(taxon)
  contrastlevels <- c("D vs C", "D vs C adj. Sex*Age + BMI", "BMI", "BMI adj. Sex*Age + MDD")
  daatab <- daatab %>% mutate(taxon=factor(taxon, levels=orderedlevels),
                              Contrast = factor(Contrast, levels = contrastlevels))

  vars2separate <- c("AllSig", "CondAdjusted", "CondOnly", "CondAdjustedOnly", "BMIAdjusted")
  vars2separatenames <- c("All significant",
                          "Different in D vs C after adjusting",
                          "Not significant after correction",
                          "Significant only after correction",
                          "Significant BMI after MDD correction")
  plots <- map2(vars2separate, vars2separatenames,\(x, y){
    tmptab <- daatab %>% dplyr::filter(!!sym(x))
    g1 <- ggplot(tmptab, aes(x=taxon, y=log2FoldChangeShrink, fill=color))+
      facet_grid( . ~ Contrast) +
      geom_hline(yintercept = 0, linetype=2, col=C_NS) +
      scale_fill_manual(values=c(C_CTRL, C_CASE, C_NS)) +
      geom_col() +
      coord_flip() +
      theme_minimal() +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
      theme(strip.text.y = element_text(size = 10,
                                        colour = "black", angle = 0, face = "italic")) +
      theme(axis.text.x = element_text(size = 8,
                                       colour = "black", angle = 0,
                                       face = "plain"))+
      theme(axis.text.y = element_text(size = 8,
                                       colour = "black", angle = 0,
                                       vjust = 0.5,
                                       face = "italic")) +
      theme(legend.position = 'none')+
      #theme(axis.text.y = element_blank())+
      xlab(y)+
      ylab("") +
      #ylim(c(-9, 9))+
      thin_barplot_lines
    return(list(plot=g1, size=length(unique(tmptab$taxon))))
  })
  names(plots) <- vars2separate
  cw <- cowplot::plot_grid(plotlist = list(plots$CondAdjusted$plot,
                                           plots$CondOnly$plot + theme(strip.text.x = element_blank()),
                                           plots$CondAdjustedOnly$plot+ theme(strip.text.x = element_blank()) + ylab("LFC")
  ),
  rel_heights = c(plots$CondAdjusted$size,
                  plots$CondOnly$size,
                  plots$CondAdjustedOnly$size),
  ncol=1)

  filename = paste0(outdir, "/", name, "_Barplot_LFC_all1.pdf")
  pdf(filename, w=8, h=18)
  print(cw)
  dev.off()

  ggsave(filename = paste0(outdir, "/", name, "_Barplot_LFC_all_BMI.pdf"), plots$BMIAdjusted$plot, w=8, h=8)
  write_tsv(daatab, file =  paste0(outdir, "/", name, "_Barplot_LFC_all1.tsv"))

  ## Cowplot including BMI
  cw <- cowplot::plot_grid(plotlist = list(plots$CondAdjusted$plot,
                                           plots$CondOnly$plot + theme(strip.text.x = element_blank()),
                                           plots$CondAdjustedOnly$plot+ theme(strip.text.x = element_blank()),
                                           plots$BMIAdjusted$plot + theme(strip.text.x = element_blank()) + ylab("LFC")
  ),
  rel_heights = c(plots$CondAdjusted$size, plots$CondOnly$size, plots$CondAdjustedOnly$size, plots$BMIAdjusted$size),
  ncol=1)

  filename = paste0(outdir, "/", name, "_Barplot_LFC_all2.pdf")
  pdf(filename, w=8, h=20)
  print(cw)
  dev.off()

  #cw <- cowplot::plot_grid(plotlist = list(plots$AllSig$plot,
  #                                         plots$CondAdjusted$plot + theme(strip.text.x = element_blank()),
  #                                         plots$CondOnly$plot + theme(strip.text.x = element_blank()),
  #                                         plots$CondAdjustedOnly$plot+ theme(strip.text.x = element_blank()),
  #                                         plots$BMIAdjusted$plot + theme(strip.text.x = element_blank()) + ylab("LFC")
  #),
  #rel_heights = c(plots$AllSig$size,
  #                plots$CondAdjusted$size,
  #                plots$CondOnly$size,
  #                plots$CondAdjustedOnly$size,
  #                plots$BMIAdjusted$size),
  #ncol=1)
  #
  #filename = paste0(outdir, "/", name, "_Barplot_LFC_all3.pdf")
  #pdf(filename, w=8, h=20)
  #print(cw)
  #dev.off()
  return(list(plots=plots, tab=daatab))
}

makeVenn <- function(vars2venn, name="VennDiagram", opt, w=5, h=5){
  gv <- ggvenn(
    vars2venn, columns = names(vars2venn),
    stroke_size = 0.5,
    stroke_color = C_NS,
    fill_color = c(C_CASE, C_CASE2, C_CTRL2, C_CTRL),show_elements = F
  )
  ggsave(filename = paste0(opt$out, name, ".pdf"), gv, width = w, height = h)
  return(gv)
}

makeFullMediationAnalysisIMC <- function(vstdf, df_all, opt, getVarsFunction, mediator_name="IMC", y_name="Condition_bin",
                                         plim=0.05, plim_plot=0.05, name="analysis_IMC_separateModel_vjust",
                                         wnet=14, hnet=12, wbars=8, hbars=10, wbars2=10, hbars2=12, use_color_scale=FALSE,
                                         fix_barplot_limits=FALSE, custom_colors=NULL, make_boxplots=TRUE, list2merge=NULL,
                                         make_power_test = FALSE, min_power=0.9){
  print(paste0("OUTPUT: ", opt$out))
  summary_df <- data.frame(variable = gsub("-", ".", vstdf$gene) %>%
                             gsub("\\[|\\]", "", .) %>%
                             gsub("\\._", "_", .) %>%
                             gsub("\\/|\\(|\\)", ".", .) %>%
                             gsub("\\.$", "", .)
  )

  assertthat::assert_that(all(summary_df$variable %in% names(df_all)))

  if(is.null(list2merge)){
    list2merge <- list(
      depr_only_padj = dea2contrasts$firstContrast$resdf,
      imc_only_padj = dea2contrasts$contrastlist2$BMI_alone$resdf,
      depr_adjimc_padj = dea2contrasts$contrastlist2$Condition_corrBMI$resdf,
      #depr_adjimc_padj = dea2contrasts$contrastlist2$Condition_corrIMC$resdf,
      imc_adjdepr_padj = dea2contrasts$contrastlist2$BMI_corrCond$resdf
    )
  }

  merged_pvals <- list2merge %>%
    lapply(\(x){
      x$taxon <- x$taxon %>%
        gsub("-", ".", .) %>%
        gsub("\\[|\\]", "", .) %>%
        gsub("\\._", "_", .) %>%
        gsub("\\/|\\(|\\)", ".", .) %>%
        gsub("\\.$", "", .)
      x <- x[match(summary_df$variable, x$taxon), ]
      #x$padj[is.na(x$padj)] <- 1
      return(x$padj)
    }) %>% bind_cols()
  summary_df <- cbind(summary_df, merged_pvals)

  vars2test <- getVarsFunction(summary_df)
  vars2venn <- list(
    "D vs C" = summary_df %>% dplyr::filter(depr_only_padj < plim & !is.na(depr_only_padj)) %>% pull(variable),
    "D vs C adj. BMI" = summary_df %>% dplyr::filter(depr_adjimc_padj < plim & ! is.na(depr_adjimc_padj)) %>% pull(variable),
    "BMI" = summary_df %>% dplyr::filter(imc_only_padj < plim & ! is.na(imc_only_padj)) %>% pull(variable),
    "BMI adj. Depr" = summary_df %>% dplyr::filter(imc_adjdepr_padj < plim & !is.na(imc_adjdepr_padj)) %>% pull(variable)
  )

  gv <- makeVenn(vars2venn, "VennDiagram_Cond_CondAdj_BMI_BMIAdj", opt)
  gv <- makeVenn(vars2venn[1:2], "VennDiagram_Cond_CondAdj", opt)
  gv <- makeVenn(vars2venn[c(1,3,4)], "VennDiagram_Cond_BMI_BMIAdj", opt)

  medresults <- lapply(vars2test, \(x, df_all){
    df <- df_all %>% dplyr::select(all_of(c(x, y_name, mediator_name)))
    with_nas <-  df %>% apply(MAR=1, \(x)any(is.na(x)))
    df <- df[!with_nas, ]
    param_names <- c("b", "cp", "a", "a*b", "cp+a*b")
    res <- tryCatch({makeMediationSimpleLavaan(df, x, y_name, mediator_name)}, error=\(x){
      xx <- data.frame(Estimate=rep(NA, 5),
                       S.E. = rep(NA, 5),
                       `Z-score`=rep(NA, 5),
                       p.value=rep(NA, 5))
      rownames(xx) <- c(param_names)
      return(list(estimates=xx))
    })
    res2 <- res$estimates %>% rownames_to_column("param") %>%
      dplyr::select(param, Estimate, p.value) %>%
      gather(key="var", value="value", Estimate, p.value) %>%
      filter(param %in% param_names) %>%
      tidyr::unite("tmp", param, var, sep="_") %>%
      spread(tmp, value) %>%
      mutate(Xvar = x, Yvar = y_name, Mediator=mediator_name, fullmodel=list(res$estimates))
    return(res2)
  }, df_all) %>% bind_rows()
  #all(medresults$Xvar %in% summary_df$variable)
  medresults_merged <- merge(medresults %>% dplyr::select(-fullmodel), summary_df, by.x="Xvar", by.y="variable")
  write_tsv(medresults_merged, file=paste0(opt$out, "mediation_analysis_IMC.tsv"))
  save(medresults_merged, file=paste0(opt$out, "mediation_analysis_IMC.RData"))

  ## Make power analysis
  if(make_power_test){
    cat("--> Entering Power Analysis\n")
    medPowertest <- map2(medresults$Xvar, medresults$fullmodel, \(x, medres_x){
      #print(class(medres_x))
      df <- df_all %>% dplyr::select(all_of(c(x, y_name, mediator_name)))
      with_nas <-  df %>% apply(MAR=1, \(x)any(is.na(x)))
      df <- df[!with_nas, ]
      #medres_x <- medresults_merged %>% filter(Xvar == x)
      #param_ests <- c(a = medres_x$a_Estimate, b=medres_x$b_Estimate, cp=medres_x$cp_Estimate)
      cat("--> Calculating power for ", x, "\n")
      res <- tryCatch({makeMediationSimplePowerCurveLavaan(df, x, y_name,
                                                     mediator_name,
                                                     med_res =medres_x,
                                                     nrep=100,
                                                     min_n=100, max_n=1000, inc_n=50, error=0.1) %>%
          dplyr::mutate(taxa = x) %>%
          dplyr::select(taxa, everything())

        }, error=\(x){
        warning(x)
        xx <- data.frame(taxa = character(0),
                         param = character(0),
                         ssize=numeric(0),
                         mean_est = character(0),
                         mean_se=numeric(0),
                         mean_zscore=character(0),
                         n = numeric(0),
                         power=numeric(0))
        return(list(pwres = NULL, pwsum=NULL, pwcurve = xx))
      })
      return(res)
    }) %>% bind_rows()
    write_tsv(medPowertest, file = paste0(opt$out, "mediation_analysis_IMC_PowerAnalysis.tsv"))

  }else{
    medPowertest <- list()
  }

  ## Transform to plot
  def_effects <- c('a', 'b', 'cp','a*b', 'cp+a*b')
  res_trans <- map(1:nrow(medresults), \(i){
    res_i <- medresults$fullmodel[[i]]
    oldnames2 <- rownames(res_i) %>%
      sapply(\(x)strsplit(x, "\\*|\\+",perl = TRUE)[[1]][1])
    resdf_i <- res_i %>%
      rownames_to_column("param") %>%
      mutate(x_labels = ifelse(param %in% def_effects,
                               medresults$Xvar[i],
                               gsub("V\\[|\\]", "", param)),
             param = ifelse(param %in% def_effects,
                            gsub("(b|a|cp)", paste0("\\1", as.character(i)), param, perl=T), param)
      )
  }) %>% bind_rows()
  bs_only <- res_trans %>% dplyr::filter(grepl("^b[0-9]+$", param))
  g_bs <- ggplot(bs_only,aes(x=Estimate))+geom_histogram()+mytheme+ggtitle("Histogram of b (IMC->Depr)")
  ggsave(filename = paste0(opt$out, "histogram_bs_IMC_separate.pdf"), g_bs)

  res_final <- rbind(
    bs_only %>%
      mutate_if(is.character, \(x)"b") %>%
      group_by(param, x_labels) %>%
      dplyr::summarise_all(mean) %>%
      dplyr::select(all_of(names(res_trans))),
    res_trans %>% dplyr::filter(!grepl("^b[0-9]+$", param)) %>%
      dplyr::mutate(param=gsub("b[0-9]+$", "b", param))
  )
  if(AJUST_PVALS){
    res_final$p_raw <- res_final$p.value
    res_final$p.value <- p.adjust(res_final$p.value, method = "BH")
  }
  write_tsv(res_final, file=paste0(opt$out, "mediation_analysis_IMC_separate_reformat.tsv"))
  write_tsv(bs_only, file=paste0(opt$out, "mediation_analysis_IMC_separate_bs.tsv"))

  plotres <- plotIMCMediationSimple(res_final, vars2test, opt$out, name, plim_plot = plim_plot,
                                    use_color_scale = use_color_scale, w=wnet, h=hnet, custom_colors=custom_colors)

  ## Make also boxplot
  bacnames <- plotres$bacorder$x_labels
  barplot_list <- list()
  if(make_boxplots){
    barplot_list$plots1 <- makePlotBySpecies(bacnames, df_all, opt$out,
                                             paste0("IMC_separate_BySpeciesPearson", as.character(plim_plot)),
                                             quantvar="IMC_log",
                                             quantvar_name = "log(IMC)",
                                             corrmethod = "pearson", w=wbars, h=hbars, plim_plot = plim_plot)
    barplot_list$plots2 <- makePlotBySpecies(bacnames,
                                             df_all,
                                             opt$out,
                                             paste0("IMC_separate_BySpeciesSpearman", as.character(plim_plot)),
                                             quantvar="IMC_log",
                                             quantvar_name = "log(IMC)", corrmethod = "spearman", w=wbars, h=hbars, plim_plot = plim_plot)
    barplot_list$plots3 <- makePlotBySpecies(bacnames, df_all, opt$out, paste0("IMC_separate_BySpeciesKendall", as.character(plim_plot)), quantvar="IMC_log",
                              quantvar_name = "log(IMC)", corrmethod = "kendall", w=wbars, h=hbars, plim_plot = plim_plot)

    barplot_list$plots4 <- makePlotBySpeciesEffects_BMI(bacnames, df_all, res_final,
                                                        opt$out,
                                                        paste0("IMC_separate_BySpeciesEffects_", as.character(plim_plot)),
                                         w=wbars2, h=hbars2, plim_plot=plim_plot, fix_limits = fix_barplot_limits)
    barplot_list$plots5 <- makePlotBySpeciesEffects_BMI2(bacnames, df_all, res_final,
                                                        opt$out,
                                                        paste0("IMC_separate_BySpeciesEffects2_", as.character(plim_plot)),
                                                        w=wbars2, h=hbars2, plim_plot=plim_plot, fix_limits = fix_barplot_limits)
  }
  return(list(plotNet=plotres, barplots = barplot_list, res_final=res_final, power_analysis=medPowertest))
}



makeFullMediationAnalysisIMC_Aggregate <- function(opt, getVarsFunctions, mediator_name="IMC", y_name="Condition_bin",
                                         plim=0.05, plim_plot=0.05, name="analysis_IMC_separateModel_vjust",
                                         wnet=14, hnet=12, wbars=8, hbars=10, wbars2=10, hbars2=12, use_color_scale=FALSE,
                                         fix_barplot_limits=FALSE, custom_colors=NULL){
  cat(opt$out)
  summary_df <- data.frame(variable = gsub("-", ".", vstdf$gene) %>%
                             gsub("\\[|\\]", "", .) %>%
                             gsub("\\._", "_", .) %>%
                             gsub("\\/|\\(|\\)", ".", .) %>%
                             gsub("\\.$", "", .)
  )

  assertthat::assert_that(all(summary_df$variable %in% names(df_all)))

  list2merge <- list(
    depr_only_padj = dea2contrasts$firstContrast$resdf,
    imc_only_padj = dea2contrasts$contrastlist2$BMI_alone$resdf,
    depr_adjimc_padj = dea2contrasts$contrastlist2$Condition_corrIMC$resdf,
    imc_adjdepr_padj = dea2contrasts$contrastlist2$BMI_corrCond$resdf
  )

  merged_pvals <- list2merge %>%
    lapply(\(x){
      x$taxon <- x$taxon %>%
        gsub("-", ".", .) %>%
        gsub("\\[|\\]", "", .) %>%
        gsub("\\._", "_", .) %>%
        gsub("\\/|\\(|\\)", ".", .) %>%
        gsub("\\.$", "", .)
      x <- x[match(summary_df$variable, x$taxon), ]
      #x$padj[is.na(x$padj)] <- 1
      return(x$padj)
    }) %>% bind_cols()
  summary_df <- cbind(summary_df, merged_pvals)

  vars2test <- getVarsFunction(summary_df)

  medresults <- lapply(vars2test, \(x, df_all){
    df <- df_all %>% dplyr::select(all_of(c(x, y_name, mediator_name)))
    with_nas <-  df %>% apply(MAR=1, \(x)any(is.na(x)))
    df <- df[!with_nas, ]
    param_names <- c("b", "cp", "a", "a*b", "cp+a*b")
    res <- tryCatch({makeMediationSimple(df, x, y_name, mediator_name)}, error=\(x){
      xx <- data.frame(Estimate=rep(NA, 5),
                       S.E. = rep(NA, 5),
                       `Z-score`=rep(NA, 5),
                       p.value=rep(NA, 5))
      rownames(xx) <- c(param_names)
      return(list(estimates=xx))
    })
    res2 <- res$estimates %>% rownames_to_column("param") %>%
      dplyr::select(param, Estimate, p.value) %>%
      gather(key="var", value="value", Estimate, p.value) %>%
      filter(param %in% param_names) %>%
      tidyr::unite("tmp", param, var, sep="_") %>%
      spread(tmp, value) %>%
      mutate(Xvar = x, Yvar = y_name, Mediator=mediator_name, fullmodel=list(res$estimates))
    return(res2)
  }, df_all) %>% bind_rows()
  #all(medresults$Xvar %in% summary_df$variable)
  medresults_merged <- merge(medresults %>% dplyr::select(-fullmodel), summary_df, by.x="Xvar", by.y="variable")
  write_tsv(medresults_merged, file=paste0(opt$out, "mediation_analysis_IMC.tsv"))
  save(medresults_merged, file=paste0(opt$out, "mediation_analysis_IMC.RData"))

  ## Transform to plot
  def_effects <- c('a', 'b', 'cp','a*b', 'cp+a*b')
  res_trans <- map(1:nrow(medresults), \(i){
    res_i <- medresults$fullmodel[[i]]
    oldnames2 <- rownames(res_i) %>%
      sapply(\(x)strsplit(x, "\\*|\\+",perl = TRUE)[[1]][1])
    resdf_i <- res_i %>%
      rownames_to_column("param") %>%
      mutate(x_labels = ifelse(param %in% def_effects,
                               medresults$Xvar[i],
                               gsub("V\\[|\\]", "", param)),
             param = ifelse(param %in% def_effects,
                            gsub("(b|a|cp)", paste0("\\1", as.character(i)), param, perl=T), param)
      )
  }) %>% bind_rows()
  bs_only <- res_trans %>% dplyr::filter(grepl("^b[0-9]+$", param))
  g_bs <- ggplot(bs_only,aes(x=Estimate))+geom_histogram()+mytheme+ggtitle("Histogram of b (IMC->Depr)")
  ggsave(filename = paste0(opt$out, "histogram_bs_IMC_separate.pdf"), g_bs)

  res_final <- rbind(
    bs_only %>%
      mutate_if(is.character, \(x)"b") %>%
      group_by(param, x_labels) %>%
      dplyr::summarise_all(mean) %>%
      dplyr::select(all_of(names(res_trans))),
    res_trans %>% dplyr::filter(!grepl("^b[0-9]+$", param)) %>%
      dplyr::mutate(param=gsub("b[0-9]+$", "b", param))
  )
  if(AJUST_PVALS){
    res_final$p_raw <- res_final$p.value
    res_final$p.value <- p.adjust(res_final$p.value, method = "BH")
  }
  write_tsv(res_final, file=paste0(opt$out, "mediation_analysis_IMC_separate_reformat.tsv"))
  write_tsv(bs_only, file=paste0(opt$out, "mediation_analysis_IMC_separate_bs.tsv"))

  plotres <- plotIMCMediationSimple(res_final, vars2test, opt$out, name, plim_plot = plim_plot,
                                    use_color_scale = use_color_scale, w=wnet, h=hnet, custom_colors=custom_colors)

  ## Make also boxplot
  bacnames <- plotres$bacorder$x_labels
  plots1 <- makePlotBySpecies(bacnames, df_all, opt$out, paste0("IMC_separate_BySpeciesPearson", as.character(plim_plot)), quantvar="IMC_log",
                              quantvar_name = "log(IMC)", corrmethod = "pearson", w=wbars, h=hbars, plim_plot = plim_plot)
  plots2 <- makePlotBySpecies(bacnames, df_all, opt$out, paste0("IMC_separate_BySpeciesSpearman", as.character(plim_plot)), quantvar="IMC_log",
                              quantvar_name = "log(IMC)", corrmethod = "spearman", w=wbars, h=hbars, plim_plot = plim_plot)
  plots3 <- makePlotBySpecies(bacnames, df_all, opt$out, paste0("IMC_separate_BySpeciesKendall", as.character(plim_plot)), quantvar="IMC_log",
                              quantvar_name = "log(IMC)", corrmethod = "kendall", w=wbars, h=hbars, plim_plot = plim_plot)

  plots4 <- makePlotBySpeciesEffects_BMI(bacnames, df_all, res_final, opt$out, paste0("IMC_separate_BySpeciesEffects", as.character(plim_plot)),
                                         w=wbars2, h=hbars2, plim_plot=plim_plot, fix_limits = fix_barplot_limits)

  return(list(plotNet=plotres, barplots = list(plots1, plots2, plots3, plots4), res_final=res_final))
}


makeFullMediationAnalysisOb <- function(opt, getVarsFunction, mediator_name="IMC", y_name="Condition_bin",
                                         plim=0.05, plim_plot=0.05, name="analysis_IMC_separateModel_vjust",
                                         wnet=14, hnet=12, wbars=8, hbars=10, wbars2=10, hbars2=12, use_color_scale=FALSE,
                                         fix_barplot_limits=FALSE, custom_colors=NULL, make_boxplots=TRUE){
  print(paste0("OUTPUT: ", opt$out))
  summary_df <- data.frame(variable = gsub("-", ".", vstdf$gene) %>%
                             gsub("\\[|\\]", "", .) %>%
                             gsub("\\._", "_", .) %>%
                             gsub("\\/|\\(|\\)", ".", .) %>%
                             gsub("\\.$", "", .)
  )

  assertthat::assert_that(all(summary_df$variable %in% names(df_all)))

  list2merge <- list(
    depr_only_padj = dea2contrasts$firstContrast$resdf,
    imc_only_padj = dea2contrasts$contrastlist2$Overweight_alone$resdf,
    depr_adjimc_padj = dea2contrasts$contrastlist2$Condition_corrOverweight$resdf,
    imc_adjdepr_padj = dea2contrasts$contrastlist2$Overweight_corrCond$resdf
  )

  merged_pvals <- list2merge %>%
    lapply(\(x){
      x$taxon <- x$taxon %>%
        gsub("-", ".", .) %>%
        gsub("\\[|\\]", "", .) %>%
        gsub("\\._", "_", .) %>%
        gsub("\\/|\\(|\\)", ".", .) %>%
        gsub("\\.$", "", .)
      x <- x[match(summary_df$variable, x$taxon), ]
      #x$padj[is.na(x$padj)] <- 1
      return(x$padj)
    }) %>% bind_cols()
  summary_df <- cbind(summary_df, merged_pvals)

  vars2test <- getVarsFunction(summary_df)
  vars2venn <- list(
    "D vs C" = summary_df %>% dplyr::filter(depr_only_padj < plim & !is.na(depr_only_padj)) %>% pull(variable),
    "D vs C adj. Ob" = summary_df %>% dplyr::filter(depr_adjimc_padj < plim & ! is.na(depr_adjimc_padj)) %>% pull(variable),
    "Ob" = summary_df %>% dplyr::filter(imc_only_padj < plim & ! is.na(imc_only_padj)) %>% pull(variable),
    "Ob adj. Depr" = summary_df %>% dplyr::filter(imc_adjdepr_padj < plim & !is.na(imc_adjdepr_padj)) %>% pull(variable)
  )

  gv <- makeVenn(vars2venn, "VennDiagram_Cond_CondAdj_Ob_ObAdj", opt)
  gv <- makeVenn(vars2venn[1:2], "VennDiagram_Cond_CondAdj", opt)
  gv <- makeVenn(vars2venn[c(1,3,4)], "VennDiagram_Cond_Ob_ObAdj", opt)

  medresults <- lapply(vars2test, \(x, df_all){
    df <- df_all %>% dplyr::select(all_of(c(x, y_name, mediator_name)))
    with_nas <-  df %>% apply(MAR=1, \(x)any(is.na(x)))
    df <- df[!with_nas, ]
    param_names <- c("b", "cp", "a", "a*b", "cp+a*b")
    res <- tryCatch({makeMediationSimple(df, x, y_name, mediator_name)}, error=\(x){
      xx <- data.frame(Estimate=rep(NA, 5),
                       S.E. = rep(NA, 5),
                       `Z-score`=rep(NA, 5),
                       p.value=rep(NA, 5))
      rownames(xx) <- c(param_names)
      return(list(estimates=xx))
    })
    res2 <- res$estimates %>% rownames_to_column("param") %>%
      dplyr::select(param, Estimate, p.value) %>%
      gather(key="var", value="value", Estimate, p.value) %>%
      filter(param %in% param_names) %>%
      tidyr::unite("tmp", param, var, sep="_") %>%
      spread(tmp, value) %>%
      mutate(Xvar = x, Yvar = y_name, Mediator=mediator_name, fullmodel=list(res$estimates))
    return(res2)
  }, df_all) %>% bind_rows()
  #all(medresults$Xvar %in% summary_df$variable)
  medresults_merged <- merge(medresults %>% dplyr::select(-fullmodel), summary_df, by.x="Xvar", by.y="variable")
  write_tsv(medresults_merged, file=paste0(opt$out, "mediation_analysis_IMC.tsv"))
  save(medresults_merged, file=paste0(opt$out, "mediation_analysis_IMC.RData"))

  ## Transform to plot
  def_effects <- c('a', 'b', 'cp','a*b', 'cp+a*b')
  res_trans <- map(1:nrow(medresults), \(i){
    res_i <- medresults$fullmodel[[i]]
    oldnames2 <- rownames(res_i) %>%
      sapply(\(x)strsplit(x, "\\*|\\+",perl = TRUE)[[1]][1])
    resdf_i <- res_i %>%
      rownames_to_column("param") %>%
      mutate(x_labels = ifelse(param %in% def_effects,
                               medresults$Xvar[i],
                               gsub("V\\[|\\]", "", param)),
             param = ifelse(param %in% def_effects,
                            gsub("(b|a|cp)", paste0("\\1", as.character(i)), param, perl=T), param)
      )
  }) %>% bind_rows()
  bs_only <- res_trans %>% dplyr::filter(grepl("^b[0-9]+$", param))
  g_bs <- ggplot(bs_only,aes(x=Estimate))+geom_histogram()+mytheme+ggtitle("Histogram of b (IMC->Depr)")
  ggsave(filename = paste0(opt$out, "histogram_bs_IMC_separate.pdf"), g_bs)

  res_final <- rbind(
    bs_only %>%
      mutate_if(is.character, \(x)"b") %>%
      group_by(param, x_labels) %>%
      dplyr::summarise_all(mean) %>%
      dplyr::select(all_of(names(res_trans))),
    res_trans %>% dplyr::filter(!grepl("^b[0-9]+$", param)) %>%
      dplyr::mutate(param=gsub("b[0-9]+$", "b", param))
  )
  if(AJUST_PVALS){
    res_final$p_raw <- res_final$p.value
    res_final$p.value <- p.adjust(res_final$p.value, method = "BH")
  }
  write_tsv(res_final, file=paste0(opt$out, "mediation_analysis_Ob_separate_reformat.tsv"))
  write_tsv(bs_only, file=paste0(opt$out, "mediation_analysis_Ob_separate_bs.tsv"))

  plotres <- plotIMCMediationSimple(res_final, vars2test, opt$out, name, plim_plot = plim_plot,
                                    use_color_scale = use_color_scale, w=wnet, h=hnet, custom_colors=custom_colors)

  ## Make also boxplot
  bacnames <- plotres$bacorder$x_labels
  barplot_list <- list()
  if(make_boxplots){
    barplot_list$plots1 <- makePlotBySpecies(bacnames, df_all, opt$out, paste0("Ob_separate_BySpeciesPearson", as.character(plim_plot)),
                                             quantvar="IMC_log",
                                             quantvar_name = "log(IMC)", corrmethod = "pearson", w=wbars, h=hbars, plim_plot = plim_plot)
    barplot_list$plots2 <- makePlotBySpecies(bacnames, df_all, opt$out, paste0("Ob_separate_BySpeciesSpearman", as.character(plim_plot)),
                                             quantvar="IMC_log",
                                             quantvar_name = "log(IMC)", corrmethod = "spearman", w=wbars, h=hbars, plim_plot = plim_plot)
    barplot_list$plots3 <- makePlotBySpecies(bacnames, df_all, opt$out, paste0("Ob_separate_BySpeciesKendall", as.character(plim_plot)),
                                             quantvar="IMC_log",
                                             quantvar_name = "log(IMC)", corrmethod = "kendall", w=wbars, h=hbars, plim_plot = plim_plot)

    barplot_list$plots4 <- makePlotBySpeciesEffects_BMI(bacnames, df_all,
                                                        res_final, opt$out,
                                                        paste0("IMC_separate_BySpeciesEffects", as.character(plim_plot)),
                                                        w=wbars2, h=hbars2, plim_plot=plim_plot, fix_limits = fix_barplot_limits)
  }
  return(list(plotNet=plotres, barplots = barplot_list, res_final=res_final))
}
