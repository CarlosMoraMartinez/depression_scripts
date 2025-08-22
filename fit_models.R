library(tidyverse)
library(ggpmisc)
library(mgcv)

editnames <- function(nn){
  gsub("-", ".", nn) %>% 
    gsub("[\\[\\]\\/\\(\\)]", "", ., perl=TRUE)
}


gam2tab <- function(mm, depvarname){
  mmsum <- summary(mm)
  part1 <- mmsum$s.table %>% as.data.frame %>% rownames_to_column("s_var") %>% 
    gather(key = "param", value="value", -s_var) %>% 
    mutate(s_var = gsub("s\\(|\\)", "", s_var, perl=TRUE)) %>% 
    mutate(s_var = paste("s_", s_var, sep="")) %>% 
    mutate(param = gsub("\\.", "", param)) %>% 
    unite("var_param", s_var, param, sep=".", remove=T) %>% 
    spread(key=var_param, value = value)
  
  part2 <-list(mmsum$p.coeff, mmsum$p.pv, mmsum$p.t) %>% bind_rows() %>% 
    mutate(Estimate = c("estimate", "pval", "zvalue") ) %>% 
    gather("variable", "value", -Estimate) %>% 
    mutate(variable = gsub("[\\(\\)]", "", variable, perl=T)) %>% 
    unite("var_estimate",  Estimate, variable, sep="_") %>% 
    spread(var_estimate, value)
  
  both <- bind_cols(part1, part2) %>% 
    mutate(DependentVariable = depvarname) %>% 
    select(DependentVariable, everything()) %>% 
    mutate(deviance = mm$deviance, 
           min_edf = mm$min.edf
           )
  return(both)
}

path_phyloseq <- paste0(opt$out, "/phyloseq")
allphyloseqlist_fname <- paste0(path_phyloseq, "/phyloseq_all_list.RData")
load(allphyloseqlist_fname)

load( paste0(opt$out, "DeSEQ2/DESEQ2_all_edad00.RData"))

functdir <-  paste0(opt$out, "Functional_age/")
functdf <- read_tsv(paste0(functdir, "functTabInput_MetaCyc_filteredByProcess.tsv"))
functdf_long <- functdf %>% 
  gather("sampleID", "Abundance", -Pathway)
functdf_t <- functdf %>% 
  column_to_rownames("Pathway") %>% 
  as.matrix %>% t %>% 
  as.data.frame %>% 
  rownames_to_column("sampleID")%>% 
  merge(sample_data(all_phyloseq$remove_tanda2) %>% data.frame, by.x="sampleID", by.y="sampleID", all=T)

outdir <- paste0(opt$out, "fit_models/")
if(!dir.exists(outdir)) dir.create(outdir)

phname <- "remove_tanda2"
tempdf <- daa_all$remove_tanda2$vstds %>% t %>% 
  as.data.frame %>% 
  rownames_to_column("sampleID") %>% 
  merge(sample_data(all_phyloseq$remove_tanda2) %>% data.frame, by.x="sampleID", by.y="sampleID", all=T)
names(tempdf) <- editnames(names(tempdf))

tempdf_norm <- daa_all$remove_tanda2$norm_counts %>% t %>% 
  as.data.frame %>% 
  rownames_to_column("sampleID") %>% 
  merge(sample_data(all_phyloseq$remove_tanda2) %>% data.frame, by.x="sampleID", by.y="sampleID", all=T)
names(tempdf_norm) <- editnames(names(tempdf_norm))

tempdf_raw <- daa_all$remove_tanda2$raw_counts %>% t %>% 
  as.data.frame %>% 
  rownames_to_column("sampleID") %>% 
  merge(sample_data(all_phyloseq$remove_tanda2) %>% data.frame, by.x="sampleID", by.y="sampleID", all=T)
names(tempdf_raw) <- editnames( names(tempdf_raw))



bacnames <- daa_all$remove_tanda2$norm_counts %>% rownames %>% 
  editnames() ##VST FAILED
table(bacnames %in% names(tempdf))

glist <- map(bacnames, \(bb){
  ggplot(tempdf, aes(x=edad_00, y=!!sym(bb)))+ 
    geom_point() +
    stat_poly_eq(use_label(c("eq", "R2", "p")), #c("eq", "R2", "f", "p", "n")
                 method="lm", small.p=T, small.r=F, label.y=0.99) +
    geom_smooth(alpha=0.4)+
    theme_bw()
  
})

pdf(paste0(outdir,  "/edad_vs_bact.pdf"), width=8, height=6)
for(gg in glist){
  print(gg)
}
dev.off()






# Better numeric
#tempdf_raw <- tempdf_raw %>% mutate(
#  edu_m_00 = as.factor(edu_m_00)
#)
all(bacnames %in% names(tempdf_raw ))
fixedvars <- c( "hospital", "Sex", "edu_m_00")
svars <- c("edad_00", "nreads_filt", "z_imc_00")
extformula <- paste( map(svars, \(x) paste0("s(", x, ")")) %>% unlist %>% paste(sep=" + ", collapse=" + " ),
                     " + ",
                     paste(fixedvars,  sep=" + ", collapse=" + " )
                     )

FILTER_QUANTILE <- 0.95
gamlist <- map(bacnames, \(bb){
  cat(bb, "\n")
  form <- as.formula(paste0(bb, " ~ ", extformula))
  q80 <- quantile(tempdf_raw[, bb] %>% unlist, 0.8)
  if(FILTER_QUANTILE < 1 & q80 > 0){
    quant <- quantile(tempdf_raw[, bb] %>% unlist, FILTER_QUANTILE) # tempdf_raw[, bb] > 0
    auxdf <- tempdf_raw %>% filter(!!sym(bb) < quant)
  }else{
    auxdf <- tempdf_raw
  }
  gam1<-mgcv::gam(form,  data = auxdf, family=nb()) 
  return(gam1)
})

names(gamlist) <- bacnames
save(gamlist, file = paste0(outdir, "GAMs_NB_filtq95.RData"))
load( paste0(outdir, "GAMs_NB_filtq95.RData"))

gamtab <- map2(gamlist, bacnames, gam2tab) %>% bind_rows
write_tsv(gamtab, file = paste0(outdir, "GAMs_NB_filtq95.tsv"))
gamtab <- read_tsv(paste0(outdir, "GAMs_NB_filtq95.tsv"))

daage <- gamtab %>% filter(`s_edad_00.p-value` < 0.05) 
deseqage <- daa_all$remove_tanda2$all_contrasts$edad_00$resdf %>% filter(padj < 0.05)

pdf(paste0(outdir, "GAMs_NB_plots_filtq95.pdf"), height = 12, width=8)
for(nn in names(gamlist)){
  
  par(mfrow=c(4,2))
  
  aux <- gamtab %>% filter(DependentVariable == nn)
  plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
  text(x = 0.5, y = 0.5, paste(nn, as.character(round(aux$`s_edad_00.p-value`, 4))), 
       cex = 1.6, col = ifelse(aux$`s_edad_00.p-value` < 0.05, "red", "blue"))
  
  
  plot(gamlist[[nn]], main=nn)
  gam.check(gamlist[[nn]])

}
dev.off()

#gam1<-mgcv::gam(list(form, form), data = tempdf_raw, family=ziplss())

## Food

foodnames <- names(tempdf_raw)[c(493:518, 524)]
fhists <- map(foodnames, \(x){ggplot(tempdf_raw, aes(x= !!sym(x))) + geom_histogram()+theme_bw()+ggtitle(x) })
pdf(paste0(outdir, "food_histograms.pdf", width=8, height=6))
for(g in fhists) print(g)
dev.off()

foodageg <- map(foodnames, \(x){
  ggplot(tempdf_raw, aes(x=edad_00, y= !!sym(x))) + 
    geom_point()+
    stat_poly_eq(use_label(c("eq", "R2", "p")), #c("eq", "R2", "f", "p", "n")
                 method="lm", small.p=T, small.r=F, label.y=0.99) +
    geom_smooth() +
    theme_bw() +
    ggtitle(x) 
  })

pdf(paste0(outdir, "food_vs_age.pdf", width=8, height=6))
for(g in foodageg) print(g)
dev.off()

foodageg <- map(foodnames, \(x){
  ggplot(tempdf_raw, aes(x=edad_00, y= log(!!sym(x)))) + 
    geom_point()+
    stat_poly_eq(use_label(c("eq", "R2", "p")), #c("eq", "R2", "f", "p", "n")
                 method="lm", small.p=T, small.r=F, label.y=0.99) +
    geom_smooth() +
    theme_bw() +
    ggtitle(x) 
})

pdf(paste0(outdir, "food_vs_age_log.pdf", width=8, height=6))
for(g in foodageg) print(g)
dev.off()




## food with ratios

foodnames <- names(tempdf_raw)[c(494:518)]
kcal <- names(tempdf_raw)[493]
tempdf_raw_foodNorm <- tempdf_raw %>% dplyr::mutate(across(all_of(foodnames), ~ 100*.x / .data[[kcal]]))

foodageg <- map(foodnames, \(x){
  ggplot(tempdf_raw_foodNorm, aes(x=edad_00, y= !!sym(x))) + 
    geom_point()+
    stat_poly_eq(use_label(c("eq", "R2", "p")), #c("eq", "R2", "f", "p", "n")
                 method="lm", small.p=T, small.r=F, label.y=0.99) +
    geom_smooth() +
    theme_bw() +
    ggtitle(x) 
})

pdf(paste0(outdir, "food_vs_age_KcalRatio.pdf", width=8, height=6))
for(g in foodageg) print(g)
dev.off()

## Model food 

FILTER_QUANTILE <- 1
foodnames <- names(tempdf_raw)[c(494:518)]

fixedvars <- c( "hospital", "Sex", "edu_m_00")
svars <- c("edad_00", "ffq_energia_00", "z_imc_00")
extformula <- paste( map(svars, \(x) paste0("s(", x, ")")) %>% unlist %>% paste(sep=" + ", collapse=" + " ),
                     " + ",
                     paste(fixedvars,  sep=" + ", collapse=" + " )
)

all(fixedvars %in% names(tempdf_raw))
all(svars %in% names(tempdf_raw))

gamlist <- map(foodnames, \(bb){
  cat(bb, "\n")
  form <- as.formula(paste0(bb, " ~ ", extformula))
  q80 <- quantile(tempdf_raw[!is.na(tempdf_raw[, bb]), bb] %>% unlist, 0.8)
  if(FILTER_QUANTILE < 1 & q80 > 0){
    quant <- quantile(tempdf_raw[!is.na(tempdf_raw[, bb]), bb] %>% unlist, FILTER_QUANTILE) # functdf_t[, bb] > 0
    auxdf <- tempdf_raw %>% filter(! is.na(!!sym(bb))) %>% filter(!!sym(bb) < quant)
  }else{
    auxdf <- tempdf_raw %>% filter(! is.na(!!sym(bb))) 
  }
  gam1<-mgcv::gam(form,  data = auxdf, family=gaussian()) 
  return(gam1)
})

names(gamlist) <- foodnames
save(gamlist, file = paste0(outdir, "GAMs_food_gaussian.RData"))

gamtab <- map2(gamlist, names(gamlist), gam2tab) %>% bind_rows
write_tsv(gamtab, file = paste0(outdir, "GAMs_food_gaussian.tsv"))
#gamtab <- read_tsv(paste0(outdir, "GAMs_NB_filtq95.tsv"))


pdf(paste0(outdir, "GAMs_food_gaussian_plots.pdf"), height = 12, width=8)
for(nn in names(gamlist)){
  
  par(mfrow=c(4,2))
  
  aux <- gamtab %>% filter(DependentVariable == nn)
  plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
  text(x = 0.5, y = 0.5, paste(gsub("_", " ", nn), '\n', as.character(round(aux$`s_edad_00.p-value`, 4))), 
       cex = 1, col = ifelse(aux$`s_edad_00.p-value` < 0.05, "red", "blue"))
  
  
  plot(gamlist[[nn]], main=nn)
  gam.check(gamlist[[nn]])
  
}
dev.off()








### Functional
fpaths <- functdf$Pathway
fplots <- map

fplots <- map(fpaths, \(bb){
  ggplot(functdf_t, aes(x=edad_00, y=!!sym(bb)))+ 
    geom_point() +
    stat_poly_eq(use_label(c("eq", "R2", "p")), #c("eq", "R2", "f", "p", "n")
                 method="lm", small.p=T, small.r=F, label.y=0.99) +
    geom_smooth(alpha=0.4)+
    theme_bw()
  
})
pdf(paste0(outdir, "functional.pdf", width=8, height=6))
for(g in fplots) print(g)
dev.off()

fplots <- map(fpaths, \(bb){
  ggplot(functdf_t, aes(x=edad_00, y=log10(!!sym(bb) + 1)))+ 
    geom_point() +
    stat_poly_eq(use_label(c("eq", "R2", "p")), #c("eq", "R2", "f", "p", "n")
                 method="lm", small.p=T, small.r=F, label.y=0.99) +
    geom_smooth(alpha=0.4)+
    theme_bw()
  
})
pdf(paste0(outdir, "functional_logs.pdf", width=8, height=6))
for(g in fplots) print(g)
dev.off()

FILTER_QUANTILE <- 0.95

functdf_t <- clean_names(functdf_t)
clean_fpaths <- make_clean_names(fpaths)

fixedvars <- c( "hospital", "sex", "edu_m_00")
svars <- c("edad_00", "nreads_filt", "z_imc_00")
extformula <- paste( map(svars, \(x) paste0("s(", x, ")")) %>% unlist %>% paste(sep=" + ", collapse=" + " ),
                     " + ",
                     paste(fixedvars,  sep=" + ", collapse=" + " )
)

all(svars %in% names(functdf_t))
all(fixedvars %in% names(functdf_t))

gamlist <- map(clean_fpaths, \(bb){
  cat(bb, "\n")
  form <- as.formula(paste0(bb, " ~ ", extformula))
  q80 <- quantile(functdf_t[, bb] %>% unlist, 0.8)
  if(FILTER_QUANTILE < 1 & q80 > 0){
    quant <- quantile(functdf_t[, bb] %>% unlist, FILTER_QUANTILE) # functdf_t[, bb] > 0
    auxdf <- functdf_t %>% filter(!!sym(bb) < quant)
  }else{
    auxdf <- functdf_t
  }
  gam1<-mgcv::gam(form,  data = auxdf, family=gaussian()) 
  return(gam1)
})

names(gamlist) <- clean_fpaths
save(gamlist, file = paste0(outdir, "GAMs_functional_gamma_filtQ95.RData"))

gamtab <- map2(gamlist, names(gamlist), gam2tab) %>% bind_rows
write_tsv(gamtab, file = paste0(outdir, "GAMs_functional_gamma_filtQ95.tsv"))
#gamtab <- read_tsv(paste0(outdir, "GAMs_NB_filtq95.tsv"))


pdf(paste0(outdir, "GAMs_functional_gamma_filtQ95_plots.pdf"), height = 12, width=8)
for(nn in names(gamlist)){
  
  par(mfrow=c(4,2))
  
  aux <- gamtab %>% filter(DependentVariable == nn)
  plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
  text(x = 0.5, y = 0.5, paste(gsub("_", " ", nn), '\n', as.character(round(aux$`s_edad_00.p-value`, 4))), 
       cex = 1, col = ifelse(aux$`s_edad_00.p-value` < 0.05, "red", "blue"))
  
  
  plot(gamlist[[nn]], main=nn)
  gam.check(gamlist[[nn]])
  
}
dev.off()
##############################
## With norm data
bacnames_norm <- daa_all$remove_tanda2$norm_counts %>% rownames %>% 
  editnames()
table(bacnames_norm %in% names(tempdf_norm))

glist <- map(bacnames_norm, \(bb){
  ggplot(tempdf_norm, aes(x=edad_00, y=!!sym(bb))) + 
    geom_point() +
    stat_poly_eq(use_label(c("eq", "R2", "p")), #c("eq", "R2", "f", "p", "n")
                 method="lm", small.p=T, small.r=F, label.y=0.99)
  geom_smooth(alpha=0.4)+
    theme_bw()
  
})

pdf(paste0(opt$out, "DeSEQ2/edad_vs_bact_Norm.pdf"), width=8, height=6)
for(gg in glist){
  print(gg)
}
dev.off()




####################### Other stuff
### test bimodality of Alistipes_sp._dk3624
gglist <- list()
for(v in names(tempdf)){
  if(v %in% c("sampleID")) next
  if(class(tempdf %>% pull(!!sym(v))) == "numeric"){
    g0 <- ggplot(tempdf, aes(x=!!sym(v), y = Alistipes_sp._dk3624)) +
      geom_point()+
      stat_poly_eq(use_label(c("eq", "R2", "p")), #c("eq", "R2", "f", "p", "n")
                   method="lm", small.p=T, small.r=F, label.y=0.99) +
      geom_smooth(alpha=0.4) +
      theme_bw()
      
  }else{
    g0 <- ggplot(tempdf, aes(x=!!sym(v), y = Alistipes_sp._dk3624, col=!!sym(v))) +
      geom_violin() +
      geom_boxplot(width=0.2, fill="darkgray") +
      geom_jitter(alpha=0.5) +
      theme_bw()
  }
  gglist[[v]] <- g0
  
  
}
pdf(paste0(outdir,  "/Alistipes_sp._dk3624_vs_all.pdf"), width=8, height=6)
for(gg in gglist){
  print(gg)
}
dev.off()




#~/projects/Dani_transposon/make_primers/ncbi_dataset_Bacteroides/ncbi_dataset/data$ mmseqs createdb */*.fna bacteroidesDB