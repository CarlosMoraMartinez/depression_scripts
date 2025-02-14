library(tidyverse)
library(ggpmisc)
library(mgcv)

path_phyloseq <- paste0(opt$out, "/phyloseq")
allphyloseqlist_fname <- paste0(path_phyloseq, "/phyloseq_all_list.RData")
load(allphyloseqlist_fname)

load( paste0(opt$out, "DeSEQ2/DESEQ2_all_edad00.RData"))

outdir <- paste0(opt$out, "fit_models/")
if(!dir.exists(outdir)) dir.create(outdir)

tempdf <- daa_all$remove_tanda2$vstds %>% t %>% 
  as.data.frame %>% 
  rownames_to_column("sampleID") %>% 
  merge(sample_data(all_phyloseq$remove_tanda2) %>% data.frame, by.x="sampleID", by.y="sampleID", all=T)

tempdf_norm <- daa_all$remove_tanda2$norm_counts %>% t %>% 
  as.data.frame %>% 
  rownames_to_column("sampleID") %>% 
  merge(sample_data(all_phyloseq$remove_tanda2) %>% data.frame, by.x="sampleID", by.y="sampleID", all=T)

tempdf_raw <- daa_all$remove_tanda2$raw_counts %>% t %>% 
  as.data.frame %>% 
  rownames_to_column("sampleID") %>% 
  merge(sample_data(all_phyloseq$remove_tanda2) %>% data.frame, by.x="sampleID", by.y="sampleID", all=T)



bacnames <- daa_all$remove_tanda2$vstds %>% rownames
table(bacnames %in% names(tempdf))

glist <- map(bacnames, \(bb){
  ggplot(tempdf, aes(x=edad_00, y=!!sym(bb)))+ 
    geom_point() +
    stat_poly_eq(use_label(c("eq", "R2", "p")), #c("eq", "R2", "f", "p", "n")
                 method="lm", small.p=T, small.r=F, label.y=0.99) +
    geom_smooth(alpha=0.4)+
    theme_bw()
  
})

pdf(paste0(opt$out, "DeSEQ2/edad_vs_bact.pdf"), width=8, height=6)
for(gg in glist){
  print(gg)
}
dev.off()

## With norm data
bacnames_norm <- daa_all$remove_tanda2$norm_counts %>% rownames
table(bacnames_norm %in% names(tempdf_norm))

glist <- map(bacnames_norm, \(bb){
  ggplot(tempdf_norm, aes(x=edad_00, y=!!sym(bb)))+ 
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

library(mgcv)

fixedvars <- c("nreads_filt", "hospital", "Sex", "edu_m_00")
svars <- c("edad_00", "nreads_filt")

extformula <- paste( map(svars, \(x) paste0("s(", x, ")")) %>% unlist %>% paste(sep=" + ", collapse=" + " ),
                     fixedvars )

gamlist <- map(bacnames, \(bb){
  form <- as.formula(paste0(bb, 
                            "~ s(edad_00)",
                            "s(nreads_filt)"))
  gam1<-mgcv::gam(list(form, form), data = tempdf_raw, family=ziplss())
  
  form <- as.formula(paste0(bb, "~ s(edad_00)"))
  gam1<-mgcv::gam(list(Cloacibacillus_porcorum ~ s(edad_00), ~ s(edad_00)), data = tempdf_raw, family=ziplss())
  
  gam1<-mgcv::gam(Cloacibacillus_porcorum ~ s(edad_00) + s(bmi_t1) + edu_m_00 + log(nreads_filt) + hospital + Sex, 
                  data = tempdf_raw, family=nb())
  
  gam1<-mgcv::gam(Cloacibacillus_porcorum ~ s(edad_00, bs = "cs")  + nreads_filt + hospital + Sex + edu_m_00, 
                  data = tempdf_raw%>% filter(Cloacibacillus_porcorum < 560), family=nb()) #Komagataeibacter_oboediens
  
  
  summary(gam1)
  par(mfrow=c(3, 2))
  plot(gam1)
  gam.check(gam1)
  
  summary(gam2)
  par(mfrow=c(3, 2))
  plot(gam2)
  gam.check(gam2)
  
})


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