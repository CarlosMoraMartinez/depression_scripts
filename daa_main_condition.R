
###########################
# Make combinations:
treatment_order <- c("REG3", "REG2", "REG1", "NO ABS", "ABS") %>% rev
region_order <- c("REG3", "REG2", "REG1") %>% rev
sd_order <- c("SD", "CONTROL") %>% rev

metad <- meta3 %>% filter(Treatment != "TRANSFER") #%>% 
  #dplyr::mutate(Treatment = factor(Treatment, levels = c("REG3", "REG2", "REG1", "NO ABS", "ABS")))

# Cominations treatment
combins_treatment <- combn(unique(metad$Treatment), 2, simplify = F) %>% 
  lapply(\(x)x[order(match(x, treatment_order))])
combins_treatment2 <- append(combins_treatment, lapply(combins_treatment, \(x) paste("SD", x, sep=" ")))

combins_treatment_all <- map(paste0("REG", 1:3), \(x) lapply(combins_treatment2, \(cc) paste(cc, x, sep=":"))) %>% flatten()

# Combinations Region
combins_region <- combn(unique(metad$Region_sequenced), 2, simplify = F) %>% 
  lapply(\(x)x[order(match(x, region_order))])
combins_region_all <- map(unique(unlist(combins_treatment2)), \(x) lapply(combins_region, \(cc) paste(x, cc, sep=":"))) %>% flatten()

# Combinations Stress
combins_stress_all <- metad %>% filter(Stress != "SD") %>% pull(Treatment_region) %>% unique %>% 
  lapply(\(x) c(paste0("SD ", x), x))

# Combinations to same level:

base_cond <- "ABS:REG1"
combins_tobaseline <- metad %>% filter(Treatment_region != base_cond) %>% 
  pull(Treatment_region) %>% unique %>% 
  lapply(\(x) c(base_cond, x))


ALL_COMBINS <- append(combins_treatment_all, combins_region_all) %>% 
  append(combins_stress_all) %>% 
  append(combins_tobaseline) %>% 
  lapply(\(x) c("Treatment_region", x))

ALL_COMBINS <- lapply(ALL_COMBINS, \(x) gsub(" |:", ".", x, perl=TRUE))
###########################

daa_all <- list()
vars2deseq <- c("Treatment_region")
vars2heatmap <- c("Treatment", "Region_sequenced", "Stress", "flowcell")
opt$mincount <- 10
opt$minsampleswithcount <- 3
phseq_to_use <- names(all_phyloseq)[1:2]


for(phname in phseq_to_use){
  # cat("Doing DESeq2 Analysys for: ", phname, "\n")
  phobj <- all_phyloseq[[phname]] 
  sdata <- sample_data(phobj) %>% data.frame
  
  sample_data(phobj)$is_transfer <- ifelse(sdata$Treatment == "TRANSFER", TRUE, FALSE)
  sample_data(phobj)$Treatment <- ifelse(sdata$Treatment == "TRANSFER", "NO ABS", sdata$Treatment)
  sample_data(phobj)$Treatment_region <- gsub("TRANSFER", "NO ABS", sdata$Treatment_region)
  sample_data(phobj)$Treatment_region <- gsub(" |:", ".", sample_data(phobj)$Treatment_region, perl=TRUE)
  
  daa_all[[phname]] <-deseq_full_pipeline(phobj, 
                                           phname, vars2deseq, opt, 
                                           doPoscounts=FALSE, 
                                           all_combins=ALL_COMBINS,
                                          plot_all = TRUE,
                                          deseqname = "DeSEQ2_v2/",
                                          vars2heatmap = vars2heatmap
                                          )
  
}
opt <- restaurar(opt)
save(daa_all, file = paste0(opt$out, "DeSEQ2/DESEQ2_all.RData"))


