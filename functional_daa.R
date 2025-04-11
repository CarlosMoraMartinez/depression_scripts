
source(opt$functional_functions)

opt <- restaurar(opt)
load(paste0(opt$out, "DeSEQ2_v9/DESEQ2_all.RData"))
opt$out <- paste0(opt$out, "Functional/")
if(!dir.exists(opt$out)) dir.create(opt$out)


s_meta <- sample_data(all_phyloseq$MGBC_plus_NRA_raw) %>% data.frame
samples2keep <- s_meta %>% pull(sampleID)
samples2keep_control <- s_meta %>% filter(Treatment == "NO ABS" & Stress == "Control") %>%  pull(sampleID)
samples2keep_stress_r1 <- s_meta %>% filter(Treatment == "NO ABS" & Region_sequenced == "REG1") %>%  pull(sampleID)
samples2keep_stress_r2 <- s_meta %>% filter(Treatment == "NO ABS" & Region_sequenced == "REG2") %>%  pull(sampleID)
samples2keep_stress_r3 <- s_meta %>% filter(Treatment == "NO ABS" & Region_sequenced == "REG3") %>%  pull(sampleID)


metacyc_ab <-  readFunctionalMatrix(opt, "humann3_merged_abundances_CPM_named.tsv", sample_substring_index=2) #%>% dplyr::select(all_of(c("Pathway",samples2keep)))
#metacyc_rxn <- readFunctionalMatrix(opt, "humann3_merged_genetables_RXN_CPM_named.tsv", sample_substring_index=2) #%>% dplyr::select(all_of(c("Pathway",samples2keep)))
#ko <- readFunctionalMatrix(opt, "humann3_merged_genetables_KO_CPM_named.tsv", sample_substring_index=2) # %>% dplyr::select(all_of(c("Pathway",samples2keep)))
#go <- readFunctionalMatrix(opt, "humann3_merged_genetables_GO_CPM_named.tsv", sample_substring_index=2) #%>% dplyr::select(all_of(c("Pathway",samples2keep)))
#cazy <- readFunctionalMatrix(opt, "humann3_merged_genetables_CAZY_CPM.tsv", sample_substring_index=2) #%>% dplyr::select(all_of(c("Pathway",samples2keep)))



tabs2test <- list("MetaCyc"=metacyc_ab) #, "MetaCyc Reactions"=metacyc_rxn, "KEGG"=ko, "Gene Ontology"=go, "CAZy" = cazy, "GOmixerModules"=modules_nosp2 )
save(tabs2test, file = paste0(opt$out, "functional_tabs.RData"))
tabs2annot <- list("MetaCyc"=NULL) #, "MetaCyc Reactions"=metacyc_mapnames, "KEGG"=ko_mapnames, "Gene Ontology"=go_mapnames, "CAZy" = cazy_mapnames, "GOmixerModules"=modules_mapnames_nosp )

walk2(tabs2test, paste0(opt$out, "/functTabInput_",names(tabs2test), ".tsv"), \(df, name)if(!is.null(df))write_tsv(df, file=name))


#### Limma
#met2use <- sample_data(all_phyloseq$filt) %>% data.frame %>%  
#  dplyr::filter(status_c2 != "Insufficient gain") %>% 
#  dplyr::filter(status_c2 != "initially_overweight")

met2use <- sample_data(all_phyloseq$MGBC_plus_416_raw) %>% data.frame 

opt$minsampleswithcount <- opt$minfreq*nrow(met2use)
filtered_samples <- lapply(tabs2test, \(x) x %>% dplyr::select(Pathway, all_of(met2use$sampleID)))
filtered_byproc <- lapply(filtered_samples, FUN=filterWholeProcessesAndFreq, opt)
walk2(filtered_byproc, paste0(opt$out, "/functTabInput_",names(tabs2test), "_filteredByProcess.tsv"), \(df, name)if(!is.null(df))write_tsv(df, file=name))

#limmares_byproc <- lapply(filtered_byproc, FUN=limma4functional, met2use, "Stress") #"status_c2" "edad_00"
limmares_byproc_r2r1 <- limma4functional(filtered_byproc$MetaCyc %>% select(Pathway, all_of(samples2keep_control)),
                 met2use %>% filter(sampleID %in% samples2keep_control),
                 "Region_sequenced", 
                 covars=c(),
                 levs2compare = c("REG1", "REG2")
)
limmares_byproc_r3r1 <- limma4functional(filtered_byproc$MetaCyc %>% select(Pathway, all_of(samples2keep_control)),
                                         met2use %>% filter(sampleID %in% samples2keep_control),
                                         "Region_sequenced", 
                                         covars=c(),
                                         levs2compare = c("REG1", "REG3")
)

limmares_byproc_r3r2 <- limma4functional(filtered_byproc$MetaCyc %>% select(Pathway, all_of(samples2keep_control)),
                                         met2use %>% filter(sampleID %in% samples2keep_control),
                                         "Region_sequenced", 
                                         covars=c(),
                                         levs2compare = c("REG2", "REG3")
)

limmares_byproc_stressr1 <- limma4functional(filtered_byproc$MetaCyc %>% select(Pathway, all_of(samples2keep_stress_r1)),
                                         met2use %>% filter(sampleID %in% samples2keep_stress_r1),
                                         "Stress", 
                                         covars=c(),
                                         levs2compare = c("Control", "SD")
)
limmares_byproc_stressr2 <- limma4functional(filtered_byproc$MetaCyc %>% select(Pathway, all_of(samples2keep_stress_r2)),
                                             met2use %>% filter(sampleID %in% samples2keep_stress_r2),
                                             "Stress", 
                                             covars=c(),
                                             levs2compare = c("Control", "SD")
)
limmares_byproc_stressr3 <- limma4functional(filtered_byproc$MetaCyc %>% select(Pathway, all_of(samples2keep_stress_r3)),
                                             met2use %>% filter(sampleID %in% samples2keep_stress_r3),
                                             "Stress", 
                                             covars=c(),
                                             levs2compare = c("Control", "SD")
)

limmares_byproc<- list(
  REG2_vs_REG1 = limmares_byproc_r2r1,
  REG3_vs_REG1 = limmares_byproc_r3r1,
  REG3_vs_REG2 = limmares_byproc_r3r2,
  REG1__Stress_vs_Control = limmares_byproc_stressr1,
  REG2__Stress_vs_Control = limmares_byproc_stressr2, 
  REG3__Stress_vs_Control = limmares_byproc_stressr3
)

limmares_byproc_annot <- mapply(limmares_byproc, tabs2annot, FUN=nameProcesses, SIMPLIFY = FALSE)
walk2(limmares_byproc_annot, paste0(opt$out, "/DAAlimma_process_",names(limmares_byproc_annot), ".tsv"), \(df, name)write_tsv(df%>% rownames_to_column("Process"), file=name))

restabs_limma <- lapply(limmares_byproc, getSummaryTablesDeseq, opt)
volcanos_limma <- lapply(names(limmares_byproc), 
                         FUN=function(name){
                           fname <- paste0("limma_volcano_SumProcesses_praw_", gsub(" ", "", name), sep="")
                           make_volcano(res = limmares_byproc[[name]], opt=opt, name=fname, pcol="pvalue")}
)
volcanos_limma_padj <- lapply(names(limmares_byproc), 
                              FUN=function(name){
                                fname <- paste0("limma_volcano_SumProcesses_padj_", gsub(" ", "", name), sep="")
                                make_volcano(res = limmares_byproc[[name]], opt=opt, name=fname, pcol="padj")}
)

plims <- rep(0.05, length(limmares_byproc)) #, 0.001, 0.001, 0.001, 0.05, 0.05)
names(plims) <- names(limmares_byproc)
barplots_limma <- lapply(names(limmares_byproc_annot), \(x){
  tab <- limmares_byproc_annot[[x]]
  include_longnames <- TRUE
  if(x == "CAZy"){
    include_longnames <- FALSE
    rownames(tab) <- paste(getCazyClass(rownames(tab)), rownames(tab), sep=" - ")
  }
  makeBarplotFunctional(tab, 
                        plim=plims[x], 
                        paste0("BarplotProcess_", x), opt$out, 
                        include_longnames = include_longnames, 
                        w=12, h=8)
})
# heatmaps_byproc <- lapply(1:length(limmares_byproc), FUN=function(i){
#             makeHeatmapFunctional(limmares_byproc[[i]], met2use, filtered_byproc[[i]],
#                         variable = "Condition",
#                         opt, 
#                         name = paste0("heatmap_byProcess_padj_", gsub(" ", "", names(limmares_byproc)[i]), sep=""), 
#                         logscale=FALSE, 
#                         ptype = "padj", w=20, h=14)
# })


### Lollipop 

# Read MetaCyc ref names
metacyc_namerefs <- read.delim("/home/carlos/projects/EvaMouse_Otcubre2024/Metadata/MetaCyc_info/MetaCyc_pasted_all.txt", 
                               head=F, sep="\t")
names(metacyc_namerefs) <- c("PathID", "Name")

metacyc_namerefs <- metacyc_namerefs %>% 
  mutate(Name = gsub("<i>", "", Name),
         Name = gsub("</i>", "", Name))
foundnames <- rownames(limmares_byproc$REG2_vs_REG1) %>% gsub(": NO_NAME", "", .)
table(foundnames %in% metacyc_namerefs$PathID)
foundnames[!foundnames %in% metacyc_namerefs$PathID]

### diff between regions
PLIM = 0.01
PLIM_COL = 0.05

region_comp <- limmares_byproc[c("REG2_vs_REG1", "REG3_vs_REG1", "REG3_vs_REG2")]
regioncomp_df <- map2(region_comp, names(region_comp), 
                      .f = \(x, y) x %>% rownames_to_column("Pathway") %>% 
                        mutate(Condition=gsub("_", " ", y)) ) %>% 
  bind_rows() %>% 
  mutate(Pathway = gsub(": NO_NAME", "", Pathway)) %>% 
  mutate(Pathway = paste(Pathway, metacyc_namerefs$Name[match(Pathway, metacyc_namerefs$PathID)], sep=":")) %>% 
  mutate(Condition = gsub("REG", "R", Condition)) %>% 
  mutate(Condition = factor(Condition, 
                              levels = c("R3 vs R2", "R3 vs R1", "R2 vs R1")))

grouptab <- regioncomp_df %>% select(Pathway, padj, Condition) %>% 
  mutate(Condition = gsub(" ", "_", Condition)) %>% 
  spread(Condition, padj) %>% 
  mutate(group = ifelse(R3_vs_R1 <= PLIM_COL & R3_vs_R2 <= PLIM_COL, "R3_vs_all", 
                        ifelse(R3_vs_R1 <= PLIM_COL, "R3_vs_R1", 
                               ifelse(R3_vs_R2 <= PLIM_COL, "R3_vs_R2", "NS"))))

group_names <- c("R3_vs_all", "R3_vs_R2", "R3_vs_R1", "NS") %>% rev
regioncomp_df <- regioncomp_df %>% 
  mutate(Condition_sig = ifelse(padj <= PLIM_COL, as.character(Condition), "NS")) %>% 
  mutate(Condition_sig = factor(Condition_sig, levels=c("NS", "R3 vs R2", "R3 vs R1", "R2 vs R1"))) %>% 
  mutate(group = grouptab$group[match(Pathway, grouptab$Pathway)]) %>% 
  mutate(group = factor(group, levels = group_names)) 

ordered_func <- regioncomp_df %>% group_by(group, Pathway) %>% 
  filter(Condition != "R2 vs R1") %>% 
  summarise(mean_LFC=mean(log2FoldChange)) %>% 
  arrange(group, mean_LFC)

regioncomp_df <- regioncomp_df %>% mutate(Pathway = factor(Pathway, levels = ordered_func$Pathway)) 

paths2show <- regioncomp_df %>% filter(padj < PLIM) %>% 
  pull(Pathway) %>% unique()


aux2 <- regioncomp_df %>% filter(Pathway %in% paths2show) 

group_breaks <- ordered_func %>% filter(Pathway %in% paths2show) |>
  dplyr::mutate(Pathway = factor(Pathway, levels = unique(Pathway))) |>
  dplyr::group_by(group) |>
  dplyr::summarise(last = max(as.numeric(Pathway)), .groups = "drop") |>
  dplyr::pull(last)
group_breaks <- group_breaks[1:(length(group_breaks)-1)]

colscsig <- pal_d3()(3)
colscsig <- c("gray", colscsig)

(glol <- ggplot(aux2, aes(x=Pathway,y=log2FoldChange, col=Condition_sig, fill=Condition_sig))+
  facet_wrap(~ Condition, nrow=1) + 
  geom_segment( aes(x=Pathway, xend=Pathway, y=0, yend=log2FoldChange)) +
  geom_point()+
  scale_color_manual(values=colscsig)+
  scale_fill_manual(values=colscsig)+
  theme_bw() +
  theme(axis.text.x = element_text(size = 12))+
  theme(strip.text.x = element_text(size = 14))+
  theme(strip.text.y = element_text(size = 12))+
  theme(axis.title.y = element_text(size = 12))+
  theme(axis.title.x = element_text(size = 12))+
  theme(axis.text.y = element_text( size = 11, face="italic")) +
  coord_flip() +
  ylab("LFC between regions ") +
    geom_vline(xintercept = group_breaks + 0.5, linetype = "dashed", color = "gray40")
)
ggsave(paste0(opt$out, "_functional_lollipop_regions.pdf"), glol, width = 10, height = 10,
       limitsize = F)


## diff between Stress and controls in different regions

PLIM = 0.05
PLIM_COL = 0.05


stress_comp <- limmares_byproc[c("REG1__Stress_vs_Control", "REG2__Stress_vs_Control", "REG3__Stress_vs_Control")]
stress_df <- map2(stress_comp, names(stress_comp), 
                      .f = \(x, y) x %>% rownames_to_column("Pathway") %>% 
                        mutate(Condition=strsplit(y, "__")[[1]][1]) ) %>% 
  bind_rows() %>% 
  mutate(Pathway = gsub(": NO_NAME", "", Pathway)) %>% 
  mutate(Pathway = paste(Pathway, metacyc_namerefs$Name[match(Pathway, metacyc_namerefs$PathID)], sep=":")) %>% 
  #mutate(Condition = gsub("REG", "R", Condition)) %>% 
  mutate(Condition = factor(Condition, 
                            levels = c("REG1", "REG2", "REG3")))

grouptab_stress <- stress_df %>% select(Pathway, padj, Condition) %>% 
  #mutate(Condition = gsub(" ", "_", Condition)) %>% 
  spread(Condition, padj) %>% 
  rowwise() %>%
  mutate(group = paste(c("R1", "R2", "R3")[c_across(REG1:REG3) <= PLIM_COL], collapse = "_") ) %>% 
  ungroup() %>% 
  mutate(group = ifelse(group == "", "NS", group))

group_names <- c("R1_R2_R3", "R1_R2", "R2_R3", "R1", "R2", "R3", "NS") %>% rev
stress_df <- stress_df %>% 
  mutate(Condition_sig = ifelse(padj <= PLIM_COL, as.character(Condition), "NS")) %>% 
  mutate(Condition_sig = factor(Condition_sig, levels=c("NS", "REG3", "REG2", "REG1"))) %>% 
  mutate(group = grouptab_stress$group[match(Pathway, grouptab_stress$Pathway)]) %>% 
  mutate(group = factor(group, levels = group_names)) 

ordered_func <- stress_df %>% group_by(group, Pathway) %>% 
  summarise(mean_LFC=mean(log2FoldChange)) %>% 
  arrange(group, mean_LFC)

stress_df <- stress_df %>% mutate(Pathway = factor(Pathway, levels = ordered_func$Pathway)) 

paths2show <- stress_df %>% filter(padj < PLIM) %>% 
  pull(Pathway) %>% unique()

aux3 <- stress_df %>% filter(Pathway %in% paths2show) 

group_breaks <- ordered_func %>% filter(Pathway %in% paths2show) |>
  dplyr::mutate(Pathway = factor(Pathway, levels = unique(Pathway))) |>
  dplyr::group_by(group) |>
  dplyr::summarise(last = max(as.numeric(Pathway)), .groups = "drop") |>
  dplyr::pull(last)
group_breaks <- group_breaks[1:(length(group_breaks)-1)]

colscsig <- pal_npg()(3) %>% rev
colscsig <- c("gray", colscsig)

(glol <- ggplot(aux3, aes(x=Pathway,y=log2FoldChange, col=Condition_sig, fill=Condition_sig))+
    facet_wrap(~ Condition, nrow=1) + 
    geom_segment( aes(x=Pathway, xend=Pathway, y=0, yend=log2FoldChange)) +
    geom_point()+
    scale_color_manual(values=colscsig)+
    scale_fill_manual(values=colscsig)+
    theme_bw() +
    theme(axis.text.x = element_text(size = 12))+
    theme(strip.text.x = element_text(size = 14))+
    theme(strip.text.y = element_text(size = 12))+
    theme(axis.title.y = element_text(size = 12))+
    theme(axis.title.x = element_text(size = 12))+
    theme(axis.text.y = element_text( size = 11)) + #, face="italic"
    coord_flip() +
    ylab("LFC between regions ") +
    geom_vline(xintercept = group_breaks + 0.5, linetype = "dashed", color = "gray40")
)
ggsave(paste0(opt$out, "_functional_lollipop_stress.pdf"), glol, width = 12, height = 8,
       limitsize = F)
