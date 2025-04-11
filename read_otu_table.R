########################################
# Read OTU table
########################################

input_tabs_dir <- paste0(opt$out, "inputdata/")
if(!dir.exists(input_tabs_dir)) dir.create(input_tabs_dir)

######

s_abund <- read_tsv(paste0(opt$indir, "best_tax_merged_freq_tax.tsv")) %>% 
  filter(!grepl("^#", id, perl=T)) %>% 
  select(-`CN_AFL25-043`) %>% # REMOVE CONTROL
  mutate(across(matches("Sib|Trb", perl=T), ~ as.numeric(.x)))

s_abund %>% sapply(class)
  
asv_tab  <- s_abund %>% select(id, Sequence, Taxon) %>% 
  mutate(ASV= paste("ASV", 1:nrow(.), sep="")) %>% 
  select(ASV, everything())

classification <- map(s_abund$Taxon, get_classif, classnames=classnames, splitchar="; ") %>% 
  bind_rows() %>% 
  mutate(ASV = asv_tab$ASV, id=asv_tab$id) %>% 
  mutate(Species = ifelse(Species == "Unclassified", ifelse(Genus=="Unclassified", 
                                                            paste("Unclassified ", ASV),
                                                            paste(Genus, " sp.", sep="")), Species)
         )%>% 
  select(ASV, id, everything())

write_tsv(classification, file = paste0(input_tabs_dir, "/classification.tsv"))
classification <- classification %>% select(-id) %>% column_to_rownames("ASV")


s_otu_tab <- s_abund %>%
  mutate(ASV = asv_tab$ASV) %>% 
  select(-id, -Sequence, - Taxon, -Confidence) %>% 
  tibble::column_to_rownames(var = "ASV") 

colnames(s_otu_tab) <- gsub("-", "_", colnames(s_otu_tab))
assertthat::assert_that(all(rownames(s_otu_tab) %in% rownames(classification)))



write_tsv(s_otu_tab %>% as.data.frame %>% rownames_to_column("taxon"), file = paste0(input_tabs_dir, "/raw_counts_table.tsv"))
write_tsv(asv_tab, file = paste0(input_tabs_dir, "/ASV_name_tab.tsv"))
