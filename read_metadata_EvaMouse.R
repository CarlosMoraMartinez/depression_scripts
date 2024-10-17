

num_reads <- s_otu_tab %>% colSums()

meta1 <- read_xlsx(opt$metadata) 


meta2 <- meta1 %>% 
  mutate(sample = gsub("A_", "", `Sample ID`),
         Mouse = gsub("IC[0-9a-zA-Z]+", "", sample, perl=TRUE),
         Mouse = ifelse(Mouse == "", "transfer", Mouse),
         Treatment_full = Treatment,
         Stress = ifelse(grepl("SD", Treatment), "SD", "Control"),
         Treatment = gsub("SD ", "", Treatment),
         num_reads = num_reads[sample],
         Treatment_region = paste(Treatment_full, Region_sequenced, sep=":")) 

assertthat::assert_that(all(meta2$sample %in% names(num_reads)))
assertthat::assert_that(all(names(num_reads) %in% meta2$sample))

meta3 <- meta2 %>% 
  inner_join(full_sample_names, by="sample") %>% 
  mutate(sampleID = sample) %>% 
  select(Num, sample, sampleID, everything())

write_tsv(meta3, paste0(input_tabs_dir, "/full_metadata.tsv"))
s_meta <- meta3 %>% column_to_rownames("sample")
