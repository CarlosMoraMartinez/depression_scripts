
meta1 <- read_xlsx(opt$metadata) 


meta2 <- meta1 %>% 
  mutate(sample = gsub("A_", "", `Sample ID`),
         Mouse = gsub("IC[0-9a-zA-Z]+", "", sample, perl=TRUE),
         Mouse = ifelse(Mouse == "", "transfer", Mouse),
         Treatment_full = Treatment,
         Stress = ifelse(grepl("SD", Treatment), "SD", "Control"),
         Treatment = gsub("SD ", "", Treatment),
         #num_reads = num_reads[sample],
         Treatment_region = paste(Treatment_full, Region_sequenced, sep=":")) 

#assertthat::assert_that(all(meta2$sample %in% names(num_reads)))
#assertthat::assert_that(all(names(num_reads) %in% meta2$sample))


all_mpas<- all_mpas %>% mutate(
  metadata= map2(reads_per_sample, full_sample_names, \( rps, fsn){
    assertthat::assert_that( all(names(rps) %in% meta2$sample) )
    assertthat::assert_that( all( meta2$sample %in%  names(rps)) )
    aux <- meta2 %>% mutate(sample = factor(sample, levels=names(rps)))
    aux <- aux[order(aux$sample), ]
    assertthat::assert_that(all(aux$sample == names(rps)))
    assertthat::assert_that(all(aux$sample == fsn$sample))
    aux <- aux %>% 
      dplyr::mutate(num_reads = rps) %>% 
      inner_join(fsn, by="sample") %>% 
      mutate(sampleID = sample) %>% 
      select(Num, sample, sampleID, everything()) %>% 
      column_to_rownames("sample")
   return(aux)
    
  })
)

write_tsv(meta3, paste0(input_tabs_dir, "/full_metadata.tsv"))
s_meta <- meta3 %>% column_to_rownames("sample")

walk2(all_mpas$Condition, all_mpas$metadata, \(condname, tab) write_tsv(tab, 
                                                                                 paste0(input_tabs_dir, 
                                                                                        "/full_metadata_", 
                                                                                        janitor::make_clean_names(condname), 
                                                                                        ".tsv")))
save(all_mpas, file=paste0(input_tabs_dir, "/all_otu_tables.RData"))
