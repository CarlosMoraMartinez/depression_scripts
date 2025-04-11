########################################
# Read MetaData
########################################

s_meta <- read_tsv((opt$metadata)) %>% 
  filter(sampleID != "CN_AFL25-043") %>% 
  mutate(sampleID =  gsub("-", "_", sampleID)) %>% 
  data.frame
rownames(s_meta) <- s_meta$sampleID

all(s_meta$sampleID %in% names(s_otu_tab) )
all(names(s_otu_tab)  %in% s_meta$sampleID)

s_otu_tab <- s_otu_tab[, s_meta$sampleID]

all(rownames(s_meta) == colnames(s_otu_tab))

write_tsv(s_meta, paste0(input_tabs_dir, "/metadata.tsv"))

