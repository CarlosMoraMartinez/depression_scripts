########################################
# Read OTU table
########################################

s_abund <- read_tsv(paste0(opt$indir, "species.mpa.combined.clean2.txt"))
s_tax_tab <- s_abund %>%
  dplyr::rename("taxonomy" = "#Classification") %>%
  dplyr::select(taxonomy) %>%
  dplyr::mutate(Species = sub('.*\\|', '', taxonomy),
                Species = gsub("s__", "", Species),
                spec_row = Species) %>%
  dplyr::select(-taxonomy) %>%
  tibble::column_to_rownames(var = "spec_row")

##Parsing Kraken's taxonomic lineage strings
classification <- gsub("[a-z]__", "", s_abund$`#Classification`)
classification <- strsplit(classification, split = "\\|")
classification <- plyr::ldply(classification, rbind)
colnames(classification) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Strain")
rownames(classification) <- s_tax_tab$Species
#rownames(classification) <- rownames(s_tax_tab)

## otu table
s_otu_tab <- s_abund %>%
  dplyr::rename("taxonomy" = "#Classification") %>%
  dplyr::mutate(taxonomy = sub('.*\\|', '', taxonomy),
                taxonomy = gsub("s__", "", taxonomy)) %>%
  tibble::column_to_rownames(var = "taxonomy")
names(s_otu_tab) <- sapply(names(s_otu_tab), FUN=function(x) strsplit(x, '_')[[1]][1]) %>% 
  gsub("G4M", "", .) %>% gsub("^0", "", .)
