########################################
# Read OTU table
########################################

riga_tandas45  <- read_tsv(opt$metadata_riga_45)
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
write_tsv(classification, file = paste0(outdir, "classification.tsv"))
## otu table
s_otu_tab <- s_abund %>%
  dplyr::rename("taxonomy" = "#Classification") %>%
  dplyr::mutate(taxonomy = sub('.*\\|', '', taxonomy),
                taxonomy = gsub("s__", "", taxonomy)) %>%
  tibble::column_to_rownames(var = "taxonomy")

s_otu_tab_full <- s_otu_tab
otus_newnames <- ifelse(colnames(s_otu_tab) %in% riga_tandas45$RigaID, paste("C", riga_tandas45$codk2[match(colnames(s_otu_tab), riga_tandas45$RigaID)], sep=""), colnames(s_otu_tab))
xx =  data.frame(newnames = otus_newnames, oldnames = colnames(s_otu_tab)) # Check assignment

#s_otu_tab <- s_otu_tab[, !grepl("CBZ", colnames(s_otu_tab))]
names(s_otu_tab) <- otus_newnames

write_tsv(s_otu_tab %>% rownames_to_column("taxon") %>% select(taxon, everything()), file = paste0(outdir, "otu_tab_names_recoded.tsv"))
write_tsv(s_otu_tab_full %>% rownames_to_column("taxon") %>% select(taxon, everything()), file = paste0(outdir, "otu_tab_original.tsv"))
