########################################
# Read OTU table
########################################

riga_tandas45  <- read_tsv(opt$metadata_riga_45)
# s_abund <- read_tsv(paste0(opt$indir, "species.mpa.combined.clean2.txt"))
# 
# s_tax_tab <- s_abund %>%
#   dplyr::rename("taxonomy" = "#Classification") %>%
#   dplyr::select(taxonomy) %>%
#   dplyr::mutate(Species = sub('.*\\|', '', taxonomy),
#                 Species = gsub("s__", "", Species),
#                 spec_row = Species) %>%
#   dplyr::select(-taxonomy) %>%
#   tibble::column_to_rownames(var = "spec_row")

classnames <- list(d="Kingdom", k="Kingdom2", p="Phylum", c="Class", o="Order", f="Family", g="Genus", s="Species", xx="Strain")

get_classif <- function(classtring, classnames=classnames){
  classvec <- strsplit(classtring, "\\|")[[1]] %>% 
    strsplit("__")
  classlist <- map(classvec, \(x)x[2]) %>% unlist
  class_init <- map(classvec, \(x)x[1]) %>% unlist
  if(length(which(class_init == "k")) > 1){
    class_init[which(class_init == "k")[1]] <- "d"
  }
  names(classlist) <- class_init
  #cat(classtring, "\n") 
  
  aux <- data.frame(matrix(ncol=length(classnames), nrow=1, dimnames = list(NULL, unlist(classnames))))
  aux[1, unlist(classnames[class_init])] <- classlist[class_init]
  return(aux)
}

s_abund <- read_tsv(paste0(opt$indir, "species.mpa.combined.clean2.changednames.txt")) # "species.mpa.combined.clean2.txt"
s_abund <- s_abund %>% select(-C630900, -C632000) # Estas dos se hicieron con otra versión de la DB de Kraken por error
## Eliminar especies que solo aparecian en esas dos especies (en su mayoria, eran especies repetidas)
spsum <- s_abund %>% select(- `#Classification`) %>% rowSums
table(spsum == 0)
s_abund <- s_abund[spsum>0, ]


s_tax_tab <- s_abund %>%
  dplyr::rename("taxonomy" = "#Classification") 
classification <- map(s_tax_tab$taxonomy, get_classif, classnames) %>% bind_rows()
write_tsv(classification, file = paste0(outdir, "classification.tsv"))

##Parsing Kraken's taxonomic lineage strings
#classification <- gsub("[a-z]__", "", s_abund$`#Classification`)
#classification <- strsplit(classification, split = "\\|")
#classification <- plyr::ldply(classification, rbind)
#colnames(classification) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Strain")
#rownames(classification) <- s_tax_tab$Species
##rownames(classification) <- rownames(s_tax_tab)
#write_tsv(classification, file = paste0(outdir, "classification.tsv"))

## otu table
s_otu_tab <- s_abund %>%
  dplyr::rename("taxonomy" = "#Classification") %>%
  #filter(taxonomy != "k__Eukaryota|k__Metazoa|p__Chordata|c__Mammalia|o__Primates|f__Hominidae|g__Homo|s__Homo_sapiens") %>% 
  dplyr::mutate(taxonomy = sub('.*\\|', '', taxonomy),
                taxonomy = gsub("s__", "", taxonomy)) %>%
  tibble::column_to_rownames(var = "taxonomy")

s_otu_tab_full <- s_otu_tab
rownames(classification) <- rownames(s_otu_tab_full)

otus_newnames <- ifelse(colnames(s_otu_tab) %in% riga_tandas45$RigaID, paste("C", riga_tandas45$codk2[match(colnames(s_otu_tab), riga_tandas45$RigaID)], sep=""), colnames(s_otu_tab))
xx =  data.frame(newnames = otus_newnames, oldnames = colnames(s_otu_tab)) # Check assignment

#s_otu_tab <- s_otu_tab[, !grepl("CBZ", colnames(s_otu_tab))]
names(s_otu_tab) <- otus_newnames

write_tsv(s_otu_tab %>% rownames_to_column("taxon") %>% select(taxon, everything()), file = paste0(outdir, "otu_tab_names_recoded.tsv"))
write_tsv(s_otu_tab_full %>% rownames_to_column("taxon") %>% select(taxon, everything()), file = paste0(outdir, "otu_tab_original.tsv"))
