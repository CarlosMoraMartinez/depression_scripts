########################################
# Read OTU table
########################################

classnames <- list(d="Kingdom", k="Kingdom2", p="Phylum", c="Class", o="Order", f="Family", g="Genus", s="Species", t="Strain")

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

s_abund <- read_tsv(paste0(opt$indir, "metaphlan_mpa_vJun23_CHOCOPhlAnSGB_202403.renamed.txt")) %>%  # "species.mpa.combined.clean2.txt" %>% 
   filter(grepl("s__", clade_name))

s_tax_tab <- s_abund %>%
  dplyr::rename("taxonomy" = "clade_name") 
classification <- map(s_tax_tab$taxonomy, get_classif, classnames) %>% bind_rows()
write_tsv(classification, file = paste0(outdir, "classification_with_strain.tsv"))

yy <- s_abund %>% column_to_rownames("clade_name") %>% 
  as.matrix %>% t %>% 
  data.frame
names(yy) <- paste(classification$Species, classification$Strain, sep="::")
pdf(paste0(outdir, "/count_histograms.pdf"))
names(yy) %>%  map( \(x) hist(log10(yy[, x]+1),  breaks=30, main=x))
dev.off()

classification_sp <- classification %>% select(-Strain) %>% 
  distinct()
write_tsv(classification_sp, file = paste0(outdir, "classification_onlySpecies.tsv"))
  
classification_reserva <- classification
classification <- classification %>% mutate(
  Species = ifelse(is.na(Strain), Species, paste(Species, Strain, sep="_"))
)
write_tsv(classification, file = paste0(outdir, "classification_with_strain_pasted.tsv"))

## otu table
s_otu_tab <- s_abund %>%
  dplyr::rename("taxonomy" = "clade_name") %>%
  #filter(taxonomy != "k__Eukaryota|k__Metazoa|p__Chordata|c__Mammalia|o__Primates|f__Hominidae|g__Homo|s__Homo_sapiens") %>% 
  dplyr::mutate(taxonomy = classification$Species) %>%
  tibble::column_to_rownames(var = "taxonomy")

s_otu_tab_full <- s_otu_tab
rownames(classification) <- rownames(s_otu_tab_full)


write_tsv(s_otu_tab_full %>% rownames_to_column("taxon"), file = paste0(outdir, "otu_tab.tsv"))


## OtuTable by Species only
s_otu_tab_sp <- s_abund %>%
  dplyr::rename("taxonomy" = "clade_name") %>%
  dplyr::mutate(taxonomy = classification_reserva$Species) %>%
  group_by(taxonomy) %>% 
  summarise_all(sum) %>% 
  tibble::column_to_rownames(var = "taxonomy")

s_otu_tab_full_sp <- s_otu_tab_sp
rownames(classification_sp) <- classification_sp$Species
classification_sp <- classification_sp[rownames(s_otu_tab_full_sp), ]

write_tsv(s_otu_tab_full_sp %>% rownames_to_column("taxon"), file = paste0(outdir, "otu_tab_summedSpecies.tsv"))


