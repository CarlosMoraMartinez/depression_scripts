########################################
# Read OTU table
########################################

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

s_abund <- read_tsv(paste0(opt$indir, "species.mpa.combined.clean2.renamed.txt")) # "species.mpa.combined.clean2.txt"

s_tax_tab <- s_abund %>%
  dplyr::rename("taxonomy" = "#Classification") 
classification <- map(s_tax_tab$taxonomy, get_classif, classnames) %>% bind_rows()
write_tsv(classification, file = paste0(outdir, "classification.tsv"))

## otu table
s_otu_tab <- s_abund %>%
  dplyr::rename("taxonomy" = "#Classification") %>%
  #filter(taxonomy != "k__Eukaryota|k__Metazoa|p__Chordata|c__Mammalia|o__Primates|f__Hominidae|g__Homo|s__Homo_sapiens") %>% 
  dplyr::mutate(taxonomy = sub('.*\\|', '', taxonomy),
                taxonomy = gsub("s__", "", taxonomy)) %>%
  tibble::column_to_rownames(var = "taxonomy")

s_otu_tab_full <- s_otu_tab
rownames(classification) <- rownames(s_otu_tab_full)


write_tsv(s_otu_tab_full %>% rownames_to_column("taxon"), file = paste0(outdir, "otu_tab.tsv"))
