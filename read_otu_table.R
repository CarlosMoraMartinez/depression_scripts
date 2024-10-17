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

input_tabs_dir <- paste0(opt$out, "inputdata/")
if(!dir.exists(input_tabs_dir)) dir.create(input_tabs_dir)

######

s_abund <- read_tsv(paste0(opt$indir, "species.mpa.combined.clean2.txt"))
#s_tax_tab <- s_abund %>%
#  dplyr::rename("taxonomy" = "#Classification") %>%
#  dplyr::select(taxonomy) %>%
#  dplyr::mutate(Species = sub('.*\\|', '', taxonomy),
#                Species = gsub("s__", "", Species),
#                spec_row = Species) %>%
#  dplyr::select(-taxonomy) %>%
#  tibble::column_to_rownames(var = "spec_row")

s_tax_tab <- s_abund %>%
  dplyr::rename("taxonomy" = "#Classification") 
classification <- map(s_tax_tab$taxonomy, get_classif, classnames) %>% bind_rows()

write_tsv(classification, file = paste0(input_tabs_dir, "/classification.tsv"))
rownames(classification) <- classification$Species

##Parsing Kraken's taxonomic lineage strings
#classification <- gsub("[a-z]__", "", s_abund$`#Classification`)
#classification <- strsplit(classification, split = "\\|")
#classification <- plyr::ldply(classification, rbind)
#colnames(classification) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Strain")
#rownames(classification) <- s_tax_tab$Species
#rownames(classification) <- rownames(s_tax_tab)

## otu table
s_otu_tab <- s_abund %>%
  dplyr::rename("taxonomy" = "#Classification") %>%
  dplyr::mutate(taxonomy = sub('.*\\|', '', taxonomy),
                taxonomy = gsub("s__", "", taxonomy)) %>%
  tibble::column_to_rownames(var = "taxonomy")

assertthat::assert_that(all(rownames(s_otu_tab) %in% rownames(classification)))

full_sample_names <- data.frame(full_name = names(s_otu_tab),
                                sample = sapply(names(s_otu_tab), FUN=function(x) strsplit(x, '_')[[1]][2]) ,
                                flowcell = sapply(names(s_otu_tab), FUN=function(x) strsplit(x, '_')[[1]][4]) 
)
names(s_otu_tab) <- full_sample_names$sample


write_tsv(s_otu_tab %>% as.data.frame %>% rownames_to_column("taxon"), file = paste0(input_tabs_dir, "/raw_counts_table.tsv"))
write_tsv(full_sample_names, file = paste0(input_tabs_dir, "/full_sample_names.tsv"))
