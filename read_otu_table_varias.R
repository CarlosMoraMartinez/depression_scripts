########################################
# Read OTU table
########################################
classnames <- list(d="Kingdom", k="Kingdom2", p="Phylum", c="Class", o="Order", f="Family", g="Genus", s="Species", t="Strain")

get_classif <- function(classtring, classnames=classnames, splitchar="\\|"){
  classvec <- strsplit(classtring, splitchar)[[1]] %>% 
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

get_classif_all <- function(otus_df, classnames=classnames, splitchar="\\|"){
  classification <- otus_df %>% 
    dplyr::pull(taxonomy) %>% 
    purrr::map(get_classif, 
               classnames=classnames, splitchar=splitchar) %>% 
    bind_rows() %>% 
    dplyr::select(-Kingdom2) 
  return(classification)
}

read_otu_table <- function(fname){
  s_abund <- read_tsv(fname) %>%
    dplyr::rename("taxonomy" = "#Classification")
  names(s_abund)[2:ncol(s_abund)] <- sapply(names(s_abund)[2:ncol(s_abund)], \(x)strsplit(x, "_")[[1]][2])
  return(s_abund)
}

get_samplenames_table <- function(fname){
  s_otu_tab <- read_tsv(fname) %>%
    dplyr::rename("taxonomy" = "#Classification") %>% 
    column_to_rownames("taxonomy") 
  full_sample_names <- data.frame(full_name = names(s_otu_tab),
                                  sample = sapply(names(s_otu_tab), FUN=function(x) strsplit(x, '_')[[1]][2]) ,
                                  flowcell = sapply(names(s_otu_tab), FUN=function(x) strsplit(x, '_')[[1]][4]) 
  )
  return(full_sample_names)
}

all_mpas <- read_tsv(opt$otu_tables_list) %>% 
  filter(Level == "species") %>% 
  mutate(otutable = map(File, read_otu_table),
         full_sample_names = map(File, get_samplenames_table),
         Classification = map(otutable, get_classif_all, 
                              classnames=classnames, splitchar="\\|"),
         otu_mat = map2(otutable, Classification, \(df, classif){
           df %>% dplyr::mutate(taxonomy = classif$Species) %>% 
             column_to_rownames("taxonomy") %>% 
             as.matrix
         }), 
         reads_per_sample = map(otu_mat, colSums)
         ) 

# check species can be used as rownames:
#map(all_mpas$Classification, \(x) nrow(x) == length(unique(x$Species)))

input_tabs_dir <- paste0(opt$out, "inputdata/")
if(!dir.exists(input_tabs_dir)) dir.create(input_tabs_dir)

walk2(all_mpas$Condition, all_mpas$otutable, \(condname, tab) write_tsv(tab, 
                                                                       paste0(input_tabs_dir, 
                                                                              "/raw_counts_table_", 
                                                                              janitor::make_clean_names(condname), 
                                                                              ".tsv")))
walk2(all_mpas$Condition, all_mpas$Classification, \(condname, tab) write_tsv(tab, 
                                                                        paste0(input_tabs_dir, 
                                                                               "/classification_", 
                                                                               janitor::make_clean_names(condname), 
                                                                               ".tsv")))
walk2(all_mpas$Condition, all_mpas$full_sample_names, \(condname, tab) write_tsv(tab, 
                                                                        paste0(input_tabs_dir, 
                                                                               "/full_sample_names_", 
                                                                               janitor::make_clean_names(condname), 
                                                                               ".tsv")))
save(all_mpas, file=paste0(input_tabs_dir, "/all_otu_tables.RData"))
