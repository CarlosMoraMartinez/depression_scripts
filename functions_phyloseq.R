
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



tax_glom_custom <- function(phobj, taxrank="Genus"){
  ttab <- tax_table(phobj) %>% data.frame %>% 
    rownames_to_column("taxname")
  otab <- otu_table(phobj) %>% data.frame %>% 
    rownames_to_column("taxname")
  snames <- names(otab)[2:ncol(otab)]
  
  all_levs <- names(ttab)[2:(which(names(ttab)==taxrank))]
  
  assertthat::assert_that(all(ttab$taxname == otab$taxname))
  
  mtab <- cbind(ttab, otab %>% select(-taxname)) %>% 
    mutate(!!taxrank := ifelse(is.na(.data[[taxrank]]),
                               paste("Unclassified ", taxname),
                               .data[[taxrank]] )
    ) %>% 
    #select(all_of(c(all_levs, snames))) %>% 
    group_by(!!sym(taxrank)) %>% 
    summarise(across(everything(), ~ ifelse(is.character(.x), 
                                            ifelse( cur_column() %in% all_levs, 
                                                    unique(.x)[1],
                                                    paste(unique(.x), sep=";", collapse=";")
                                            ), 
                                            sum(.x) 
    )
    ))
  ttab2 <- mtab %>% select(all_of(names(ttab))) %>% as.data.frame
  otab2 <- mtab %>% select(all_of(snames)) %>% as.data.frame
  rownames(ttab2) <- ttab2 %>% pull(!!sym(taxrank))
  rownames(otab2) <- rownames(ttab2)
  
  new_ps <-  phyloseq(sample_data(phobj),
                      otu_table(otab2, taxa_are_rows = TRUE),
                      tax_table(as.matrix(ttab2)))
  return(new_ps)
}