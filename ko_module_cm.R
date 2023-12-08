#' @title omixer-rpmrR Celiacos
#' @author Sergio Romera

library(omixerRpm)
library(stringr)

dat <- read.table("/home/ccarlos/Documentos/202309_DEPRESION/results_cluster/mg13_humann3_all/merge2/humann3_merged_genetables_KO_CPM.tsv",
                  sep = "\t",
                  header = TRUE,
                  comment.char = "")

#' Quitamos los no-agrupados y no filtrados (sería interesante hacer unos
#' gráficos explicando cuántos tenemos de cada uno)
patterns <- c("UNMAPPED", "UNGROUPED")
dat_filt <- dat[(grep(paste0(patterns, collapse="|"), dat$X..Gene.Family, invert = TRUE )),]

# Vamos a generar una columna con la anotación, que tendremos que colocar justo antes de la anotación KO
gene_family <- dat_filt$X..Gene.Family

ko_terms <- sapply(gene_family, function(x){substr(x, 1, 6)})
names(ko_terms) <- ""
taxonomy <- sapply(gene_family, function(x){substr(x, 8,10000)})
names(taxonomy) <- ""


colnames(dat_filt)[1] <- "Gene_Families"
dat_filt$Gene_Families <- ko_terms
dat_filt$Taxonomy <- taxonomy
dat_filt <- dat_filt[,c(ncol(dat_filt),1:(ncol(dat_filt))-1)]
dat_filt$Taxonomy[dat_filt$Taxonomy==""] <- "Combined_taxa"


#' Lanzando omixR
mods <- rpm(dat_filt, minimum.coverage=0.3, annotation = 2, java.mem = 6)
mods_nosp <- rpm(dat_filt[dat_filt$Taxonomy=="Combined_taxa", 2:ncol(dat_filt)], annotation = 1, java.mem = 2)

# Load the default mapping database
db <- loadDB("GMMs.v1.07")
getNames(db, mods@annotation[2,2])

# # get the abundance|coverage as a data.frame with module id and description
coverage <- asDataFrame(mods, "coverage")
abundance <-  asDataFrame(mods, "abundance")

coverage_nosp <- asDataFrame(mods_nosp, "coverage")
abundance_nosp <-  asDataFrame(mods_nosp, "abundance")

#' Generar archivo de path_abundance como el de Humman3

mf_info <- paste0(abundance$Description,": ", abundance$`theObject@db@module.names[as.character(theObject@annotation$Module), `)
# taxonomy_clean <- gsub("Combined_taxa","", abundance$Taxon)
pathway <- paste0(mf_info,"|", abundance$Taxon)
pathway <- gsub("(\\|Combined_taxa)", "", pathway)


humann_coverage <- coverage
humann_abundance <- abundance

humann_coverage$`# Pathway` <- pathway
humann_coverage <- humann_coverage[,4:ncol(humann_coverage)]
humann_coverage <- humann_coverage[,c(ncol(humann_coverage),1:(ncol(humann_coverage))-1)]
colnames(humann_coverage) <- gsub("_merged_Abundance.CPM", "", colnames(humann_coverage))

humann_abundance$`# Pathway` <- pathway
humann_abundance <- humann_abundance[,4:ncol(humann_abundance)]
humann_abundance <- humann_abundance[,c(ncol(humann_abundance),1:(ncol(humann_abundance))-1)]
colnames(humann_abundance) <- gsub("_merged_Abundance.CPM", "", colnames(humann_abundance))

#' Separando los géneros y las especies para más agrupaciones a nivel taxonómico

genus_regex <-"g__([a-zA-Z]+)"
species_regex <- ".s__(\\w+)"

genus_names <- str_extract(coverage$Taxon, genus_regex)
genus_names <- gsub("g__","",genus_names)

species_names <- str_extract(coverage$Taxon, species_regex)
species_names <- gsub(".s__","",species_names)

abundance$Genus <- coverage$Genus <- genus_names
abundance$Species <- coverage$Species <- species_names

abundance <- abundance[,2:ncol(abundance)] %>% dplyr::select(Genus, Species, everything())

#abundance <- abundance[,c(132:133 ,1:(ncol(abundance))-1)]
#abundance <- abundance[,1:(ncol(abundance))-1]
colnames(abundance) <- gsub("_merged_Abundance.CPM", "", colnames(abundance))
colnames(abundance)[3] <- "Module"
colnames(abundance)[4] <- "Description"

coverage <- coverage[,2:ncol(coverage)]  %>% dplyr::select(Genus, Species, everything())
#coverage <- coverage[,c(132:133 ,1:(ncol(coverage))-1)]
#coverage <- coverage[,1:(ncol(coverage))-1]
colnames(coverage) <- gsub("_merged_Abundance.CPM", "", colnames(coverage))
colnames(coverage)[3] <- "Module"
colnames(coverage)[4] <- "Description"

#' Guardamos los datos

outdir <- "/home/ccarlos/Documentos/202309_DEPRESION/funcional1/kegg_modules/"
save(mods,file = paste0(outdir, "omixr_obj.Rda"))
save(abundance, file = paste0(outdir, "path_abundance.Rda"))
write_tsv(abundance, file = paste0(outdir, "path_abundance.tsv"))
save(coverage, file = paste0(outdir, "path_coverage.Rda"))
write_tsv(coverage, file = paste0(outdir, "path_coverage.tsv"))
save(humann_abundance, file = paste0(outdir, "path_abundance_humann.Rda"))
save(humann_coverage, file = paste0(outdir, "path_coverage_humann.Rda"))

save(abundance_nosp, file = paste0(outdir, "path_abundance_nosp.Rda"))
write_tsv(abundance_nosp, file = paste0(outdir, "path_abundance_nosp.tsv"))
save(coverage_nosp, file = paste0(outdir, "path_coverage_nosp.Rda"))
write_tsv(coverage_nosp, file = paste0(outdir, "path_coverage_nosp.tsv"))
