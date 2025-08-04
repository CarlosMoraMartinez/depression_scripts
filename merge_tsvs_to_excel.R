library(openxlsx)
library(tidyverse)

folder <- "/home/carlos/Escritorio/202311_DEPRESION/ReviewJune2025/Results_rstudio/results2/DAA_tables_2_save/"
files <- list.files(folder, pattern = "\\.tsv$", full.names = TRUE)

wb <- createWorkbook()

for (file in files) {
  sheet_name <- tools::file_path_sans_ext(basename(file))
  sheet_name <- gsub("\\*", "Int", sheet_name)
  sheet_name <- gsub("\\+", "Plus", sheet_name)
  if(nchar(sheet_name) > 31) sheet_name <- paste(strsplit(sheet_name, "")[[1]][1:31], sep="", collapse="")

  data <- read_tsv(file) %>% 
    dplyr::mutate(padj=ifelse(is.na(padj), 1, padj))
  addWorksheet(wb, sheet_name)
  writeData(wb, sheet_name, data)
}

saveWorkbook(wb, file = file.path(folder, "DAA_results_2.xlsx"), overwrite = TRUE)


## Linda DAA

folder <- "/home/carlos/Escritorio/202311_DEPRESION/ReviewJune2025/Results_rstudio/results2/DAA_linda/"
files <- list.files(folder, pattern = "\\.tsv$", full.names = TRUE)
files <- files[grep("^LinDA_", basename(files), perl=T)]

wb <- createWorkbook()

for (file in files) {
  sheet_name <- tools::file_path_sans_ext(basename(file))
  sheet_name <- gsub("\\*", "Int", sheet_name)
  sheet_name <- gsub("\\+", "Plus", sheet_name)
  if(nchar(sheet_name) > 31) sheet_name <- paste(strsplit(sheet_name, "")[[1]][1:31], sep="", collapse="")
  
  data <- read_tsv(file) %>% 
    dplyr::mutate(padj=ifelse(is.na(padj), 1, padj))
  addWorksheet(wb, sheet_name)
  writeData(wb, sheet_name, data)
}

saveWorkbook(wb, file = file.path(folder, "DAA_results_LinDA.xlsx"), overwrite = TRUE)



## Mediation

folder <- "/home/carlos/Escritorio/202311_DEPRESION/ReviewJune2025/Results_rstudio/results2/use/mediation/"
files <- list.files(folder, pattern = "\\.tsv$", full.names = TRUE)

wb <- createWorkbook()

for (file in files) {
  sheet_name <- tools::file_path_sans_ext(basename(file))
  sheet_name <- gsub("\\*", "Int", sheet_name)
  sheet_name <- gsub("\\+", "Plus", sheet_name)
  if(nchar(sheet_name) > 31) sheet_name <- paste(strsplit(sheet_name, "")[[1]][1:31], sep="", collapse="")
  
  if(! grepl("power_analysis_curve", file)){
  data <- read_tsv(file) %>% 
    dplyr::mutate(padj = p.adjust(p_raw, method = "BH")) %>% 
    dplyr::select(taxon, Effect, Estimate, S.E., `z-score`, p_raw, padj, color)

  }else{
    data <- read_tsv(file)
  }
  addWorksheet(wb, sheet_name)
  writeData(wb, sheet_name, data)
}

saveWorkbook(wb, file = file.path(folder, "mediation_results.xlsx"), overwrite = TRUE)

