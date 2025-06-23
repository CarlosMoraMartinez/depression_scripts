library(tidyverse)
library(wesanderson)
library(LinDA)
library(phyloseq)
library(ggvenn)
library(ggVennDiagram)

pal <- wes_palette("AsteroidCity1", 2)
source("/home/carlos/Escritorio/202311_DEPRESION/ReviewJune2025/scripts/depression_scripts/metagenomics_core_functions.R")
source("/home/carlos/Escritorio/202311_DEPRESION/ReviewJune2025/scripts/depression_scripts/mediation_functions.R")


outdir <- "/home/carlos/Escritorio/202311_DEPRESION/ReviewJune2025/Results_rstudio/results2/DAA_linda/"
if(!dir.exists(outdir)) dir.create(outdir)

#load("/home/carlos/Escritorio/202311_DEPRESION/ReviewJune2025/Results_rstudio/results2/phyloseq_original//phyloseq_all_list.RData")
load("/home/carlos/Escritorio/202311_DEPRESION/ReviewJune2025/Results_rstudio/results2/DESeq2_AgeSexInteraction/Integrate1/phyloseq_used_remove_tanda2.RData")


otu.tab <- otu_table(phobj) %>% data.frame
meta <- sample_data(phobj) %>% data.frame

linda.obj_cond <- linda(otu.tab, meta, formula = '~Condition', 
                   alpha = 0.05,
                   prev.cut = 0.05, 
                   lib.cut = 1000, 
                   winsor.quan = 0.97)

linda.obj_bmi <- linda(otu.tab, meta, formula = '~BMI_log', 
                   alpha = 0.05,
                   prev.cut = 0.05, 
                   lib.cut = 1000, 
                   winsor.quan = 0.97)

linda.obj_bmiAdj <- linda(otu.tab, meta, formula = '~BMI_log + Sex*Age_log', 
                   alpha = 0.05,
                   prev.cut = 0.05, 
                   lib.cut = 1000, 
                   winsor.quan = 0.97)

linda.obj_condAdj <- linda(otu.tab, meta, formula = '~Condition + Sex*Age_log', 
                   alpha = 0.05,
                   prev.cut = 0.05, 
                   lib.cut = 1000, 
                   winsor.quan = 0.97)

linda.obj_adjAll <- linda(otu.tab, meta, formula = '~Condition + BMI_log + Sex*Age_log', 
                   alpha = 0.05,
                   prev.cut = 0.05, 
                   lib.cut = 1000, 
                   winsor.quan = 0.97)


linda.obj_adjAll_nolog <- linda(otu.tab, meta, formula = '~Condition + BMI + Sex*Age', 
                          alpha = 0.05,
                          prev.cut = 0.05, 
                          lib.cut = 1000, 
                          winsor.quan = 0.97)

linda.obj_age_nolog <- linda(otu.tab, meta, formula = '~Age', 
                                alpha = 0.05,
                                prev.cut = 0.05, 
                                lib.cut = 1000, 
                                winsor.quan = 0.97)

linda.obj_age_log <- linda(otu.tab, meta, formula = '~Age_log', 
                                alpha = 0.05,
                                prev.cut = 0.05, 
                                lib.cut = 1000, 
                                winsor.quan = 0.97)

pp <- linda.plot(linda.obj_adjAll, c('Condition'), 
           titles = c('Control', 'MDD'), alpha = 0.05, lfc.cut = 1,
           legend = TRUE, directory = NULL, width = 11, height = 8)

#L <- matrix(c(0, 1, 0, 0, 0, 1), nrow = 2, byrow = TRUE)
#L <- matrix(c(0, 1, 0), nrow = 1, byrow = TRUE)
#linda.wald.test(linda.obj, L, 'LM', alpha = 0.1, p.adj.method="BH")

lindalist <- list(
  "D_vs_C_alone"=linda.obj_cond$output$ConditionDepression,
  "D_vs_C_adj.Age*Sex"=linda.obj_condAdj$output$ConditionDepression,
  "D_vs_C_adj.Age*Sex+BMI"=linda.obj_adjAll$output$ConditionDepression,
  "BMI_alone"=linda.obj_bmi$output$BMI_log,
  "BMI_adj.Age*Sex"=linda.obj_bmiAdj$output$BMI_log,
  "BMI_adj.Age*Sex+Depr"=linda.obj_adjAll$output$BMI_log
) %>% 
  map(\(x){
    aux <- x %>% 
      dplyr::mutate(
        log2FoldChangeShrink = log2FoldChange,
        lfcSE_Shrink = lfcSE
      ) %>% 
      rownames_to_column("taxon")
    return(list(resdf=aux))
  })
lindalist <- list(
  firstContrast = lindalist[[1]],
  contrastlist2 = lindalist[2:6]
)
save(lindalist, file = paste0(outdir, "lindalist.RData"))

load(paste0("/home/carlos/Escritorio/202311_DEPRESION/ReviewJune2025/Results_rstudio/results2/", "DESeq2_AgeSexInteraction_log/Integrate1/LFC_Comparison_AgeISex_BMI_allCombos.RData"))
load(paste0("/home/carlos/Escritorio/202311_DEPRESION/ReviewJune2025/Results_rstudio/results2/", "DESeq2_AgeSexInteraction_log/all_DESeq2.RData"))


######

## Venn diagram of MDD related bacteria
plim <- 0.05 
vars2venn <- list(
  "D vs C" = lindalist$firstContrast$resdf %>% dplyr::filter(padj < plim) %>% pull(taxon),
  "D vs C\nadj Sex*Age" = lindalist$contrastlist2$`D_vs_C_adj.Age*Sex`$resdf %>% dplyr::filter(padj < plim) %>% pull(taxon),
  "D vs C\nadj Sex*Age+BMI" = lindalist$contrastlist2$`D_vs_C_adj.Age*Sex+BMI`$resdf %>% dplyr::filter(padj < plim) %>% pull(taxon)
)
(gv <- ggvenn(
  vars2venn, columns = names(vars2venn),
  stroke_size = 0.5,
  stroke_color = C_NS,
  fill_color = c(C_CASE, C_CTRL2, C_CTRL, C_CASE2),show_elements = F
)+
  theme(
    plot.margin = margin(t = 20, r = 20, b = 20, l = 20),  # top, right, bottom, left
    text = element_text(size = 12)  # optional: tweak font size
  )
)
ggsave(filename = paste0(outdir, "VennDiagram_CvsD_controlSexIntAgeBMI_p", as.character(plim), ".pdf"), gv, width = 8, height = 8)


(gv <- ggVennDiagram(vars2venn) +
    scale_fill_gradient(low = "white", high = "#FD8B2F") +
    theme(
      plot.margin = margin(100, 100, 100, 100),
      text = element_text(size=12),
      legend.position = "none"
    )
)
#venn_data <- process_data(Venn(vars2venn))

ggsave(filename = paste0(outdir, "VennDiagram_CvsD_controlSexIntAgeBMI_p", as.character(plim), "_2.pdf"), gv, width = 8, height = 8)

## Venn diagram of BMI related bacteria
vars2venn <- list(
  "BMI" = lindalist$contrastlist2$BMI_alone$resdf %>% dplyr::filter(padj < plim) %>% pull(taxon),
  "BMI\nadj Sex*Age" = lindalist$contrastlist2$`BMI_adj.Age*Sex`$resdf %>% dplyr::filter(padj < plim) %>% pull(taxon),
  "BMI\nadj Sex*Age+MDD" = lindalist$contrastlist2$`BMI_adj.Age*Sex+Depr`$resdf %>% dplyr::filter(padj < plim) %>% pull(taxon)
)
gv <- ggvenn(
  vars2venn, columns = names(vars2venn),
  stroke_size = 0.5,
  stroke_color = C_NS,
  fill_color = c(C_CASE, C_CTRL2, C_CTRL, C_CASE2),show_elements = F
)+
  theme(
    plot.margin = margin(t = 20, r = 20, b = 20, l = 20),  # top, right, bottom, left
    text = element_text(size = 12)  # optional: tweak font size
  )
ggsave(filename = paste0(outdir, "VennDiagram_BMI_controlSexIntAgeCond_p", as.character(plim), ".pdf"), gv, width = 8, height = 8)

(gv <- ggVennDiagram(vars2venn) +
    scale_fill_gradient(low = "white", high = "#FD8B2F") +
    theme(
      plot.margin = margin(100, 100, 100, 100),
      text = element_text(size=12),
      legend.position = "none"
    )
)
#venn_data <- process_data(Venn(vars2venn))

ggsave(filename = paste0(outdir, "VennDiagram_BMI_controlSexIntAgeCond_p", as.character(plim), "_2.pdf"), gv, width = 8, height = 8)

#####

## Venn diagram of MDD vs BMI related bacteria, uncorrected
vars2venn <- list(
  "D vs C" = lindalist$firstContrast$resdf %>% dplyr::filter(padj < plim) %>% pull(taxon),
  #"D vs C\nadj Sex*Age" = lindalist$contrastlist2$`D_vs_C_adj.Age*Sex`$resdf %>% dplyr::filter(padj < plim) %>% pull(taxon),
  "D vs C\nadj Sex*Age+BMI" = lindalist$contrastlist2$`D_vs_C_adj.Age*Sex+BMI`$resdf %>% dplyr::filter(padj <plim) %>% pull(taxon),
  
  "BMI" = lindalist$contrastlist2$BMI_alone$resdf %>% dplyr::filter(padj < plim) %>% pull(taxon),
  #"BMI\nadj Sex*Age" = lindalist$contrastlist2$`BMI_adj.Age*Sex`$resdf %>% dplyr::filter(padj < plim) %>% pull(taxon),
  "BMI\nadj Sex*Age+MDD" = lindalist$contrastlist2$`BMI_adj.Age*Sex+Depr`$resdf %>% dplyr::filter(padj < plim) %>% pull(taxon)
)
gv <- ggvenn(
  vars2venn, columns = names(vars2venn),
  stroke_size = 0.5,
  stroke_color = C_NS,
  fill_color = c(C_CASE, C_CTRL2, C_CTRL, C_CASE2),show_elements = F
)+
  theme(
    plot.margin = margin(t = 20, r = 20, b = 20, l = 20),  # top, right, bottom, left
    text = element_text(size = 12)  # optional: tweak font size
  )
ggsave(filename = paste0(outdir, "VennDiagram_BMICond_alone_and_all_p", as.character(plim), ".pdf"), gv, width = 8, height = 8)

(gv <- ggVennDiagram(vars2venn) +
    scale_fill_gradient(low = "white", high = "#FD8B2F") +
    theme(
      plot.margin = margin(100, 100, 100, 100),
      text = element_text(size=12),
      legend.position = "none"
    )
)
#venn_data <- process_data(Venn(vars2venn))
ggsave(filename = paste0(outdir, "VennDiagram_BMICond_alone_and_all_2_p", as.character(plim), ".pdf"), gv, width = 8, height = 8)

######

## Venn diagram of MDD vs BMI related bacteria, uncorrected
vars2venn <- list(
  #"D vs C" = lindalist$firstContrast$resdf %>% dplyr::filter(padj < 0.05) %>% pull(taxon),
  "D vs C\nadj Sex*Age" = lindalist$contrastlist2$`D_vs_C_adj.Age*Sex`$resdf %>% dplyr::filter(padj < plim) %>% pull(taxon),
  "D vs C\nadj Sex*Age+BMI" = lindalist$contrastlist2$`D_vs_C_adj.Age*Sex+BMI`$resdf %>% dplyr::filter(padj < plim) %>% pull(taxon),
  
  #"BMI" = lindalist$contrastlist2$BMI_alone$resdf %>% dplyr::filter(padj < 0.05) %>% pull(taxon),
  "BMI\nadj Sex*Age" = lindalist$contrastlist2$`BMI_adj.Age*Sex`$resdf %>% dplyr::filter(padj < plim) %>% pull(taxon),
  "BMI\nadj Sex*Age+MDD" = lindalist$contrastlist2$`BMI_adj.Age*Sex+Depr`$resdf %>% dplyr::filter(padj < plim) %>% pull(taxon)
)
gv <- ggvenn(
  vars2venn, columns = names(vars2venn),
  stroke_size = 0.5,
  stroke_color = C_NS,
  fill_color = c(C_CASE, C_CTRL2, C_CTRL, C_CASE2),show_elements = F
)+
  theme(
    plot.margin = margin(t = 20, r = 20, b = 20, l = 20),  # top, right, bottom, left
    text = element_text(size = 12)  # optional: tweak font size
  )
ggsave(filename = paste0(outdir, "VennDiagram_BMICond_adjInt_and_all_p", as.character(plim), ".pdf"), gv, width = 8, height = 8)

(gv <- ggVennDiagram(vars2venn) +
    scale_fill_gradient(low = "white", high = "#FD8B2F") +
    theme(
      plot.margin = margin(100, 100, 100, 100),
      text = element_text(size=12),
      legend.position = "none"
    )
)
#venn_data <- process_data(Venn(vars2venn))
ggsave(filename = paste0(outdir, "VennDiagram_BMICond_adjInt_and_all_p", as.character(plim), "_2.pdf"), gv, width = 8, height = 8)

######


### Plot bars of DEA contrasts
daalist <- list(
  "D_vs_C" = lindalist$firstContrast$resdf,
  #"D_vs_C_adj_SexIAge" = lindalist$contrastlist2$`D_vs_C_adj.Age*Sex`$resdf,
  "D_vs_C_adj_SexInterAge_plus_BMI" = lindalist$contrastlist2$`D_vs_C_adj.Age*Sex+BMI`$resdf,
  
  "BMI" = lindalist$contrastlist2$BMI_alone$resdf,
  #"BMI_adj_SexIAge" = lindalist$contrastlist2$`BMI_adj.Age*Sex`$resdf,
  "BMI_adj_SexInterAge_plus_MDD" = lindalist$contrastlist2$`BMI_adj.Age*Sex+Depr`$resdf
  
)

batplotsdaa <- makeBarplotDAA3_Int(daalist, outdir, plim=0.05, name="corrSexIAgeBMI_b")
batplotsdaa <- makeBarplotDAA3_Int(daalist, outdir, plim=0.01, name="corrSexIAgeBMIp01_b")

daalist <- list(
  "D_vs_C" =lindalist$contrastlist2$`D_vs_C_adj.Age*Sex`$resdf,
  #"D_vs_C_adj_SexIAge" = lindalist$contrastlist2$`D_vs_C_adj.Age*Sex`$resdf,
  "D_vs_C_adj_SexInterAge_plus_BMI" = lindalist$contrastlist2$`D_vs_C_adj.Age*Sex+BMI`$resdf,
  
  "BMI" = lindalist$contrastlist2$`BMI_adj.Age*Sex`$resdf,
  #"BMI_adj_SexIAge" = lindalist$contrastlist2$`BMI_adj.Age*Sex`$resdf,
  "BMI_adj_SexInterAge_plus_MDD" = lindalist$contrastlist2$`BMI_adj.Age*Sex+Depr`$resdf
  
)

batplotsdaa <- makeBarplotDAA3_Int(daalist, outdir, plim=0.05, name="corrSexIAgeBMIAllAdj")
batplotsdaa <- makeBarplotDAA3_Int(daalist, outdir, plim=0.01, name="corrSexIAgeBMIp01AllAdj")



#### Compare with DESeq2

## Condition
vars2venn <- list(
  "D vs C\nDESeq2" = dea2contrasts$firstContrast$resdf %>% dplyr::filter(padj < 0.05) %>% pull(taxon),
  "D vs C\nLinDA" = lindalist$firstContrast$resdf %>% dplyr::filter(padj < 0.05) %>% pull(taxon)
  #"D vs C\nadj Sex*Age" = lindalist$contrastlist2$`D_vs_C_adj.Age*Sex`$resdf %>% dplyr::filter(padj < plim) %>% pull(taxon),
  #"D vs C\nadj Sex*Age+BMI" = lindalist$contrastlist2$`D_vs_C_adj.Age*Sex+BMI`$resdf %>% dplyr::filter(padj < plim) %>% pull(taxon),
  #"BMI" = lindalist$contrastlist2$BMI_alone$resdf %>% dplyr::filter(padj < 0.05) %>% pull(taxon),
  #"BMI\nadj Sex*Age" = lindalist$contrastlist2$`BMI_adj.Age*Sex`$resdf %>% dplyr::filter(padj < plim) %>% pull(taxon),
  #"BMI\nadj Sex*Age+MDD" = lindalist$contrastlist2$`BMI_adj.Age*Sex+Depr`$resdf %>% dplyr::filter(padj < plim) %>% pull(taxon)
  
  )
gv <- ggvenn(
  vars2venn, columns = names(vars2venn),
  stroke_size = 0.5,
  stroke_color = C_NS,
  fill_color = c(C_CASE, C_CTRL2, C_CTRL, C_CASE2),show_elements = F
)+
  theme(
    plot.margin = margin(t = 20, r = 20, b = 20, l = 20),  # top, right, bottom, left
    text = element_text(size = 12)  # optional: tweak font size
  )
ggsave(filename = paste0(outdir, "VennDiagram_CompLindaDeseq_Condition_p", as.character(plim), ".pdf"), gv, width = 8, height = 8)

(gv <- ggVennDiagram(vars2venn) +
    scale_fill_gradient(low = "white", high = "#FD8B2F") +
    theme(
      plot.margin = margin(100, 100, 100, 100),
      text = element_text(size=12),
      legend.position = "none"
    )
)
ggsave(filename = paste0(outdir, "VennDiagram_CompLindaDeseq_Condition_p", as.character(plim), "_2.pdf"), gv, width = 8, height = 8)


## Condition Adj all
vars2venn <- list(
  "D vs C adj.\nDESeq2" = dea2contrasts$contrastlist2$`D_vs_C_adj.Age*Sex+BMI`$resdf %>% dplyr::filter(padj < 0.05) %>% pull(taxon),
  "D vs C adj.\nLinDA" = lindalist$contrastlist2$`D_vs_C_adj.Age*Sex+BMI`$resdf %>% dplyr::filter(padj < 0.05) %>% pull(taxon)
  #"D vs C\nadj Sex*Age" = lindalist$contrastlist2$`D_vs_C_adj.Age*Sex`$resdf %>% dplyr::filter(padj < plim) %>% pull(taxon),
  #"D vs C\nadj Sex*Age+BMI" = lindalist$contrastlist2$`D_vs_C_adj.Age*Sex+BMI`$resdf %>% dplyr::filter(padj < plim) %>% pull(taxon),
  #"BMI" = lindalist$contrastlist2$BMI_alone$resdf %>% dplyr::filter(padj < 0.05) %>% pull(taxon),
  #"BMI\nadj Sex*Age" = lindalist$contrastlist2$`BMI_adj.Age*Sex`$resdf %>% dplyr::filter(padj < plim) %>% pull(taxon),
  #"BMI\nadj Sex*Age+MDD" = lindalist$contrastlist2$`BMI_adj.Age*Sex+Depr`$resdf %>% dplyr::filter(padj < plim) %>% pull(taxon)
  
)
gv <- ggvenn(
  vars2venn, columns = names(vars2venn),
  stroke_size = 0.5,
  stroke_color = C_NS,
  fill_color = c(C_CASE, C_CTRL2, C_CTRL, C_CASE2),show_elements = F
)+
  theme(
    plot.margin = margin(t = 20, r = 20, b = 20, l = 20),  # top, right, bottom, left
    text = element_text(size = 12)  # optional: tweak font size
  )
ggsave(filename = paste0(outdir, "VennDiagram_CompLindaDeseq_ConditionAdjSexIAgeBMI_p", as.character(plim), ".pdf"), gv, width = 8, height = 8)

(gv <- ggVennDiagram(vars2venn) +
    scale_fill_gradient(low = "white", high = "#FD8B2F") +
    theme(
      plot.margin = margin(100, 100, 100, 100),
      text = element_text(size=12),
      legend.position = "none"
    )
)
ggsave(filename = paste0(outdir, "VennDiagram_CompLindaDeseq_ConditionAdjSexIAgeBMI_p", as.character(plim), "_2.pdf"), gv, width = 8, height = 8)

## BMI
vars2venn <- list(
  "BMI\nDESeq2" = dea2contrasts$contrastlist2$BMI_alone$resdf %>% dplyr::filter(padj < 0.05) %>% pull(taxon),
  "BMI\nLinDA" = lindalist$contrastlist2$BMI_alone$resdf %>% dplyr::filter(padj < 0.05) %>% pull(taxon)
  #"D vs C\nadj Sex*Age" = lindalist$contrastlist2$`D_vs_C_adj.Age*Sex`$resdf %>% dplyr::filter(padj < plim) %>% pull(taxon),
  #"D vs C\nadj Sex*Age+BMI" = lindalist$contrastlist2$`D_vs_C_adj.Age*Sex+BMI`$resdf %>% dplyr::filter(padj < plim) %>% pull(taxon),
  #"BMI" = lindalist$contrastlist2$BMI_alone$resdf %>% dplyr::filter(padj < 0.05) %>% pull(taxon),
  #"BMI\nadj Sex*Age" = lindalist$contrastlist2$`BMI_adj.Age*Sex`$resdf %>% dplyr::filter(padj < plim) %>% pull(taxon),
  #"BMI\nadj Sex*Age+MDD" = lindalist$contrastlist2$`BMI_adj.Age*Sex+Depr`$resdf %>% dplyr::filter(padj < plim) %>% pull(taxon)
  
)
gv <- ggvenn(
  vars2venn, columns = names(vars2venn),
  stroke_size = 0.5,
  stroke_color = C_NS,
  fill_color = c(C_CASE, C_CTRL2, C_CTRL, C_CASE2),show_elements = F
)+
  theme(
    plot.margin = margin(t = 20, r = 20, b = 20, l = 20),  # top, right, bottom, left
    text = element_text(size = 12)  # optional: tweak font size
  )
ggsave(filename = paste0(outdir, "VennDiagram_CompLindaDeseq_BMI_p", as.character(plim), ".pdf"), gv, width = 8, height = 8)

(gv <- ggVennDiagram(vars2venn) +
    scale_fill_gradient(low = "white", high = "#FD8B2F") +
    theme(
      plot.margin = margin(100, 100, 100, 100),
      text = element_text(size=12),
      legend.position = "none"
    )
)
ggsave(filename = paste0(outdir, "VennDiagram_CompLindaDeseq_BMI_p", as.character(plim), "_2.pdf"), gv, width = 8, height = 8)


## BMI Adj all
vars2venn <- list(
  "BMI adj.\nDESeq2" = dea2contrasts$contrastlist2$`BMI_adj.Age*Sex+Depr`$resdf %>% dplyr::filter(padj < 0.05) %>% pull(taxon),
  "BMI adj.\nLinDA" = lindalist$contrastlist2$`BMI_adj.Age*Sex+Depr`$resdf %>% dplyr::filter(padj < 0.05) %>% pull(taxon)
  #"D vs C\nadj Sex*Age" = lindalist$contrastlist2$`D_vs_C_adj.Age*Sex`$resdf %>% dplyr::filter(padj < plim) %>% pull(taxon),
  #"D vs C\nadj Sex*Age+BMI" = lindalist$contrastlist2$`D_vs_C_adj.Age*Sex+BMI`$resdf %>% dplyr::filter(padj < plim) %>% pull(taxon),
  #"BMI" = lindalist$contrastlist2$BMI_alone$resdf %>% dplyr::filter(padj < 0.05) %>% pull(taxon),
  #"BMI\nadj Sex*Age" = lindalist$contrastlist2$`BMI_adj.Age*Sex`$resdf %>% dplyr::filter(padj < plim) %>% pull(taxon),
  #"BMI\nadj Sex*Age+MDD" = lindalist$contrastlist2$`BMI_adj.Age*Sex+Depr`$resdf %>% dplyr::filter(padj < plim) %>% pull(taxon)
  
)
gv <- ggvenn(
  vars2venn, columns = names(vars2venn),
  stroke_size = 0.5,
  stroke_color = C_NS,
  fill_color = c(C_CASE, C_CTRL2, C_CTRL, C_CASE2),show_elements = F
)+
  theme(
    plot.margin = margin(t = 20, r = 20, b = 20, l = 20),  # top, right, bottom, left
    text = element_text(size = 12)  # optional: tweak font size
  )
ggsave(filename = paste0(outdir, "VennDiagram_CompLindaDeseq_BMIAdjSexIAgeCond_p", as.character(plim), ".pdf"), gv, width = 8, height = 8)

(gv <- ggVennDiagram(vars2venn) +
    scale_fill_gradient(low = "white", high = "#FD8B2F") +
    theme(
      plot.margin = margin(100, 100, 100, 100),
      text = element_text(size=12),
      legend.position = "none"
    )
)
ggsave(filename = paste0(outdir, "VennDiagram_CompLindaDeseq_BMIAdjSexIAgeCond_p", as.character(plim), "_2.pdf"), gv, width = 8, height = 8)


## Cond alone and Adj
vars2venn <- list(
  "D vs C\nDESeq2" = dea2contrasts$firstContrast$resdf %>% dplyr::filter(padj < 0.05) %>% pull(taxon),
  "D vs C\nLinDA" = lindalist$firstContrast$resdf %>% dplyr::filter(padj < 0.05) %>% pull(taxon),
  "D vs C adj.\nDESeq2" = dea2contrasts$contrastlist2$`D_vs_C_adj.Age*Sex+BMI`$resdf %>% dplyr::filter(padj < 0.05) %>% pull(taxon),
  "D vs C adj.\nLinDA" = lindalist$contrastlist2$`D_vs_C_adj.Age*Sex+BMI`$resdf %>% dplyr::filter(padj < 0.05) %>% pull(taxon)
)
gv <- ggvenn(
  vars2venn, columns = names(vars2venn),
  stroke_size = 0.5,
  stroke_color = C_NS,
  fill_color = c(C_CASE, C_CTRL2, C_CTRL, C_CASE2),show_elements = F
)+
  theme(
    plot.margin = margin(t = 20, r = 20, b = 20, l = 20),  # top, right, bottom, left
    text = element_text(size = 12)  # optional: tweak font size
  )
ggsave(filename = paste0(outdir, "VennDiagram_CompLindaDeseq_CondAloneAndAdj_p", as.character(plim), ".pdf"), gv, width = 8, height = 8)

(gv <- ggVennDiagram(vars2venn) +
    scale_fill_gradient(low = "white", high = "#FD8B2F") +
    theme(
      plot.margin = margin(100, 100, 100, 100),
      text = element_text(size=12),
      legend.position = "none"
    )
)
ggsave(filename = paste0(outdir, "VennDiagram_CompLindaDeseq_CondAloneAndAdj_p", as.character(plim), "_2.pdf"), gv, width = 8, height = 8)

### bars with all

allist <- lindalist$contrastlist2
allist[["D_vs_C_alone"]] <- lindalist$firstContrast

allist2 <- dea2contrasts$contrastlist2
allist2[["D_vs_C_alone"]] <- dea2contrasts$firstContrast

for(nn in names(allist)){
  fname <- gsub("\\+", "", nn) %>% 
    gsub("\\*", "I", .)
  write_tsv(allist[[nn]]$resdf, paste0(outdir, "LinDA_" , fname, ".tsv"))
}

resall <- data.frame()
for(nn in names(allist)){
  aux <- allist[[nn]]$resdf %>% 
    dplyr::mutate(Contrast = nn, Software = "LinDA")
  
  aux2 <- allist2[[nn]]$resdf %>% 
    dplyr::mutate(Contrast = nn, Software = "DESeq2")
  
  resall <- resall %>% 
    bind_rows(aux2) %>% 
    bind_rows(aux %>% dplyr::select(all_of(names(aux2))))
  
}

PLIM = 0.05
resall <- resall %>% 
  dplyr::mutate(padj = ifelse(is.na(padj), 1, padj)) %>% 
  dplyr::mutate(Sig = ifelse(padj <PLIM, 
                             ifelse(log2FoldChange < 0 , "Down", "Up"),
                             "NS"))
write_tsv(resall, paste0(outdir, "all_contrasts_long_LinDA_and_DESeq2.tsv"))


#
plotdf <- resall %>% 
  dplyr::mutate(taxon = gsub("_", " ", taxon), 
                taxon = gsub("[\\[\\]]", "", taxon),
                Sig = factor(Sig, levels=c("Down", "Up", "NS"))
  ) %>% 
  filter(Contrast %in% c("D_vs_C_alone", "D_vs_C_adj.Age*Sex+BMI")) 

tax2plot <- plotdf %>% filter(padj<0.01) %>% pull(taxon) %>% unique

taxorder <- plotdf %>% 
  filter(Contrast == "D_vs_C_alone" & Software == "DESeq2") %>% 
  filter(taxon %in% tax2plot) %>% 
  arrange(log2FoldChangeShrink) %>% 
  pull(taxon)

taxorder <- plotdf %>% 
  filter(taxon %in% tax2plot) %>% 
  tidyr::unite(col = "tmp", Software, Contrast, sep="__") %>%  
  dplyr::select(log2FoldChangeShrink, padj, taxon, tmp) %>% 
  gather("var", "val", log2FoldChangeShrink, padj) %>% 
  tidyr::unite(col = "tmp2", tmp, var, sep="__") %>% 
  spread(tmp2, val) %>% 
  dplyr::mutate(
    siggroup = ifelse(
      `DESeq2__D_vs_C_alone__padj` <= 0.05 &
      `DESeq2__D_vs_C_adj.Age*Sex+BMI__padj` <= 0.05 &
        `LinDA__D_vs_C_alone__padj` <= 0.05 &
        `LinDA__D_vs_C_adj.Age*Sex+BMI__padj` <= 0.05, 1, 
            ifelse(
              `DESeq2__D_vs_C_alone__padj` <= 0.05 &
                `DESeq2__D_vs_C_adj.Age*Sex+BMI__padj` <= 0.05 &
                `LinDA__D_vs_C_alone__padj` <= 0.05, 2, 
              ifelse(
                `DESeq2__D_vs_C_alone__padj` <= 0.05 &
                  `DESeq2__D_vs_C_adj.Age*Sex+BMI__padj` <= 0.05, 3, 
                ifelse(
                  `DESeq2__D_vs_C_alone__padj` <= 0.05 & 
                    `LinDA__D_vs_C_alone__padj` <= 0.05, 4, 
                  ifelse(`DESeq2__D_vs_C_alone__padj` <= 0.05, 5, 
                         ifelse(`LinDA__D_vs_C_alone__padj` <= 0.05, 6, 7))))))
            ) %>% 
  dplyr::mutate(siggroup = factor(siggroup, ordered = TRUE)) %>% 
  group_by(siggroup) %>% 
  dplyr::arrange(desc(siggroup), desc(DESeq2__D_vs_C_alone__log2FoldChangeShrink))


plotdf2 <- plotdf %>% 
  filter(taxon %in% tax2plot) %>% 
  dplyr::mutate(
    taxon = factor(taxon, levels =taxorder$taxon),
    Contrast = gsub("_", " ", Contrast),
    Contrast = gsub(" alone", "", Contrast),
    Contrast = gsub("adj.", "adj. ", Contrast),
    Contrast = factor(Contrast, levels= c("D vs C", "D vs C adj. Age*Sex+BMI")),
    Software = factor(Software, levels = c("DESeq2", "LinDA"))
  )

hlines <- taxorder %>% 
  ungroup() %>% 
  dplyr::mutate(pos = 1:n()) %>% 
  group_by(siggroup) %>% 
  dplyr::summarise(maxpos = max(pos) + 0.5) %>% 
  tail(nrow(.)-1)


(gb <- ggplot(plotdf2, aes(x=taxon, y=log2FoldChangeShrink, fill=Sig))+
    facet_grid(~ Contrast + Software) +
    geom_col()+
    coord_flip()+
    theme_minimal() +
    geom_hline(yintercept = 0, linetype=2, col=C_NS) +
    scale_fill_manual(values=c(C_CTRL, C_CASE, C_NS)) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    theme(strip.text.y = element_text(size = 10,
                                      colour = "black", angle = 0, face = "italic")) +
    theme(axis.text.x = element_text(size = 8,
                                     colour = "black", angle = 0,
                                     face = "plain"))+
    theme(axis.text.y = element_text(size = 8,
                                     colour = "black", angle = 0,
                                     vjust = 0.5,
                                     face = "italic")) +
    theme(legend.position = 'none')+
    #theme(axis.text.y = element_blank())+
    xlab("Species")+
    ylab("LFC") +
    geom_vline(xintercept = hlines$maxpos, linetype=2)
    #ylim(c(-9, 9))+
   # thin_barplot_lines
  
)
ggsave(filename = paste0(outdir, "CompareDESeqLinda_Condition_Barplot_show01plot05.pdf"), width = 10, height = 14)

## Check the ones with BMI and MDD effects

#
plotdf <- resall %>% 
  dplyr::mutate(taxon = gsub("_", " ", taxon), 
                taxon = gsub("[\\[\\]]", "", taxon, perl=TRUE),
                Sig = factor(Sig, levels=c("Down", "Up", "NS"))
  ) %>% 
  filter(Contrast %in% c("D_vs_C_alone", "D_vs_C_adj.Age*Sex+BMI", "BMI_alone", "BMI_adj.Age*Sex+Depr")) 

tax2plot <- plotdf %>% filter(padj<0.05 & Software == "DESeq2") %>% 
  group_by(taxon) %>% 
  dplyr::summarise(nsig = sum(padj <= 0.05)) %>% 
  filter(nsig == 4) %>% 
  pull(taxon) # n = 16

taxorder <- plotdf %>% 
  filter(Contrast == "D_vs_C_alone" & Software == "DESeq2") %>% 
  filter(taxon %in% tax2plot) %>% 
  arrange(log2FoldChangeShrink)

plotdf2 <- plotdf %>% 
  filter(taxon %in% tax2plot) %>% 
  dplyr::mutate(
    taxon = factor(taxon, levels =taxorder$taxon),
    Contrast = gsub("_", " ", Contrast),
    Contrast = gsub(" alone", "", Contrast),
    Contrast = gsub("adj.", "adj. ", Contrast),
    Contrast = factor(Contrast, levels= c("D vs C", "D vs C adj. Age*Sex+BMI",  "BMI", "BMI adj. Age*Sex+Depr")),
    Software = factor(Software, levels = c("DESeq2", "LinDA"))
  )

(gb <- ggplot(plotdf2, aes(x=taxon, y=log2FoldChangeShrink, fill=Sig))+
    facet_grid(~ Contrast + Software) +
    geom_col()+
    coord_flip()+
    theme_minimal() +
    geom_hline(yintercept = 0, linetype=2, col=C_NS) +
    scale_fill_manual(values=c(C_CTRL, C_CASE, C_NS)) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    theme(strip.text.y = element_text(size = 10,
                                      colour = "black", angle = 0, face = "italic")) +
    theme(axis.text.x = element_text(size = 8,
                                     colour = "black", angle = 0,
                                     face = "plain"))+
    theme(axis.text.y = element_text(size = 8,
                                     colour = "black", angle = 0,
                                     vjust = 0.5,
                                     face = "italic")) +
    theme(legend.position = 'none')+
    #theme(axis.text.y = element_blank())+
    xlab("Species")+
    ylab("LFC") #+
   # geom_vline(xintercept = hlines$maxpos, linetype=2)
  #ylim(c(-9, 9))+
  # thin_barplot_lines
  
)
ggsave(filename = paste0(outdir, "CompareDESeqLinda_Condition_Barplot_16SigFromMediation.pdf"), width = 14, height = 4)

####### study F praus

linda.obj_adjAll$output$Age_log %>% rownames_to_column("taxon") %>% filter(grepl("praus", taxon))
linda.obj_adjAll$output$SexFemale %>% rownames_to_column("taxon") %>% filter(grepl("praus", taxon))
linda.obj_adjAll$output$`SexFemale:Age_log` %>% rownames_to_column("taxon") %>% filter(grepl("praus", taxon))


linda.obj_adjAll_nolog$output$Age %>% rownames_to_column("taxon") %>% filter(grepl("praus", taxon))
linda.obj_adjAll_nolog$output$SexFemale %>% rownames_to_column("taxon") %>% filter(grepl("praus", taxon))
linda.obj_adjAll_nolog$output$`SexFemale:Age` %>% rownames_to_column("taxon") %>% filter(grepl("praus", taxon))

linda.obj_age_nolog$output$Age %>% rownames_to_column("taxon") %>% filter(grepl("praus", taxon))
linda.obj_age_log$output$Age %>% rownames_to_column("taxon") %>% filter(grepl("praus", taxon))

# Age Not sig in any
# Sex sig when adjusting with no log


samples <- sample_data(phobj)$sampleID[sample_data(phobj)$Sex == "Female"]
phobj_women <- phyloseq::prune_samples(samples, phobj)

otu.tab_women <- otu_table(phobj_women) %>% data.frame
meta_women <- sample_data(phobj_women) %>% data.frame

linda.obj_age_nolog_wmn <- linda(otu.tab_women, meta_women, formula = '~Age', 
                             alpha = 0.05,
                             prev.cut = 0.05, 
                             lib.cut = 1000, 
                             winsor.quan = 0.97)

linda.obj_age_log_wmn <- linda(otu.tab_women, meta_women, formula = '~Age_log', 
                           alpha = 0.05,
                           prev.cut = 0.05, 
                           lib.cut = 1000, 
                           winsor.quan = 0.97)

linda.obj_age_nolog_wmn$output$Age %>% rownames_to_column("taxon") %>% filter(grepl("praus", taxon)) # 0.63
linda.obj_age_log_wmn$output$Age %>% rownames_to_column("taxon") %>% filter(grepl("praus", taxon)) # 0.56
