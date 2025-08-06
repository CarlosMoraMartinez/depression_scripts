library(tidyverse)
library(phyloseq)
library(vegan)
library(ggsci)
library(scico)
library(ggnewscale)
library(patchwork)

options(ggplot2.discrete.fill = c("#A1C6EA","#FD8B2F", "#00AA5A",
                                  "#8E7BFF","#00D1EE", "#00E6BB",
                                  "#F9F871", "#F45680", "#A5ABBD",
                                  "#B60E50"))
options(ggplot2.discrete.colour = c("#A1C6EA","#FD8B2F","#00AA5A",
                                    "#8E7BFF","#00D1EE", "#00E6BB",
                                    "#F9F871", "#F45680", "#A5ABBD",
                                    "#B60E50"))


load("/home/carlos/Escritorio/202311_DEPRESION/ReviewJune2025/Results_rstudio/results2/phyloseq/pre_phyloseq_filt5.RData")

outdir <- paste0("/home/carlos/Escritorio/202311_DEPRESION/ReviewJune2025/Results_rstudio/results2/why_alpha_diversity/")
if(!dir.exists(outdir)) dir.create(outdir)

s_meta <- sample_data(pre_phyloseq_filt) %>% data.frame

otus <- pre_phyloseq_filt %>% otu_table() %>% data.frame %>%
  mutate_all( ~ 100*.x/sum(.x)) %>%
  rownames_to_column("taxon") %>%
  gather("sample", "prop", -taxon) %>%
  group_by(sample) %>%
  dplyr::mutate(order = rank(-prop, ties.method = "first")) %>%
  dplyr::mutate(prop_norm = prop/max(prop)) %>%
  dplyr::mutate(sample=gsub("X", "", sample)) %>%
  merge(s_meta, by.x="sample", by.y="sampleID", all.x=TRUE) %>%
  filter(prop>0) %>%
  group_by(sample)

divcalc <- estimate_richness(
  pre_phyloseq_filt,
  measures = c("Observed", "Chao1", "ACE", "Shannon", "Simpson", "InvSimpson", "Fisher")) %>%
  rownames_to_column("sample") %>%
  dplyr::mutate(sample=gsub("X", "", sample))

s_meta <- s_meta %>% merge(divcalc, by.x="sampleID", by.y="sample")

# test fit
B = 0.5

aux <- otus %>% filter(sample=="1")
nls(prop_norm ~ exp(- (a * order)^B),
    data = aux,
    start = list(a = 0.1))


otus_exp <-  otus %>%
  group_by(sample) %>%
  group_modify( .f=~ .x %>%
                  dplyr::mutate(fit_exp = summary(nls(prop_norm ~ exp(- (a * order)^B),
                     data = .x,
                     start = list(a = 0.1)))$coefficients[1, 1])
  ) %>%
  dplyr::mutate(predicted = exp(- (fit_exp*order)^B))

aux <- otus_exp %>% dplyr::select(sample, taxon, Condition, fit_exp, order, predicted)
write_tsv(aux, file = paste0(outdir, "exponential_predicted_by_sample.tsv"))
aux <- aux %>% dplyr::select(-taxon)
exp_sum <- otus_exp %>%
  group_by(sample) %>%
  dplyr::summarise(
            fit_exp = unique(fit_exp),
            Cond=unique(Condition)
            ) %>%
  group_by(Cond) %>%
  dplyr::summarise(
    mean=mean(fit_exp),
    median=median(fit_exp),
    q05=quantile(fit_exp, 0.05),
    q25=quantile(fit_exp, 0.25),
    q75=quantile(fit_exp, 0.75),
    q95=quantile(fit_exp, 0.95),
    sd=sd(fit_exp),
  ) %>%
  gather("measure", "value", -Cond)

write_tsv(exp_sum, file=paste0(outdir, "exponents_by_group.tsv"))

levs <- 1:max(otus_exp$order)

all_preds <- exp_sum %>%
  filter(measure != "sd") %>%
  group_by(Cond, measure) %>%
  group_modify(.f =~ data.frame(
                            order=levs,
                            prop=exp(-.x$value*levs)) %>%
                 dplyr::mutate(prop_norm = prop/max(prop)))

preds2use <- all_preds %>%  
  dplyr::filter(measure %in% c("q25", "q75", "median")) %>%
  dplyr::select(-prop) %>%
  #pivot_wider(names_from = measure, values_from = prop_norm) %>%
  spread(measure, prop_norm) %>%
  dplyr::select(Cond, order, q25, median, q75)


max_order <- max(otus_exp$order)
a_control=exp_sum %>% dplyr::filter(Cond == "Control" & measure == "median") %>% pull(value)%>% round(3)
a_depr=exp_sum %>% dplyr::filter(Cond == "Depression" & measure == "median") %>% pull(value) %>% round(3)


(ggp <- ggplot(otus_exp, aes(x=order, y=prop_norm, group=sample, col=Condition)) +
    geom_ribbon(data = preds2use,
                aes(ymin = q25, ymax = q75, x=order,
                    fill = Cond),
                alpha = 0.5,
                inherit.aes = FALSE) +
  geom_line(size=0.5, linetype=3) +
  geom_line(data=preds2use,
            aes(x=order, y=median, group=Cond, col=Cond),
            linetype=1, size=1.2) +
  #geom_point() +
  theme_bw() +
    scale_x_log10() +
    annotation_logticks(sides = "b")
 +
  ylab("Proportion normalized to 1") +
  xlab("Order")+
  labs(fill = "Condition", color = "Condition") +
    geom_segment(aes(x = max_order * 0.3, xend = max_order * 0.45,  # adjust coords
                     y = 0.9, yend = 0.9),
                 color = "#A1C6EA", size = 1, linetype = 1, inherit.aes = FALSE) +
    annotate("text", x = max_order * 0.47, y = 0.9, label = bquote(italic("a=") * .(a_control)), hjust = 0) +

    # Add orange dashed line segment for Treatment
    geom_segment(aes(x = max_order * 0.3, xend = max_order * 0.45,
                     y = 0.85, yend = 0.85),
                 color = "#FD8B2F", size = 1, linetype = 1, inherit.aes = FALSE) +
    annotate("text", x = max_order * 0.47, y = 0.85, label = bquote(italic("a=") * .(a_depr)), hjust = 0)
)
ggsave(filename = paste0(outdir, "exp_curve2.pdf"), ggp, width = 6, height = 4)

## make simulations

getDistr <- function(num_species, a, b=B){
  dd <- exp(-(a*1:num_species)^b)
  dd/sum(dd)
}

min_exp <- 0 # min(otus_exp$fit_exp)
max_exp <- max(otus_exp$fit_exp)
min_num  <- 11 # otus_exp %>% group_by(sample) %>% group_map( ~ max(.x$order)) %>% unlist %>% min
max_num  <- otus_exp %>% group_by(sample) %>% group_map( ~ max(.x$order)) %>% unlist %>% max

orders <- (min_num - 10):(max_num + 10)
exponents <- seq(round(min_exp*0.8, 2), round(max_exp*1.2, 2), by=0.01)

simdf <- expand.grid(orders, exponents)
names(simdf) <- c("max_order", "exponent")

simdf <- simdf %>%
  dplyr::mutate(Shannon = map2_vec(max_order, exponent,
                               \(x, y) vegan::diversity(getDistr(x, y), index="shannon")))

(gsim <- ggplot(simdf, aes(x =exponent, y = max_order, z = Shannon)) +
  geom_contour_filled()+
  theme_minimal() +
    #scale_fill_uchicago() +
  ylab("Number of species")
  #geom_contour_filled()
)


exp_by_sample <- otus_exp %>%
  group_by(sample) %>%
  dplyr::summarise(
    fit_exp = unique(fit_exp),
    Cond=unique(Condition)
  )

pointstab <- s_meta %>% filter(Tanda == 1) %>%
  merge(exp_by_sample, by.x="sampleID", by.y="sample", all.x=T, all.y=F)

(gsim <- ggplot(simdf, aes(x =exponent, y = max_order, fill = Shannon)) +
    geom_tile()+

    geom_contour(aes(z = Shannon), color = "black", size = 0.3) +
    scale_fill_scico(palette = "oslo") +
    ggnewscale::new_scale_fill() +
    geom_point(data= pointstab, aes(x=fit_exp, y=Observed, fill=Cond),
                                    col="black",
                                    shape = 21,
                                    stroke = 0.8) +

    theme_minimal() +
    stat_ellipse(data = pointstab, aes(x = fit_exp,
                                       y = Observed,
                                       group = Cond, color = Cond),
                 type = "norm",
                 linetype = "dashed", size = 0.8) +
    ylab("Number of species") +
    xlab(expression("value of " * italic(a))) +
    labs(fill = "Condition", color = "Condition")
  #geom_contour_filled()
)

ggsave(filename = paste0(outdir, "shannon_map.pdf"), gsim, width = 6, height = 4)


#gsim <- gsim +
#  theme(plot.margin = margin(0, 100, 0, 0))

library(cowplot)
cowplot::plot_grid(plotlist=list(ggp, gsim), ncol=1,
                   align = "v", rel_widths = c(1, 1))

panelplot <- ggp + gsim +
  plot_layout(widths = c(2, 1.2)) + # 2, 1.4
  plot_annotation(tag_levels = 'a')
panelplot

ggsave(filename = paste0(outdir, "exponential_panel_min.pdf"), panelplot, width = 12, height = 4) # w=10

# finally, get most abundant taxa:

ranked_species <- otus_exp %>% group_by(taxon) %>%
  summarise(median_order = median(order),
            mean_order = mean(order)) %>%
  arrange(median_order)

# top10:
#1 Faecalibacterium_prausnitzii            1       2.88
#2 Bacteroides_uniformis                   6       9.40
#3 Phocaeicola_vulgatus                   11      17.0
#4 Roseburia_intestinalis                 12      13.3
#5 Dorea_longicatena                      15      17.8
#6 Blautia_massiliensis                   16      19.7
#7 Collinsella_aerofaciens                17      34.7
#8 Vescimonas_coprocola                   17      19.6
#9 Anaerobutyricum_hallii                 18      22.6
#10 Coprococcus_comes                      19      22.8
