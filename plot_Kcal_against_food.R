
load("/home/carlos/Documentos/CORALS/results_rstudio/results_Agosto25_1//foodPCA/phyloseq_list_foodPCA_withNMF_withInc.RData")

s_meta <- sample_data(all_phyloseq$remove_tanda2) %>% data.frame

food_names <- names(s_meta)[grep("_clr", names(s_meta))]
food_names <- gsub("_clr", "", food_names)

slong <- s_meta %>% gather("food_type", "amount", all_of(food_names))
g1 <- ggplot(slong, aes(x=energy_kcal, y=amount)) + 
  facet_wrap(. ~ food_type, scales="free") + 
  geom_point(size=0.2, alpha=0.2, col="tomato") + 
  geom_smooth() + 
  theme_classic() 

ggsave("/home/carlos/Documentos/CORALS/results_rstudio/results_Agosto25_1/energy_vs_food_regr_gam.pdf", g1, width = 12, height = 8)

g2 <-   ggplot(slong, aes(x=energy_kcal, y=amount)) + 
  facet_wrap(. ~ food_type, scales="free") + 
  geom_point(size=0.2, alpha=0.2, col="tomato") + 
  geom_smooth(method="lm") +
  ggpmisc::stat_poly_eq() +
  theme_classic() 

ggsave("/home/carlos/Documentos/CORALS/results_rstudio/results_Agosto25_1/energy_vs_food_regr_lm.pdf", g2, width = 12, height = 8)
