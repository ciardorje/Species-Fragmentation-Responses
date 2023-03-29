rm(list=ls()); gc()

library(pacman)
p_load(tidyverse, cowplot, ggpubr)

setwd()

load('AF_DBs_Species_Habitat_Data.RData')

#####Calculate Habitat Preferences for each observed species#####

#Calculate matrix captures and standardize by no. of traps
matrix <- matrix[matrix$Fragment != 'HFA',]
matrix_abun <- matrix %>% 
  select(Binomial) %>%
  group_by(Binomial) %>% 
  mutate(Matrix = n()) %>%
  unique() %>%
  mutate(Matrix = Matrix/(23*3)) #Standardise to abundance/trap - 23 fragments, 3 matrix traps for each

#Calculate no. of traps in all fragments
trap_no <- dbs %>% 
  select(Fragment, Core_Traps, Edge_Traps) %>%
  unique() %>%
  mutate(Core_Traps = sum(Core_Traps)) %>%
  mutate(Edge_Traps = sum(Edge_Traps)) %>%
  select(-Fragment) %>%
  unique()

#Calculate centre and edge captures and standardize by number of traps
frag_abun <- dbs %>% 
  select(Habitat, Binomial, Count) %>%
  group_by(Habitat, Binomial) %>%
  mutate(Count = sum(Count)) %>%
  ungroup() %>%
  unique() %>%
  pivot_wider(names_from = Habitat,
              values_from = Count) %>%
  unique() %>%
  mutate(Fragment = (Core + Edge)/(sum(trap_no))) %>%
  mutate(Core = Core/trap_no$Core_Traps) %>%
  mutate(Edge = Edge/trap_no$Edge_Traps)

colnames(frag_abun)[2] <- 'Centre'

#Merge matrix and forest counts
habitat_abun <- merge(frag_abun, matrix_abun, all = T)
habitat_abun[is.na(habitat_abun)] <- 0
habitat_abun$Total <- habitat_abun$Centre + habitat_abun$Edge + habitat_abun$Matrix

#Calculate habitat preference
habitat_abun <- 
  pivot_longer(habitat_abun,
               cols = c(Centre, Edge, Matrix),
               names_to = 'Habitat',
               values_to = 'Std_Abun') %>%
  mutate(Proportional_Abun = Std_Abun/Total) %>% 
  select(-Std_Abun) %>%
  pivot_wider(names_from = Habitat,
              values_from = Proportional_Abun) %>%
  mutate(FHP = ((Edge - Matrix) + (Centre - Edge) + 1) / 2) %>%
  mutate(FE = ifelse(Matrix == 0, 1, 0)) %>%
  pivot_longer(cols = c(Centre, Edge, Matrix),
               names_to = 'Habitat',
               values_to = 'Proportional_Abun')

#Save
FHP_table <- habitat_abun %>% select(Binomial, FHP, FE) %>% unique()
write.csv(FHP_table, 'AF_DBs_Fragment-Related_Habitat_Preferences.csv')

#Plot

#Split into individual species
sp_abun <- split(habitat_abun, f = habitat_abun$Binomial)

#Order by FHP index
FHPs <- unique(data.frame(Binomial = habitat_abun$Binomial, FHP = habitat_abun$FHP))
FHPs <- FHPs[order(FHPs$FHP, decreasing = T),]
sp_abun <- sp_abun[FHPs$Binomial]

abun_plots <- list()

#Plot for each species
for(i in 1:length(sp_abun)){
  abun_plots[[i]] <- 
    ggplot(sp_abun[[i]], aes(x = factor(Habitat, level = c('Matrix', 'Edge', 'Centre')), 
                             y = Proportional_Abun, colour = Habitat)) +
    geom_path(aes(group = Binomial), color = 'black', size = 0.8, lty = 1) +
    geom_point(size = 2.75) +
    theme_bw() +
    scale_y_continuous(limits = c(0, 1)) +
    scale_x_discrete() +
    labs(title = paste0(sp_abun[[i]]$Binomial), 
         subtitle = paste0('FSI Score = ', round(sp_abun[[i]]$FHP, 3)),
         y = 'Proportional Abundance') +
    theme(panel.grid = element_blank(),
          axis.text = element_text(face = 'bold', color = 'black'),
          legend.position = 'none',
          plot.title = element_text(size = 10, face = 'bold.italic', vjust = 3),
          plot.subtitle = element_text(size = 9, face = 'bold', vjust = 3),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = 9, face = 'bold', vjust = 3),
          axis.ticks = element_line(size = 0.8),
          axis.line = element_line(size = 0.7),
          plot.margin = unit(c(.5,.5,.5,.5), 'cm'),
          panel.border = element_blank())
  
  if(sum(sp_abun[[i]][["FE"]]) > 0){
    abun_plots[[i]] <- abun_plots[[i]] + 
      geom_text(label = '*', x = Inf, y = Inf,
                hjust = 1, vjust = 1, colour = 'black', size = 12) 
  }
  
  if(i != 1 & ((i-1)/3) %% 1 != 0){
      abun_plots[[i]]$labels$y <- NULL
  }
  
}

#Split plots into groups of 12
sps <- 1:86
sps_chunks <- split(sps, ceiling(seq_along(sps) / 12))
sps_grouped <- list()
for(i in 1:length(sps_chunks)){
  sps_grouped[[i]] <- abun_plots[sps_chunks[[i]]]}

#Plot in groups on 1 page
for(i in 1:(length(sps_grouped)-1)){
  multi_sps <- plot_grid(sps_grouped[[i]][[1]], sps_grouped[[i]][[2]],
                         sps_grouped[[i]][[3]], sps_grouped[[i]][[4]], 
                         sps_grouped[[i]][[5]], sps_grouped[[i]][[6]],
                         sps_grouped[[i]][[7]], sps_grouped[[i]][[8]],
                         sps_grouped[[i]][[9]], sps_grouped[[i]][[10]], 
                         sps_grouped[[i]][[11]], sps_grouped[[i]][[12]],
                         ncol = 3, nrow = 4, align = 'hv')
  ggsave(filename = paste('./Figures/Habitat_Preference_Plots/Species_Plots_', i, '.png'),
         height = 297, width = 210, dpi = 500, device = 'png', units = 'mm', limitsize = F)
}

multi_sps <- plot_grid(sps_grouped[[8]][[1]], sps_grouped[[8]][[2]],
                       ncol = 3, nrow = 4, align = 'hv')
ggsave(filename = paste('./Figures/Habitat_Preference_Plots/Species_Plots_', 15, '.png'),
       height = 297, width = 210, dpi = 500, device = 'png', units = 'mm', limitsize = F)

#Plot for all species together
ggplot(habitat_abun, aes(x = Habitat, y = Proportional_Abun, colour = Binomial)) +
  geom_path(aes(group = Binomial), size = 0.1, lty = 1) +
  geom_point(size = 2.75) +
  theme_bw() +
  labs(y = '% of Total Capture Prevalance (per Trap)') +
  theme(panel.grid = element_blank(),
        axis.text = element_text(face = 'bold', color = 'black'),
        legend.position = 'none',
        plot.title = element_text(face = 'bold.italic', vjust = 3),
        plot.subtitle = element_text(face = 'bold', vjust = 3),
        axis.title.x = element_blank(),
        axis.title.y = element_text(face = 'bold', vjust = 3),
        axis.ticks = element_line(size = 0.8),
        axis.line = element_line(size = 0.7),
        plot.margin = unit(c(.5,.5,.5,.5), 'cm'),
        panel.border = element_blank())

#Calculate forest-dependancy as per Palmeirim et al. (2020)
FD_abun <- merge(frag_abun, matrix_abun, all = T) %>% 
  select(Binomial, Fragment, Matrix)
FD_abun[is.na(FD_abun)] <- 0
FD_abun$Fragment <- FD_abun$Fragment + 0.01
FD_abun$Matrix <- FD_abun$Matrix + 0.01
FD_abun$FD <- log10(FD_abun$Fragment/FD_abun$Matrix)

metrics <- habitat_abun %>% 
  select(Binomial, FHP) %>% 
  merge(select(FD_abun, Binomial, FD)) %>% 
  unique()

metric_lm <- lm(FHP ~ FD, data = metrics)

ggplot(data = metrics, aes(x = FD, y = FHP, colour = Binomial)) +
  geom_abline(slope = coef(metric_lm)[['FD']],
              intercept = coef(metric_lm)[['(Intercept)']],
              size = 1, lty = 2) +
  geom_point(size = 3, pch = 16) +
  theme_bw() +
  labs(x = 'Forest Dependency Index', 
       y = 'Forest Specificity Index') +
  scale_x_continuous(limits = c(-3, 3), labels = seq(-3, 3, 1), breaks =  seq(-3, 3, 1)) +
  theme(#panel.grid = element_blank(),
        axis.text = element_text(face = 'bold', color = 'black'),
        legend.position = 'none',
        axis.title.x = element_text(face = 'bold', vjust = 0),
        axis.title.y = element_text(face = 'bold', vjust = 3),
        axis.ticks = element_line(size = 0.8),
        axis.line = element_line(size = 0.7),
        plot.margin = unit(c(.5,.5,.5,.5), 'cm'),
        panel.border = element_blank())
