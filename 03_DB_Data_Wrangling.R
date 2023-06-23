rm(list=ls()); gc()

library(tidyverse)

setwd('C:/Users/ciarn/Desktop/PhD/Dung_Beetles/Raw_Data/')

dbs <- read.csv('scarabs_clean.csv')

#Translate and Correct fragment/habitat codes
colnames(dbs) <- c('Fragment', 'Habitat', 'Trap_No.', 'Date', 
                   'Trap_Type', 'Genus', 'Species', 'Binomial', 'Notes')
dbs$Habitat[dbs$Habitat == 'Borda'] <- 'Edge'
dbs$Habitat[dbs$Habitat == 'Pasto'] <- 'Matrix'
dbs$Habitat[dbs$Habitat == 'Centro'] <- 'Core'
dbs$Fragment[dbs$Fragment == 'Hotel Floresta Amazonica'] <- 'HFA'
dbs$Fragment <- sub('Fragmento ', '', dbs$Fragment) 
dbs$Fragment[dbs$Fragment == '153a'] <- '157' #153 and 153a are actually separate, re-code

#Drop FIT traps as only one per site, can't model detection
dbs <- filter(dbs, Trap_Type == 'Pitfall') %>% select(-Trap_Type)

#Seperate fragment codes into fragment and trap array
dbs <- dbs %>% 
  select(-Notes) %>% 
  mutate(Trap_Plot = Fragment) %>%
  mutate(Fragment = sub("(\\d)[^0-9]+$", "\\1", Fragment)) #Remove trap codes from fragment code

#Seperate out records without trap number
woTrapCode <- filter(dbs, Trap_No. == '' | Trap_No. == 'FIT')
dbs <- filter(dbs, Trap_No. != '' & Trap_No. != 'FIT')

#Seperate out matrix samples
matrix <- filter(dbs, Habitat == 'Matrix')
dbs <- filter(dbs, Habitat != 'Matrix')

#Unite trap plot and number to give unique trap ID
dbs <- dbs %>% unite(Trap, c(Trap_Plot, Trap_No.), remove = T) 

#Count traps/habitat and remove areas with only 1 trap
dbs <- dbs %>% 
  group_by(Fragment) %>%
  mutate(Fragment_Samples = n_distinct(Trap)) %>%
  filter(Fragment_Samples > 1)

#Count traps / habitat within each fragment
Core_Traps <- dbs %>% select(Fragment, Habitat, Trap) %>%
  filter(Habitat == 'Core') %>% select(-Habitat) %>% group_by(Fragment) %>%
  mutate(Core_Traps = length(unique(Trap))) %>% select(-Trap) %>% unique()
Edge_Traps <- dbs %>% select(Fragment, Habitat, Trap) %>%
  filter(Habitat == 'Edge') %>% select(-Habitat) %>% group_by(Fragment) %>%
  mutate(Edge_Traps = length(unique(Trap))) %>% select(-Trap) %>% unique()

dbs <- merge(dbs, Core_Traps) %>% merge(Edge_Traps)

#Summarise to count data
dbs <- dbs %>% 
  group_by(Fragment, Date, Habitat, Fragment_Samples, 
           Core_Traps, Edge_Traps, Trap, Binomial) %>%
  summarise(Count = n(), .groups = 'keep') #count sps. obs/fragment

#Include all sps with 0 occurrences within sites
all_sites <- as.data.frame(unique(dbs[,1:7])) 
sps <- as.data.frame(unique(dbs$Binomial))
site_sps_combos <- crossing(all_sites, sps)
colnames(site_sps_combos)[8] <- 'Binomial'
dbs <- merge(dbs, site_sps_combos, all = T)
dbs$Count[is.na(dbs$Count)] <- 0

#Convert to detection/non-detection
dbs$Trap_Detection <- ifelse(dbs$Count > 0, 1, 0)
dbs <- dbs %>% 
  group_by(Fragment, Binomial) %>%
  mutate(Fragment_Detection = sum(Trap_Detection)) %>%
  ungroup()

#Spread to fragment x n_detection matrix
dbs_wide_incidence <- dbs %>% 
  select(Fragment, Fragment_Samples, Core_Traps, 
         Edge_Traps, Binomial, Fragment_Detection) %>%
  unique() %>%
  mutate(Binomial = sub(' ', '_', Binomial)) %>%
  pivot_wider(names_from = Binomial,
              values_from = Fragment_Detection)

#Spread to fragment x abundance matrix
dbs_wide_abundance <- dbs %>% 
  select(Fragment, Fragment_Samples, Core_Traps, 
         Edge_Traps, Binomial, Count) %>%
  unique() %>%
  mutate(Binomial = sub(' ', '_', Binomial)) %>%
  pivot_wider(names_from = Binomial,
              values_from = Count,
              values_fn = sum)

#Create Trapping Event x Site x Species Observation Array
dbs_array <- array(NA, dim = c(nrow(dbs_wide_incidence), length(unique(dbs$Binomial)),
                               max(dbs$Fragment_Samples)))
Trap_Habitat <- matrix(NA, ncol = max(dbs$Fragment_Samples), nrow = nrow(dbs_wide_incidence))

dbs <- dbs[order(dbs$Fragment, dbs$Habitat, dbs$Trap, dbs$Binomial),]

Traps <- dbs %>% select(Fragment, Trap, Habitat) %>% unique() 
Fragments <- unique(Traps$Fragment)
Species <- length(unique(dbs$Binomial))

for(i in 1:length(Fragments)){
  frag_traps <- filter(Traps, Fragment %in% Fragments[i])
  for(n in 1:nrow(frag_traps)){
    Trap_Habitat[i, n] <- frag_traps$Habitat[n]
    trap_capts <- filter(dbs, Trap %in% frag_traps$Trap[n])
    for(x in 1:Species){
      dbs_array[i, x, n] <- trap_capts$Trap_Detection[x]
    }
  }
}

rownames(Trap_Habitat) <- Fragments

#Read Fragment variables
frag_vars <- read.csv('../AF_DB_Fragment_Metrics.csv', row.names = 1) 
dbs_wide_incidence <- merge(frag_vars, dbs_wide_incidence)
dbs_wide_incidence <- dbs_wide_incidence[order(dbs_wide_incidence$Fragment),]
dbs_wide_abundance <- merge(frag_vars, dbs_wide_abundance)
dbs_wide_abundance <- dbs_wide_abundance[order(dbs_wide_abundance$Fragment),]

save(dbs_wide_incidence, dbs_wide_abundance, Trap_Habitat, dbs_array, 
     file = '../AF_DB_HMSOM_Data.RData')
save(matrix, dbs, file = '../AF_DBs_Species_Habitat_Data.RData')
write.csv(dbs_wide_incidence, '../AF_DB_Detections_and_Covariates.csv')


sps <- unique(dbs[,5:7])
matrix <- filter(dbs, Habitat == 'Matrix')
core <- filter(dbs, Habitat == 'Core')
edge <- filter(dbs, Habitat == 'Edge')

core <- core[,5:7]
core_count <- as.data.frame(table(core$Binomial))
colnames(core_count) <- c('Binomial', 'Observed_Abundance_Core')
core_count[,1] <- as.character(core_count[,1])

edge <- edge[,5:7]
edge_count <- as.data.frame(table(edge$Binomial))
colnames(edge_count) <- c('Binomial', 'Observed_Abundance_Edge')
edge_count[,1] <- as.character(edge_count[,1])

matrix <- matrix[,5:7]
matrix_count <- as.data.frame(table(matrix$Binomial))
colnames(matrix_count) <- c('Binomial', 'Observed_Abundance_Matrix')
matrix_count[,1] <- as.character(matrix_count[,1])

counts <- merge(sps, core_count, all = T)
counts <- merge(counts, edge_count, all = T)
counts <- merge(counts, matrix_count, all = T)
counts[is.na(counts)] <- 0

write.csv(counts, 'DBs_Species_List.csv')




