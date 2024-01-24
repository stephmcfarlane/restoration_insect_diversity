# January 27, 2023
# by: Stephanie McFarlane
# Cleaning code for insect analyses

#Load packages####
library(tidyverse)

# Site History Data ####
# Upload easement treatment (aka restoration category) data & clean ##
rest_history <- read_csv("../vegetation_outcomes_of_restoration/clean/rest_history.csv") %>% rename(EasementID = SiteID)

##Fire History####
# file appears to be missing fire data for a few remnants
fire_years  <- read_csv("raw/restoration_history.csv") %>% 
  select(EasementID, Fire_Years) %>%
  mutate(EasementID = recode(EasementID, '792' = "00792")) %>% 
  separate(Fire_Years, c("Fire1", "Fire2", "Fire3")) 

#Insect Data ####
  all_data <- read_csv("raw/All_Insect_Data.csv") %>% 
  mutate(EasementID = replace(EasementID, EasementID == "792", "00792")) %>% 
  left_join(rest_history) %>% 
  mutate(RestorationCategory = factor(RestorationCategory, levels = c("No Seed", "Seed", "Seed + Fire", "Remnant"))) %>% #####this code reorders the treatments --the default is alphabetical
  select(-'...1') %>% #remove random column with nothing in it
  relocate(., RestorationCategory, .before = "Date" ) %>% 
  separate(., Date, c('Month', 'Day', 'Year'), sep="/") %>% 
  filter(EasementID != "Borah Creek") %>% 
  filter(EasementID != "Lulu lake") %>%
  filter(EasementID != "Scuppernog") %>% 
  filter(EasementID != "Hauser Road") %>% 
  filter(EasementID != "UW-Arboretum") %>% 
  filter(EasementID != "Snapper Prairie") %>% 
  filter(EasementID !="00MDP"|Year!="2020") %>%  ##Need to filter out 00MDP 2020 bc we stopped data collection due to COVID, but not 00MDP 2019
  select(EasementID, Year, Month, Day, Sample, Total, Family) %>% 
  group_by(EasementID,  Month, Day,Sample, Year, Family)%>% 
  summarise(Total = sum(Total))
#write_csv(all_data, "clean/all_data.csv")

all_data <- read_csv("../restoration_insect_diversity/clean/all_data.csv")

# Considered removing data sweeps that were dominated by an insect family that tends to aggragate; however ultimately decided against this#### 
# all_data_rm_agg <- all_data %>% 
#  group_by(EasementID, Year, Family) %>% 
#  summarise(Total = sum(Total)) %>% 
#  ungroup() %>% 
#  group_by(EasementID, Year) %>% 
#  mutate(site_total = sum(Total)) %>% 
#  mutate(percent_total_ind = Total/site_total) %>% 
#  filter(percent_total_ind < 0.3)



## Natural History ####
top_families <-  all_data %>% 
  group_by(Family) %>% 
  summarise(total_per_fam = sum(Total))

top_families_cat <- all_data %>% 
  left_join(rest_history) %>% 
  group_by(Family, RestorationCategory) %>% 
  summarise(total_per_fam = sum(Total))
