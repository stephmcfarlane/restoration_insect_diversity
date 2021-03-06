library(tidyverse)

#upload data
TPE.AllData <- read_csv("raw/3.6.21_insect_data.csv") %>% 
  left_join(treatment) %>% 
  relocate(., RestorationCategory, .before = "Date" ) %>% 
  filter(RestorationCategory == "Seeded + Fire")


###########Reading in Treatment data and cleaning#############
treatment <- read_csv("raw/TreatmentNov20.csv") %>% 
  filter(!is.na(RestorationCategory)) %>% 
  mutate(EasementID = replace(EasementID, EasementID == "792", "00792")) %>% 
  mutate(RestorationCategory = replace(RestorationCategory, RestorationCategory =="Seed+ Fire", "Seeded + Fire"),
         RestorationCategory = replace(RestorationCategory, RestorationCategory == "Seed", "Seeded Only"),
         RestorationCategory = replace(RestorationCategory, RestorationCategory == "No Seed", "Not Seeded")) %>% 
  mutate(RestorationCategory = factor(RestorationCategory, levels = c("Not Seeded", "Seeded Only", "Seeded + Fire")))  #####this code reorders the treatments --the default is alphabetical

insects <- read_csv("raw/3.6.21_insect_data.csv") %>% 
  select(EasementID, Date, Sample, Total, Order, Family) %>% 
  filter(!is.na(EasementID)) %>% 
  group_by(EasementID) %>% 
  filter(!duplicated(Family))%>% 
  summarise(fam_rich = n()) %>% 
  left_join(treatment)
  
Coleoptera_rich <- read_csv("raw/3.6.21_insect_data.csv") %>% 
  select(EasementID, Date, Sample, Total, Order, Family) %>% 
  filter(Order == "Coleoptera") %>% 
  filter(!is.na(EasementID)) %>% 
  group_by(EasementID) %>% 
  filter(!duplicated(Family))%>% 
  summarise(fam_rich = n()) %>% 
  left_join(treatment)

Hymenoptera_rich <- read_csv("raw/3.6.21_insect_data.csv") %>% 
  select(EasementID, Date, Sample, Total, Order, Family) %>% 
  filter(Order == "Hymenoptera") %>% 
  filter(!is.na(EasementID)) %>% 
  group_by(EasementID) %>% 
  filter(!duplicated(Family))%>% 
  summarise(fam_rich = n()) %>% 
  left_join(treatment)

Hemiptera_rich <- read_csv("raw/3.6.21_insect_data.csv") %>% 
  select(EasementID, Date, Sample, Total, Order, Family) %>% 
  filter(Order == "Hemiptera") %>% 
  filter(!is.na(EasementID)) %>% 
  group_by(EasementID) %>% 
  filter(!duplicated(Family))%>% 
  summarise(fam_rich = n()) %>% 
  left_join(treatment)

Orthoptera_rich <- read_csv("raw/3.6.21_insect_data.csv") %>% 
  select(EasementID, Date, Sample, Total, Order, Family) %>% 
  filter(Order == "Orthoptera") %>% 
  filter(!is.na(EasementID)) %>% 
  group_by(EasementID) %>% 
  filter(!duplicated(Family))%>% 
  summarise(fam_rich = n()) %>% 
  left_join(treatment)



fam_rich <- read_csv("raw/3.6.21_insect_data.csv") %>% 
  select(EasementID, Date, Sample, Total, Order, Family) %>% 
  filter(!is.na(EasementID)) %>% 
  group_by(EasementID) %>% 
  filter(!duplicated(Family))%>% 
  summarise(fam_rich = n()) %>% 
  left_join(treatment)

Orthoptera_rich %>% 
  ggplot(aes(RestorationCategory, fam_rich, fill=RestorationCategory)) +
  geom_boxplot(outlier.alpha = 0) +
  theme_classic()+
  scale_fill_manual(values = c("slategray3", "darkolivegreen3", "thistle3"),
                    name = "Site Categories")+
  labs(title = "Insect Family Richness \n per Restoration Category")+
  xlab("\n Restoration Category") +
  ylab("Average Family Richness")+
  theme(plot.title = element_text(family = "Palatino", size = 18, hjust = 0.5, vjust = 2),
        axis.title.x = element_text(family = "Palatino", size = 15, vjust = 4),
        axis.title.y = element_text(family = "Palatino", size = 15, hjust = 0.5, vjust =4),
        axis.text.x =  element_text(family = "Palatino", size = 12),
        axis.text.y = element_text(family = "Palatino", size = 12),
        legend.title = element_blank(),
        legend.position = "none",
        plot.margin = unit(c(1,1,2,1), "lines")) +
  geom_jitter(alpha=.2, width = .1, size = 3)
  
  