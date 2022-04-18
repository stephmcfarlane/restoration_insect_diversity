##Beginning Visualizations for Insect Data####

# Calculate the family richness for each Order ####
coleoptera_rel_rich <- all_data %>% 
  select(EasementID, RestorationCategory, Date, Sample, Total, Order, Family) %>% 
  filter(Order == "Coleoptera") %>% 
  filter(!is.na(EasementID)) %>% 
  group_by(EasementID, RestorationCategory, Date, Sample) %>% 
  filter(!duplicated(Family))%>% 
  summarise(fam_rich = n()) %>% 
  group_by(EasementID, RestorationCategory) %>% 
  summarise(ave_rich = mean(fam_rich)) 

hymenoptera_rel_rich <- read_csv("raw/All_Insect_Data.csv") %>% 
  select(EasementID, Date, Sample, Total, Order, Family) %>% 
  filter(Order == "Hymenoptera") %>% 
  filter(!is.na(EasementID)) %>% 
  group_by(EasementID, Date, Sample) %>% 
  filter(!duplicated(Family))%>% 
  summarise(fam_rich = n()) %>% 
  group_by(EasementID) %>% 
  summarise(ave_rich = mean(fam_rich)) %>% 
  left_join(treatment)

hemiptera_rel_rich <- read_csv("raw/All_Insect_Data.csv") %>% 
  select(EasementID, Date, Sample, Total, Order, Family) %>% 
  filter(Order == "Hemiptera") %>% 
  filter(!is.na(EasementID)) %>% 
  group_by(EasementID, Date, Sample) %>% 
  filter(!duplicated(Family))%>% 
  summarise(fam_rich = n()) %>% 
  group_by(EasementID) %>% 
  summarise(ave_rich = mean(fam_rich)) %>% 
  left_join(treatment)

orthoptera_rel_rich <- read_csv("raw/All_Insect_Data.csv") %>% 
  select(EasementID, Date, Sample, Total, Order, Family) %>% 
  filter(Order == "Orthoptera") %>% 
  filter(!is.na(EasementID)) %>% 
  group_by(EasementID, Date, Sample) %>% 
  filter(!duplicated(Family))%>% 
  summarise(fam_rich = n()) %>% 
  group_by(EasementID) %>% 
  summarise(ave_rich = mean(fam_rich)) %>% 
  left_join(treatment)

diptera_rel_rich <- read_csv("raw/All_Insect_Data.csv") %>% 
  select(EasementID, Date, Sample, Total, Order, Family) %>% 
  filter(Order == "Diptera") %>% 
  filter(!is.na(EasementID)) %>% 
  group_by(EasementID, Date, Sample) %>% 
  filter(!duplicated(Family))%>% 
  summarise(fam_rich = n()) %>% 
  group_by(EasementID) %>% 
  summarise(ave_rich = mean(fam_rich)) %>% 
  left_join(treatment)

thysanoptera_ave_rich <- read_csv("raw/All_Insect_Data.csv") %>% 
  select(EasementID, Date, Sample, Total, Order, Family) %>% 
  filter(Order == "Thysanoptera") %>% 
  filter(!is.na(EasementID)) %>% 
  group_by(EasementID, Date, Sample) %>% 
  filter(!duplicated(Family))%>% 
  summarise(fam_rich = n()) %>% 
  group_by(EasementID) %>% 
  summarise(ave_rich = mean(fam_rich)) %>% 
  left_join(treatment)

odonata_ave_rich <- read_csv("raw/All_Insect_Data.csv") %>% 
  select(EasementID, Date, Sample, Total, Order, Family) %>% 
  filter(Order == "Odonata") %>% 
  filter(!is.na(EasementID)) %>% 
  group_by(EasementID, Date, Sample) %>% 
  filter(!duplicated(Family))%>% 
  summarise(fam_rich = n()) %>% 
  group_by(EasementID) %>% 
  summarise(ave_rich = mean(fam_rich)) %>% 
  left_join(treatment)

neuroptera_ave_rich <- read_csv("raw/All_Insect_Data.csv") %>% 
  select(EasementID, Date, Sample, Total, Order, Family) %>% 
  filter(Order == "Neuroptera") %>% 
  filter(!is.na(EasementID)) %>% 
  group_by(EasementID, Date, Sample) %>% 
  filter(!duplicated(Family))%>% 
  summarise(fam_rich = n()) %>% 
  group_by(EasementID) %>% 
  summarise(ave_rich = mean(fam_rich)) %>% 
  left_join(treatment)

ephemeroptera_ave_rich <- read_csv("raw/All_Insect_Data.csv") %>% 
  select(EasementID, Date, Sample, Total, Order, Family) %>% 
  filter(Order == "Ephemeroptera") %>% 
  filter(!is.na(EasementID)) %>% 
  group_by(EasementID, Date, Sample) %>% 
  filter(!duplicated(Family))%>% 
  summarise(fam_rich = n()) %>% 
  group_by(EasementID) %>% 
  summarise(ave_rich = mean(fam_rich)) %>% 
  left_join(treatment)

lepidoptera_ave_rich <- read_csv("raw/All_Insect_Data.csv") %>% 
  select(EasementID, Date, Sample, Total, Order, Family) %>% 
  filter(Order == "Lepidoptera") %>% 
  filter(!is.na(EasementID)) %>% 
  group_by(EasementID, Date, Sample) %>% 
  filter(!duplicated(Family))%>% 
  summarise(fam_rich = n()) %>% 
  group_by(EasementID) %>% 
  summarise(ave_rich = mean(fam_rich)) %>% 
  left_join(treatment)

strepsiptera_ave_rich <- read_csv("raw/All_Insect_Data.csv") %>% 
  select(EasementID, Date, Sample, Total, Order, Family) %>% 
  filter(Order == "Strepsiptera") %>% 
  filter(!is.na(EasementID)) %>% 
  group_by(EasementID, Date, Sample) %>% 
  filter(!duplicated(Family))%>% 
  summarise(fam_rich = n()) %>% 
  group_by(EasementID) %>% 
  summarise(ave_rich = mean(fam_rich)) %>% 
  left_join(treatment)

psocodea_ave_rich <- read_csv("raw/All_Insect_Data.csv") %>% 
  select(EasementID, Date, Sample, Total, Order, Family) %>% 
  filter(Order == "Psocodea") %>% 
  filter(!is.na(EasementID)) %>% 
  group_by(EasementID, Date, Sample) %>% 
  filter(!duplicated(Family))%>% 
  summarise(fam_rich = n()) %>% 
  group_by(EasementID) %>% 
  summarise(ave_rich = mean(fam_rich)) %>% 
  left_join(treatment)
