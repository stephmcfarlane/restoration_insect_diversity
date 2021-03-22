#                     Load Libraries                                                 #
######################################################################################

### load libraries
library(tidyverse)
library(vegan)
library(ggplot2)
library(dplyr)
library(grid)
library(gridExtra)
library(janitor)

######################################################################################
#                     Enter/orangize insects                                         #
######################################################################################

#setwd
setwd("~/r.Studio/Bio Sem/")



### Load in 
Insect.Data <- read_csv("./Bio.Sem.Insect.Data.csv")

### Remove variables 
Insect.Data <-Insect.Data%>%
  select(Date, EasementID, Total, Sample, Date, Order, Family)

### Remove empty rows
Insect.Data <- Insect.Data %>% filter_all(all_vars(!is.na(.)))
Insect.Data <- Insect.Data %>% filter_all(all_vars(complete.cases(.)))  

### Remove remnants
Remnants <- c('Arena Hill Prairie', 'Black Earth Rettenmund', 'Borah Creek', 'Oliver Prairie', 'Pleasant Valley Conservancy', 'Westport Drumlin')
Insect.Data <- Insect.Data[ !grepl(paste(Remnants, collapse="|"), Insect.Data$EasementID),]




######################################################################################
#                     Insect Richness                                                #
######################################################################################

### Insect richness averaged for the number of transects completed. 
Insect.Ave.Richness <- Insect.Data %>% 
  select(EasementID, Date, Sample, Family) %>% 
  filter(!is.na(EasementID)) %>% 
  group_by(EasementID, Sample, Date) %>% 
  filter(!duplicated(Family)) %>% #remove duplicates of Family/easement
  summarise(fam_rich = n()) %>% #count the number of families/easement
  group_by(EasementID ) %>% 
  summarise(Ave.Richness = mean(fam_rich)) #calculate average families

### Total insect richness, not averaged for the number of transects completed. 
Insect.Total.Richness <- Insect.Data%>%
  group_by(EasementID) %>%
  filter(!duplicated(Family))%>%
  dplyr::summarise(FamRichness = n())

### Compiling
Insect.Richness <- full_join(Insect.Ave.Richness, Insect.Total.Richness)


######################################################################################
#                     Insect Abundance                                               #
######################################################################################
### Code used from MS1 analysis 

### Number of bags id'd for 2019
num_trans19 <- Insect.Data%>% 
  select(EasementID, Date, Sample) %>%
  group_by(EasementID) %>% 
  separate(., Date, c('Month', 'Day', 'Year'), sep="/") %>% 
  filter(Year == '2019') %>% 
  filter(!duplicated(Sample)) %>% 
  summarise(Trans2019 = n())

### Find number of bags id'd for 2020
num_trans20 <- Insect.Data%>% 
  select(EasementID, Date, Sample) %>%
  group_by(EasementID) %>% 
  separate(., Date, c('Month', 'Day', 'Year'), sep="/") %>% 
  filter(Year == '2020') %>% 
  filter(!duplicated(Sample)) %>% 
  summarise(Trans2020 = n())

### combine dataframes into one. This gives us the total number of samples (bags or transects) that we currently have for each easement
num_trans <- num_trans19 %>% 
  full_join(num_trans20) %>% 
  replace(is.na(.), 0) %>% 
  mutate(total_trans = Trans2019 + Trans2020) %>% #calculate total samples by adding 2019 and 2020
  select(EasementID, total_trans) #simply df

### Calculate average abundance per Family per easement
Insect.Ave.Abundance <- Insect.Data %>% 
  select(EasementID, Date, Sample, Total, Family) %>% 
  filter(!is.na(EasementID)) %>% 
  group_by(EasementID) %>% 
  summarise(abundance= sum(Total)) %>% 
  full_join(num_trans) %>% 
  mutate(Ave.Abundance = abundance/total_trans) %>%
  select(EasementID, Ave.Abundance)

### Compiling
Insect.Summarized <-  full_join(Insect.Ave.Abundance,Insect.Richness)

######################################################################################
#                     Family abundance matrix/ rank abudnance                                       #
######################################################################################

### Gets the total average abundance (total, and then divided for the number of transects)
Total.Abundance <- Insect.Data %>% 
  select(EasementID, Date, Sample, Total, Family) %>% 
  filter(!is.na(EasementID)) %>% 
  group_by(EasementID) %>% 
  summarise(TotalAbundance= sum(Total))%>%
  full_join(num_trans)%>%
  mutate(TotalAveAbundance = TotalAbundance/total_trans)

### This gets the average abundance per family
Average.Abundance <- Insect.Data %>% 
  select(EasementID, Date, Sample, Total, Family) %>% 
  filter(!is.na(EasementID)) %>% 
  group_by(EasementID, Family) %>% 
  summarise(abundance= sum(Total)) %>% 
  full_join(num_trans) %>% 
  mutate(AveAbund = abundance/total_trans) 

### This combines the abundances for family, and total, and then creating the proportional abundances for the families. 
Abundance <-left_join(Average.Abundance, Total.Abundance)%>%
  mutate(PropAbundance = AveAbund/TotalAveAbundance)

### This creates a proportion table of the proportional abundances of the families in a site. 
Prop.Abundance <- Abundance %>% 
  select(EasementID, Family, PropAbundance) %>%
  pivot_wider(names_from = Family, values_from = PropAbundance) %>%
  replace(is.na(.), 0) %>%
  ungroup()

### Switches rows and columns
Prop.Abundance <- t(Prop.Abundance)

##### Need to get the 1st row as the title of the columns
##### then, need to sort from highest to lowest abundance at a site. 


######################################################################################
#                     Enter/Organize Plants                                          #
######################################################################################

### Load in
Plant.Data <- read_csv("./CVSPlantCommunityData_.csv")

### Only year needed    
#Plant.Data <- filter(Plant.Data, Year == 2019)

### Remove one weird case
#Plant.Data <- Plant.Data[-c(1),]

### Fix 00792
Plant.Data <- Plant.Data%>%
  mutate(EasementID = recode(EasementID, '792' = "00792"))


######################################################################################
#                     Plant CC Scores                                                #
######################################################################################

### Load all plant CC values
Plant.All.CC <- read_csv("./cc_trait_20190319 (Updated).csv")

### Cleaning -- Selecting CC scores for WI spp, Adding column with only SppCode and removing unneeded data
Plant.WI.CC <- Plant.All.CC %>% 
  select(WorZ_SppCode, WI_C) %>% 
  mutate(SppCode=substr(WorZ_SppCode, 2,7)) %>% 
  filter(!is.na(WI_C))

### Attaching C value to each plant at each site
Plant.Data <- Plant.Data %>% 
  left_join(Plant.WI.CC, "SppCode")

### Condesing df
Plant.Easement.CC <- Plant.Data%>%
  select(EasementID, SppCode, WI_C,PctCover1x1)

### Getting the average c value at each site.
Plant.Ave.CC <- Plant.Easement.CC %>%
  group_by(EasementID) %>% 
  dplyr::summarise(AveCC = mean(WI_C,na.rm=T)) 

### Merging for info at each site
Summarized.Plant.Insect <- left_join(Insect.Summarized, Plant.Ave.CC)


######################################################################################
#                     Plant native % cover                                           #
######################################################################################

### #Read in Native/NonNative data and select data needed
Plant.Native <- read_csv("./ExternalTraits (Updated).csv")%>%
  rename(SppCode = WSppCode) %>% 
  select(SppCode, USDANativeStatus)

### join Native/NonNative to the plant data 
Plant.Data <-Plant.Data %>% 
  left_join(Plant.Native, by = "SppCode") 

### Condensing df
Plant.Easement.Native <- Plant.Data%>%
  select(EasementID, PctCover1x1, Module, Corner, SppCode, USDANativeStatus, WI_C)%>% 
  na.omit(PctCover1x1)

### Summarize data to get mean native cover and mean nonnative cover per easem
Plant.Percent.Native<- Plant.Easement.Native %>% 
  group_by(EasementID, Module, Corner, USDANativeStatus) %>% 
  summarize(cover=sum(PctCover1x1)) %>% 
  group_by(EasementID, USDANativeStatus) %>% 
  summarise(avecover=mean(cover))

### Summarize data to get mean c by cover
Plant.Percent.CC<- Plant.Easement.Native %>% 
  group_by(EasementID, Module, Corner, WI_C) %>% 
  summarize(cover=sum(PctCover1x1)) %>% 
  group_by(EasementID, WI_C) %>% 
  summarise(avecover=mean(cover))

### get I/N own column, % by easement
Plant.Percent.I <- Plant.Percent.Native %>% 
  group_by(EasementID) %>%
  filter(USDANativeStatus == "I")
Plant.Percent.I <- rename(Plant.Percent.I, I = avecover)%>%
  select(EasementID, I)

Plant.Percent.N <- Plant.Percent.Native %>% 
  group_by(EasementID) %>%
  filter(USDANativeStatus == "N")
Plant.Percent.N <- rename(Plant.Percent.N, N = avecover)%>%
  select(EasementID, N)

Plant.Percent.NI <- Plant.Percent.Native %>% 
  group_by(EasementID) %>%
  filter(USDANativeStatus == "NI")
Plant.Percent.NI <- rename(Plant.Percent.NI, NI = avecover)%>%
  select(EasementID, NI)


### Merging to all summarized data by easement. 
Summarized.Plant.Insect <- left_join(Summarized.Plant.Insect, Plant.Percent.I)
Summarized.Plant.Insect <- left_join(Summarized.Plant.Insect, Plant.Percent.N)
Summarized.Plant.Insect <- left_join(Summarized.Plant.Insect, Plant.Percent.NI)%>%
  replace(is.na(.), 0)


Ctesting <- Summarized.Plant.Insect %>%
  select(-EasementID,)%>%
  replace(is.na(.), 0)


### Create proportion of native plants at a site. 
Summarized.Plant.Insect <- mutate(Summarized.Plant.Insect, PropN= ((N)/(N+I+NI)))
testing <- mutate(testing, PropN= ((N)/(N+I+NI)))

### Create proportion of invasive plants at a site. 
Summarized.Plant.Insect <- mutate(Summarized.Plant.Insect, PropI= ((I)/(N+I+NI)))
testing <- mutate(testing, PropI= ((I)/(N+I+NI)))

######################################################################################
#                     Plant richness                                                 #
######################################################################################

### Get the plant richness at each easement 
Plant.Richness <- Plant.Easement.CC%>%
  group_by(EasementID) %>%
  filter(!duplicated(SppCode))%>%
  dplyr::summarise(SppCode = n())

Summarized.Plant.Insect <- left_join(Summarized.Plant.Insect, Plant.Richness)

######################################################################################
#                    plant % gram or forb cover                                        #
######################################################################################

### create df for isolating form
Plant.Form <- Plant.All.CC %>% 
  select(WorZ_SppCode, USDAHabit) %>% 
  mutate(SppCode=substr(WorZ_SppCode, 2,7)) %>% 
  filter(!is.na(USDAHabit))

### merge the form to the plant data
Plant.Data <- Plant.Data%>% 
  left_join(Plant.Form, "SppCode")

  ### Condensing df
  Plant.Form.Percent <- Plant.Data%>%
    select(EasementID, PctCover1x1, Module, Corner, SppCode, USDAHabit, )%>% 
    na.omit(PctCover1x1)
  
  ### Summarize data to get mean gram and forb per easem
  Plant.Form.Percent<- Plant.Form.Percent %>% 
    group_by(EasementID, Module, Corner, USDAHabit) %>% 
    summarize(cover=sum(PctCover1x1)) %>% 
    group_by(EasementID, USDAHabit) %>% 
    summarise(avecover=mean(cover))
  
  ### get forb/gram own column, % by easement
  Plant.Form.Percent.Gram <- Plant.Form.Percent%>% 
    group_by(EasementID) %>%
    filter(USDAHabit == "Graminoid")
  Plant.Form.Percent.Gram <- rename(Plant.Form.Percent.Gram, Gram = avecover)%>%
    select(EasementID, Gram)
  
  Plant.Form.Percent.Forb <- Plant.Form.Percent%>% 
    group_by(EasementID) %>%
    filter(USDAHabit == "Forb_herb")
  Plant.Form.Percent.Forb <- rename(Plant.Form.Percent.Forb, Forb = avecover)%>%
    select(EasementID, Forb)
  

  ### Merging to all summarized data by easement. 
  Summarized.Plant.Form <- full_join(Plant.Form.Percent.Gram, Plant.Form.Percent.Forb)%>%
    replace(is.na(.), 0)
  
  Summarized.Plant.Insect <- left_join(Summarized.Plant.Insect, Summarized.Plant.Form)
 

  testing <- Summarized.Plant.Insect %>%
    select(-EasementID,)%>%
    replace(is.na(.), 0)
  
  ### Create proportion of form plants at a site. 
  Summarized.Plant.Insect <- mutate(Summarized.Plant.Insect, PropGram= ((Gram)/(Gram + Forb)))
  testing <- mutate(testing, PropGram= ((N)/(Gram + Forb)))
  
  ### Create proportion of form plants at a site. 
  Summarized.Plant.Insect <- mutate(Summarized.Plant.Insect, PropForb= ((Forb)/(Gram + Forb)))
  testing <- mutate(testing, PropForb= ((Forb)/(Gram + Forb)))
  
  
######################################################################################
#                     ETC                                                            #
######################################################################################

### Change the NA to zero
#Summarized.Plant.Insect[is.na(Summarized.Plant.Insect)] <- 0 

### Change fam and ord richness to numeric values from integers 
#Summarized.Plant.Insect[2:3] <- lapply(Summarized.Plant.Insect[2:3], as.numeric)

######################################################################################
#                     Correlations                                                   #
######################################################################################

### Gets the pearson value for all variables 
pearson <- cor(testing, method = c("pearson"))

### Gets the p value and condenses to one df
pv1 <- cor.test(testing$Ave.Richness, testing$AveCC, method = "pearson")
pv1$p.value
pv2<- cor.test(testing$Ave.Abundance, testing$AveCC, method = "pearson")
pv2$p.value
#pv3<- cor.test(testing$Abundance, testing$AveCC, method = "pearson")
#pv3$p.value
pv4 <- cor.test(testing$Ave.Abundance, testing$PropN, method = "pearson")
pv4$p.value
pv5<- cor.test(testing$Ave.Richness, testing$PropN, method = "pearson")
pv5$p.value
#pv6<- cor.test(testing$Abundance, testing$PropN, method = "pearson")
#pv6$p.value
pv7<- cor.test(testing$Ave.Abundance, testing$SppCode, method = "pearson")
pv7$p.value
pv8<- cor.test(testing$Ave.Richness, testing$SppCode, method = "pearson")
pv8$p.value

pv9<- cor.test(testing$Ave.Richness, testing$PropGram, method = "pearson")
pv9$p.value
pv10<- cor.test(testing$Ave.Abundance, testing$PropGram, method = "pearson")
pv10$p.value
pv11<- cor.test(testing$Ave.Richness, testing$PropForb, method = "pearson")
pv11$p.value
pv12<- cor.test(testing$Ave.Abundance, testing$PropForb, method = "pearson")
pv12$p.value

#working on creating a df that easily shows all the correlations and p values

#Pearson <- rbind(pv1$p.value, pv2$p.value, pv4$p.value, pv5$p.value, pv7$p.value, pv8$p.value,pv9$p.value,pv10$p.value,pv11$p.value,pv12$p.value)


#pv <- combine(pv1, pv2, pv3, pv4, pv5, pv6, pv7, pv8)

#pearson <- pearson[c(3, 8, 12, 17, 21, 26, 30, 35, 39, 44, 48, 53, 57, 62, 66, 71)]
#remove(pv1, pv2, pv3, pv4, pv5, pv6, pv7, pv8)


######################################################################################
#                     Graphing                                                       #
######################################################################################

p1 <- ggplot(data = testing, aes(x = AveCC, y = Ave.Richness)) + 
  geom_point(aes(size=0.7),show.legend = F,color='black') +
  geom_smooth(aes(colour="green"),show.legend = F,method= "lm", se=FALSE)+
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  theme(text = element_text(size=22))+
  labs(y = "
       Insect Family Richness") +
  theme(axis.text.y = element_text(angle=90))+ 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

p2 <- ggplot(data = testing, aes(x = SppCode, y = Ave.Richness)) + 
  geom_point(aes(size=0.7),show.legend = F,color='black') +
  geom_smooth(aes(colour="green"),show.legend = F,method= "lm", se=FALSE)+
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  theme(axis.title.y=element_blank(), axis.text.y=element_blank(),
        axis.ticks.y=element_blank())+ 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

p3 <- ggplot(data = testing, aes(x = AveCC, y = Ave.Abundance)) + 
  geom_point(aes(size=0.7),show.legend = F,color='black') +
  geom_smooth(aes(colour="green"),show.legend = F,method= "lm", se=FALSE)+
  theme(text = element_text(size=22))+
  labs(x= "Average Plant C Score 
  (Coefficient of Conservatism)", y = "
       Insect Abundance") +
  theme(axis.text.y = element_text(angle=90))+ 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())


p4 <- ggplot(data = testing, aes(x = SppCode, y = Ave.Abundance)) + 
  geom_point(aes(size=0.7),show.legend = F,color='black') +
  geom_smooth(aes(colour="green"),show.legend = F,method= "lm", se=FALSE)+
  theme(axis.title.y=element_blank(), axis.text.y=element_blank(),
        axis.ticks.y=element_blank())+
  theme(text = element_text(size=22))+
  labs(x= "Plant Species Richness") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())


graphs1 <- grid.arrange(p1, p2, p3, p4, nrow=2, 
                       top = textGrob("Relationship between insect communities and prairie quality", gp=gpar(fontsize=25,font=1)), 
                       bottom = textGrob("Prairie Quality
                               ", gp=gpar(fontsize=20,font=1)))

### more graphs to test


p5 <- ggplot(data = testing, aes(x = PropN, y = Ave.Abundance)) + 
  geom_point(aes(size=0.7),show.legend = F,color='black') +
  geom_smooth(aes(colour="green"),show.legend = F,method= "lm", se=FALSE)+
  theme(text = element_text(size=22))+
  labs(x= "propN", y= "ave abundance") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

p6 <- ggplot(data = testing, aes(x = PropN, y = Ave.Richness)) + 
  geom_point(aes(size=0.7),show.legend = F,color='black') +
  geom_smooth(aes(colour="green"),show.legend = F,method= "lm", se=FALSE)+
  theme(text = element_text(size=22))+
  labs(x= "propN", y= "ave richness") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

graphs2 <- grid.arrange(p5, p6, nrow=1)


### more graphs to test


p7 <- ggplot(data = testing, aes(x = PropGram, y = Ave.Abundance)) + 
  geom_point(aes(size=0.7),show.legend = F,color='black') +
  geom_smooth(aes(colour="green"),show.legend = F,method= "lm", se=FALSE)+
  theme(text = element_text(size=22))+
  labs(x= "PropGram", y= "ave abundance") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

p8 <- ggplot(data = testing, aes(x = PropGram, y = Ave.Richness)) + 
  geom_point(aes(size=0.7),show.legend = F,color='black') +
  geom_smooth(aes(colour="green"),show.legend = F,method= "lm", se=FALSE)+
  theme(text = element_text(size=22))+
  labs(x= "PropGram", y= "ave richness") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p9 <- ggplot(data = testing, aes(x = PropForb, y = Ave.Abundance)) + 
  geom_point(aes(size=0.7),show.legend = F,color='black') +
  geom_smooth(aes(colour="green"),show.legend = F,method= "lm", se=FALSE)+
  theme(text = element_text(size=22))+
  labs(x= "PropForb", y= "ave abundance") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

p10 <- ggplot(data = testing, aes(x = PropForb, y = Ave.Richness)) + 
  geom_point(aes(size=0.7),show.legend = F,color='black') +
  geom_smooth(aes(colour="green"),show.legend = F,method= "lm", se=FALSE)+
  theme(text = element_text(size=22))+
  labs(x= "PropForb", y= "ave richness") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

graphs3 <- grid.arrange(p7, p8, p9,p10,nrow=2)



###When exporting to clipboard, set height to 1500, preview, then copy. 


######################################################################################
#                     Notes/Extras                                                   #
######################################################################################