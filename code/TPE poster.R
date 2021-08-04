library(tidyverse)





treatment <- read_csv("../PrairieRestoration/datasets/TreatmentNov20.csv") %>% 
  mutate(EasementID = replace(EasementID, EasementID == "792", "00792")) %>% 
  mutate(RestorationCategory = replace(RestorationCategory, RestorationCategory =="Seed+ Fire", "Seeded + Fire"),
         RestorationCategory = replace(RestorationCategory, RestorationCategory == "Seed", "Seeded Only"),
         RestorationCategory = replace(RestorationCategory, RestorationCategory == "No Seed", "Not Seeded")) %>% 
  mutate(RestorationCategory = factor(RestorationCategory, levels = c("Not Seeded", "Seeded Only", "Seeded + Fire")))  #####this code reorders the treatments --the default is alphabetical

insect <- read_csv("raw/Insects_data_02_21.csv") %>% 
  select(EasementID, Date, Sample, Total, Order, Family) %>% 
  left_join(treatment) %>% 
  replace_na(RestorationCategory, "Remnant")


######################################################################################
#                     Enter data: insects,sites,grass                                #
######################################################################################

#upload data
TPE.AllData <- read_csv("2.6.21.TPE.AllData.csv")


#make sure that the 1st empty line is deleted in the collection sites csv file
TPE.Sites <- read_csv("2.6.21.TPE.CollectionSites.csv")
TPE.Sites = subset(TPE.Sites, select = (Category:EasementID))


Grass <- read_csv("grassland1.5km.csv")

Summary <- full_join(TPE.Sites, Grass, by = "EasementID")


######################################################################################
#                    other anova code                                    #
######################################################################################
richness_anova <- aov(richness ~ RestorationCategory, data = spprichness)
summary(richness_anova)
TukeyHSD(richness_anova)

######################################################################################
#                     Insect Abundance & Graphs                                      #
######################################################################################

#get abundance numbers
TPE.sum <- TPE.AllData %>%
  group_by(EasementID)%>%
  summarize(Abundance=sum(Total))

Summary <- full_join(TPE.sum, Summary, by = "EasementID")%>%
  slice(1:25)



# Percent grassland and total insect abundance
AbundancexGrass <- ggplot(Summary, aes(x = PercentGrass, y = log(Abundance), color = Category)) +
  geom_point(size = 4) + theme_classic() + 
  xlab("Percent grassland") +
  ylab("Log Total\nInsect Abundance") +
  theme(axis.title = element_text(size = 25))  + 
  geom_smooth(method = "lm", se = FALSE, linetype = 2)
AbundancexGrass + scale_color_manual(values=c("burlywood4", "darkolivegreen4", "black", "yellow")) + 
  theme(legend.position="none", axis.text.x = element_text(face="bold",
                                                           size=18),
        axis.text.y = element_text(face="bold", 
                                   size=18))


# Restoration Category and Total insect abundance
AbundancexRestoration <- ggplot(Summary, aes(x = Category, y = log(Abundance), fill = Category)) +
  geom_boxplot(size = 1, outlier.shape = NA) + theme_classic() + 
  ggtitle("") +
  xlab("Restoration Category") +
  ylab("Log Insect Abundance") +
  theme(axis.title = element_text(size = 25))  + 
  stat_summary(fun.y = "mean", shape = 17, size = 3, geom = "point") + 
  geom_jitter(width = 0.1, size = 3)
AbundancexRestoration + scale_fill_manual(values=c("burlywood4", "darkolivegreen4","black", "yellow")) +
  theme(legend.position="none", axis.text.x = element_text(face="bold",
                                                           size=18),
        axis.text.y = element_text(face="bold", 
                                   size=18))


######################################################################################
#                     Family x Site Matrix    = richness/diversity                   #
######################################################################################

# create the site x family matrix
FamilyMatrix <- TPE.AllData %>% 
  pivot_wider(
    names_from = Family,
    values_from = Total)

# site x species matrix with one row per site
FamilyMatrix <- FamilyMatrix %>% 
  mutate_at(vars(Formicidae:Psychidae), ~replace(., is.na(.), 0)) %>%
  group_by(EasementID) %>%
  summarise_at(vars(Formicidae:Psychidae), sum)%>%
  slice(1:25)

FamilyMatrix$Richness <- specnumber(FamilyMatrix[2:166], MARGIN = 1) 

FamilyMatrix$Simpson <- diversity(FamilyMatrix[2:166], index = "simpson", MARGIN = 1, base = exp(1))

FamilyMatrix = subset(FamilyMatrix, select = c(EasementID,Richness,Simpson))

Summary <- full_join(FamilyMatrix, Summary, by = "EasementID")



######################################################################################
#                     Insect Richness/Diveristy Graphs                               #
######################################################################################

#percent grass and insect richness
ggplot(Summary, aes(x = PercentGrass, y = log(Richness), color = Category)) +
  geom_point() + theme_classic() + 
  ggtitle("") +
  xlab("% Area Grassland within 1.5 km") +
  ylab("Log Family Richness") +
  geom_smooth(method = "lm", se = FALSE)

# Restoration Category and family richness
ggplot(Summary, aes(x = Category, y = log(Richness), fill = Category)) +
  geom_boxplot() + theme_classic() + 
  ggtitle("") +
  ylab("Log Family Richness") +
  xlab("Restoration Category")
theme(plot.title = element_text(size = 25), axis.title = element_text(size = 15))  + 
  stat_summary(fun.y = "mean", colour = "red", size = 2, geom = "point") + 
  geom_jitter(width = 0.1) 

# Percent grassland and Simpson's Diversity
DiversityxGrass<- ggplot(Summary, aes(x = PercentGrass, y = Simpson, color = Category)) +
  geom_point(size = 4) + theme_classic() + 
  xlab("% Area Grassland within 1.5 km") +
  ylab("Simpson's Diversity") +
  theme(axis.title = element_text(size = 25)) +
  geom_smooth(method = "lm", se = FALSE, linetype = 2)
DiversityxGrass + scale_color_manual(values=c("burlywood4", "darkolivegreen4", "black", "yellow")) +
  theme(axis.text.x = element_text(face="bold",
                                   size=18),
        axis.text.y = element_text(face="bold", 
                                   size=18),
        legend.position = c(.8,.3),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 20),
        legend.background = element_rect(linetype = "solid"))

# Restoration Category and Simpson's Diversity
DiversityxRestoration <- ggplot(Summary, aes(x = Category, y = Simpson, fill = Category)) +
  geom_boxplot(size = 1) + theme_classic() + 
  xlab("Restoration Category") +
  ylab("Simpsons Diversity") +
  theme(axis.title = element_text(size = 25))  + 
  stat_summary(fun.y = "mean", shape = 17, size = 3, geom = "point") + 
  geom_jitter(width = 0.1, size = 3) 
DiversityxRestoration + scale_fill_manual(values = c("burlywood4", "darkolivegreen4", "black", "yellow")) +
  theme(axis.text.x = element_text(face="bold",
                                   size=18),
        axis.text.y = element_text(face="bold", 
                                   size=18), 
        legend.position = "none")
