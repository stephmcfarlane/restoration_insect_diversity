

#### data testing on 10/17/2021


library(tidyverse)
library(vegan)
# NOTE: adding four # after a note makes it into a section heading!

## Load data #####
# Upload easement treatment (aka restoration category) data & clean ##
rest_history <- read_csv("raw/Seeding and Fire checklist.csv") %>% 
  select(EasementID, RestorationCategory)

# Upload raw insect sweep identification data
# This csv was downloaded from the Google Doc data entry spreadsheet on 6 March 2021

all_data <- read_csv("raw/All_Insect_Data.csv") %>% 
  left_join(rest_history) %>% 
  mutate(EasementID = replace(EasementID, EasementID == "792", "00792")) %>% 
  mutate(RestorationCategory = factor(RestorationCategory, levels = c("No Seed", "Seed", "Seed + Fire", "Remnant"))) %>% #####this code reorders the treatments --the default is alphabetical
  select(-X1) %>% #remove random column with nothing in it
  relocate(., RestorationCategory, .before = "Date" ) %>% 
  separate(., Date, c('Month', 'Day', 'Year'), sep="/") %>% 
  filter(EasementID != "Borah Creek") %>% 
  filter(EasementID != "Lulu lake") %>%
  filter(EasementID != "Scuppernog") %>% 
  filter(EasementID !="00MDP"|Year!="2020") ##Need to filter out 00MDP 2020, but not 00MDP 2019

rest_year <- read_csv("clean/enroll_rest.csv")


# Color Palettes #
palette1<-  c("#767171", "#A9D18E", "#548235", "#A49988")
palette2 <-  c("tomato3", "#A49988","#5F9EA0", "#006887")



##### test for normality with residuals in abundance

## Get mean abundance for all insects, get the difference from the mean for all samples

# get total abundance per site, per year
abundance_family <- all_data %>%
  select(EasementID, RestorationCategory, Year, Sample, Total, Family) %>% 
  filter(!is.na(EasementID)) %>% 
  group_by(EasementID, Year) %>% 
  summarise(abundance= sum(Total))%>%
  select(EasementID, abundance)

mean1<-mean(abundance_family$abundance)

abundance_family$residual <- abundance_family$abundance-mean1

shapiro.test(abundance_family$residual)

ggplot(data = abundance_family) + geom_histogram(mapping = aes(x = residual), bins = 6)


## Get mean richness for all insects, get the difference from the mean for all samples

# Calculate Family richness per easement, per year
insects_rich_2019 <- all_data %>% 
  select(EasementID, RestorationCategory, Year, Sample, Family) %>% 
  filter(!is.na(EasementID)) %>% 
  filter(Year == '2019') %>% 
  group_by(EasementID, Year) %>% 
  filter(!duplicated(Family)) %>% #remove duplicates of Family/easement
  summarise(fam_rich = n()) #count the number of families/easement
insects_rich_2020 <- all_data %>% 
  select(EasementID, RestorationCategory, Year, Sample, Family) %>% 
  filter(!is.na(EasementID)) %>% 
  filter(Year == '2020') %>% 
  group_by(EasementID, Year) %>% 
  filter(!duplicated(Family)) %>% #remove duplicates of Family/easement
  summarise(fam_rich = n())  #count the number of families/easement
##combine the richness per year per easment

insect_rich <- insects_rich_2019%>%
  full_join(insects_rich_2020)


mean2<-mean(insect_rich$fam_rich)

insect_rich$residual <- insect_rich$fam_rich-mean2

shapiro.test(insect_rich$residual)

ggplot(data = insect_rich) + geom_histogram(mapping = aes(x = residual))



#diversity 
sitebyfam <- all_data %>% 
  select(EasementID, RestorationCategory, Year, Sample, Total, Family) %>% 
  filter(!is.na(EasementID)) %>% 
  group_by(EasementID, RestorationCategory, Year, Family) %>% 
  summarise(abundance= sum(Total)) %>% 
  select(EasementID, RestorationCategory, Family, Year, abundance) %>%
  pivot_wider(names_from = Family, values_from = abundance) %>%
  replace(is.na(.), 0) %>%
  ungroup()

##Shannon div per easement per year

sitebyfam$Shannon <- diversity(sitebyfam[4:182], index = "shannon", MARGIN = 1, base = exp(1))

insect_div <- sitebyfam%>%
  select(EasementID, Year, Shannon)

mean3<-mean(insect_div$Shannon)

insect_div$residual <- insect_div$Shannon-mean3

shapiro.test(insect_div$residual)

ggplot(data = insect_div) + geom_histogram(mapping = aes(x = residual), bins = 6)




# Raw data visualizations -------------------------------------------------



## Calculate average Family richness per easement
insects_ave_rich <- all_data %>% 
  select(EasementID, RestorationCategory, Date, Sample, Family) %>% 
  filter(!is.na(EasementID)) %>% 
  group_by(EasementID, RestorationCategory, Date, Sample) %>% 
  filter(!duplicated(Family)) %>% #remove duplicates of Family/easement
  summarise(fam_rich = n()) %>% #count the number of families/easement
  group_by(EasementID, RestorationCategory) %>% 
  summarise(ave_rich = mean(fam_rich)) #calculate average families

insects_rich_trans <- all_data %>% 
  select(EasementID, RestorationCategory, Date, Sample, Family) %>% 
  filter(!is.na(EasementID)) %>% 
  group_by(EasementID, RestorationCategory, Date, Sample) %>% 
  filter(!duplicated(Family)) %>% #remove duplicates of Family/easement
  summarise(fam_rich = n()) 

# visualize family richment per restoration age
rich_year <-  rest_year %>% 
  inner_join(insects_ave_rich) %>% 
  mutate(time_since = (2021- rest_year))


rich_year %>% 
  ggplot(aes(time_since, ave_rich, color = RestorationCategory)) +
  geom_point(outlier.alpha = 0) +
  geom_smooth(aes(color = RestorationCategory), method = "lm") +
  theme_classic()+
  scale_color_manual(values = palette2,
                     name = "Site Categories")+
  labs(title = "Average Family Richness per Restoration Year")+
  xlab("\n Years Since Restoration") +
  ylab("Ave Family Richness")+
  theme(plot.title = element_text(family = "Palatino", size = 18, hjust = 0.5, vjust = 2),
        axis.title.x = element_text(family = "Palatino", size = 15, vjust = 4),
        axis.title.y = element_text(family = "Palatino", size = 15, hjust = 0.5, vjust =4),
        axis.text.x =  element_text(family = "Palatino", size = 12),
        axis.text.y = element_text(family = "Palatino", size = 12),
        legend.title = element_blank(),
        legend.position = "none",
        plot.margin = unit(c(1,1,2,1), "lines")) +
  geom_jitter(alpha=.2, width = .1, size = 3)

#Lydia ####
# Calculate Family richness per easement, per year
insects_rich_2019 <- all_data %>% 
  select(EasementID, RestorationCategory, Year, Sample, Family) %>% 
  filter(!is.na(EasementID)) %>% 
  filter(Year == '2019') %>% 
  group_by(EasementID, Year) %>% 
  filter(!duplicated(Family)) %>% #remove duplicates of Family/easement
  summarise(fam_rich = n()) #count the number of families/easement


insects_rich_2020 <- all_data %>% 
  select(EasementID, RestorationCategory, Year, Sample, Family) %>% 
  filter(!is.na(EasementID)) %>% 
  filter(Year == '2020') %>% 
  group_by(EasementID, Year) %>% 
  filter(!duplicated(Family)) %>% #remove duplicates of Family/easement
  summarise(fam_rich = n())  #count the number of families/easement
##combine the richness per year per easment

insect_rich <- insects_rich_2019%>%
  full_join(insects_rich_2020)

## Boxplot of insect family richness by restoration category ###
ggplot(insect_rich, aes(x=Year, y=fam_rich)) +
  geom_boxplot()


# Jade ####
# 10 Oct 2021
# Calculate abundance per Family per easement
abundance_family <- all_data %>%
  select(EasementID, RestorationCategory, Year, Sample, Total, Family) %>% 
  filter(!is.na(EasementID)) %>% 
  group_by(EasementID, RestorationCategory, Year, Family) %>% 
  summarise(abundance= sum(Total))

# Calculate abundance per easement
abundance_total<- all_data %>% 
  select(EasementID, RestorationCategory, Year, Sample, Total) %>% 
  filter(!is.na(EasementID)) %>% 
  group_by(EasementID, RestorationCategory, Year) %>% 
  summarise(abundance= sum(Total))

# plot abundance by year
abundance_total %>% 
  ggplot(aes(Year, abundance, fill = "EasementID")) +
  geom_boxplot() +
  theme_classic()+
  scale_fill_manual(values = c(palette2),
                    name = "EasementID")+
  labs(title = "Insect Abundance \n per Year")+
  xlab("\n Year") +
  ylab("Total Abundance")+
  geom_jitter(alpha=.7, width = .1, size = 3, aes(colour = EasementID))

# see if 2020 is significantly lower abundance
fit1 <- lm(abundance ~ Year + RestorationCategory, data = abundance_total)
summary(fit1) # Year explains a significant amount of variation in the data
# Need to include year as a variable in models, at least w/abundance



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

# Total family richness per easement
fam_rich <- read_csv("raw/All_Insect_Data.csv") %>% 
  select(EasementID, Date, Sample, Total, Order, Family) %>% 
  filter(!is.na(EasementID)) %>% 
  group_by(EasementID) %>% 
  filter(!duplicated(Family))%>% 
  summarise(fam_rich = n()) %>% 
  left_join(treatment)

## Boxplot of insect family richness by restoration category ####
insects_ave_rich %>% 
  ggplot(aes(RestorationCategory, ave_rich, fill=RestorationCategory)) +
  geom_boxplot(outlier.alpha = 0) +
  theme_classic()+
  scale_fill_manual(values = c(palette2),
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

# Community composition analyses ------------------------------------------
## Data Wrangling ####
# This section is where we are wrangling and cleaning the raw data in order to reformat the structure of the dataframe

##List of ALL families observed on NRCS easements##
insect_fam <- read_csv("raw/All_Insect_Data.csv") %>% 
  select(EasementID, Date, Sample, Family) %>% 
  filter(!is.na(EasementID)) %>% 
  filter(!duplicated(Family)) %>% 
  select(Family)



## Calculate number of transects ID'd per year, per easement in order to get total # of bags id'd.  We needed this number to get the average abundance of Family per site##
# Find number of bags id'd for 2019

num_trans19 <- read_csv("raw/All_Insect_Data.csv") %>% 
  select(EasementID, Date, Sample) %>%
  group_by(EasementID) %>% 
  separate(., Date, c('Month', 'Day', 'Year'), sep="/") %>% 
  filter(Year == '2019') %>% 
  filter(!duplicated(Sample)) %>% 
  summarise(Trans2019 = n())

# Find number of bags id'd for 2020
num_trans20 <- read_csv("raw/All_Insect_Data.csv") %>% 
  select(EasementID, Date, Sample) %>%
  group_by(EasementID) %>% 
  separate(., Date, c('Month', 'Day', 'Year'), sep="/") %>% 
  filter(Year == '2020') %>% 
  filter(!duplicated(Sample)) %>% 
  summarise(Trans2020 = n())

# combine dataframes into one
# this gives us the total number of samples (bags or transects) that we currently have for each easement
num_trans <- num_trans19 %>% 
  full_join(num_trans20) %>% 
  replace(is.na(.), 0) %>% 
  mutate(total_trans = Trans2019 + Trans2020) %>% #calculate total samples by adding 2019 and 2020
  select(EasementID, total_trans) #simply df

# Calculate average abundance per Family per easement
insect_ave_abundance <- all_data %>% 
  select(EasementID, RestorationCategory, Date, Sample, Total, Family) %>% 
  filter(!is.na(EasementID)) %>% 
  group_by(EasementID, RestorationCategory, Family) %>% 
  summarise(abundance= sum(Total)) %>% 
  full_join(num_trans) %>% 
  mutate(AveAbund = abundance/total_trans) %>%
  select(EasementID, RestorationCategory, Family, AveAbund)

### Site x Family Matrix ####
# Make a new dataframe where each row is a site (easement) and each column is an insect Family. 
# The number in each cell will be the average abundance for that Family

# I need to combine insect_fam + num_trans dataframes
# I need EasementID + RestorationCategory + transpose insect_fam into column. 
sitebyfam <- all_data %>% 
  select(EasementID, RestorationCategory, Date, Sample, Total, Family) %>% 
  filter(!is.na(EasementID)) %>% 
  group_by(EasementID, RestorationCategory, Family) %>% 
  summarise(abundance= sum(Total)) %>% 
  full_join(num_trans) %>% 
  mutate(AveAbund = abundance/total_trans) %>%
  select(EasementID, RestorationCategory, Family, AveAbund) %>%
  pivot_wider(names_from = Family, values_from = AveAbund) %>%
  replace(is.na(.), 0) %>%
  ungroup()

ShannonDiv <- sitebyfam %>%
  select(-EasementID, -RestorationCategory) %>%
  diversity()

#### Lydia, 10/10/21 #####
## getting shannon div per easement per year. 

sitebyfam <- all_data %>% 
  select(EasementID, RestorationCategory, Date, Sample, Total, Family) %>% 
  filter(!is.na(EasementID)) %>% 
  separate(., Date, c('Month', 'Day', 'Year'), sep="/") %>% 
  group_by(EasementID, RestorationCategory, Year, Family) %>% 
  summarise(abundance= sum(Total)) %>% 
  select(EasementID, RestorationCategory, Family, Year, abundance) %>%
  pivot_wider(names_from = Family, values_from = abundance) %>%
  replace(is.na(.), 0) %>%
  ungroup()

##Shannon div per easement per year

sitebyfam$Shannon <- diversity(sitebyfam[4:183], index = "shannon", MARGIN = 1, base = exp(1))

insect_div <- sitebyfam%>%
  select(EasementID, Year, Shannon)

ggplot(insect_div, aes(x=Year, y=Shannon)) +
  geom_boxplot()

### Vegan Code ####

#take the square root to give less to abundant spp
# NEED TO DISCUSS W/GROUP
#community.matrix <- sqrt(community.matrix)

# create a matrix with only the family and average abundance
# these need to be whole numbers for PERMANOVA to work
fam_matrix <- sitebyfam %>%
  select(-EasementID, -RestorationCategory) %>%
  mutate_all(round, 0) #round up to nearest whole number

# site (aka environmental-this is any independent variable) variables in separate df
env_matrix <- sitebyfam %>%
  dplyr::select(EasementID, RestorationCategory) ##because vegan is also loading, must indicate that we are using select from the dplyr package

# STUFF FROM JADE'S
# NOT SURE IF RELEVANT YET
#community.env <- inner_join(community.env, bombus)
#community.env <- community.env %>%
# select(c(EasementID, RestorationCategory, percent.natural, Year))

#spp.sum <- rbind(community.matrix, colSums(community.matrix))
#spp.sum$N <- rowSums(spp.sum)

# remove species that made up less than 5% of total 
#community.matrix <- community.matrix %>%
#  dplyr::select(-c(BOMAFF, BOMBOR, BOMFER, BOMCIT))

### PERMANOVA ####
perm <- adonis2(fam_matrix ~ RestorationCategory, data = env_matrix, method="bray", by="margin", permutations=9999, model = "reduced")
#model = "reduced" determines the method of permuations. 
#This permutes the residuals under a reduced model

# view results, comparing centroid of each category
perm


### Pairwise: package doesn't exist????
# pairwise permanova comparisons between restoration categories
#library(pairwiseAdonis)
#pairwise.adonis(community.matrix, community.env$RestorationCategory,sim.function = "vegdist", sim.method = "bray")

### Calculate Dispersion 
# Dispersion is basically the community composition variance
bray2 <- vegdist(fam_matrix, method = "bray")

dispersion<-betadisper(bray2, group=env_matrix$RestorationCategory)
dispersion
permutest(dispersion, pairwise = T) #p=0.554 so groups do not have different variance (aka dispersions). this model assumption is met!

## Ordination Figures: NMDS ####
## create the basic nmds ordination
nmds <- metaMDS(fam_matrix, distance = "bray")
plot(nmds) #basic plot

# fit the site data (environmental info) onto the ordination
community.envfit <- env_matrix[-c(1)]
(fit <- envfit(nmds, community.envfit, perm = 999))
head(fit)
scores(fit, "vectors")  #extracting XYs

# extract site data to plot it
env.scores <- as.data.frame(scores(fit, display = "vectors"))
env.scores <- cbind(env.scores, Species = rownames(env.scores))

plot(nmds) #basic plot again
plot(fit, col="black") #annotate the plot with restoration categories

## Extract Scores
# extract NMDS scores (x and y coordinates)
# need scores in a dateframe in order to make a nmds plot using ggplot
nmds.scores <- as.data.frame(scores(nmds))

#add identifying columns to dataframe
nmds.scores$EasementID <- env_matrix$EasementID
nmds.scores <- inner_join(nmds.scores, env_matrix)

#### Basic Plot
# create hulls
grp.control <- nmds.scores[nmds.scores$RestorationCategory == "Not Seeded", ][chull(nmds.scores[nmds.scores$RestorationCategory == 
                                                                                                  "Not Seeded", c("NMDS1", "NMDS2")]), ]  # hull values for grp A

grp.seeded <- nmds.scores[nmds.scores$RestorationCategory == "Seeded Only", ][chull(nmds.scores[nmds.scores$RestorationCategory == 
                                                                                                  "Seeded Only", c("NMDS1", "NMDS2")]), ]  # hull values for grp B

grp.burned <- nmds.scores[nmds.scores$RestorationCategory == "Seeded + Fire", ][chull(nmds.scores[nmds.scores$RestorationCategory == 
                                                                                                    "Seeded + Fire", c("NMDS1", "NMDS2")]), ]  # hull values for grp C

grp.remnant <- nmds.scores[nmds.scores$RestorationCategory == "Remnant", ][chull(nmds.scores[nmds.scores$RestorationCategory == 
                                                                                               "Remnant", c("NMDS1", "NMDS2")]), ]  # hull values for grp D

hull.data <- rbind(grp.control, grp.seeded, grp.burned, grp.remnant)  #combine groups
hull.data

#specify order of categories 
color.nmds <- factor(nmds.scores$RestorationCategory, levels = c("Not Seeded", "Seeded Only", "Seeded + Fire", "Remnant"))
color.hull <- factor(hull.data$RestorationCategory, levels = c("Not Seeded", "Seeded Only", "Seeded + Fire", "Remnant"))

# Make NMDS PLot in ggplot
p <- ggplot() + 
  geom_polygon(data=hull.data,aes(x=NMDS1,y=NMDS2, group=color.hull),alpha=0.2) + 
  geom_point(data = nmds.scores, aes(x=NMDS1,y=NMDS2, colour=color.nmds, shape = color.nmds), size = 2) + # add the point markers
  scale_fill_manual(values = palette2) +
  scale_color_manual(values = palette2) +
  scale_shape_manual(values = c(19,19,19,19)) +
  theme_bw() + 
  theme(axis.text.x = element_blank(),  # remove x-axis text
        axis.text.y = element_blank(), # remove y-axis text
        axis.ticks = element_blank(),  # remove axis ticks
        axis.title.x = element_text(size=12), # remove x-axis labels
        axis.title.y = element_text(size=12), # remove y-axis labels
        panel.background = element_blank(), 
        panel.grid.major = element_blank(),  #remove major-grid labels
        panel.grid.minor = element_blank(),  #remove minor-grid labels
        plot.background = element_blank(),
        legend.title = element_blank(),
        legend.background = element_rect(color = "black",linetype = "solid"),
        legend.position = "right")

p #view the plot


## Barchart for Relative Abundance ####
# This provides a closer look at the community composition
# TOO MANY FAMILIES TO PROPERLY SEE BARS! 
# NEED TO DISCUSS HOW TO VISUALIZE

# this example is using Order so there are fewer bars
relative_abund <- all_data %>% 
  select(EasementID, RestorationCategory, Date, Sample, Total, Order) %>% 
  filter(!is.na(EasementID)) %>% 
  group_by(EasementID, RestorationCategory, Order) %>% 
  summarise(abundance= sum(Total)) %>% 
  full_join(num_trans) %>% 
  mutate(AveAbund = abundance/total_trans) %>%
  select(EasementID, RestorationCategory, Order, AveAbund) %>%
  mutate(RestorationCategory = factor(RestorationCategory, levels = c("Not Seeded", "Seeded Only", "Seeded + Fire", "Remnant")))  #####this code reorders the treatments --the default is alphabetical

relative_abund %>% ggplot(aes(x = RestorationCategory, fill = Order)) + 
  geom_bar(position = "fill", width = 0.75, colour = "black") + 
  theme_bw() + 
  ylab("Relative Abundance") +
  theme(legend.position = "right",
        axis.text.x = element_text(size = 12),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 12)) 
















# Below is work from 10/10/2021


# Run this section every time! --------------------------------------------

library(tidyverse)
library(vegan)
# NOTE: adding four # after a note makes it into a section heading!

## Load data #####
# Upload easement treatment (aka restoration category) data & clean ##
treatment <- read_csv("raw/All_Insect_Data.csv") %>% 
  left_join(read_csv("raw/TreatmentNov20.csv")) %>% 
  select(EasementID, RestorationCategory) %>% 
  filter(!duplicated(EasementID)) %>% 
  replace(is.na(.), "Remnant") %>% 
  mutate(EasementID = replace(EasementID, EasementID == "792", "00792")) %>% 
  mutate(RestorationCategory = replace(RestorationCategory, RestorationCategory =="Seed+ Fire", "Seeded + Fire"),
         RestorationCategory = replace(RestorationCategory, RestorationCategory == "Seed", "Seeded Only"),
         RestorationCategory = replace(RestorationCategory, RestorationCategory == "No Seed", "Not Seeded")) %>% 
  
  mutate(RestorationCategory = factor(RestorationCategory, levels = c("Not Seeded", "Seeded Only", "Seeded + Fire", "Remnant")))  #####this code reorders the treatments --the default is alphabetical

# Upload raw insect sweep identification data
# This csv was downloaded from the Google Doc data entry spreadsheet on 6 March 2021
all_data <- read_csv("raw/All_Insect_Data.csv") %>% 
  select(-X1) %>% #remove random column with nothing in it
  left_join(treatment) %>% #join with treatment df
  relocate(., RestorationCategory, .before = "Date" ) 

rest_year <- read_csv("clean/enroll_rest.csv")


# Color Palettes #
palette1<-  c("#767171", "#A9D18E", "#548235", "#A49988")
palette2 <-  c("tomato3", "#A49988","#5F9EA0", "#006887")


# Raw data visualizations -------------------------------------------------



##### LYDIA EDITS

# Calculate Family richness per easement, per year
insects_rich_2019 <- all_data %>% 
  select(EasementID, RestorationCategory, Date, Sample, Family) %>% 
  filter(!is.na(EasementID)) %>% 
  separate(., Date, c('Month', 'Day', 'Year'), sep="/") %>%
  filter(Year == '2019') %>% 
  group_by(EasementID, Year) %>% 
  filter(!duplicated(Family)) %>% #remove duplicates of Family/easement
  summarise(fam_rich = n()) #count the number of families/easement

  
insects_rich_2020 <- all_data %>% 
  select(EasementID, RestorationCategory, Date, Sample, Family) %>% 
  filter(!is.na(EasementID)) %>% 
  separate(., Date, c('Month', 'Day', 'Year'), sep="/") %>%
  filter(Year == '2020') %>% 
  group_by(EasementID, Year) %>% 
  filter(!duplicated(Family)) %>% #remove duplicates of Family/easement
  summarise(fam_rich = n())  #count the number of families/easement
##combine the richness per year per easment

insect_rich <- insects_rich_2019%>%
  full_join(insects_rich_2020)

## Boxplot of insect family richness by restoration category ####
ggplot(insect_rich, aes(x=Year, y=fam_rich)) +
  geom_boxplot()







# Calculate average Family richness per easement
insects_ave_rich <- all_data %>% 
  select(EasementID, RestorationCategory, Date, Sample, Family) %>% 
  filter(!is.na(EasementID)) %>% 
  group_by(EasementID, RestorationCategory, Date, Sample) %>% 
  filter(!duplicated(Family)) %>% #remove duplicates of Family/easement
  summarise(fam_rich = n()) %>% #count the number of families/easement
  group_by(EasementID, RestorationCategory) %>% 
  summarise(ave_rich = mean(fam_rich)) #calculate average families

insects_rich_trans <- all_data %>% 
  select(EasementID, RestorationCategory, Date, Sample, Family) %>% 
  filter(!is.na(EasementID)) %>% 
  group_by(EasementID, RestorationCategory, Date, Sample) %>% 
  filter(!duplicated(Family)) %>% #remove duplicates of Family/easement
  summarise(fam_rich = n()) 

# visualize family richment per restoration age
rich_year <-  rest_year %>% 
  inner_join(insects_ave_rich) %>% 
  mutate(time_since = (2021- rest_year))


rich_year %>% 
  ggplot(aes(time_since, ave_rich, color = RestorationCategory)) +
  geom_point(outlier.alpha = 0) +
  geom_smooth(aes(color = RestorationCategory), method = "lm") +
  theme_classic()+
  scale_color_manual(values = palette2,
                     name = "Site Categories")+
  labs(title = "Average Family Richness per Restoration Year")+
  xlab("\n Years Since Restoration") +
  ylab("Ave Family Richness")+
  theme(plot.title = element_text(family = "Palatino", size = 18, hjust = 0.5, vjust = 2),
        axis.title.x = element_text(family = "Palatino", size = 15, vjust = 4),
        axis.title.y = element_text(family = "Palatino", size = 15, hjust = 0.5, vjust =4),
        axis.text.x =  element_text(family = "Palatino", size = 12),
        axis.text.y = element_text(family = "Palatino", size = 12),
        legend.title = element_blank(),
        legend.position = "none",
        plot.margin = unit(c(1,1,2,1), "lines")) +
  geom_jitter(alpha=.2, width = .1, size = 3)

# Calculate the family richness for each Order
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

# Total family richness per easement
fam_rich <- read_csv("raw/All_Insect_Data.csv") %>% 
  select(EasementID, Date, Sample, Total, Order, Family) %>% 
  filter(!is.na(EasementID)) %>% 
  group_by(EasementID) %>% 
  filter(!duplicated(Family))%>% 
  summarise(fam_rich = n()) %>% 
  left_join(treatment)

## Boxplot of insect family richness by restoration category ####
insects_ave_rich %>% 
  ggplot(aes(RestorationCategory, ave_rich, fill=RestorationCategory)) +
  geom_boxplot(outlier.alpha = 0) +
  theme_classic()+
  scale_fill_manual(values = c(palette2),
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

# Community composition analyses ------------------------------------------
## Data Wrangling ####
# This section is where we are wrangling and cleaning the raw data in order to reformat the structure of the dataframe

##List of ALL families observed on NRCS easements##
insect_fam <- read_csv("raw/All_Insect_Data.csv") %>% 
  select(EasementID, Date, Sample, Family) %>% 
  filter(!is.na(EasementID)) %>% 
  filter(!duplicated(Family)) %>% 
  select(Family)



## Calculate number of transects ID'd per year, per easement in order to get total # of bags id'd.  We needed this number to get the average abundance of Family per site##
# Find number of bags id'd for 2019

num_trans19 <- read_csv("raw/All_Insect_Data.csv") %>% 
  select(EasementID, Date, Sample) %>%
  group_by(EasementID) %>% 
  separate(., Date, c('Month', 'Day', 'Year'), sep="/") %>% 
  filter(Year == '2019') %>% 
  filter(!duplicated(Sample)) %>% 
  summarise(Trans2019 = n())

# Find number of bags id'd for 2020
num_trans20 <- read_csv("raw/All_Insect_Data.csv") %>% 
  select(EasementID, Date, Sample) %>%
  group_by(EasementID) %>% 
  separate(., Date, c('Month', 'Day', 'Year'), sep="/") %>% 
  filter(Year == '2020') %>% 
  filter(!duplicated(Sample)) %>% 
  summarise(Trans2020 = n())

# combine dataframes into one
# this gives us the total number of samples (bags or transects) that we currently have for each easement
num_trans <- num_trans19 %>% 
  full_join(num_trans20) %>% 
  replace(is.na(.), 0) %>% 
  mutate(total_trans = Trans2019 + Trans2020) %>% #calculate total samples by adding 2019 and 2020
  select(EasementID, total_trans) #simply df

# Calculate average abundance per Family per easement
insect_ave_abundance <- all_data %>% 
  select(EasementID, RestorationCategory, Date, Sample, Total, Family) %>% 
  filter(!is.na(EasementID)) %>% 
  group_by(EasementID, RestorationCategory, Family) %>% 
  summarise(abundance= sum(Total)) %>% 
  full_join(num_trans) %>% 
  mutate(AveAbund = abundance/total_trans) %>%
  select(EasementID, RestorationCategory, Family, AveAbund)




### Site x Family Matrix ####
# Make a new dataframe where each row is a site (easement) and each column is an insect Family. 
# The number in each cell will be the average abundance for that Family

# I need to combine insect_fam + num_trans dataframes
# I need EasementID + RestorationCategory + transpose insect_fam into column. 
sitebyfam <- all_data %>% 
  select(EasementID, RestorationCategory, Date, Sample, Total, Family) %>% 
  filter(!is.na(EasementID)) %>% 
  separate(., Date, c('Month', 'Day', 'Year'), sep="/") %>% 
  group_by(EasementID, RestorationCategory, Year, Family) %>% 
  summarise(abundance= sum(Total)) %>% 
  select(EasementID, RestorationCategory, Family, Year, abundance) %>%
  pivot_wider(names_from = Family, values_from = abundance) %>%
  replace(is.na(.), 0) %>%
  ungroup()

##Shannon div per easement per year

sitebyfam$Shannon <- diversity(sitebyfam[4:183], index = "shannon", MARGIN = 1, base = exp(1))

insect_div <- sitebyfam%>%
  select(EasementID, Year, Shannon)

ggplot(insect_div, aes(x=Year, y=Shannon)) +
  geom_boxplot()

#### Lydia, 10/10/21 #####
## getting shannon div per easement per year. 

sitebyfam <- all_data %>% 
  select(EasementID, RestorationCategory, Date, Sample, Total, Family) %>% 
  filter(!is.na(EasementID)) %>% 
  separate(., Date, c('Month', 'Day', 'Year'), sep="/") %>% 
  group_by(EasementID, RestorationCategory, Year, Family) %>% 
  summarise(abundance= sum(Total)) %>% 
  select(EasementID, RestorationCategory, Family, Year, abundance) %>%
  pivot_wider(names_from = Family, values_from = abundance) %>%
  replace(is.na(.), 0) %>%
  ungroup()

##Shannon div per easement per year

sitebyfam$Shannon <- diversity(sitebyfam[4:183], index = "shannon", MARGIN = 1, base = exp(1))

insect_div <- sitebyfam%>%
  select(EasementID, Year, Shannon)

ggplot(insect_div, aes(x=Year, y=Shannon)) +
  geom_boxplot()










### Vegan Code ####

#take the square root to give less to abundant spp
# NEED TO DISCUSS W/GROUP
#community.matrix <- sqrt(community.matrix)

# create a matrix with only the family and average abundance
# these need to be whole numbers for PERMANOVA to work
fam_matrix <- sitebyfam %>%
  select(-EasementID, -RestorationCategory) %>%
  mutate_all(round, 0) #round up to nearest whole number

# site (aka environmental-this is any independent variable) variables in separate df
env_matrix <- sitebyfam %>%
  dplyr::select(EasementID, RestorationCategory) ##because vegan is also loading, must indicate that we are using select from the dplyr package

# STUFF FROM JADE'S
# NOT SURE IF RELEVANT YET
#community.env <- inner_join(community.env, bombus)
#community.env <- community.env %>%
# select(c(EasementID, RestorationCategory, percent.natural, Year))

#spp.sum <- rbind(community.matrix, colSums(community.matrix))
#spp.sum$N <- rowSums(spp.sum)

# remove species that made up less than 5% of total 
#community.matrix <- community.matrix %>%
#  dplyr::select(-c(BOMAFF, BOMBOR, BOMFER, BOMCIT))

### PERMANOVA ####
perm <- adonis2(fam_matrix ~ RestorationCategory, data = env_matrix, method="bray", by="margin", permutations=9999, model = "reduced")
#model = "reduced" determines the method of permuations. 
#This permutes the residuals under a reduced model

# view results, comparing centroid of each category
perm


### Pairwise: package doesn't exist????
# pairwise permanova comparisons between restoration categories
#library(pairwiseAdonis)
#pairwise.adonis(community.matrix, community.env$RestorationCategory,sim.function = "vegdist", sim.method = "bray")

### Calculate Dispersion 
# Dispersion is basically the community composition variance
bray2 <- vegdist(fam_matrix, method = "bray")

dispersion<-betadisper(bray2, group=env_matrix$RestorationCategory)
dispersion
permutest(dispersion, pairwise = T) #p=0.554 so groups do not have different variance (aka dispersions). this model assumption is met!

## Ordination Figures: NMDS ####
## create the basic nmds ordination
nmds <- metaMDS(fam_matrix, distance = "bray")
plot(nmds) #basic plot

# fit the site data (environmental info) onto the ordination
community.envfit <- env_matrix[-c(1)]
(fit <- envfit(nmds, community.envfit, perm = 999))
head(fit)
scores(fit, "vectors")  #extracting XYs

# extract site data to plot it
env.scores <- as.data.frame(scores(fit, display = "vectors"))
env.scores <- cbind(env.scores, Species = rownames(env.scores))

plot(nmds) #basic plot again
plot(fit, col="black") #annotate the plot with restoration categories

## Extract Scores
# extract NMDS scores (x and y coordinates)
# need scores in a dateframe in order to make a nmds plot using ggplot
nmds.scores <- as.data.frame(scores(nmds))

#add identifying columns to dataframe
nmds.scores$EasementID <- env_matrix$EasementID
nmds.scores <- inner_join(nmds.scores, env_matrix)

#### Basic Plot
# create hulls
grp.control <- nmds.scores[nmds.scores$RestorationCategory == "Not Seeded", ][chull(nmds.scores[nmds.scores$RestorationCategory == 
                                                                                                  "Not Seeded", c("NMDS1", "NMDS2")]), ]  # hull values for grp A

grp.seeded <- nmds.scores[nmds.scores$RestorationCategory == "Seeded Only", ][chull(nmds.scores[nmds.scores$RestorationCategory == 
                                                                                                  "Seeded Only", c("NMDS1", "NMDS2")]), ]  # hull values for grp B

grp.burned <- nmds.scores[nmds.scores$RestorationCategory == "Seeded + Fire", ][chull(nmds.scores[nmds.scores$RestorationCategory == 
                                                                                                    "Seeded + Fire", c("NMDS1", "NMDS2")]), ]  # hull values for grp C

grp.remnant <- nmds.scores[nmds.scores$RestorationCategory == "Remnant", ][chull(nmds.scores[nmds.scores$RestorationCategory == 
                                                                                               "Remnant", c("NMDS1", "NMDS2")]), ]  # hull values for grp D

hull.data <- rbind(grp.control, grp.seeded, grp.burned, grp.remnant)  #combine groups
hull.data

#specify order of categories 
color.nmds <- factor(nmds.scores$RestorationCategory, levels = c("Not Seeded", "Seeded Only", "Seeded + Fire", "Remnant"))
color.hull <- factor(hull.data$RestorationCategory, levels = c("Not Seeded", "Seeded Only", "Seeded + Fire", "Remnant"))

# Make NMDS PLot in ggplot
p <- ggplot() + 
  geom_polygon(data=hull.data,aes(x=NMDS1,y=NMDS2, group=color.hull),alpha=0.2) + 
  geom_point(data = nmds.scores, aes(x=NMDS1,y=NMDS2, colour=color.nmds, shape = color.nmds), size = 2) + # add the point markers
  scale_fill_manual(values = palette2) +
  scale_color_manual(values = palette2) +
  scale_shape_manual(values = c(19,19,19,19)) +
  theme_bw() + 
  theme(axis.text.x = element_blank(),  # remove x-axis text
        axis.text.y = element_blank(), # remove y-axis text
        axis.ticks = element_blank(),  # remove axis ticks
        axis.title.x = element_text(size=12), # remove x-axis labels
        axis.title.y = element_text(size=12), # remove y-axis labels
        panel.background = element_blank(), 
        panel.grid.major = element_blank(),  #remove major-grid labels
        panel.grid.minor = element_blank(),  #remove minor-grid labels
        plot.background = element_blank(),
        legend.title = element_blank(),
        legend.background = element_rect(color = "black",linetype = "solid"),
        legend.position = "right")

p #view the plot


## Barchart for Relative Abundance ####
# This provides a closer look at the community composition
# TOO MANY FAMILIES TO PROPERLY SEE BARS! 
# NEED TO DISCUSS HOW TO VISUALIZE

# this example is using Order so there are fewer bars
relative_abund <- all_data %>% 
  select(EasementID, RestorationCategory, Date, Sample, Total, Order) %>% 
  filter(!is.na(EasementID)) %>% 
  group_by(EasementID, RestorationCategory, Order) %>% 
  summarise(abundance= sum(Total)) %>% 
  full_join(num_trans) %>% 
  mutate(AveAbund = abundance/total_trans) %>%
  select(EasementID, RestorationCategory, Order, AveAbund) %>%
  mutate(RestorationCategory = factor(RestorationCategory, levels = c("Not Seeded", "Seeded Only", "Seeded + Fire", "Remnant")))  #####this code reorders the treatments --the default is alphabetical

relative_abund %>% ggplot(aes(x = RestorationCategory, fill = Order)) + 
  geom_bar(position = "fill", width = 0.75, colour = "black") + 
  theme_bw() + 
  ylab("Relative Abundance") +
  theme(legend.position = "right",
        axis.text.x = element_text(size = 12),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 12)) 

















#playing with NMDS - Lydia

stressplot(nmds)
stressplot(nmds2)

#adding a 3rd dimension?

nmds2 <- metaMDS(fam_matrix, distance = "bray", k=3)
plot(nmds2)

community.envfit <- env_matrix[-c(1)]
(fit <- envfit(nmds2, community.envfit, perm = 999))
head(fit)
scores(fit, "vectors")  #extracting XYs

# extract site data to plot it
env.scores <- as.data.frame(scores(fit, display = "vectors"))
env.scores <- cbind(env.scores, Species = rownames(env.scores))

plot(nmds2) #basic plot again
plot(fit, col="black") #annotate the plot with restoration categories

nmds.scores2 <- as.data.frame(scores(nmds2))


#add identifying columns to dataframe
nmds.scores2$EasementID <- env_matrix$EasementID
nmds.scores2 <- inner_join(nmds.scores2, env_matrix)

grp.control <- nmds.scores2[nmds.scores2$RestorationCategory == "Not Seeded", ][chull(nmds.scores2[nmds.scores2$RestorationCategory == 
                                                                                                     "Not Seeded", c("NMDS1", "NMDS2")]), ]  # hull values for grp A

grp.seeded <- nmds.scores2[nmds.scores2$RestorationCategory == "Seeded Only", ][chull(nmds.scores2[nmds.scores2$RestorationCategory == 
                                                                                                     "Seeded Only", c("NMDS1", "NMDS2")]), ]  # hull values for grp B

grp.burned <- nmds.scores2[nmds.scores2$RestorationCategory == "Seeded + Fire", ][chull(nmds.scores2[nmds.scores2$RestorationCategory == 
                                                                                                       "Seeded + Fire", c("NMDS1", "NMDS2")]), ]  # hull values for grp C

grp.remnant <- nmds.scores2[nmds.scores2$RestorationCategory == "Remnant", ][chull(nmds.scores2[nmds.scores2$RestorationCategory == 
                                                                                                  "Remnant", c("NMDS1", "NMDS2")]), ]  # hull values for grp D

hull.data <- rbind(grp.control, grp.seeded, grp.burned, grp.remnant)  #combine groups
hull.data

color.nmds <- factor(nmds.scores2$RestorationCategory, levels = c("Not Seeded", "Seeded Only", "Seeded + Fire", "Remnant"))
color.hull <- factor(hull.data$RestorationCategory, levels = c("Not Seeded", "Seeded Only", "Seeded + Fire", "Remnant"))

# Make NMDS PLot in ggplot
p2 <- ggplot() + 
  geom_polygon(data=hull.data,aes(x=NMDS1,y=NMDS2, group=color.hull),alpha=0.2) + 
  geom_point(data = nmds.scores2, aes(x=NMDS1,y=NMDS2, colour=color.nmds, shape = color.nmds), size = 2) + # add the point markers
  scale_fill_manual(values = palette2) +
  scale_color_manual(values = palette2) +
  scale_shape_manual(values = c(19,19,19,19)) +
  theme_bw() + 
  theme(axis.text.x = element_blank(),  # remove x-axis text
        axis.text.y = element_blank(), # remove y-axis text
        axis.ticks = element_blank(),  # remove axis ticks
        axis.title.x = element_text(size=12), # remove x-axis labels
        axis.title.y = element_text(size=12), # remove y-axis labels
        panel.background = element_blank(), 
        panel.grid.major = element_blank(),  #remove major-grid labels
        panel.grid.minor = element_blank(),  #remove minor-grid labels
        plot.background = element_blank(),
        legend.title = element_blank(),
        legend.background = element_rect(color = "black",linetype = "solid"),
        legend.position = "right")

p2 #view the plot
