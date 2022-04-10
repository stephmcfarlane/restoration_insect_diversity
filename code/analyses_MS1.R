#Load packages####
library(iNEXT)
library(tidyverse)
library(hablar)
library(vegan)
library(devtools)

# Color palettes ####
palette3<-  c("#767171", "#A9D18E", "#548235", "#2A4117") #3B5422

# Site History Data ####
# Upload easement treatment (aka restoration category) data & clean ##

### Restoration History ####
rest_history <- read_csv("raw/restoration_history.csv") %>% 
  select(EasementID, RestorationCategory = RestorationCategory_NEW, Size = Contiguous_Seeding_Size) %>% 
  mutate(EasementID = recode(EasementID, '792' = "00792")) %>%
  mutate(RestorationCategory = factor(RestorationCategory, levels = c("No Seed", "Seed", "Seed + Fire", "Remnant")))

## Restoration Age
age <- read_csv("clean/rest_age.csv") %>% 
  select(EasementID = SiteID, rest_year)

##Fire History####
fire_years  <- read_csv("raw/restoration_history.csv") %>% 
  select(SiteID = EasementID, RestorationCategory = RestorationCategory_NEW, Fire_Years) %>%
  mutate(SiteID = recode(SiteID, '792' = "00792")) %>%
  filter(RestorationCategory != "Seed") %>% 
  filter(RestorationCategory != "No Seed") %>% 
  separate(Fire_Years, c("Fire1", "Fire2", "Fire3")) %>%
  left_join(age) %>% 
  mutate(sample_year = '2020') %>% 
  select(-rest_year, -rest_age)  %>% 
  mutate(Fire3 = ifelse(Fire3 > sample_year, NA, Fire3)) %>% 
  mutate(Fire2 = ifelse(Fire2 > sample_year, NA, Fire2)) %>% 
  mutate(Fire1 = ifelse(Fire1 > sample_year, NA, Fire1))

###REMOVE FIRE FREQ####
fire_frequency <- fire_years %>% 
  pivot_longer(cols = c("Fire1", "Fire2", "Fire3"), values_to = "Fires") %>% 
  filter(!is.na(Fires)) %>% 
  group_by(SiteID) %>% 
  summarise(number_fires = n()) %>% 
  left_join(age) %>% 
  mutate(fire_frequency = (rest_age/number_fires)) 


time_since_fire <- fire_years %>% 
  mutate_at(vars(c(Fire1, Fire2, Fire3)), ~replace(., is.na(.), "")) %>% 
  mutate(last_fire = ifelse(Fire3 > 0, Fire3, Fire2)) %>% 
  mutate(last_fire = ifelse(last_fire > 0, last_fire, Fire1)) %>% 
  convert(int(last_fire)) %>% 
  mutate(years_since_fire = (sample_year - last_fire))

#Insect Data ####
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
  filter(EasementID != "Hauser Road") %>% 
  filter(EasementID != "UW-Arboretum") %>% 
  filter(EasementID != "Snapper Prairie") %>% 
  filter(EasementID !="00MDP"|Year!="2020") ##Need to filter out 00MDP 2020, but not 00MDP 2019

#site 00VTR may be an outlier in the data...

rest_year <- read_csv("clean/enroll_rest.csv") 




# Normality ----
###### LD,GB,GW checking for normality
######10/17/2021

# Abundance ####
# get total abundance per site, per year
abundance_by_year <- all_data %>%
  select(EasementID, RestorationCategory, Year, Sample, Total, Family) %>% 
  filter(!is.na(EasementID)) %>% 
  group_by(EasementID, Year) %>% 
  summarise(abundance= sum(Total))

mean1<-mean(abundance_by_year$abundance)
abundance_by_year$residual <- abundance_by_year$abundance-mean1
shapiro.test(abundance_by_year$residual)
ggplot(data = abundance_by_year) + geom_histogram(mapping = aes(x = residual), bins = 6)

#abundance has normally distributed residuals

#Diversity metrics ####
### Species Richness ####
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
  filter(Year == '2020') %>% 
  group_by(EasementID, Year) %>% 
  filter(!duplicated(Family)) %>% #remove duplicates of Family/easement
  summarise(fam_rich = n())  #count the number of families/easement
##combine the richness per year per easment

##This is not averaging family richness--it is just removing duplicates which
##which will result in higher family richness at sites visited two years in a row
insect_rich<- all_data %>% 
 select(EasementID, RestorationCategory, Year, Sample, Family) %>% 
 group_by(EasementID) %>% 
 filter(!duplicated(Family)) %>% #remove duplicates of Family/easement
 summarise(fam_rich = n()) %>% 
  left_join(rest_history)

## can't do this...need to remove on year randomly
insect_rich<- all_data %>% 
  select(EasementID, RestorationCategory, Year, Sample, Family) %>% 
  group_by(EasementID, Year) %>% 
  filter(!duplicated(Family)) %>% #remove duplicates of Family/easement
  summarise(fam_rich = n()) %>% 
  ungroup() %>% 
  group_by(EasementID) %>% 
  sum(ave_rich = mean(fam_rich))

####Plotting richness ####
insect_rich %>% 
  ggplot(aes(RestorationCategory, fam_rich, fill=RestorationCategory)) +
  geom_boxplot(outlier.alpha = 0) +
  theme_classic()+
  scale_fill_manual(values = palette3,
                    name = "Site Categories")+
  labs(title = "Insect Family Richness \n per Site Category")+
  xlab("\n Site category") +
  ylab("Average family richness")+
  theme(plot.title = element_text(family = "Palatino", size = 18, hjust = 0.5, vjust = 2),
        axis.title.x = element_text(family = "Palatino", size = 15, vjust = 4),
        axis.title.y = element_text(family = "Palatino", size = 15, hjust = 0.5, vjust =4),
        axis.text.x =  element_text(family = "Palatino", size = 12),
        axis.text.y = element_text(family = "Palatino", size = 12),
        legend.title = element_blank(),
        legend.position = "none",
        plot.margin = unit(c(1,1,2,1), "lines")) +
  geom_jitter(alpha=.2, width = .1, size = 3)

## 8 sites have species richness for two years...how should we deal with that?
insect_rich_combined <- insects_rich_2019%>%
  full_join(insects_rich_2020) %>% 
  

mean2<-mean(insect_rich$fam_rich)
insect_rich$residual <- insect_rich$fam_rich-mean2
shapiro.test(insect_rich$residual)
ggplot(data = insect_rich) + geom_histogram(mapping = aes(x = residual), bins = 6)

#richness just barely has normally distributed residuals

### Shannon's Diversity Index ####
sitebyfam <- all_data %>% 
  select(EasementID, Year, Sample, Total, Family) %>% 
  group_by(EasementID, Year, Family) %>% 
  summarise(abundance= sum(Total)) %>% 
  select(EasementID, Family, Year, abundance) %>%
  pivot_wider(names_from = Family, values_from = abundance) %>%
  replace(is.na(.), 0) %>%
  ungroup() %>% 
  unite(., "x", EasementID, Year, sep = "_") %>% ## joining columns and 
  column_to_rownames(var = "x")

ShannonDiv <- sitebyfam %>% 
  unite(., "x", EasementID, Year, sep = "_") %>% ## joining columns and 
  column_to_rownames(var = "x") %>% ##making column rownames so matrix is just numbers
  diversity() %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "rowname") %>% 
  separate(., rowname, c("EasementID", "Year"), sep = "_") %>% 
  rename(Shannon = ".")

##Shannon div per easement per year
#sitebyfam$Shannon <- diversity(sitebyfam[4:182], index = "shannon", MARGIN = 1, base = exp(1))
#insect_div <- sitebyfam%>%
 # select(EasementID, Year, Shannon) %>% 
 # left_join(rest_history) %>% 
 # filter(EasementID != "00VTR")

#### Plotting Diversity ####
insect_div %>% 
  ggplot(aes(RestorationCategory, Shannon, fill=RestorationCategory)) +
  geom_boxplot(outlier.alpha = 0) +
  theme_classic()+
  scale_fill_manual(values = palette3,
                    name = "Site Categories")+
  labs(title = "Insect Family Diversity \n per Site Category")+
  xlab("\n Site Category") +
  ylab("Average family diversity")+
  theme(plot.title = element_text(family = "Palatino", size = 18, hjust = 0.5, vjust = 2),
        axis.title.x = element_text(family = "Palatino", size = 15, vjust = 4),
        axis.title.y = element_text(family = "Palatino", size = 15, hjust = 0.5, vjust =4),
        axis.text.x =  element_text(family = "Palatino", size = 12),
        axis.text.y = element_text(family = "Palatino", size = 12),
        legend.title = element_blank(),
        legend.position = "none",
        plot.margin = unit(c(1,1,2,1), "lines")) +
  geom_jitter(alpha=.2, width = .1, size = 3)

mean3<-mean(insect_div$Shannon)
insect_div$residual <- insect_div$Shannon-mean3
shapiro.test(insect_div$residual)
ggplot(data = insect_div) + geom_histogram(mapping = aes(x = residual), bins = 6)

#diversity does not have normally distributed residuals

## iNEXT: Interpolation Exterpolation####

#### Formatting Data into a List of Lists ####
plot_ids <- unique(all_data[c("EasementID", "Year")])  ##creates unique list of plot-year combo

# output list of lists, start with blank list
abundance_list <-  list()

list_names <- mutate(plot_ids, paste0(EasementID, "_", Year)) %>%
  select(names = 3) %>% ##naming 3rd column names
  pull() ##changing to a vector of names which is needed to rename lists in abundance list

## for loop -- converting species' abundances into a list for each site/date combo and adding to blank list, abundance_list
for(i in 1:nrow(plot_ids)) 
  {
  df.out <- all_data %>% 
    select(EasementID, RestorationCategory, Year, Sample, Total, Family) %>% 
    group_by(EasementID, RestorationCategory, Year, Family) %>% 
    summarise(abundance = sum(Total)) %>%
    filter(EasementID==plot_ids[i,1]) %>% ##filtering easement by the site name (column 1)
    filter(Year==plot_ids[i,2]) %>%  #& by year (column 2)
    ungroup() %>%
    select(abundance) %>%
    as.list()
  
  abundance_list <- append(abundance_list, df.out)
  }

names(abundance_list) <- list_names ## assigning names to lists

#### iNEXT calculations ####
##setting sample sizes for R/E estimations
insect <- iNEXT(abundance_list, q=c(0, 1), datatype = "abundance", nboot = 100) ##can set endpoint too, referring to number of indivduals 

## output into separate data frames
insect_abundance <- as.data.frame(insect[["DataInfo"]])  ##site name, sample size, Family richness, number of singletons, doubletons, etc.
#insect_div_est <- as.data.frame(insect[["iNextEst"]]) ## diversity estimates with rarefied and extrapolated samples for q = 0,1,2
insect_div_asy <- as.data.frame(insect[["AsyEst"]]) ## asympototic diversity estimates with bootstrap s.e. and confidence intervals for q = 0,1,2

write_csv(insect$DataInfo, "clean/insect_datainfo.csv")
write_csv(insect$AsyEst, "clean/insect_div_asym.csv")
write_csv(insect_div_est, "clean/insect_div_est.csv")

## abundance
insect_abund <- insect_datainfo %>% 
  separate(., site, c("EasementID", "Year"), sep = "_")
  
##determining estimated species richness at N=450
insect_rich <- iNEXT(abundance_list, q=0, datatype = "abundance", endpoint = 450, knots= 100, nboot = 100) 

insect_rich_est <- as.data.frame(insect_rich[["iNextEst"]])

##determining estimated species diversity with enpoint = 2500
insect_shannon <- iNEXT(abundance_list, q=1, datatype = "abundance", knots = 100) 
insect_shannon_est <- as.data.frame(insect_shannon[["iNextEst"]]) 

## Cleaning spp richness estimates ####
rich_est <- insect_rich_est %>% 
  select(ends_with(".qD")) %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "ID") %>% 
  #separate(., rowname, c("EasementID", "Year"), sep = "_") %>% 
  mutate_at("ID", str_replace, ".qD", "") %>% 
  mutate_at("ID", str_replace_all, '\\.', " ") %>% 
  mutate(ID = if_else(str_starts(ID, "X"), str_replace(ID,"X", ""), ID)) %>% 
  select(ID,Rich_est=V100)

shannon_est <- insect_div_asy %>% 
  filter(Diversity == 'Shannon diversity') %>% 
  select(Site, Shannon_iNext = Observed, Shannon_asy = Estimator)
#%>% 
#  separate(., Site, c("EasementID", "Year"), sep = "_")
## Cleaning up diveristy estimates
#shannon_est <- insect_shannon_est %>% 
#  select(ends_with(".qD")) %>% 
#  t() %>% 
#  as.data.frame() %>% 
#  rownames_to_column(var = "rowname") %>% 
#  separate(., rowname, c("EasementID", "Year"), sep = "_") %>% 
#  mutate_at("Year", str_replace, ".qD", "") %>% 
#  mutate_at("EasementID", str_replace_all, '\\.', " ") %>% 
#  mutate(EasementID = if_else(str_starts(EasementID, "X"), str_replace(EasementID,"X", ""), EasementID)) %>% 
#  rename(Shannon_est=V1)

#All Variables ####
## Must remove sites with abundance less than 225 when modeling estimated richness #filter(Abund > 200)
insect_response <- insect_abundance %>% 
  select(ID = site, Abund = n, Rich_obs = S.obs) %>% 
  left_join(rich_est) %>% 
  left_join(shannon_est %>% rename(ID=Site)) %>% 
  separate(., ID, c("EasementID", "Sample_year"), sep = "_") %>%  
  left_join(rest_history) %>% 
  left_join(age) %>% 
  convert(int(Sample_year, rest_year)) %>% 
  mutate(rest_age =  (Sample_year - rest_year)) %>% 
  select(-rest_year)

write_csv(insect_response, "clean/insect_response.csv")

#Removing random year####  
insect_response_YR <- insect_response %>% 
  unite(., "x", EasementID, Sample_year, sep = "_") %>% 
  filter(x != "0078S_2020") %>% 
  filter(x !="007C6_2020") %>% 
  filter(x !="007XH_2019") %>% 
  filter(x !="00MBG_2020") %>% 
  filter(x !="00MJC_2020") %>% 
  filter(x !="00NCQ_2019") %>% 
  filter(x !="00WV0_2020") %>% 
  separate(., x, c("EasementID", "Sample_year"), sep = "_")

##Random number to remove year, duplicate site with largest random number kept, smaller number removed   
year_ran_num <- insect_response %>% 
  mutate(random_num = sample(c(1:37))) %>% 
  relocate(., random_num, .before = "EasementID") %>% 
  select(EasementID, Sample_year, random_num, Rich_est,Shannon_asy)
write_csv(year_ran_num, "clean/year_random_number.csv")


##time since fire code
  

##Plotting iNEXT example with remanant sites...not really useful with 40 sites
remnant <- all_data %>% 
  filter(RestorationCategory == "Remnant")

rem_plot_ids <- unique(remnant[c("EasementID", "Year")])

rem_list_names <- mutate(rem_plot_ids, paste0(EasementID, "_", Year)) %>%
  select(names = 3) %>%
  pull()
rem_abundance_list <-  list()

for(i in 1:nrow(rem_plot_ids)) {
  
  df.out <- remnant %>% 
    select(EasementID, RestorationCategory, Year, Sample, Total, Family) %>% 
    group_by(EasementID, RestorationCategory, Year, Family) %>% 
    summarise(abundance= sum(Total)) %>%
    filter(EasementID==rem_plot_ids[i,1]) %>% 
    filter(Year==rem_plot_ids[i,2]) %>%
    ungroup() %>%
    select(abundance) %>%
    as.list()
  #pivot_wider(names_from = Family, values_from = abundance) %>%
  #split(1:nrow(.))  #split(df.out, 1:nrow(df.out))
  
  rem_abundance_list <- append(rem_abundance_list, df.out)
}

names(rem_abundance_list) <- rem_list_names

rem_abundance_list <- rem_abundance_list[c(-3, -5)]  ##removing Faville and Oliver prairies, due to extremely low abundance

remnant_insect <- iNEXT(rem_abundance_list, q=c(0, 1, 2), datatype = "abundance", endpoint = 1500) 

write_csv(remnant_insect$DataInfo, "clean/remnant_insect_data_info.csv")
write_csv(remnant_insect$AsyEst, "clean/remnant_insect_AsyEst.csv")

ggiNEXT(remnant_insect, type=1, facet.var="site")  
# Sample-size-based R/E curves, separating plots by "order" 
ggiNEXT(remnant_insect, type=1, facet.var="order", color.var="site") +
  theme_bw() +
  scale_color_manual(values = c("#A9D18E", "#548235", "#767171", "#A9D18E", "#548235", "#2A4117", "#2A4117")) +
  scale_fill_manual(values = c("#A9D18E", "#548235", "#767171", "#A9D18E", "#548235", "#2A4117", "#2A4117"))





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
 
 
# Jade ###
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
insect_response %>% 
   ggplot(aes(Year, Abund, fill = "EasementID")) +
   geom_boxplot() +
   theme_classic()+
   scale_fill_manual(values = c(palette3),
                     name = "EasementID")+
   labs(title = "Insect Abundance \n per Year")+
   xlab("\n Year") +
   ylab("Total Abundance")+
   geom_jitter(alpha=.7, width = .1, size = 3, aes(colour = EasementID))
 
# see if 2020 is significantly lower abundance
fit1 <- lm(abundance ~ Year + RestorationCategory, data = abundance_total)
summary(fit1) # Year explains a significant amount of variation in the data
# Need to include year as a variable in models, at least w/abundance
 
 
 

# Total family richness per easement
fam_rich <- read_csv("raw/All_Insect_Data.csv") %>% 
  select(EasementID, Date, Sample, Total, Order, Family) %>% 
  filter(!is.na(EasementID)) %>% 
  group_by(EasementID) %>% 
  filter(!duplicated(Family))%>% 
  summarise(fam_rich = n()) %>% 
  left_join(treatment)

# Plotting diversity by restoration category ####
## Estimated richness ####
plot_rich <- insect_response_YR %>% filter(Abund > 200) 
plot_rich %>% 
  ggplot(aes(RestorationCategory, Rich_est, fill=RestorationCategory)) +
  geom_boxplot(outlier.alpha = 0) +
  theme_classic()+
  scale_fill_manual(values = c(palette3),
                    name = "Site Categories")+
  labs(title = "Interpolated and Extrapolated Insect\n Richness per Restoration Category")+
  xlab("\n Site Category") +
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

## Estimated diversity ####
plot_div <- insect_response_YR 
plot_div %>% 
  ggplot(aes(RestorationCategory, Shannon_asy, fill=RestorationCategory)) +
  geom_boxplot(outlier.alpha = 0) +
  theme_classic()+
  scale_fill_manual(values = c(palette3),
                    name = "Site Categories")+
  labs(title = "Extrapolated Insect Diversity\n per Restoration Category")+
  xlab("\n Site Category") +
  ylab("Average Family Diversity")+
  theme(plot.title = element_text(family = "Palatino", size = 18, hjust = 0.5, vjust = 2),
        axis.title.x = element_text(family = "Palatino", size = 15, vjust = 4),
        axis.title.y = element_text(family = "Palatino", size = 15, hjust = 0.5, vjust =4),
        axis.text.x =  element_text(family = "Palatino", size = 12),
        axis.text.y = element_text(family = "Palatino", size = 12),
        legend.title = element_blank(),
        legend.position = "none",
        plot.margin = unit(c(1,1,2,1), "lines")) +
  geom_jitter(alpha=.2, width = .1, size = 3)



# Model  ANOVA & Tukey----------------------------------------
insect_rest <-  insect_response %>% 
  filter(RestorationCategory != "No Seed")

anova <- aov(Rich_est ~ RestorationCategory, data=plot_div)
summary(anova)

anova <- aov(Shannon_asy ~ RestorationCategory + Size, data = plot_rich)
#there is no signifcant variation of insect diversity between restoration categories, nor collection year
#Restoration category = 0.0572
#year = 0.2471

TukeyHSD(anova)
#this gives intersting results for the differences in the means

sitebyfam %>% #visualizing restoration category and insect diversity
  ggplot(aes(RestorationCategory, Shannon, fill = "EasementID")) +
  geom_boxplot() +
  theme_classic()+
  xlab("Restoration Category") +
  ylab("Insect Diversity")


# Ordination & PERMANOVA --------------------------------------------------
## PCA ####


insect_pca <-  princomp(sitebyfam[,-1], cor=TRUE)
sitebyfam2 <- all_data %>% #matrix with only numeric values
  select(EasementID, RestorationCategory, Year, Sample, Total, Family) %>% #note that this is including year, so a site can have up to 2 entries
  filter(!is.na(EasementID)) %>% 
  group_by(EasementID, RestorationCategory, Year, Family) %>% 
  summarise(abundance= sum(Total)) %>% 
  select(EasementID, RestorationCategory, Family, Year, abundance) %>%
  pivot_wider(names_from = Family, values_from = abundance) %>%
  replace(is.na(.), 0) %>%
  ungroup()%>% 
  select(-Year, -RestorationCategory, -EasementID)


#the model statement for the capscale
  #model.cap <-capscale(dataonly ~ group, groupings, dist="bray")

div_can <- capscale(sitebyfam2 ~ RestorationCategory, sitebyfam, dist="bray")
div_can
  #creating the ordination

#plotting the ordination
ordiplot(div_can,type="n", main = "Insect Community")
orditorp(div_can, display="sites", priority, labels= sitebyfam$RestorationCategory,)
ordihull(div_can, groups= sitebyfam$RestorationCategory, label= TRUE, draw= "polygon")


div_permanova <- adonis(formula = sitebyfam2 ~ RestorationCategory, data = sitebyfam, permutations = 999, method ="bray")
div_permanova
  #permanova on the ordination: testing if there are differences in the dispersion/spread of groups
    #does not assume normality (non-parametric)


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

#### Lydia, 10/10/21 ###

#shannon div per easement per year. ----

sitebyfam <- all_data %>% 
  select(EasementID, Year, Sample, Total, Family) %>% 
  filter(!is.na(EasementID)) %>% 
  group_by(EasementID, Year, Family) %>% 
  summarise(abundance= sum(Total)) %>% 
  select(EasementID,Family, Year, abundance) %>%
  pivot_wider(names_from = Family, values_from = abundance) %>%
  ungroup() %>% 
  replace(is.na(.), 0) %>%
  unite(., "x", EasementID, Year, sep = "_") %>% 
  column_to_rownames(var = "x")

##Shannon div per easement per year

sitebyfam$Shannon <- diversity(sitebyfam[3:181], index = "shannon", MARGIN = 1, base = exp(1))

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



# Lydia's Rarefied richness  Code---------------------------------
# rarefaction standardizes the sample sizes for the richness based on the smallest sample size

#2020
sitebyfam20 <- all_data %>% 
  select(EasementID, RestorationCategory, Year, Sample, Total, Family) %>% 
  filter(!is.na(EasementID)) %>% 
  group_by(EasementID, RestorationCategory, Year, Family) %>% 
  summarise(abundance= sum(Total)) %>% 
  filter(Year == '2020') %>% 
  select(EasementID, RestorationCategory, Family, Year, abundance) %>%
  pivot_wider(names_from = Family, values_from = abundance) %>%
  replace(is.na(.), 0) %>%
  ungroup()
sitebyfam20$Shannon <- diversity(sitebyfam20[4:182], index = "shannon", MARGIN = 1, base = exp(1))

library(janitor)
fambysite20 <- data.frame(t(sitebyfam20)) %>% 
  row_to_names(row_number = 1)
fambysite20 <-   fambysite20[-c(1,2),]
rowremove <- c("Shannon")
fambysite20 <- fambysite20[!(row.names(fambysite20)%in%rowremove),]
fambysite20 <- mutate_all(fambysite20,function(x)as.numeric(as.character(x)))

# the lowest abundance for 2020 was 60
rare_richness20 <- rarefy(fambysite20, 60, MARGIN = 2)
rare_richness20
rarerich20 <- data.frame(rare_richness20)

#2019
sitebyfam19 <- all_data %>% 
  select(EasementID, RestorationCategory, Year, Sample, Total, Family) %>% 
  filter(!is.na(EasementID)) %>% 
  group_by(EasementID, RestorationCategory, Year, Family) %>% 
  summarise(abundance= sum(Total)) %>% 
  filter(Year == '2019') %>% 
  select(EasementID, RestorationCategory, Family, Year, abundance) %>%
  pivot_wider(names_from = Family, values_from = abundance) %>%
  replace(is.na(.), 0) %>%
  ungroup()
sitebyfam19$Shannon <- diversity(sitebyfam19[4:182], index = "shannon", MARGIN = 1, base = exp(1))


fambysite19 <- data.frame(t(sitebyfam19)) %>% 
  row_to_names(row_number = 1)
fambysite19 <-   fambysite19[-c(1,2),]
rowremove <- c("Shannon")
fambysite19 <- fambysite19[!(row.names(fambysite19)%in%rowremove),]
fambysite19 <- mutate_all(fambysite19,function(x)as.numeric(as.character(x)))

# the lowest abundance in 2019 was 366
rare_richness19 <- rarefy(fambysite19, 366, MARGIN = 2)
rare_richness19
rarerich19 <- data.frame(rare_richness19)

library(dplyr)
rarerich19 <- tibble::rownames_to_column(rarerich19, "VALUE")
colnames(rarerich19)[colnames(rarerich19) == 'VALUE'] <- 'EasementID'
rarerich19 <-  rarerich19%>% 
  left_join(insects_rich_2019)
attach(rarerich19)
plot(rare_richness19, fam_rich, xlab = "fam_rich", ylab = "rare_richness19")
abline(0, 1)
detach(rarerich19)

rarecurve(rarerich19, step = 20, sample = 366, col = "blue", cex = 0.6)

#When comparing, must do within year, not across years, because they are standardized differently. 


#checking the residuals for each year's rarefied richness

#2019
mean4<-mean(rarerich19$rare_richness19)
rarerich19$residual <- rarerich19$rare_richness19-mean4
shapiro.test(rarerich19$residual)
ggplot(data = rarerich19) + geom_histogram(mapping = aes(x = residual))
#not normally distributed p value 0.03392

#2020
mean5<-mean(rarerich20$rare_richness20)
rarerich20$residual <- rarerich20$rare_richness20-mean5
shapiro.test(rarerich20$residual)
ggplot(data = rarerich20) + geom_histogram(mapping = aes(x = residual))
#normally distributed p value 0.9187
