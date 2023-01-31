# January 27, 2023
# by: Stephanie McFarlane
# Final analyses for manuscript
#This code does not remove aggregating insects

#Load packages####
library(tidyverse)
library(hablar)
library(vegan)
library(devtools)
#library(BiodiversityR)


# Color palettes ####
palette2 <- c("#548235", "#2A4117")
palette2b<-  c("#A9D18E", "#548235")
palette3<-  c("#767171", "#A9D18E", "#548235") #3B5422
palette4<-  c("#767171", "#A9D18E", "#548235", "#2A4117") #3B5422

#  Load data ####
##Site History Data ####
# Upload easement treatment (aka restoration category) data & clean ##

### Restoration History ####
rest_history <- read_csv("clean/restoration_history.csv") %>% 
  mutate(EasementID = recode(EasementID, '792' = "00792")) 

### Restoration Age ####
age <- read_csv("clean/rest_age.csv") %>% 
  select(EasementID = SiteID, rest_year)

rest_year <- read_csv("clean/enroll_rest.csv") 

###Fire History####
fire_years  <- read_csv("raw/restoration_history.csv") %>% 
  select(EasementID, Fire_Years) %>%
  mutate(EasementID = recode(EasementID, '792' = "00792")) %>% 
  separate(Fire_Years, c("Fire1", "Fire2", "Fire3")) 

## Vegetation Summaries ####
nonnative_metrics <- read_csv("../vegetation_outcomes_of_restoration/clean/nonnative_metrics.csv") %>% 
  select(EasementID = SiteID, nonnativerich = nativerichness, nonnativecover = avenativecover)

veg_metrics <- read_csv("../vegetation_outcomes_of_restoration/clean/veg_metrics.csv") %>% 
  select(EasementID = SiteID,no_seeded_spp, richness, nativerichness, avenativecover, AveCC) %>% 
  left_join(nonnative_metrics)



# #Insect iNext data ####
insect_abundance <- read_csv("clean/insect_abundance.csv")
#write_csv(insect_rich$DataInfo, "clean/abundance_withAgg.csv")
estimates <- read_csv("clean/estimates.csv")
estimates_120 <- read_csv("clean/estimates_120.csv")

## abundance
insect_abund <- insect_abundance %>% 
  separate(., Assemblage, c("EasementID", "Year"), sep = "_")

rich_est <- estimates %>% 
  filter(Diversity == "Species richness") %>% 
  select(Assemblage, Rich_obs = Observed, Rich_est = Estimator) %>% 
  mutate(Rich_est = round(Rich_est, digits = 0))

shannon_est <- estimates %>% 
  filter(Diversity == "Shannon diversity") %>% 
  select(Assemblage, Shannon_obs = Observed, Shannon_est = Estimator) %>% 
  mutate(Shannon_est = round(Shannon_est, digits = 0))

rich_est_120 <- estimates_120 %>% 
  filter(Diversity == "Species richness") %>% 
  select(Assemblage, Rich_obs = Observed, Rich_est_120 = Estimator) %>% 
  mutate(Rich_est_120 = round(Rich_est_120, digits = 0))
  


#All Variables ####
## Must remove sites with abundance less than 225 when modeling estimated richness #filter(Abund > 200)
insect_response <- insect_abundance %>% 
  select(Assemblage, Abund = n, Rich_obs = S.obs) %>% 
  left_join(rich_est) %>% 
  left_join(rich_est_120) %>% 
  left_join(shannon_est)  %>% #rename(ID=Site)) %>% 
  separate(., Assemblage, c("EasementID", "Sample_year"), sep = "_") %>%  
  left_join(rest_history) %>% 
  left_join(age) %>% 
  convert(int(Sample_year, rest_year)) %>% 
  mutate(rest_age =  (Sample_year - rest_year)) %>% 
  select(-rest_year) %>% 
  filter(EasementID != "007XH")  %>% 
  mutate(RestorationCategory = factor(RestorationCategory, levels = c("No Seed", "Seed", "Seed + Fire", "Remnant"))) 
 

##Removing random year####  
insect_response_YR <- insect_response %>% 
  unite(., "x", EasementID, Sample_year, sep = "_") %>% 
  filter(x != "0078S_2020") %>% 
  filter(x !="007C6_2020") %>% 
  filter(x !="00MBG_2020") %>% 
  filter(x !="00MJC_2020") %>% 
  filter(x !="00NCQ_2019") %>% 
  filter(x !="00WV0_2020") %>% 
  separate(., x, c("EasementID", "Sample_year"), sep = "_") %>% 
  #mutate(Rich_est = ifelse(Abund < 200, "NA", Rich_est)) ## Must remove sites with abundance less than 225 when modeling estimated richness
  left_join(veg_metrics)%>% 
  mutate(RestorationCategory = recode(RestorationCategory, 'No Seed' = "Old field")) %>% 
  mutate(RestorationCategory = factor(RestorationCategory, levels = c("Old field", "Seed", "Seed + Fire", "Remnant")))

## NRCS Sites ####
NRCS_sites <- insect_response_YR %>% 
  filter(RestorationCategory != "Remnant") %>% 
  filter(nativerichness != "NA")

## Time since fire -burned sites ####
burn_insect <- insect_response_YR %>% 
  filter(RestorationCategory %in% c("Seed + Fire", "Remnant")) %>% 
  left_join(fire_years) %>% 
  mutate(Fire3 = ifelse(Fire3 > Sample_year, NA, Fire3)) %>% 
  mutate(Fire2 = ifelse(Fire2 > Sample_year, NA, Fire2)) %>% 
  mutate(Fire1 = ifelse(Fire1 > Sample_year, NA, Fire1)) %>% 
  mutate_at(vars(c(Fire1, Fire2, Fire3)), ~replace(., is.na(.), "")) %>% 
  mutate(last_fire = ifelse(Fire3 > 0, Fire3, Fire2)) %>% 
  mutate(last_fire = ifelse(last_fire > 0, last_fire, Fire1)) %>% 
  convert(int(last_fire)) %>% 
  convert(int(Sample_year)) %>% 
  mutate(years_since_fire = (Sample_year - last_fire)) %>% 
  select(-Fire3, -Fire2, -Fire1) %>% 
  left_join(veg_metrics)

##Seeded sites ####
seeded_insect <- insect_response_YR %>%
  filter(RestorationCategory != "No Seed") %>% 
  filter(RestorationCategory != "Remnant") %>% 
  left_join(veg_metrics)

# Normality & heteroscedaticity#### 

## Abundance Normality####
ggplot(abund_norm)+ geom_histogram(mapping = aes(x = Abund), bins = 6)
shapiro.test(insect_response_YR$Abund) 
#W = 0.96663, p-value = 0.3932

## Estimated Richness Normality####
shapiro.test(insect_response_YR$Rich_est) 
# W = 0.938, p-value = 0.05943

## Shannon's diversity Normality####
shapiro.test(insect_response_YR$Shannon_est)
# W = 0.95535, p-value = 0.1903
# Plotting diversity by restoration category ####
## Abundance ####
insect_response_YR %>% 
  ggplot(aes(RestorationCategory, Abund, fill=RestorationCategory)) +
  geom_boxplot(outlier.alpha = 0) +
  theme_classic()+
  scale_fill_manual(values = c(palette4),
                    name = "Site Categories")+
  labs(title = "")+
  xlab("\n Site Category") +
  ylab("Average Abundance")+
  theme(plot.title = element_text( size = 18, hjust = 0.5, vjust = 2),
        axis.title.x = element_text(size = 15, vjust = 4),
        axis.title.y = element_text(size = 15, hjust = 0.5, vjust = 2),
        axis.text.x =  element_text( size = 12),
        axis.text.y = element_text( size = 12),
        legend.title = element_blank(),
        legend.position = "none",
        plot.margin = unit(c(1,1,2,1), "lines")) +
  geom_jitter( width = .1, size = 1.4)

## Estimated richness ####
plot_rich <- insect_response_YR #%>% filter(Abund > 200) 
plot_rich %>% 
  ggplot(aes(RestorationCategory, Rich_est, fill=RestorationCategory)) +
  geom_boxplot(outlier.alpha = 0) +
  theme_classic()+
  scale_fill_manual(values = c(palette4),
                    name = "Site categories")+
  labs(title = "")+
  xlab("\n Site category") +
  ylab("Insect family richness \n (hill number, q = 0)")+
  theme(plot.title = element_text(size = 18, hjust = 0.5, vjust = 2),
        axis.title.x = element_text(size = 15, vjust = 4),
        axis.title.y = element_text(size = 15, hjust = 0.5, vjust =2),
        axis.text.x =  element_text(size = 12),
        axis.text.y = element_text(size = 12),
        legend.title = element_blank(),
        legend.position = "none",
        plot.margin = unit(c(1,1,2,1), "lines")) +
  geom_jitter(width = .1, size = 1.4)

## Estimated diversity ####
plot_div <- insect_response_YR 
plot_div %>% 
  ggplot(aes(RestorationCategory, Shannon_est, fill=RestorationCategory)) +
  geom_boxplot(outlier.alpha = 0) +
  theme_classic()+
  scale_fill_manual(values = c(palette4),
                    name = "Site Categories")+
  labs(title = "")+
  xlab("\n Site category") +
  ylab("Insect family diversity \n (hill number, q = 1)")+
  theme(plot.title = element_text( size = 18, hjust = 0.5, vjust = 2),
        axis.title.x = element_text( size = 15, vjust = 4),
        axis.title.y = element_text( size = 15, hjust = 0.5, vjust =2),
        axis.text.x =  element_text( size = 12),
        axis.text.y = element_text( size = 12),
        legend.title = element_blank(),
        legend.position = "none",
        plot.margin = unit(c(1,1,2,1), "lines")) +
  geom_jitter( width = .1, size = 1.4)

## Plant species richness####
## Diversity
NRCS_sites %>% 
  ggplot(aes(richness, Shannon_est, color = RestorationCategory)) +
  geom_point() +
  #geom_smooth(aes(color = RestorationCategory), method = "lm") + 
  geom_smooth(aes(), color = "black", method = "lm", se=FALSE) +
  theme_classic()+
  scale_color_manual(values = palette3,
                     name = "Site Categories")+
  labs(title = "")+
  xlab("Plant species richness") +
  ylab("Insect family diversity \n (hill number, q = 1)")+
  theme(plot.title = element_text( size = 18, hjust = 1, vjust = 1),
        axis.title.x = element_text( size = 15, vjust = 1),
        axis.title.y = element_text( size = 15, hjust = 0.5, vjust =2),
        axis.text.x =  element_text( size = 12),
        axis.text.y = element_text( size = 12),
        legend.title = element_blank(),
        legend.position = "bottom",
        plot.margin = unit(c(1,1,2,1), "lines")) 

## Richness
NRCS_sites %>% 
  ggplot(aes(richness, Rich_est, color = RestorationCategory)) +
  geom_point() +
  #geom_smooth(aes(color = RestorationCategory), method = "lm") + 
  geom_smooth(aes(), color = "black", method = "lm", se=FALSE) +
  theme_classic()+
  scale_color_manual(values = palette3,
                     name = "Site Categories")+
  labs(title = "")+
  xlab("Plant species richness") +
  ylab("Insect family richness \n (hill number, q = 0)")+
  theme(plot.title = element_text( size = 18, hjust = 0.1, vjust = 1),
        axis.title.x = element_text( size = 15, vjust = 1),
        axis.title.y = element_text( size = 15, hjust = 0.5, vjust =2),
        axis.text.x =  element_text( size = 12),
        axis.text.y = element_text( size = 12),
        legend.title = element_blank(),
        legend.position = "bottom",
        plot.margin = unit(c(1,1,2,1), "lines")) 

## Abundance
NRCS_sites %>% 
  ggplot(aes(richness, Abund, color = RestorationCategory)) +
  geom_point() +
  #geom_smooth(aes(color = RestorationCategory), method = "lm") + 
  geom_smooth(aes(), color = "black", method = "lm", se=FALSE) +
  theme_classic()+
  scale_color_manual(values = palette3,
                     name = "Site Categories")+
  labs(title = "")+
  xlab("Plant species richness") +
  ylab("Insect abundance")+
  theme(plot.title = element_text( size = 18, hjust = 0.1, vjust = 1),
        axis.title.x = element_text( size = 15, vjust = 1),
        axis.title.y = element_text( size = 15, hjust = 0.5, vjust =2),
        axis.text.x =  element_text( size = 12),
        axis.text.y = element_text( size = 12),
        legend.title = element_blank(),
        legend.position = "bottom",
        plot.margin = unit(c(1,1,2,1), "lines")) 

## Non-native plant species richness####
## Diversity
NRCS_sites %>% 
  ggplot(aes(nonnativerich, Shannon_est, color = RestorationCategory)) +
  geom_point() +
  #geom_smooth(aes(color = RestorationCategory), method = "lm") + 
  geom_smooth(aes(), color = "black", method = "lm", se=FALSE) +
  theme_classic()+
  scale_color_manual(values = palette3,
                     name = "Site Categories")+
  labs(title = "")+
  xlab("Non-native plant richness") +
  ylab("Insect family diversity \n (hill number, q = 1)")+
  theme(plot.title = element_text( size = 18, hjust = 1, vjust = 1),
        axis.title.x = element_text( size = 15, vjust = 1),
        axis.title.y = element_text( size = 15, hjust = 0.5, vjust =2),
        axis.text.x =  element_text( size = 12),
        axis.text.y = element_text( size = 12),
        legend.title = element_blank(),
        legend.position = "bottom",
        plot.margin = unit(c(1,1,2,1), "lines")) 

## Richness
NRCS_sites %>% 
  ggplot(aes(nonnativerich, Rich_est, color = RestorationCategory)) +
  geom_point() +
  #geom_smooth(aes(color = RestorationCategory), method = "lm") + 
  geom_smooth(aes(), color = "black", method = "lm", se=FALSE) +
  theme_classic()+
  scale_color_manual(values = palette3,
                     name = "Site Categories")+
  labs(title = "")+
  xlab("Non-native plant richness") +
  ylab("Insect family richness \n (hill number, q = 0)")+
  theme(plot.title = element_text( size = 18, hjust = 0.1, vjust = 1),
        axis.title.x = element_text( size = 15, vjust = 1),
        axis.title.y = element_text( size = 15, hjust = 0.5, vjust =2),
        axis.text.x =  element_text( size = 12),
        axis.text.y = element_text( size = 12),
        legend.title = element_blank(),
        legend.position = "bottom",
        plot.margin = unit(c(1,1,2,1), "lines")) 

## Abundance
NRCS_sites %>% 
  ggplot(aes(nonnativerich, Abund, color = RestorationCategory)) +
  geom_point() +
  #geom_smooth(aes(color = RestorationCategory), method = "lm") + 
  geom_smooth(aes(), color = "black", method = "lm", se=FALSE) +
  theme_classic()+
  scale_color_manual(values = palette3,
                     name = "Site Categories")+
  labs(title = "")+
  xlab("Non-native plant richness") +
  ylab("Insect abundance")+
  theme(plot.title = element_text( size = 18, hjust = 0.1, vjust = 1),
        axis.title.x = element_text( size = 15, vjust = 1),
        axis.title.y = element_text( size = 15, hjust = 0.5, vjust =2),
        axis.text.x =  element_text( size = 12),
        axis.text.y = element_text( size = 12),
        legend.title = element_blank(),
        legend.position = "bottom",
        plot.margin = unit(c(1,1,2,1), "lines")) 

## Time since fire ####
## Diversity
burn_insect%>% 
  ggplot(aes(years_since_fire, Shannon_est, color = RestorationCategory)) +
  geom_point(outlier.alpha = 0) +
  #geom_smooth(aes(color = RestorationCategory), method = "lm") + 
  geom_smooth(aes(), color = "black", method = "lm", se=FALSE) +
  theme_classic()+
  #scale_color_manual(values = c("Seed" = "#A9D18E", "Seed + Fire" = "#548235")) +
  scale_color_manual(values = palette2,
                     name = "Site Categories")+
  labs(title = "")+
  xlab("\n Years since last fire") +
  ylab("Insect family diversity \n (hill number, q = 1)")+
  theme(plot.title = element_text( size = 18, hjust = 0.1, vjust = 1),
        axis.title.x = element_text( size = 15, vjust = 4),
        axis.title.y = element_text( size = 15, hjust = 0.5, vjust =2),
        axis.text.x =  element_text( size = 12),
        axis.text.y = element_text( size = 12),
        legend.title = element_blank(),
        legend.position = "bottom",
        plot.margin = unit(c(1,1,2,1), "lines")) 

## Richness
burn_insect%>% 
  ggplot(aes(years_since_fire, Rich_est, color = RestorationCategory)) +
  geom_point(outlier.alpha = 0) +
  #geom_smooth(aes(color = RestorationCategory), method = "lm") + 
  geom_smooth(aes(), color = "black", method = "lm", se=FALSE) +
  theme_classic()+
  scale_color_manual(values = palette2,
                     name = "Site Categories")+
  labs(title = "")+
  xlab("\n Years since last fire") +
  ylab("Insect family richness \n (hill number, q = 0)")+
  theme(plot.title = element_text( size = 18, hjust = 0.1, vjust = 1),
        axis.title.x = element_text( size = 15, vjust = 4),
        axis.title.y = element_text( size = 15, hjust = 0.5, vjust =2),
        axis.text.x =  element_text( size = 12),
        axis.text.y = element_text( size = 12),
        legend.title = element_blank(),
        legend.position = "bottom",
        plot.margin = unit(c(1,1,2,1), "lines")) 

## Abundance
burn_insect%>% 
  ggplot(aes(years_since_fire, Abund, color = RestorationCategory)) +
  geom_point(outlier.alpha = 0) +
  #geom_smooth(aes(color = RestorationCategory), method = "lm") + 
  geom_smooth(aes(), color = "black", method = "lm", se=FALSE) +
  theme_classic()+
  scale_color_manual(values = palette2,
                     name = "Site categories")+
  labs(title = "")+
  xlab("\n Years since last fire") +
  ylab("Insect abundance")+
  theme(plot.title = element_text( size = 18, hjust = 0.1, vjust = 1),
        axis.title.x = element_text( size = 15, vjust = 4),
        axis.title.y = element_text( size = 15, hjust = 0.5, vjust =2),
        axis.text.x =  element_text( size = 12),
        axis.text.y = element_text( size = 12),
        legend.title = element_blank(),
        legend.position = "bottom",
        plot.margin = unit(c(1,1,2,1), "lines")) 

## Seedmix richness ####
## Diversity
seeded_insect%>% filter(no_seeded_spp != "NA") %>% 
  ggplot(aes(no_seeded_spp, Shannon_est, color = RestorationCategory)) +
  geom_point(outlier.alpha = 0) +
  #geom_smooth(aes(color = RestorationCategory), method = "lm") + 
  geom_smooth(aes(), color = "black", method = "lm", se=FALSE) +
  theme_classic()+
  scale_color_manual(values = palette2b,
                     name = "Site Categories")+
  labs(title = "")+
  xlab("\n Seedmix richness") +
  ylab("Insect family diversity \n (hill number, q = 1)")+
  theme(plot.title = element_text( size = 18, hjust = 0.1, vjust = 1),
        axis.title.x = element_text( size = 15, vjust = 4),
        axis.title.y = element_text( size = 15, hjust = 0.5, vjust =2),
        axis.text.x =  element_text( size = 12),
        axis.text.y = element_text( size = 12),
        legend.title = element_blank(),
        legend.position = "bottom",
        plot.margin = unit(c(1,1,2,1), "lines")) 

## Richness
seeded_insect%>% filter(no_seeded_spp != "NA") %>% 
  ggplot(aes(no_seeded_spp, Rich_est, color = RestorationCategory)) +
  geom_point(outlier.alpha = 0) +
  #geom_smooth(aes(color = RestorationCategory), method = "lm") + 
  geom_smooth(aes(), color = "black", method = "lm", se=FALSE) +
  theme_classic()+
  scale_color_manual(values = palette2b,
                     name = "Site Categories")+
  labs(title = "")+
  xlab("\n Seedmix richness") +
  ylab("Insect family richness \n (hill number, q = 0)")+
  theme(plot.title = element_text( size = 18, hjust = 0.1, vjust = 1),
        axis.title.x = element_text( size = 15, vjust = 4),
        axis.title.y = element_text( size = 15, hjust = 0.5, vjust =2),
        axis.text.x =  element_text( size = 12),
        axis.text.y = element_text( size = 12),
        legend.title = element_blank(),
        legend.position = "bottom",
        plot.margin = unit(c(1,1,2,1), "lines")) 

## Abundance
seeded_insect%>% filter(no_seeded_spp != "NA") %>% 
  ggplot(aes(no_seeded_spp, Abund, color = RestorationCategory)) +
  geom_point(outlier.alpha = 0) +
  #geom_smooth(aes(color = RestorationCategory), method = "lm") + 
  geom_smooth(aes(), color = "black", method = "lm", se=FALSE) +
  theme_classic()+
  scale_color_manual(values = palette2b,
                     name = "Site categories")+
  labs(title = "")+
  xlab("\n Seedmix richness") +
  ylab("Insect abundance")+
  theme(plot.title = element_text( size = 18, hjust = 0.1, vjust = 1),
        axis.title.x = element_text( size = 15, vjust = 4),
        axis.title.y = element_text( size = 15, hjust = 0.5, vjust =2),
        axis.text.x =  element_text( size = 12),
        axis.text.y = element_text( size = 12),
        legend.title = element_blank(),
        legend.position = "bottom",
        plot.margin = unit(c(1,1,2,1), "lines")) 


# Models ####
## Possible variables to include: Restoration Category, plant spp richness, native species richness, seeded species rich , size, and age.  Note that age can only be used for restored sites and not for any analysis that includes remnant

mod <- aov(Shannon_est ~ RestorationCategory, data = insect_response_YR)
summary(mod)
TukeyHSD(mod)

mod <- aov(Rich_est ~ RestorationCategory, data = insect_response_YR)
summary(mod)
TukeyHSD(mod)

mod <- aov(Abund ~ RestorationCategory, data = insect_response_YR)
summary(mod)
TukeyHSD(mod)

mod <- aov(nonnativerich ~ RestorationCategory, data = insect_response_YR)
summary(mod)
TukeyHSD(mod)

#### Abundance ####
mod_abund <- lm(Abund ~ RestorationCategory + nativerichness, data = insect_response_YR)
car::Anova(mod_abund,type="II")
summary(mod_abund)

ggResidpanel::resid_panel(mod_abund)

#### Richness ####
mod_rich <- lm(Rich_est ~ RestorationCategory + nativerichness, data = insect_response_YR)
car::Anova(mod_rich,type="II")
summary(mod_rich)

ggResidpanel::resid_panel(mod_rich)

#### Diversity ####
mod_div <- lm(Shannon_est ~ RestorationCategory + nativerichness, data = insect_response_YR)
car::Anova(mod_div,type="II")
summary(mod_div)

ggResidpanel::resid_panel(mod_div)

#### Plant richness ####
#Abundance
abund_rich <- lm(Abund ~ richness, data = NRCS_sites)
car::Anova(abund_rich,type="II")
summary(abund_rich)

ggResidpanel::resid_panel(abund_rich)

#Diversity
rich_rich <- lm(Rich_est ~ richness, data = NRCS_sites)
car::Anova(rich_rich,type="II")
summary(rich_rich)

ggResidpanel::resid_panel(rich_rich)
#Diversity
div_rich <- lm(Shannon_est ~ richness, data = NRCS_sites)
car::Anova(div_rich,type="II")
summary(div_rich)

ggResidpanel::resid_panel(div_rich)

#### Non-native plant richness ####
#Abundance
abund_rich <- lm(Abund ~ nonnativerich, data = NRCS_sites)
car::Anova(abund_rich,type="II")
summary(abund_rich)

ggResidpanel::resid_panel(abund_rich)

#Diversity
rich_rich <- lm(Rich_est ~ nonnativerich, data = NRCS_sites)
car::Anova(rich_rich,type="II")
summary(rich_rich)

ggResidpanel::resid_panel(rich_rich)
#Diversity
div_rich <- lm(Shannon_est ~ nonnativerich, data = NRCS_sites)
car::Anova(div_rich,type="II")
summary(div_rich)

ggResidpanel::resid_panel(div_rich)

#### Non-native plant cover ####
#Abundance
abund_rich <- lm(Abund ~ nonnativecover , data = NRCS_sites)
car::Anova(abund_rich,type="II")
summary(abund_rich)

ggResidpanel::resid_panel(abund_rich)

#Diversity
rich_rich <- lm(Rich_est ~ nonnativecover , data = NRCS_sites)
car::Anova(rich_rich,type="II")
summary(rich_rich)

ggResidpanel::resid_panel(rich_rich)
#Diversity
div_rich <- lm(Shannon_est ~ nonnativecover , data = NRCS_sites)
car::Anova(div_rich,type="II")
summary(div_rich)

ggResidpanel::resid_panel(div_rich)

#### Time since fire ####
#Abundance
mod_since_fire <- lm(Abund ~ years_since_fire, data = burn_insect)
car::Anova(mod_since_fire,type="II")
summary(mod_since_fire)

ggResidpanel::resid_panel(mod_since_fire)

#Richness
mod_since_fire <- lm(Rich_est ~  years_since_fire, data = burn_insect)
car::Anova(mod_since_fire,type="II")
summary(mod_since_fire)

#Diversity
mod_since_fire <- lm(Shannon_est ~  years_since_fire, data = burn_insect)
car::Anova(mod_since_fire,type="II")
summary(mod_since_fire)

ggResidpanel::resid_panel(mod_since_fire)

#### Seed mix richness ####
#Abundance 
mod_seed <- lm(Abund ~  no_seeded_spp, data = seeded_insect)
car::Anova(mod_seed,type="II")
summary(mod_seed)

ggResidpanel::resid_panel(mod_seed)

# Richness
mod_seed <- lm(Rich_est ~  no_seeded_spp, data = seeded_insect)
car::Anova(mod_seed,type="II")
summary(mod_seed)

ggResidpanel::resid_panel(mod_seed)
#Diversity
mod_seed <- lm(Shannon_est ~ no_seeded_spp, data = seeded_insect)
car::Anova(mod_seed,type="II")
summary(mod_seed)

ggResidpanel::resid_panel(mod_seed)

### Testing for quadratic fit ####
## per John Orrock suggestion Oct 24, 2022
all_data2 <- all_data %>% 
  mutate(rest_age2 = rest_age^2)

mod_natQ <- lm(nativerichness ~ RestorationCategory + rest_age + rest_age2, data = all_data2)

## plot quadratic fit 
p <- ggplot(burn_insect, aes(x = years_since_fire, y = Shannon_iNext)) + geom_point()
p + stat_smooth(method = "lm", formula = y ~ x + I(x^2), size = 1)

models <- list(mod_natQ, mod_nat)

#specify model names
mod.names <- c('quadratic', 'linear')

anova(mod_nat)
anova(mod_natQ)

## not really quadratic in plot
## worse AIC, so not going to consider quadratic
###End Test####
