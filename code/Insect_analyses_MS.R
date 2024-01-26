# January 18, 2024
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
palette2b<- c("#A9D18E", "#548235")
palette3 <- c("#767171", "#A9D18E", "#548235") #3B5422
palette4 <- c("#767171", "#A9D18E", "#548235", "#2A4117") #3B5422

#  Load data ####
## Site histories ####
rest_history <- read_csv("../vegetation_outcomes_of_restoration/clean/rest_history.csv") %>% 
  select(EasementID = SiteID, RestorationCategory, size, rest_year) 


fire_years  <- read_csv("raw/restoration_history.csv") %>%
  select(EasementID, Fire_Years) %>%
  mutate(EasementID = recode(EasementID, '792' = "00792")) %>% 
  separate(Fire_Years, c("Fire1", "Fire2", "Fire3")) 

## Vegetation Summaries ####
# Non-native plant richness and cover
nonnative_metrics <- read_csv("../vegetation_outcomes_of_restoration/clean/nonnative_metrics.csv") %>% 
  select(EasementID = SiteID, nonnativerich = nativerichness, nonnativecover = avenativecover)

veg_metrics <- read_csv("../vegetation_outcomes_of_restoration/clean/veg_metrics.csv") %>% 
  select(EasementID = SiteID,no_seeded_spp, richness, Shannon, nativerichness, avenativecover, AveCC) %>% 
  left_join(nonnative_metrics)

# Non-iNext data, year removed ####
data <- read_csv("clean/all_data.csv") 

sitebyfam <- data %>% 
  select(EasementID, Year, Sample, Total, Family) %>% 
  group_by(EasementID, Year, Family) %>% 
  summarise(abundance= sum(Total)) %>% 
  select(EasementID, Family, Year, abundance) %>%
  pivot_wider(names_from = Family, values_from = abundance) %>%
  replace(is.na(.), 0) %>%
  ungroup() %>% 
  unite(., "x", EasementID, Year, sep = "_") %>% ## joining columns and 
  column_to_rownames(var = "x")

insect_shannon <- data.frame(diversity(sitebyfam, index = "shannon")) %>% rownames_to_column(var = "x") %>% 
  rename(shannon_div = 'diversity.sitebyfam..index....shannon..') %>% 
  separate(., x, c("EasementID", "Sample_year"), sep = "_") 
  

##Insect iNext data ####
insect_abund <-  read_csv("clean/insect_abundance.csv")# %>% separate(., Assemblage, c("EasementID", "Year"), sep = "_")

rich_est <- read_csv("clean/estimates.csv") %>% 
  filter(Diversity == "Species richness") %>% 
  select(Assemblage, Rich_obs = Observed, Rich_est = Estimator) %>% 
  mutate(Rich_est = round(Rich_est, digits = 0))

## iNext shannon's diversity estimates
shannon_est <- read_csv("clean/estimates.csv") %>% 
  filter(Diversity == "Shannon diversity") %>% 
  select(Assemblage, Shannon_obs = Observed, Shannon_est = Estimator) %>% 
  mutate(Shannon_est = round(Shannon_est, digits = 0))

## iNext Simpson's estimates
simpson_est <- read_csv("clean/estimates.csv") %>% 
  filter(Diversity == "Simpson diversity") %>% 
  select(Assemblage, Simpson_obs = Observed, Simpson_est = Estimator) %>% 
  mutate(Simpson_est = round(Simpson_est, digits = 0)) 

 
#All Variables ####
insect_response <- insect_abund %>% 
  select(Assemblage, Abund = n, Rich_obs = S.obs) %>% 
  left_join(rich_est) %>% 
  left_join(shannon_est)  %>%
  left_join(simpson_est)  %>%
  separate(., Assemblage, c("EasementID", "Sample_year"), sep = "_") %>%  
  left_join(insect_shannon)  %>%
  mutate(pielou_even = shannon_div/log(Rich_obs)) %>%    
  mutate(pielou_iNext = shannon_div/log(Rich_est)) %>% 
  left_join(rest_history) %>% 
  #left_join(age) %>% 
  #convert(int(Sample_year, rest_year)) %>% 
  #mutate_at(rest_age =  (Sample_year - rest_year)) %>% 
  #filter(EasementID != "007XH")  %>% 
  mutate(RestorationCategory = factor(RestorationCategory, levels = c("No Seed", "Seed", "Seed + Fire", "Remnant"))) 
 

##Removing random year####  
insect_response_YR <- insect_response %>% 
  unite(., "x", EasementID, Sample_year, sep = "_") %>% 
  filter(x != "0078S_2020") %>% 
  filter(x !="007C6_2020") %>% 
  filter(x !="007XH_2019") %>% 
  filter(x !="00MBG_2020") %>% 
  filter(x !="00MJC_2020") %>% 
  filter(x !="00NCQ_2019") %>% 
  filter(x !="00WV0_2020") %>% 
  filter(x !="0129G_2019") %>% 
  separate(., x, c("EasementID", "Sample_year"), sep = "_") %>% 
  left_join(veg_metrics)%>%
  mutate(RestorationCategory = recode(RestorationCategory, 'No Seed' = "Old field")) %>% 
  mutate(RestorationCategory = factor(RestorationCategory, levels = c("Old field", "Seed", "Seed + Fire", "Remnant")))

#write_csv(insect_response_YR, "clean/insect_response_YR.csv")



## NRCS Sites ####
NRCS_sites <- insect_response_YR %>% 
  filter(RestorationCategory != "Remnant") %>% 
  filter(nativerichness != "NA")

## Time since fire-burned sites ####
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



# Shapiro-Wilk - Normality####
##Shapiro-wilk test- testing for normality; testing probability that data is normal.  If p-value is low, data is non-normal

shapiro.test(insect_response_YR$Abund) 
#W = 0.96862, p-value = 0.443

shapiro.test(insect_response_YR$Rich_est) 
# W = 0.93619, p-value = 0.05269

shapiro.test(insect_response_YR$Shannon_est)
# W = 0.95758, p-value = 0.2206

## Pielou's evenness Normality####
shapiro.test(insect_response_YR$pielou_even)
#W = 0.92115, p-value = 0.01972

## All metrics are normal, except evenness  #####

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
mod_abund <- aov(Abund ~ RestorationCategory, data = insect_response_YR)
car::Anova(mod_abund,type="II")
summary(mod_abund)

ggResidpanel::resid_panel(mod_abund)

#### Richness ####
mod_rich <- aov(Rich_est ~ RestorationCategory, data = insect_response_YR)
car::Anova(mod_rich,type="II")
summary(mod_rich)
TukeyHSD(mod_rich)


ggResidpanel::resid_panel(mod_rich)

#### Diversity ####
mod_div <- aov(Shannon_est ~ RestorationCategory, data = insect_response_YR)
car::Anova(mod_div,type="II")
TukeyHSD(mod_div)

ggResidpanel::resid_panel(mod_div)

#### Evenness ####
mod_even <- aov(pielou_even ~ RestorationCategory, data = insect_response_YR)
car::Anova(mod_even,type="II")
summary(mod_even)
TukeyHSD(mod_even)

ggResidpanel::resid_panel(mod_div)

#### Plant richness ####
#Abundance
abund_rich <- lm(Abund ~ richness, data = NRCS_sites)
car::Anova(abund_rich,type="II")
summary(abund_rich)

ggResidpanel::resid_panel(abund_rich)

#Insect richness
rich_rich <- lm(Rich_est ~ richness, data = NRCS_sites)
car::Anova(rich_rich,type="II")
summary(rich_rich)

ggResidpanel::resid_panel(rich_rich)
#Diversity
div_rich <- lm(Shannon_est ~ richness, data = NRCS_sites)
car::Anova(div_rich,type="II")
summary(div_rich)

ggResidpanel::resid_panel(div_rich)

#evenness
even_rich <- lm(pielou_even ~ richness, data = NRCS_sites)
car::Anova(even_rich,type="II")
summary(even_rich)

ggResidpanel::resid_panel(even_rich)
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

#Evenness
even_rich <- lm(pielou_even ~ nonnativerich, data = NRCS_sites)
car::Anova(even_rich,type="II")
summary(even_rich)



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

mod_since_fire <- lm(pielou_even ~  years_since_fire, data = burn_insect)
car::Anova(mod_since_fire,type="II")
summary(mod_since_fire)

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

# Plotting ####
#Printing Plots ####
ggsave(path="figures", filename="rich by cat.tiff", width = 5, height = 4, dpi=700)

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
  ylab("Average insect abundance")+
  theme(plot.title = element_text( size = 18, hjust = 0.5, vjust = 2),
        axis.title.x = element_text(size = 15, vjust = 4),
        axis.title.y = element_text(size = 15, hjust = 0.5, vjust = 2),
        axis.text.x =  element_text( size = 12),
        axis.text.y = element_text( size = 12),
        legend.title = element_blank(),
        legend.position = "none",
        plot.margin = unit(c(1,1,2,1), "lines")) +
  geom_jitter(width = .1, size = 2.5)

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
  geom_jitter(width = .1, size = 2.5)

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
  ylab("Insect family Shannon's diversity \n (hill number, q = 1)")+
  theme(plot.title = element_text( size = 18, hjust = 0.5, vjust = 2),
        axis.title.x = element_text( size = 15, vjust = 4),
        axis.title.y = element_text( size = 15, hjust = 0.5, vjust =2),
        axis.text.x =  element_text( size = 12),
        axis.text.y = element_text( size = 12),
        legend.title = element_blank(),
        legend.position = "none",
        plot.margin = unit(c(1,1,2,1), "lines")) +
  geom_jitter(width = .1, size = 2.5)

## Estimated evenness ####
plot_even <- insect_response_YR 
plot_even %>% 
  ggplot(aes(RestorationCategory, pielou_even, fill=RestorationCategory)) +
  geom_boxplot(outlier.alpha = 0) +
  theme_classic()+
  scale_fill_manual(values = c(palette4),
                    name = "Site Categories")+
  labs(title = "")+
  xlab("\n Site category") +
  ylab("Insect family Pielou's evenness")+
  theme(plot.title = element_text( size = 18, hjust = 0.5, vjust = 2),
        axis.title.x = element_text( size = 15, vjust = 4),
        axis.title.y = element_text( size = 15, hjust = 0.5, vjust =2),
        axis.text.x =  element_text( size = 12),
        axis.text.y = element_text( size = 12),
        legend.title = element_blank(),
        legend.position = "none",
        plot.margin = unit(c(1,1,2,1), "lines")) +
  geom_jitter(width = .1, size = 2.5)

## Plant species richness####
## evenness
NRCS_sites %>% 
  ggplot(aes(richness, pielou_even, color = RestorationCategory)) +
  geom_point(aes(),size=2.5) +
  #geom_smooth(aes(color = RestorationCategory), method = "lm") + 
  geom_smooth(aes(), color = "black", method = "lm", se=FALSE) +
  theme_classic()+
  scale_color_manual(values = palette3,
                     name = "Site Categories")+
  labs(title = "")+
  xlab("Plant species richness") +
  ylab("Insect family \n Pielou's evenness")+
  theme(plot.title = element_text( size = 18, hjust = 1, vjust = 1),
        axis.title.x = element_text( size = 15, vjust = 1),
        axis.title.y = element_text( size = 15, hjust = 0.5, vjust =2),
        axis.text.x =  element_text( size = 12),
        axis.text.y = element_text( size = 12),
        legend.title = element_blank(),
        legend.position = "bottom",
        plot.margin = unit(c(1,1,2,1), "lines")) 

## Diversity
NRCS_sites %>% 
  ggplot(aes(richness, Shannon_est, color = RestorationCategory)) +
  geom_point(aes(),size=2.5) +
  #geom_smooth(aes(color = RestorationCategory), method = "lm") + 
  geom_smooth(aes(), color = "black", method = "lm", se=FALSE) +
  theme_classic()+
  scale_color_manual(values = palette3,
                     name = "Site Categories")+
  labs(title = "")+
  xlab("Plant species richness") +
  ylab("Insect family Shannon's \n diversity (q = 1)")+
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
  geom_point(aes(), size = 2.5) +
  #geom_smooth(aes(color = RestorationCategory), method = "lm") + 
  geom_smooth(aes(), color = "black", method = "lm", se=FALSE) +
  theme_classic()+
  scale_color_manual(values = palette3,
                     name = "Site Categories")+
  labs(title = "")+
  xlab("Plant species richness") +
  ylab("Insect family richness \n (q = 0)")+
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
  geom_point(aes(), size = 2.5) +
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
## evenness
NRCS_sites %>% 
  ggplot(aes(nonnativerich, pielou_even, color = RestorationCategory)) +
  geom_point(aes(), size = 2.5) +
  #geom_smooth(aes(color = RestorationCategory), method = "lm") + 
  geom_smooth(aes(), color = "black", method = "lm", se=FALSE) +
  theme_classic()+
  scale_color_manual(values = palette3,
                     name = "Site Categories")+
  labs(title = "")+
  xlab("Non-native plant richness") +
  ylab("Insect family \n Pielou's evenness")+
  theme(plot.title = element_text( size = 18, hjust = 1, vjust = 1),
        axis.title.x = element_text( size = 15, vjust = 1),
        axis.title.y = element_text( size = 15, hjust = 0.5, vjust =2),
        axis.text.x =  element_text( size = 12),
        axis.text.y = element_text( size = 12),
        legend.title = element_blank(),
        legend.position = "bottom",
        plot.margin = unit(c(1,1,2,1), "lines")) 

## Diversity
NRCS_sites %>% 
  ggplot(aes(nonnativerich, Shannon_est, color = RestorationCategory)) +
  geom_point(aes(), size = 2.5) +
  #geom_smooth(aes(color = RestorationCategory), method = "lm") + 
  geom_smooth(aes(), color = "black", method = "lm", se=FALSE) +
  theme_classic()+
  scale_color_manual(values = palette3,
                     name = "Site Categories")+
  labs(title = "")+
  xlab("Non-native plant richness") +
  ylab("Insect family Shannon's \n diversity (q = 1)")+
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
  geom_point(aes(), size = 2.5) +
  #geom_smooth(aes(color = RestorationCategory), method = "lm") + 
  geom_smooth(aes(), color = "black", method = "lm", se=FALSE) +
  theme_classic()+
  scale_color_manual(values = palette3,
                     name = "Site Categories")+
  labs(title = "")+
  xlab("Non-native plant richness") +
  ylab("Insect family richness \n (q = 0)")+
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
  geom_point(aes(), size = 2.5) +
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
## Evenness
burn_insect%>% 
  ggplot(aes(years_since_fire, pielou_even, color = RestorationCategory)) +
  geom_point(aes(), size = 2.5, outlier.alpha = 0) +
  #geom_smooth(aes(color = RestorationCategory), method = "lm") + 
  geom_smooth(aes(), color = "black", method = "lm", se=FALSE) +
  theme_classic()+
  #scale_color_manual(values = c("Seed" = "#A9D18E", "Seed + Fire" = "#548235")) +
  scale_color_manual(values = palette2,
                     name = "Site Categories")+
  labs(title = "")+
  xlab("\n Years since last fire") +
  ylab("Insect family \n Pielou's evenness")+
  theme(plot.title = element_text( size = 18, hjust = 0.1, vjust = 1),
        axis.title.x = element_text( size = 15, vjust = 4),
        axis.title.y = element_text( size = 15, hjust = 0.5, vjust =2),
        axis.text.x =  element_text( size = 12),
        axis.text.y = element_text( size = 12),
        legend.title = element_blank(),
        legend.position = "bottom",
        plot.margin = unit(c(1,1,2,1), "lines")) 

## Diversity
burn_insect%>% 
  ggplot(aes(years_since_fire, Shannon_est, color = RestorationCategory)) +
  geom_point(aes(), size = 2.5, outlier.alpha = 0) +
  #geom_smooth(aes(color = RestorationCategory), method = "lm") + 
  geom_smooth(aes(), color = "black", method = "lm", se=FALSE) +
  theme_classic()+
  #scale_color_manual(values = c("Seed" = "#A9D18E", "Seed + Fire" = "#548235")) +
  scale_color_manual(values = palette2,
                     name = "Site Categories")+
  labs(title = "")+
  xlab("\n Years since last fire") +
  ylab("Insect family Shannon's \n diversity (q = 1)")+
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
  geom_point(aes(), size = 2.5, outlier.alpha = 0) +
  #geom_smooth(aes(color = RestorationCategory), method = "lm") + 
  geom_smooth(aes(), color = "black", method = "lm", se=FALSE) +
  theme_classic()+
  scale_color_manual(values = palette2,
                     name = "Site Categories")+
  labs(title = "")+
  xlab("\n Years since last fire") +
  ylab("Insect family richness \n (q = 0)")+
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
  geom_point(aes(), size = 2.5, outlier.alpha = 0) +
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



