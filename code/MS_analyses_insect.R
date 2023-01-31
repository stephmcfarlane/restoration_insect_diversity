library(tidyverse)
library(hablar)
library(vegan)
library(devtools)
library(BiodiversityR)


# Color palettes ####
palette2 <- c("#548235", "#2A4117")
palette2b<-  c("#A9D18E", "#548235")
palette3<-  c("#767171", "#A9D18E", "#548235") #3B5422
palette4<-  c("#767171", "#A9D18E", "#548235", "#2A4117") #3B5422

# Site History Data ####
# Upload easement treatment (aka restoration category) data & clean ##

## Restoration History ####
rest_history <- read_csv("raw/restoration_history.csv") %>% 
  select(EasementID, RestorationCategory = RestorationCategory_NEW, Size = Contiguous_Seeding_Size) %>% 
  mutate(EasementID = recode(EasementID, '792' = "00792")) %>%
  mutate(RestorationCategory = factor(RestorationCategory, levels = c("No Seed", "Seed", "Seed + Fire", "Remnant"))) 


#write_csv(rest_history, "clean/restoration_history.csv")

## Restoration Age####
age <- read_csv("clean/rest_age.csv") %>% 
  select(EasementID = SiteID, rest_year)

##Fire History####
fire_years  <- read_csv("raw/restoration_history.csv") %>% 
  select(EasementID, Fire_Years) %>%
  mutate(EasementID = recode(EasementID, '792' = "00792")) %>% 
  separate(Fire_Years, c("Fire1", "Fire2", "Fire3")) 

# Vegetation Summaries ####
veg_metrics <- read_csv("../vegetation_outcomes_of_restoration/clean/veg_metrics.csv") %>% 
  select(EasementID = SiteID,no_seeded_spp, richness, nativerichness, avenativecover, AveCC)

# Clean insect Data - iNext ####
insect_response <- read_csv("clean/insect_response.csv") %>% 
 select(-RestorationCategory) %>% 
 left_join(rest_history) %>% 
 left_join(veg_metrics)

####Removing random year####  
insect_response_YR <- insect_response %>% 
  unite(., "x", EasementID, Sample_year, sep = "_") %>% 
  filter(x != "0078S_2020") %>% 
  filter(x !="007C6_2020") %>% 
  filter(x !="00MBG_2020") %>% 
  filter(x !="00MJC_2020") %>% 
  filter(x !="00NCQ_2019") %>% 
  filter(x !="00WV0_2020") %>% 
  separate(., x, c("EasementID", "Sample_year"), sep = "_") 

#### Time since fire -burned sites ####
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

#### Seeded sites ####
seeded_insect <- insect_response_YR %>%
  filter(RestorationCategory != "No Seed") %>% 
  filter(RestorationCategory != "Remnant") %>% 
  left_join(veg_metrics)
  

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
        axis.title.y = element_text(size = 15, hjust = 0.5, vjust =4),
        axis.text.x =  element_text( size = 12),
        axis.text.y = element_text(, size = 12),
        legend.title = element_blank(),
        legend.position = "none",
        plot.margin = unit(c(1,1,2,1), "lines")) +
  geom_jitter( width = .1, size = 3)

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
  ylab("Insect family richness")+
  theme(plot.title = element_text(size = 18, hjust = 0.5, vjust = 2),
        axis.title.x = element_text(size = 15, vjust = 4),
        axis.title.y = element_text(size = 15, hjust = 0.5, vjust =4),
        axis.text.x =  element_text( size = 12),
        axis.text.y = element_text(size = 12),
        legend.title = element_blank(),
        legend.position = "none",
        plot.margin = unit(c(1,1,2,1), "lines")) +
  geom_jitter(width = .1, size = 3)

## Estimated diversity ####
plot_div <- insect_response_YR 
plot_div %>% 
  ggplot(aes(RestorationCategory, Shannon_asy, fill=RestorationCategory)) +
  geom_boxplot(outlier.alpha = 0) +
  theme_classic()+
  scale_fill_manual(values = c(palette4),
                    name = "Site Categories")+
  labs(title =l "")+
  xlab("\n Site category") +
  ylab("Insect family diversity")+
  theme(plot.title = element_text( size = 18, hjust = 0.5, vjust = 2),
        axis.title.x = element_text( size = 15, vjust = 4),
        axis.title.y = element_text( size = 15, hjust = 0.5, vjust =4),
        axis.text.x =  element_text( size = 12),
        axis.text.y = element_text( size = 12),
        legend.title = element_blank(),
        legend.position = "none",
        plot.margin = unit(c(1,1,2,1), "lines")) +
  geom_jitter( width = .1, size = 3)


## Time since fire ####
## Diversity
burn_insect%>% 
  ggplot(aes(years_since_fire, Rich_est, color = RestorationCategory)) +
  geom_point(outlier.alpha = 0) +
  geom_smooth(aes(color = RestorationCategory), method = "lm") +
  theme_classic()+
  scale_color_manual(values = palette2b,
                     name = "Site Categories")+
  labs(title = "")+
  xlab("\n Years since last fire") +
  ylab("Insect diversity \n (hills number, q=1)")+
  theme(plot.title = element_text( size = 18, hjust = 0.1, vjust = 1),
        axis.title.x = element_text( size = 15, vjust = 4),
        axis.title.y = element_text( size = 15, hjust = 0.5, vjust =4),
        axis.text.x =  element_text( size = 12),
        axis.text.y = element_text( size = 12),
        legend.title = element_blank(),
        legend.position = "bottom",
        plot.margin = unit(c(1,1,2,1), "lines")) 

## Abundance
burn_insect%>% 
  ggplot(aes(years_since_fire, Abund, color = RestorationCategory)) +
  geom_point(outlier.alpha = 0) +
  geom_smooth(aes(color = RestorationCategory), method = "lm") +
  theme_classic()+
  scale_color_manual(values = palette2,
                     name = "Site categories")+
  labs(title = "")+
  xlab("\n Years since last fire") +
  ylab("Insect abundance")+
  theme(plot.title = element_text( size = 18, hjust = 0.1, vjust = 1),
        axis.title.x = element_text( size = 15, vjust = 4),
        axis.title.y = element_text( size = 15, hjust = 0.5, vjust =4),
        axis.text.x =  element_text( size = 12),
        axis.text.y = element_text( size = 12),
        legend.title = element_blank(),
        legend.position = "bottom",
        plot.margin = unit(c(1,1,2,1), "lines")) 

## Seedmix richness ####
## Diversity
seeded_insect%>% 
  ggplot(aes(no_seeded_spp, Shannon_asy, color = RestorationCategory)) +
  geom_point(outlier.alpha = 0) +
  geom_smooth(aes(color = RestorationCategory), method = "lm") +
  theme_classic()+
  scale_color_manual(values = palette2b,
                     name = "Site Categories")+
  labs(title = "")+
  xlab("\n Seedmix richness") +
  ylab("Insect diversity \n (hills number, q=1)")+
  theme(plot.title = element_text( size = 18, hjust = 0.1, vjust = 1),
        axis.title.x = element_text( size = 15, vjust = 4),
        axis.title.y = element_text( size = 15, hjust = 0.5, vjust =4),
        axis.text.x =  element_text( size = 12),
        axis.text.y = element_text( size = 12),
        legend.title = element_blank(),
        legend.position = "bottom",
        plot.margin = unit(c(1,1,2,1), "lines")) 

## Abundance
seeded_insect%>% 
  ggplot(aes(no_seeded_spp, Abund, color = RestorationCategory)) +
  geom_point(outlier.alpha = 0) +
  geom_smooth(aes(color = RestorationCategory), method = "lm") +
  theme_classic()+
  scale_color_manual(values = palette2b,
                     name = "Site categories")+
  labs(title = "")+
  xlab("\n Seedmix richness") +
  ylab("Insect abundance")+
  theme(plot.title = element_text( size = 18, hjust = 0.1, vjust = 1),
        axis.title.x = element_text( size = 15, vjust = 4),
        axis.title.y = element_text( size = 15, hjust = 0.5, vjust =4),
        axis.text.x =  element_text( size = 12),
        axis.text.y = element_text( size = 12),
        legend.title = element_blank(),
        legend.position = "bottom",
        plot.margin = unit(c(1,1,2,1), "lines")) 


# Models ####
## Possible variables to include: Restoration Category, plant spp richness, native species richness, seeded species rich , size, and age.  Note that age can only be used for restored sites and not for any analysis that includes remnant

mod <- aov(Shannon_iNext ~ RestorationCategory, data = insect_response_YR)
summary(mod)
TukeyHSD(mod)

mod <- aov(Abund ~ RestorationCategory, data = insect_response_YR)
summary(mod)
TukeyHSD(mod)

#### Abundance ####
mod_abund <- lm(Abund ~ RestorationCategory + nativerichness, data = insect_response_YR)
ggResidpanel::resid_panel(mod_abund)

car::Anova(mod_abund,type="II")
summary(mod_abund)


#### Richness ####
mod_rich <- lm(Rich_est ~ RestorationCategory + nativerichness+ richness, data = insect_response_YR)
ggResidpanel::resid_panel(mod_rich)

car::Anova(mod_rich,type="II")
summary(mod_rich)

#### Diversity ####
mod_div <- lm(Shannon_iNext ~ RestorationCategory + nativerichness, data = insect_response_YR)
ggResidpanel::resid_panel(mod_div)

car::Anova(mod_div,type="II")
summary(mod_div)

#### Time since fire ####
#Abundance
mod_since_fire <- lm(Abund ~ RestorationCategory + years_since_fire, data = burn_insect)
ggResidpanel::resid_panel(mod_since_fire)
car::Anova(mod_since_fire,type="II")
summary(mod_since_fire)

#Diversity
mod_since_fire <- lm(Shannon_iNext ~ RestorationCategory + years_since_fire, data = burn_insect)
ggResidpanel::resid_panel(mod_since_fire)

car::Anova(mod_since_fire,type="II")
summary(mod_since_fire)

#### Seed mix richness ####
#Abundance 
mod_seed <- lm(Abund ~ RestorationCategory + nativerichness + no_seeded_spp, data = seeded_insect)
ggResidpanel::resid_panel(mod_seed)
car::Anova(mod_seed,type="II")
summary(mod_seed)

#Diversity
mod_seed <- lm(Shannon_iNext ~ RestorationCategory + nativerichness + no_seeded_spp, data = seeded_insect)
ggResidpanel::resid_panel(mod_seed)
car::Anova(mod_seed,type="II")
summary(mod_seed)

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
