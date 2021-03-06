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
  

################################################################################
# community composition analyses from Bombus stuff. NOT cleaned or annotated so beware :)

# remove extraneous information. 
community.matrix <- sitexspp  %>%
  dplyr::select(BOMAFF, BOMAUR, BOMBIM, BOMBOR, BOMCIT, BOMFER, BOMGRI, BOMIMP, BOMRUF, VAGSAN)

#take the square root to give less to abundant spp
community.matrix <- sqrt(community.matrix)

#round up to whole numbers
community.matrix <- community.matrix %>%
  mutate(BOMAFF = ceiling(BOMAFF),
         BOMAUR = ceiling(BOMAUR),
         BOMBIM = ceiling(BOMBIM),
         BOMBOR = ceiling(BOMBOR),
         BOMCIT = ceiling(BOMCIT),
         BOMFER = ceiling(BOMFER),
         BOMGRI = ceiling(BOMGRI),
         BOMIMP = ceiling(BOMIMP),
         BOMRUF = ceiling(BOMRUF),
         VAGSAN = ceiling(VAGSAN))

# env variables in separate df
community.env <- sitexspp %>%
  dplyr::select(EasementID, RestorationCategory)

community.env <- inner_join(community.env, bombus)
community.env <- community.env %>%
  select(c(EasementID, RestorationCategory, percent.natural, Year))

spp.sum <- rbind(community.matrix, colSums(community.matrix))
spp.sum$N <- rowSums(spp.sum)

# remove species that made up less than 5% of total 
community.matrix <- community.matrix %>%
  dplyr::select(-c(BOMAFF, BOMBOR, BOMFER, BOMCIT))

## PERMANOVA
perm <- adonis2(community.matrix ~ RestorationCategory + percent.natural + Year, data = community.env, method="bray", by="margin", permutations=9999, model = "reduced")
#model = "reduced" determines the method of permuations. 
#This permutes the residuals under a reduced model

# view results
perm

### Pairwise

# pairwise permanova comparisons between restoration categories
library(pairwiseAdonis)
pairwise.adonis(community.matrix, community.env$RestorationCategory,sim.function = "vegdist", sim.method = "bray")

### Dispersion
bray2 <- vegdist(community.matrix, method = "bray")

dispersion<-betadisper(bray2, group=community.env$RestorationCategory)
dispersion
permutest(dispersion, pairwise = T) #p=0.003 so groups have different dispersions

# Figures: NMDS
## create the nmds with transformed data
nmds <- metaMDS(community.matrix, distance = "bray")
plot(nmds)

community.envfit <- community.env[-c(1)]

(fit <- envfit(nmds, community.envfit, perm = 999))
head(fit)
scores(fit, "vectors")

env.scores <- as.data.frame(scores(fit, display = "vectors"))
env.scores <- cbind(env.scores, Species = rownames(env.scores))

plot(nmds)
plot(fit, col="black")

## Extract Scores
#extract NMDS scores (x and y coordinates)
#ggplot needs these in a df
nmds.scores <- as.data.frame(scores(nmds))

#add identifying columns to dataframe
nmds.scores$EasementID <- community.env$EasementID
nmds.scores <- inner_join(nmds.scores, community.env)

#### Basic Plot
# create hulls
grp.control <- nmds.scores[nmds.scores$RestorationCategory == "Control", ][chull(nmds.scores[nmds.scores$RestorationCategory == 
                                                                                               "Control", c("NMDS1", "NMDS2")]), ]  # hull values for grp A

grp.seeded <- nmds.scores[nmds.scores$RestorationCategory == "Seeded", ][chull(nmds.scores[nmds.scores$RestorationCategory == 
                                                                                             "Seeded", c("NMDS1", "NMDS2")]), ]  # hull values for grp B

grp.burned <- nmds.scores[nmds.scores$RestorationCategory == "Burned", ][chull(nmds.scores[nmds.scores$RestorationCategory == 
                                                                                             "Burned", c("NMDS1", "NMDS2")]), ]  # hull values for grp C

hull.data <- rbind(grp.control, grp.seeded, grp.burned)  #combine groups
hull.data

#specify order of categories 
color.nmds <- factor(nmds.scores$RestorationCategory, levels = c("Control", "Seeded", "Burned"))
color.hull <- factor(hull.data$RestorationCategory, levels = c("Control", "Seeded", "Burned"))

# grey scale 
## scale_color_manual(values = c("#252525", "#969696", "#636363", "#252525"))
## scale_shape_manual(values = c(1,19,19,15))

# grey and green
## scale_color_manual(values = c("#767171", "#A9D18E", "#548235"))
## scale_shape_manual(values = c(19,19,19,15))

p <- ggplot() + 
  geom_polygon(data=hull.data,aes(x=NMDS1,y=NMDS2, group=color.hull),alpha=0.2) + 
  geom_point(data = nmds.scores, aes(x=NMDS1,y=NMDS2, colour=color.nmds, shape = color.nmds), size = 2) + # add the point markers
  scale_fill_manual(values = c("#767171", "#A9D18E", "#548235")) +
  scale_color_manual(values = c("#767171", "#A9D18E", "#548235")) +
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

p

legend.restore <- get_legend(p)

# relative abundance
# filter out nest searching
relative <- timed %>% 
  filter(SppCode != "n/a") %>%
  filter(!is.na(SppCode)) 

#order the categories so they show up specifically on the x-axis.                                 
order2 <- factor(relative$RestorationCategory, levels = c("Control", "Seeded", "Burned"))

relative.abund <- ggplot(relative, aes(x = order2, fill = SppCode)) + 
  geom_bar(position = "fill", width = 0.75) + 
  theme_bw() + 
  ylab("Relative Abundance") +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 12),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 12)) +
  scale_fill_brewer(palette = "BrBG")

relative.abund

#ggsave("./Figures/relative_abundance_ESA.jpeg",plot = relative.abund, width = 7, height = 5)