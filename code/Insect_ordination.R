library(tidyverse)
library(hablar)
library(vegan)
library(devtools)
library(BiodiversityR)

# Color palettes ####
palette2 <- c("#548235", "#2A4117")
palette3<-  c("#767171", "#A9D18E", "#548235") #3B5422
palette4<-  c("#767171", "#A9D18E", "#548235", "#2A4117") #3B5422

# Site History Data ####
# Upload easement treatment (aka restoration category) data & clean ##
rest_history <- read_csv("../vegetation_outcomes_of_restoration/clean/rest_history.csv") %>% 
  select(EasementID = SiteID, RestorationCategory) %>% 
  mutate(RestorationCategory = factor(RestorationCategory, levels = c("No Seed", "Seed", "Seed + Fire", "Remnant")))


#Insect Data ####
data <- read_csv("clean/all_data.csv") %>% 
  unite(., "x", EasementID, Year, sep = "_") %>%  #removing random year
  filter(x != "0078S_2020") %>% 
  filter(x !="007C6_2020") %>% 
  filter(x !="00MBG_2020") %>% 
  filter(x !="00MJC_2020") %>% 
  filter(x !="00NCQ_2019") %>% 
  filter(x !="00WV0_2020") %>% 
  separate(., x, c("EasementID", "Year"), sep = "_") %>% 
  left_join(rest_history)

## Site by Family matrix ####
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



## Environmental matrix ####
# Ordination --------------------------------------------------
##must use a site by family matrix
presence_site_by_family <- sitebyfam %>% 
  replace(.>1, 1)

## NMDS ####
#This next line runs the actual NMDS. There are more possible arguments, but the ones that are in it at the moment are the dataset, k, and the number of attempts R will do to try to reach a solution. 
nmds.insect<-metaMDS(sitebyfam,k=5,trymax=250,distance="bray") 
#The Bray-Curtis dissimilarity is based on occurrence data (abundance), while the Jaccard distance is based on presence/absence data (does not include abundance information)
#nmds.insect <- metaMDS(presence_site_by_family, k=4, trymax = 250, distance = "jaccard")

#Now we need to make some plots!
#Set plot window dimensions
par(mfcol=c(0.7,1.5))

#First is a stress plot to make sure that things look ok. If the points are widely scattered away from the generated line, the variation is not best squished into two dimensions.
stressplot(nmds.insect)
nmds.insect$stress

#Generic plot:the red crosses are the columns (insect families) and the black, unfilled circles are the rows (sites). Ordinarily, it is the rows that people are interested in, though this will obviously vary with whatever question you are asking.
plot(nmds.insect)

#In order to color code by restoration category, we need to assign them into groups. To do this, you need to know the order of them as they reside in the matrix.
matrix_env <- sitebyfam %>% 
  rownames_to_column(var = "rowname") %>% 
  separate(., rowname, c("EasementID", "Year"), sep = "_") %>% 
  left_join(rest_history) %>% 
  unite(., "id", EasementID, Year, sep = "_") %>% 
  column_to_rownames(var = "id")

site_groups <- matrix_env$RestorationCategory

#This next line generates the plotting space for the NMDS. The type argument specifies that it will be blank, which allows us to then come and plot things the way we want.
ordiplot(nmds.insect,type="n",main="NMDS Comparing Restored Prairies\n and Remnant Prairies")
points(nmds.insect,display="sites",pch=16,col= palette4[site_groups])

#This is the line for the complex hulls.
ordihull(nmds.insect, site_groups)

#This is simple legend.
legend("bottomright",y=NULL,expression(italic("No Seed"),italic("Seed"),italic("Seed + Fire"), "Remnant"),pch=16,col=palette4 ,cex=0.8)

##Plotting different axes ####
insect_envir <- matrix_env %>% 
  select(-Size)

nmds.insect<-metaMDS(sitebyfam,k=5,trymax=250,distance="bray")
insect.envfit <- envfit(nmds.insect, insect_envir, permutations = 999) # this fits environmental vectors
insect_fam_fit <- envfit(nmds.insect, sitebyfam, permutations = 999)# this fits spp vectors

##To plot the output from the mds using ggplot a new datasheet needs to be created which contains the x,y points for each site. You can do this by calling the scores of you mds.
site.scrs <- as.data.frame(scores(nmds.insect, display = "sites"))  #Saves nmds axes results into a dataframe
site.scrs <- cbind(site.scrs, RestorationCategory = insect_envir$RestorationCategory) #adds restoration category
#site.scrs <- cbind(site.scrs, RestorationCategory = insect_envir$SizeCategory)  # could add another grouping variable if desired

spp.scrs <- as.data.frame(scores(insect_fam_fit, display = "vectors"))  #Saves species intrinsic values into a dataframe

##A new dataset containing species data also needs to be made to look at species vectors. This is not necessary if you don’t want to show the species on the final graph.
spp.scrs <- cbind(spp.scrs, Family = rownames(spp.scrs)) #add family names to dataframe
spp.scrs <- cbind(spp.scrs, pval = insect_fam_fit$vectors$pvals) #add pvalues to dataframe so you can select species which are significant
sig.spp.scrs <- subset(spp.scrs, pval<=0.05)


nmds.plot.insect <- ggplot(site.scrs, aes(x=NMDS1, y=NMDS5))+
  geom_point(aes(NMDS1, NMDS5, colour = factor(RestorationCategory)), size = 2) + 
  coord_fixed()+
  theme_classic()+
  theme(panel.background = element_rect(fill = NA, colour = "black", size = 1, linetype = "solid")) +
  labs(colour = "RestorationCategory") +
  theme(legend.position = "right", legend.text = element_text(size = 12), legend.title = element_text(size = 12), axis.text = element_text(size = 10))

nmds.plot.insect
#ordihull(nmds.insect, site_groups)

nmds.plot.insect +
  geom_segment(data = sig.spp.scrs, aes(x=0, xend=NMDS1, y=0, yend=NMDS5), arrow= arrow(length = unit(0.25, "cm")), colour = "grey10", lwd=0.3) +
  ggrepel::geom_label_repel(data = sig.spp.scrs, aes(x=NMDS1, y=NMDS5, label= Family), cex = 3, direction = "both", segment.size = 0.25)

## PERMANOVA ####
adonis2(matrix_env[,c(1:179)]~ #Take only the numeric iris data
          matrix_env[,180],    #Test by species
        method="bray", #Method to convert numeric data into a distance matrix
        permutations=10000) #Number of permutations to run.


#Testing for significance of homogeneity of multivariate dispersions. If the result of betadisper is significant, the results are potentially biased. 
test.insect<-betadisper(vegdist(sitebyfam, method = "bray"),site_groups) 
anova(test.insect)
#the result is non-significant, so we assume not bias in our results
## if this test was significant, we can look pairwise differences in our dispersion

##permutest(test.insect, pairwise = T, permutations = 999)  ## if it was significant, we would run this code

## CCA - Constrained Correspondance Analysis ####

cca_insect <-  cca(sitebyfam ~ RestorationCategory, data = matrix_env)

ordiplot(cca_insect,type="n",main="CCA Comparing Restored\n and Remnant Prairies")
points(cca_insect,display="sites",pch=16,col= palette4[site_groups])

#This is the line for the complex hulls.
ordihull(cca_insect, site_groups)

#This is simple legend.
legend("bottomright",y=NULL,expression(italic("No Seed"),italic("Seed"),italic("Seed + Fire"), "Remnant"),pch=16,col=palette4 ,cex=0.8)
### Vegan Code ####

#take the square root to give less to abundant spp
# NEED TO DISCUSS W/GROUP
#community.matrix <- sqrt(community.matrix)

# create a matrix with only the family and average abundance
# these need to be whole numbers for PERMANOVA to work
fam_matrix <- sitebyfam %>%
  mutate_all(round, 0) #round up to nearest whole number

# site (aka environmental-this is any independent variable) variables in separate df
env_ %>% matrix <- sitebyfam %>%
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
  scale_fill_manual(values = palette3) +
  scale_color_manual(values = palette3) +
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
relative_abund <- data %>% 
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
















