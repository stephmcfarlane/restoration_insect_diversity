# January 19, 2024
# by: Stephanie McFarlane
# Final analyses for manuscript

library(tidyverse)
library(hablar)
library(vegan)
library(devtools)
library(BiodiversityR)

# Color palettes ####
palette4<-  c("#767171", "#A9D18E", "#548235", "#2A4117") 

# Site History Data ####
rest_history <- read_csv("../vegetation_outcomes_of_restoration/clean/rest_history.csv") %>% 
  select(EasementID = SiteID, RestorationCategory) %>% 
  mutate(RestorationCategory = factor(RestorationCategory, levels = c("No Seed", "Seed", "Seed + Fire", "Remnant")))


#Insect Data ####
data <- read_csv("clean/all_data.csv") %>% 
  unite(., "x", EasementID, Year, sep = "_") %>%  #removing random year
  filter(x != "0078S_2020") %>% 
  filter(x !="007C6_2020") %>% 
  filter(x !="007XH_2019") %>% 
  filter(x !="00MBG_2020") %>% 
  filter(x !="00MJC_2020") %>% 
  filter(x !="00NCQ_2019") %>% 
  filter(x !="00WV0_2020") %>% 
  filter(x !="0129G_2019") %>% 
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
# CCA - Constrained Correspondance Analysis --------------------------------------------------
##must use a site by family matrix
presence_site_by_family <- sitebyfam %>% 
  replace(.>1, 1)

#In order to color code by restoration category, we need to assign them into groups. To do this, you need to know the order of them as they reside in the matrix.
matrix_env <- sitebyfam %>% 
  rownames_to_column(var = "rowname") %>% 
  separate(., rowname, c("EasementID", "Year"), sep = "_") %>% 
  left_join(rest_history) %>% 
  unite(., "id", EasementID, Year, sep = "_") %>% 
  column_to_rownames(var = "id")

site_groups <- matrix_env$RestorationCategory


## Simple CCA ####
## defines CCA plot
cca_insect <-  cca(sitebyfam ~ RestorationCategory, data = matrix_env)

summary(cca_insect)

ordiplot(cca_insect,type="n",main="CCA Comparing Restored\n and Remnant Prairies")
points(cca_insect,display="sites",pch=16,col= palette4[site_groups])

#This is the line for the complex hulls.
ordihull(cca_insect, site_groups)

#This is simple legend.
legend("bottomright",y=NULL,expression(italic("No Seed"),italic("Seed"),italic("Seed + Fire"), "Remnant"),pch=16,col=palette4 ,cex=0.8)


## CCA with vectors showing families driving patterns####

## determines which families are driving patterns observed on the CCA
fit <- envfit(sitebyfam~RestorationCategory, data = matrix_env)
family.fit <- envfit(cca_insect, sitebyfam, permutations = 999)
head(family.fit)

## Plotting CCA of sites, with vectors of insect families driving the pattern 

plot(cca_insect, main="Insect CCA", type="n") # Creating empty scatter plot for ordination type="n" prevents automatic plotting
points(cca_insect,display="sites", col= palette4[site_groups], cex=1.2, pch=21, bg=palette4[site_groups]) #plots sites on scatterplot, with color (col=) corresponding to restoration category. cex= sets size of points, pch=2 sets assigns points to filled circles
ordihull(cca_insect, site_groups)
with(matrix_env, legend("topleft", legend=levels(RestorationCategory), col=palette4, pch=21, pt.bg=palette4)) #adds the legend
plot(family.fit, p.max = 0.05, col = "black", cex = 0.5) #adds vectors for environmental fitting to the CCA plot. p.max = sets the maximum p-value from permutation of families driving pattern











