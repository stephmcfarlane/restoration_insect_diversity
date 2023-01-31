#Insect Data ####
all_data <- read_csv("../restoration_insect_diversity/clean/all_data.csv")
#Data with aggragate families removed#### 
## Use for diversity estimate



#site 00VTR may be an outlier in the data...but not filtering out aggregating families
data <- all_data %>% 
  group_by(EasementID, Year) %>% 
  mutate(site_total = sum(Total)) %>% 
  mutate(percent_total_ind = Total/site_total)

## iNEXT: Interpolation Exterpolation####
#### Formatting Data into a List of Lists ####
plot_ids <- unique(data[c("EasementID", "Year")])  ##creates unique list of plot-year combo

# output list of lists, start with blank list
abundance_list <-  list()

list_names <- mutate(plot_ids, paste0(EasementID, "_", Year)) %>%
  select(names = 3) %>% ##naming 3rd column names
  pull() ##changing to a vector of names which is needed to rename lists in abundance list

## for loop -- converting species' abundances into a list for each site/date combo and adding to blank list, abundance_list
for(i in 1:nrow(plot_ids)) 
{
  df.out <- data %>% 
    select(EasementID, Year, Sample, Total, Family) %>% 
    group_by(EasementID, Year, Family) %>% 
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
# N = 120 #### 
# Determining estimated species richness at N=120
insect_rich_120 <- iNEXT(abundance_list, q=0, datatype = "abundance", endpoint = 120, knots= 100, nboot = 100) 
estimates_120 <- as.data.frame(insect_rich_120$AsyEst)
write_csv(insect_rich_120$AsyEst, "clean/estimates_120.csv")

insect_rich_est_120 <- as.data.frame(insect_rich[["iNextEst"]])
# output into separate data frames
# with rarefied and extrapolated samples for q = 0,1,2
insect_abundance_120 <- as.data.frame(insect_rich[["DataInfo"]])  ##site name, sample size, Family richness, number of singletons, doubletons, etc.

# Specifying Richness and Shannon's only, q = 0 & 1 ####
insect <- iNEXT(abundance_list, q=c(0, 1), datatype = "abundance", nboot = 100) ##can set endpoint too, referring to number of indivduals 

## output into separate data frames
insect_abundance <- as.data.frame(insect[["DataInfo"]])  ##site name, sample size, Family richness, number of singletons, doubletons, etc.
#insect_div_est <- as.data.frame(insect[["iNextEst"]]) ## diversity estimates with rarefied and extrapolated samples for q = 0,1,2
insect_div_asy <- as.data.frame(insect[["AsyEst"]]) ## asympototic diversity estimates with bootstrap s.e. and confidence intervals for q = 0,1,2

#write_csv(insect$DataInfo, "clean/insect_datainfo.csv")
#write_csv(insect$AsyEst, "clean/insect_div_asym.csv")
#write_csv(insect_abundance, "clean/insect_abundance.csv")

# N = 450 #### 
#determining estimated species richness at N=450
insect_rich <- iNEXT(abundance_list, q=0, datatype = "abundance", endpoint = 450, knots= 100, nboot = 100) 

insect_rich_est <- as.data.frame(insect_rich[["iNextEst"]])
## output into separate data frames
insect_div_est <- as.data.frame(insect_rich[["iNextEst"]]) ## diversity estimates with rarefied and extrapolated samples for q = 0,1,2
insect_abundance <- as.data.frame(insect_rich[["DataInfo"]])  ##site name, sample size, Family richness, number of singletons, doubletons, etc.

#write_csv(insect_rich$DataInfo, "clean/abundance_withAgg.csv")
#write_csv(insect_rich$AsyEst, "clean/estimates.csv")
estimates <- as.data.frame(insect_rich$AsyEst)



