######################################################################################
#                     Load Libraries                                                 #
######################################################################################

### load libraries
library(tidyverse)
library(vegan)
library(ggplot2)
library(dplyr)
library(grid)
library(gridExtra)


##########################################################################################


#setwd

#upload insect data
TPE.AllData <- read_csv("2.14.21.TPE.AllData.csv")


#upload for the sites and the restoration category
TPE.Sites <- read_csv("2.6.21.TPE.CollectionSites.csv")
TPE.Sites = subset(TPE.Sites, select = (Category:EasementID))

#Upload the grassland info
Grass <- read_csv("grassland1.5km.csv")

#Link all the data together
Summary <- full_join(TPE.Sites, Grass, by = "EasementID")



######################################################################################
#                     Family x Site Matrix, Proportion                               #
######################################################################################

# create the site x family matrix
FamilyMatrix <- TPE.AllData %>% 
  pivot_wider(
    names_from = Family,
    values_from = Total)

# site x species matrix with one row per site
FamilyMatrix <- FamilyMatrix %>% 
  mutate_at(vars(Formicidae:Melyridae), ~replace(., is.na(.), 0)) %>%
  group_by(EasementID) %>%
  summarise_at(vars(Formicidae:Melyridae), sum)%>%
  slice(1:25)

#family matrix now has the abundance of each family at each site. 
#Now need to get the proportions of each at a site. This is where I struggled to make a proportion table in r,
#so i exported it to excel and created the proportion table there easily. 


#this exports the current family matrix to excel
#write.csv(FamilyMatrix, "D:/OneDrive/Documents/r.Studio/TPE/FamilyMatrix.csv")
#then I created the prop table
#export back into r
AbundProp <- read_csv("FamilyMatrix.csv")







######################################################################################
#                     Sort sites from largest to smallest                            #
#                                Rank abundances                                     #    
######################################################################################

#So now I need to order each site from highest to lowest proportional abundance. Everytime I tried order it said there
#was an error in the vectors. Instead, I isolated the site into its own df, exported it, reordered it in excel, 
#and then brought it back. Then I had a df for each site that had the family proportions in descending order 
#so then I could add them to the rank abundance. Definitely more of a pain but it wasn't working for me

#Since all of these have been done, the only line that needs running are the ones titled "all.remnant or all.sf" 
#because they are compiled of all the sites, and compiled rank at the very end



##########################################################################

# all of the controls : 00MDP, 007CQ, 00NCQ, 007R5, 007XH

##########################################################################


##############################################################007R5

#I annotated this one so that you see all the steps, then repeated this for all the sites. 


#get site alone with a new df and then slicing out the other sites

Site.007R5 <- AbundProp[-(1:7),]
Site.007R5 <- Site.007R5[-(2:18),]

#change the df from horizontal to verticle
site.007R5 <- t(Site.007R5)

#export to excel to sort them into descending order
#write.csv(Site.007R5,"D:/OneDrive/Documents/r.Studio/TPE/Site.007R5.csv")



#bring back
Site.Sort.007R5 <- read.csv("Site.007R5.csv")

#remove zero abundance families so that there arn't empty cases in the rank abundance. 
Site.Sort.007R5 <- Site.Sort.007R5[-c(75:165),]


##############################################################00MDP
site.00mdp <- AbundProp[-c(1:10),]
site.00mdp <- site.00mdp[-c(2:15),]

site.00mdp <- t(site.00mdp)

write.csv(site.00mdp,"D:/OneDrive/Documents/r.Studio/TPE/site.00mdp.csv")


#only run this
Site.Sort.00MDP <- read.csv("site.00mdp.csv")

Site.Sort.00MDP <- Site.Sort.00MDP[-c(78:165),]


##############################################################0078Z
site.007CQ <- AbundProp[-c(1:5),]
site.007CQ <- site.007CQ[-c(2:20),]

site.007CQ <- t(site.007CQ)

write.csv(site.007CQ,"D:/OneDrive/Documents/r.Studio/TPE/site.007CQ.csv")

Site.Sort.007CQ <- read.csv("site.007CQ.csv")

Site.Sort.007CQ <- Site.Sort.007CQ[-c(75:165),]


##############################################################00NCQ
site.00NCQ <- AbundProp[-c(1:14),]
site.00NCQ <- site.00NCQ[-c(2:11),]

site.00NCQ <- t(site.00NCQ)

write.csv(site.00NCQ,"D:/OneDrive/Documents/r.Studio/TPE/site.00NCQ.csv")

Site.Sort.00NCQ <- read.csv("site.00NCQ.csv")

Site.Sort.00NCQ <- Site.Sort.00NCQ[-c(92:165),]


##############################################################007XH
site.007XH <- AbundProp[-c(1:9),]
site.007XH <- site.007XH[-c(2:16),]

site.007XH <- t(site.007XH)

write.csv(site.007XH,"D:/OneDrive/Documents/r.Studio/TPE/site.007XH.csv")

Site.Sort.007XH <- read.csv("site.007XH.csv")

Site.Sort.007XH <- Site.Sort.007XH[-c(81:165),]


#############################################################################################################

#now I merged each of site df into one in excel because i took out the families and just complied the sites in their
#descending proportional abudnances 

#00MDP, 007CQ, 00NCQ, 007R5, 007XH

All.Controls <- read.csv("All.Controls.csv")


# for the rank abundancd plot:

All.Controls[All.Controls == 0] <- NA

rankabunof2 <-plot (All.Controls$X007R5, type ="o", xlim= c(1,60), ylim=c(0,0.4),
                    lwd =2, col="red",pch=1, cex=2, xlab="rank", ylab="prop of abund", 
                    main="All the controls") +
  lines(All.Controls$X007XH,type="o",lwd=2,col="blue",pch=1,cex=2)+
  lines(All.Controls$X00NCQ,type="o",lwd=2,col="yellow",pch=1,cex=2)+  
  lines(All.Controls$X007CQ,type="o",lwd=2,col="green",pch=1,cex=2)+
  lines(All.Controls$X00MDP,type="o",lwd=2,col="black",pch=1,cex=2)







###########################################################################


# all of the SF : ,007V7,,00MJC,007F5,0077R,,,00WV0,00792,00MJ8,00MJ9,,,,007BY,00WY5


##############################################################007V7
site.007V7 <- AbundProp[-c(1:8),]
site.007V7 <- site.007V7[-c(2:17),]

site.007V7 <- t(site.007V7)

write.csv(site.007V7,"D:/OneDrive/Documents/r.Studio/TPE/site.007V7.csv")

Site.Sort.007V7 <- read.csv("site.007V7.csv")

Site.Sort.007V7 <- Site.Sort.007V7[-c(61:165),]


##############################################################00MJC
site.00MJC <- AbundProp[-c(1:13),]
site.00MJC <- site.00MJC[-c(2:12),]

site.00MJC <- t(site.00MJC)

write.csv(site.00MJC,"D:/OneDrive/Documents/r.Studio/TPE/site.00MJC.csv")

Site.Sort.00MJC <- read.csv("site.00MJC.csv")

Site.Sort.00MJC <- Site.Sort.00MJC[-c(59:165),]


##############################################################00WY5
site.00WY5 <- AbundProp[-c(1:16),]
site.00WY5 <- site.00WY5[-c(2:9),]

site.00WY5 <- t(site.00WY5)

write.csv(site.00WY5,"D:/OneDrive/Documents/r.Studio/TPE/site.00WY5.csv")

Site.Sort.00WY5 <- read.csv("site.00WY5.csv")

Site.Sort.00WY5 <- Site.Sort.00WY5[-c(49:165),]


##############################################################007F5
site.007F5 <- AbundProp[-c(1:6),]
site.007F5 <- site.007F5[-c(2:19),]

site.007F5 <- t(site.007F5)

write.csv(site.007F5,"D:/OneDrive/Documents/r.Studio/TPE/site.007F5.csv")

Site.Sort.007F5 <- read.csv("site.007F5.csv")

Site.Sort.007F5 <- Site.Sort.007F5[-c(54:165),]


##############################################################007BY
site.007BY<- AbundProp[-c(1:4),]
site.007BY <- site.007BY[-c(2:21),]

site.007BY<- t(site.007BY)

write.csv(site.007BY,"D:/OneDrive/Documents/r.Studio/TPE/site.007BY.csv")

Site.Sort.007BY <- read.csv("site.007BY.csv")

Site.Sort.007BY <- Site.Sort.007BY[-c(49:165),]



##############################################################00WV0
site.00WV0 <- AbundProp[-c(1:15),]
site.00WV0 <- site.00WV0[-c(2:10),]

site.00WV0 <- t(site.00WV0)

write.csv(site.00WV0,"D:/OneDrive/Documents/r.Studio/TPE/site.00WV0.csv")

Site.Sort.00WV0 <- read.csv("site.00WV0.csv")

Site.Sort.00WV0 <- Site.Sort.00WV0[-c(58:165),]

##############################################################00792
site.00792 <- AbundProp[-c(1:3),]
site.00792 <- site.00792[-c(2:22),]

site.00792 <- t(site.00792)

write.csv(site.00792,"D:/OneDrive/Documents/r.Studio/TPE/site.00792.csv")

Site.Sort.00792 <- read.csv("site.00792.csv")

Site.Sort.00792 <- Site.Sort.00792[-c(76:165),]

##############################################################00MJ8
site.00MJ8 <- AbundProp[-c(1:11),]
site.00MJ8 <- site.00MJ8[-c(2:14),]

site.00MJ8 <- t(site.00MJ8)

write.csv(site.00MJ8,"D:/OneDrive/Documents/r.Studio/TPE/site.00MJ8.csv")

Site.Sort.00MJ8 <- read.csv("site.00MJ8.csv")

Site.Sort.00MJ8 <- Site.Sort.00MJ8[-c(53:165),]

##############################################################00MJ9
site.00MJ9 <- AbundProp[-c(1:12),]
site.00MJ9 <- site.00MJ9[-c(2:13),]

site.00MJ9 <- t(site.00MJ9)

write.csv(site.00MJ9,"D:/OneDrive/Documents/r.Studio/TPE/site.00MJ9.csv")

Site.Sort.00MJ9 <- read.csv("site.00MJ9.csv")

Site.Sort.00MJ9 <- Site.Sort.00MJ9[-c(56:165),]



#############################################################################################################

#now merge all the sites to one

#00MDP, 007CQ, 00NCQ, 007R5, 007XH

All.SF <- read.csv("All.SF.csv")


# for the plot:

All.SF[All.SF == 0] <- NA

rankabunof2 <-plot (All.SF$X007V7, type ="o", xlim= c(1,35), ylim=c(0,0.4),
                    lwd =2, col="red",pch=1, cex=2, xlab="rank", ylab="prop of abund", 
                    main="All the SF") +
  lines(All.SF$X00MJC,type="o",lwd=2,col="blue",pch=1,cex=2)+
  lines(All.SF$X007F5,type="o",lwd=2,col="yellow",pch=1,cex=2)+  
  lines(All.SF$X0077R,type="o",lwd=2,col="green",pch=1,cex=2)+
  lines(All.SF$X00WV0,type="o",lwd=2,col="pink",pch=1,cex=2)+
  lines(All.SF$X00792,type="o",lwd=2,col="purple",pch=1,cex=2)+
  lines(All.SF$X00MJ8,type="o",lwd=2,col="orange",pch=1,cex=2)+
  lines(All.SF$X00MJ9,type="o",lwd=2,col="brown",pch=1,cex=2)+
  lines(All.SF$X00WY5,type="o",lwd=2,col="brown",pch=1,cex=2)+
  lines(All.SF$X007BY,type="o",lwd=2,col="brown",pch=1,cex=2)


##################################################################################






##################################################################################

# all of the remnant : oliver, westport, black earth, borah, pleasant


##############################################################oliver
site.oliver <- AbundProp[-c(1:22),]
site.oliver <- site.oliver[-c(2:3),]

site.oliver <- t(site.oliver)

write.csv(site.oliver,"D:/OneDrive/Documents/r.Studio/TPE/site.oliver.csv")

Site.Sort.oliver <- read.csv("site.oliver.csv")

Site.Sort.oliver <- Site.Sort.oliver[-c(27:165),]


##############################################################westport
site.westport <- AbundProp[-c(1:24),]

site.westport <- t(site.westport)

write.csv(site.westport,"D:/OneDrive/Documents/r.Studio/TPE/site.westport.csv")

Site.Sort.westport<- read.csv("site.westport.csv")

Site.Sort.westport <- Site.Sort.westport[-c(73:165),]


##############################################################black
site.black <- AbundProp[-c(1:20),]
site.black <- site.black[-c(2:5),]

site.black <- t(site.black)

write.csv(site.black,"D:/OneDrive/Documents/r.Studio/TPE/site.black.csv")

Site.Sort.black <- read.csv("site.black.csv")

Site.Sort.black <- Site.Sort.black[-c(69:165),]



##############################################################borah
site.borah <- AbundProp[-c(1:21),]
site.borah <- site.borah[-c(2:4),]

site.borah <- t(site.borah)

write.csv(site.borah,"D:/OneDrive/Documents/r.Studio/TPE/site.borah.csv")

Site.Sort.borah<- read.csv("site.borah.csv")

Site.Sort.borah <- Site.Sort.borah[-c(53:165),]

##############################################################pleasant
site.pleasant <- AbundProp[-c(1:23),]
site.pleasant <- site.pleasant[-c(2),]

site.pleasant <- t(site.pleasant)

write.csv(site.pleasant,"D:/OneDrive/Documents/r.Studio/TPE/site.pleasant.csv")

Site.Sort.pleasant <- read.csv("site.pleasant.csv")

Site.Sort.pleasant <- Site.Sort.pleasant[-c(57:165),]


#############################################################################################################

#now merge all the sites to one


All.remnant <- read.csv("All.remnant.csv")



# for the plot:

All.remnant[All.remnant == 0] <- NA

rankabunof2 <-plot (All.remnant$Pleasant.Valley.Conservancy, type ="o", xlim= c(1,35), ylim=c(0,0.4),
                    lwd =2, col="red",pch=1, cex=2, xlab="rank", ylab="prop of abund", 
                    main="All the remnants") +
  lines(All.remnant$Oliver.Prairie,type="o",lwd=2,col="blue",pch=1,cex=2)+
  lines(All.remnant$Borah.Creek,type="o",lwd=2,col="yellow",pch=1,cex=2)+  
  lines(All.remnant$Westport.Drumlin,type="o",lwd=2,col="green",pch=1,cex=2)+
  lines(All.remnant$Black.Earth.Rettenmund,type="o",lwd=2,col="pink",pch=1,cex=2)



##################################################################################


#############################################################################################################


# now to average all the lines within each site

#To do this, in excel, I just merged all the sites, and then avaeraged by category the abundance of each family. 

#Imported back into r
compiledrank <- read.csv("rank abun compiled.csv")




#this then creates the compiled graph
compiledrank[compiledrank == 0] <- NA

rankabunof2 <-plot (compiledrank$SF, type ="o", xlim= c(1,20), ylim=c(0,0.2),
                    lwd =4, col="#5F92C0",pch=1, cex=2, xlab="Abundance Rank", ylab="Average proptions of abundance", 
                    main="") +
  lines(compiledrank$Con,type="o",lwd=4,col="#DEEB0F",pch=1,cex=2)+
  lines(compiledrank$REM,type="o",lwd=4,col="#8FE033",pch=1,cex=2)






#the colors used 
"#DEEB0F", "#8FE033", "#5F92C0"






























































#############################################################################################################
#############################################################################################################
#############################################################################################################
#############################################################################################################
#############################################################################################################
#############################################################################################################
