## Name: Isty Rysava
## Date: 6/4/2024
## Code: Takes epi and climate variable data and pre-processes for regression analysis 

rm(list=ls())

### Libraries
library(tidycensus)
library(tidyverse)
library(stringr)
library(SpatialEpi)
library(sf)
library(harrypotter)

## Roger Bivand, Nick Bearman, Nicholas Lewin-Koh 
leglabs <- function(vec, under = "under", over = "over", between = "-", 
                    reverse = FALSE) {
  x <- vec
  lx <- length(x)
  if (lx < 3) 
    stop("vector too short")
  if (reverse) {
    x <- rev(x)
    under <- "over"
    over <- "under"
  }
  res <- character(lx - 1)
  res[1] <- paste(under, x[2])
  for (i in 2:(lx - 2)) res[i] <- paste(x[i], between, x[i + 
                                                           1])
  res[lx - 1] <- paste(over, x[lx - 1])
  res
}

### Data
fl <- read_sf("data/shp/Florida_adm1.shp", layer="Florida_adm1")
vars <- read.csv(file="output/variables/acs_dat_2009_2022.csv")
density <- read.csv("output/variables/pop_density_2009_2022.csv")
climall <- read.csv("output/variables/FL_climate_all.csv")
clim <- read.csv("output/variables/FL_summaries_a.csv")

### Plot mean alpha across counties
summary(clim$meanA)
breaks <- seq(0, 0.5, by=.05)
pal <- hp(n = length(breaks), house = "NewtScamander")
plot(1:length(breaks), 1:length(breaks), col=pal, pch=16, cex=3)
# match with the shapefile

setdiff(clim$county, fl$NAME_2)
clim$county[which(clim$county=="St. Johns")] <- "Saint Johns"
clim$county[which(clim$county=="St. Lucie")] <- "Saint Lucie"
fl$NAME_2[which(fl$NAME_2=="Miami-Dade")] <- "MiamiDade"
fl$meanA <- NA
meanA.inds=match(clim$county, fl$NAME_2)
fl$meanA[meanA.inds] <- clim$meanA

## draw map
pdf("figs/Mean_a.pdf", height=5, width=6)
plot(st_geometry(fl), col="white",border=NA)
plot(st_geometry(fl), add=T, border=NA, col=pal[findInterval(as.vector(as.character(fl$meanA)),
                                                             breaks,all.inside=T)])
plot(st_geometry(fl[which(fl$NAME_2 %in% c("Bay", "Franklin", "Okaloosa")),]), add=T, col="black")
plot(st_geometry(fl), add=T, border="grey70", lwd=.1)
legend("bottomleft", 
       legend=c(leglabs(breaks, under="<",over=">"), "no cases", "poor estimates"),
       pch=c(rep(22,(length(breaks)+1))),
       lty=c(rep(NA, 12)),
       col = c(rep("black",(length(breaks)+1))),
       pt.bg=c(pal[-length(pal)], "white", "black"), pt.cex=2, lwd=.1,
       title=c("Mean contact (a)"),
       cex=.8, bty="n")
dev.off()


### Check that all variables have unique names!
head(vars)
unique(vars$label)
unique(vars$concept[vars$label=="Estimate!!Total!!Male"])                                                                       
unique(vars$concept[vars$label=="Estimate!!Total!!Female"])                                                                     
unique(vars$concept[vars$label=="Estimate!!Total"])
vars$label[which(vars$concept=="TOTAL POPULATION" & vars$label=="Estimate!!Total")] <- "TotPop"
vars$label[which(vars$concept=="HOUSEHOLD TYPE BY RELATIVES AND NONRELATIVES FOR POPULATION IN HOUSEHOLDS" & vars$label=="Estimate!!Total")] <- "HHtyperelatives"
vars$label[which(vars$concept=="OCCUPANCY STATUS" & vars$label=="Estimate!!Total")] <- "TotOccup"
vars$label[which(vars$concept=="TOTAL POPULATION IN OCCUPIED HOUSING UNITS BY TENURE" & vars$label=="Estimate!!Total")] <- "TotPopByTenure"
vars$label[which(vars$concept=="AVERAGE HOUSEHOLD SIZE OF OCCUPIED HOUSING UNITS BY TENURE" & vars$label=="Estimate!!Total")] <- "AvrHHsizeByTenure"
vars$label[which(vars$concept=="PLUMBING FACILITIES FOR ALL HOUSING UNITS" & vars$label=="Estimate!!Total")] <- "TotPlumb"
vars$label[which(vars$concept=="TENURE BY OCCUPANTS PER ROOM" & vars$label=="Estimate!!Total")] <- "OccupPerRoom"

# remove family household
vars <- vars[-which(vars$label=="HHtyperelatives"),]
vars <- vars[-which(vars$label=="Estimate!!Total!!In nonfamily households"),]
vars <- vars[-which(vars$label=="Estimate!!Total!!In family households"),]

### Average over years
vars_avr <- tibble(vars)%>% 
  mutate(GEOID = as.factor(GEOID)) %>%
  mutate(NAME = as.factor(NAME)) %>%
  group_by(GEOID, NAME,label) %>%
  summarize(estimate = mean(estimate)) 

density_avr <- tibble(density)%>% 
  mutate(GEOID = as.factor(GEOID)) %>%
  group_by(GEOID) %>%
  summarize(density = mean(density)) 

# match
ind <- match(density_avr$GEOID, vars_avr$GEOID)
density_avr$NAME <- vars_avr$NAME[ind]

### Sort out the demo/epi names & think about how to turn census into meaningful variables 
## (e.g., % of occupied vs. not occupied, immigrant, etc.)
vars_avr$label <- gsub("Estimate!!Total!!", "", vars_avr$label)
vars_avr$label <- gsub("Estimate!!", "", vars_avr$label)
unique(vars_avr$label)
temp <- vars_avr %>% pivot_wider(names_from = label, values_from = estimate)

# Total Pop
plot(temp$Female+temp$Male, temp$TotPop)
temp <- temp %>% select(!(c(Male, Female)))

# Occupancy
temp$Poccup <- temp$Occupied/(temp$Occupied+temp$Vacant)
temp$Pvac <- temp$Vacant/(temp$Occupied+temp$Vacant)
temp$Powned <- temp$`Owner occupied`/(temp$`Owner occupied`+ temp$`Renter occupied`)
temp$Prented <- temp$`Renter occupied`/(temp$`Owner occupied`+ temp$`Renter occupied`)
# plot(temp$TotPopByTenure, (temp$Occupied+temp$Vacant))

# U.S. per origin
tot <- temp$`Not a U.S. citizen` + 
  temp$`U.S. citizen by naturalization` +
  temp$`U.S. citizen, born abroad of American parent(s)` +
  temp$`U.S. citizen, born in Puerto Rico or U.S. Island Areas` +
  temp$`U.S. citizen, born in the United States`
temp$PbornUS <- temp$`U.S. citizen, born in the United States`/tot
temp$PnatUS <- temp$`U.S. citizen by naturalization`/tot
temp$PampUS <- temp$`U.S. citizen, born abroad of American parent(s)`/tot
temp$PislUS <- temp$`U.S. citizen, born in Puerto Rico or U.S. Island Areas`/tot
temp$PnotUS <- temp$`Not a U.S. citizen`/tot

# Rooms
temp$AvrHHsizeByTenure # keep
PperR <- temp$OccupPerRoom/temp$`Aggregate number of rooms`
colnames(temp)[which(colnames(temp)=="Median number of rooms!!Total")] <- "MedianNoRooms"

# Plumbing
Pplumb <- temp$`Complete plumbing facilities`/temp$TotPlumb

# Income
colnames(temp)[which(colnames(temp)=="Median family income in the past 12 months (in 2009 inflation-adjusted dollars)")] <- "MedianFamilyInc"
colnames(temp)[which(colnames(temp)=="Median household income in the past 12 months (in 2009 inflation-adjusted dollars)")] <- "MedianHHInc"
colnames(temp)[which(colnames(temp)=="Per capita income in the past 12 months (in 2009 inflation-adjusted dollars)")] <- "PerCapitaInc"

# Remove excess
temp <- temp[,-c(4, 9:21, 23)]
vars_avr <- temp %>%
  pivot_longer(!c(GEOID, NAME),  names_to = "label", values_to = "estimate")

### Combine density with climate data
## make sure names match
setdiff(clim$county, density_avr$NAME) 
setdiff(climall$county, density_avr$NAME) 
setdiff(density_avr$NAME, climall$county)
clim$county[which(clim$county=="MiamiDade")] <- "Miami-Dade"
clim$county[which(clim$county=="Desoto")] <- "DeSoto"
clim$county[which(clim$county=="Saint Johns")] <- "St. Johns"
clim$county[which(clim$county=="Saint Lucie")] <- "St. Lucie"
climall$county[which(climall$county=="MiamiDade")] <- "Miami-Dade"
climall$county[which(climall$county=="Desoto")] <- "DeSoto"

ind <- match(climall$county, density_avr$NAME)
climall$density <- density_avr$density[ind]

## remove counties with no cases
no_inc <- setdiff(density_avr$NAME, clim$county)

density_avr <- density_avr %>%
  filter(! NAME %in% no_inc)
vars_avr <- vars_avr %>%
  filter(! NAME %in% no_inc)

## match and merge 
ind <- match(clim$county, density_avr$NAME)
clim$density <- density_avr$density[ind]

### Combine demo/epi with climate data
## split alpha!
alpha <- clim[,c(1, 6:8)]
clim <- clim[,c(1:5, 9:10)]

clim <- clim %>%
  mutate(mosq = as.factor(mosq)) %>%
  mutate(mosq = as.numeric(mosq)) %>% # aegyp=1, albo=2, both=3
  pivot_longer(!county,  names_to = "label", values_to = "estimate")
names(clim)[1] <- "NAME"
names(alpha)[1] <- "NAME"

climall <- climall %>%
  mutate(mosq = as.factor(mosq)) %>%
  mutate(mosq = as.numeric(mosq)) %>% # aegyp=1, albo=2, both=3
  pivot_longer(!county,  names_to = "label", values_to = "estimate")
names(climall)[1] <- "NAME"

## combine the files
vars_avr <- vars_avr[,-1] 
head(vars_avr); head(clim); head(climall)
vars_all <- bind_rows(clim, vars_avr)
vars_allcounties <- bind_rows(climall, vars_avr)

names <- unique(vars_all$NAME)
vars_all$meanA <- NA
vars_all$upperA <- NA
vars_all$lowerA <- NA
for(i in 1:length(names)){
  varind <- which(vars_all$NAME==names[i])
  alphind <- which(alpha$NAME==names[i])
  vars_all$meanA[varind] <- rep(alpha$meanA[alphind], length(varind))
  vars_all$upperA[varind] <- rep(alpha$upperA[alphind], length(varind))
  vars_all$lowerA[varind] <- rep(alpha$lowerA[alphind], length(varind))
}

# filter(vars_all, NAME=="Martin")
# filter(vars_all, NAME=="Orange")
# filter(vars_all, label=="density")$meanA
# filter(vars_all, label=="minT")$meanA
head(vars_all); tail(vars_all)
write.csv(vars_all, "output/variables/vars_all_a.csv", row.names=F)
# write.csv(vars_all, "output/variables/vars_all_counties.csv", row.names=F)

### make regression figs
pdf("figs/Vars_vs_a.pdf", width=9, height=7)
ggplot(vars_all, aes(x=estimate, y=meanA)) +
  geom_point() +
  geom_smooth(method=lm) +
  facet_wrap(~label, ncol = 4, scales = "free") +
  xlab(" ") + ylab("Contact (a)")
dev.off()



