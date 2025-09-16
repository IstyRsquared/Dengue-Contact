## Name: Isty Rysava
## Date: 23/4/2024
## Code: Extract epidemiological variables for FL counties - mosquito

rm(list=ls())

## Initialize
library(tidyverse)
library(sf)

## Prep for files
surv <- "Poiss"
ranges <- c("MD", "MT", "MR", "HL", "OR", "ES", "BR", "PB", "DU", "VL", "CO", "AL", "BY", "OS", "CH",
            "BV", "CI", "CL", "CB", "DE", "FL", "FR", "GU", "HE", "HR", "HG", "IR", "LA", "LE", "LO",
            "MA", "MN", "NA", "OK", "PA", "PI", "PO", "PU", "SR", "SA", "SE", "SJ", "SL", "SU", "SW",
            "WA")
names <- c("MiamiDade", "Martin", "Monroe", "Hillsborough", "Orange", "Escambia", "Broward",
           "Palm Beach", "Duval", "Volusia", "Collier", "Alachua", "Bay", "Osceola", "Charlotte",
           "Brevard", "Citrus", "Clay", "Columbia", "Desoto",
           "Flagler", "Franklin", "Gulf", "Hendry", "Hernando", "Highlands", "Indian River", "Lake", "Lee", "Leon",
           "Manatee", "Marion", "Nassau", "Okaloosa", "Pasco", "Pinellas", "Polk", "Putnam", "Santa Rosa", "Sarasota",
           "Seminole", "St. Johns", "St. Lucie", "Sumter", "Suwannee", "Washington")

NSEQ=5000 
NP= 500 

### Get epi + climate 
for (idx in 1:length(ranges)){
  range <- paste0("Biting2_", ranges[idx], "_364")
  
  ## Read in output
  y.dengueMDm <- read_rds(paste0("output/param_search_fl_", surv, "_", NSEQ, "param", NP, "np", range, ".rds"))
  
  ## Get subsets
  subMDm <- y.dengueMDm[order(y.dengueMDm$logLik, decreasing = TRUE),][1:250,] 
  
  ## Calculate summaries
  tibble(subMDm) %>%
    summarize(upperA = quantile(a, 0.975),
              lowerA = quantile(a, 0.025),
              meanA = mean(a)) -> CI_MDm
  
  ## Get R0s, temp, precipitation 
  epi <- matrix(NA, 1, 4)
  
  R=5;M=2 # medium
  R0.temp <- read.csv(paste0("data/test_EstR0perT_", names[idx], "_R", R, "M", M, ".csv"))
  epi[1] <- mean(R0.temp$meanT)
  epi[2] <- mean(R0.temp$minT)
  epi[3] <- mean(R0.temp$meanP)
  epi[4] <- mean(R0.temp$totP)
  
  ## Combine all
  summaries.temp <- cbind(names[idx], epi, CI_MDm)
  if(idx==1){
    summaries <- summaries.temp
  }else{
    summaries <- rbind(summaries, summaries.temp)
  }
}

summaries <- data.frame(summaries)
colnames(summaries)[1:5] <- c("county", "meanT", "minT", "meanP", "totP")
head(summaries); tail(summaries)

## Get mosquito species
mosq <- read.csv("data/mosquito_majority.csv")
ind <- match(mosq$county, summaries$county)
summaries$mosq[ind] <- mosq$species

write.csv(summaries, "output/variables/FL_summaries_a.csv", row.names=F)

################################################################################################################################
### Get variables for all counties
## Prep for files
# fl <- read_sf("data/shp/Florida_adm1.shp", layer="Florida_adm1")

names <- fl$NAME_2
names[which(names=="Miami-Dade")] <- "MiamiDade"
names[which(names=="Saint Johns")] <- "St. Johns"
names[which(names=="Saint Lucie")] <- "St. Lucie"

### Get epi + climate 
for (idx in 1:length(names)){
  R=5;M=2 # medium
  epi <- matrix(NA, 1, 4)
  R0.temp <- read.csv(paste0("data/test_EstR0perT_", names[idx], "_R", R, "M", M, ".csv"))
  epi[1] <- mean(R0.temp$meanT)
  epi[2] <- mean(R0.temp$minT)
  epi[3] <- mean(R0.temp$meanP)
  epi[4] <- mean(R0.temp$totP)
  
  ## Combine all
  summaries.temp <- cbind(names[idx], epi)
  if(idx==1){
    summaries <- summaries.temp
  }else{
    summaries <- rbind(summaries, summaries.temp)
  }
  print(idx)
}

summaries <- data.frame(summaries)
colnames(summaries)[1:5] <- c("county", "meanT", "minT", "meanP", "totP")
head(summaries); tail(summaries)

## Get mosquito species
mosq <- read.csv("data/mosquito_majority.csv")
ind <- match(mosq$county, summaries$county)
summaries$mosq <- "aegypt"
summaries$mosq[ind] <- mosq$species

write.csv(summaries, "output/variables/FL_climate_all.csv", row.names=F)
