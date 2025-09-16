## Name: Isty Rysava
## Date: 10/6/2024
## Code: Summarize and plot transmission chains for FL counties 

rm(list=ls())

### Libraries
library(tidyverse)

### Data
load(file="./output/county_chains_alpha.Rdata")
chains_alpha <- outbreaks
load(file="./output/county_chains_a.Rdata")
chains_a <- outbreaks
load(file="./output/county_chains_a+alpha.Rdata")
chains_a_alpha <- outbreaks
load(file="./output/county_chains_a+k.Rdata")
chains_a_k <- outbreaks
rm(outbreaks)

all <- list(chains_alpha, chains_a, chains_a_alpha, chains_a_k)

### Initialize
names <- c("MiamiDade", "Martin", "Monroe", "Hillsborough", "Orange", "Escambia", "Broward",
           "Palm Beach", "Duval", "Volusia", "Collier", "Alachua", "Bay", "Osceola", "Charlotte",
           "Brevard", "Citrus", "Clay", "Columbia", "Desoto",
           "Flagler", "Franklin", "Gulf", "Hendry", "Hernando", "Highlands", "Indian River", "Lake", "Lee", "Leon",
           "Manatee", "Marion", "Nassau", "Okaloosa", "Pasco", "Pinellas", "Polk", "Putnam", "Santa Rosa", "Sarasota",
           "Seminole", "St. Johns", "St. Lucie", "Sumter", "Suwannee", "Washington")

type <- c("alpha", "a", "a+alpha", "a+k")

### Get outbreak length by county
for(j in 1:length(type)){
  chains <- all[[j]]
  
  matWinter <- matSpring <- matSummer <- matFall <- matWet <- matDry <- data.frame(matrix(NA, nrow=length(names), ncol=4))
  colnames(matWinter) <- colnames(matSpring) <- colnames(matSummer) <- colnames(matFall) <- colnames(matWet) <- colnames(matDry) <- c("mean", "highPI", "lowPI", "county")
  
  for (idx in 1:length(names)){
    temp <- chains[[idx]]
    temp[,2:7] <- temp[,2:7] %>% 
      mutate_if(is.character,as.numeric)
    matWinter[idx, 1] <- round(temp[1,2], digits=2); matWinter[idx, 2] <- round(temp[1,4], digits=2); matWinter[idx, 3] <-  round(temp[1,5], digits=2)
    matSpring[idx, 1] <- round(temp[2,2], digits=2); matSpring[idx, 2] <- round(temp[2,4], digits=2); matSpring[idx, 3] <-  round(temp[2,5], digits=2)
    matSummer[idx, 1] <- round(temp[3,2], digits=2); matSummer[idx, 2] <- round(temp[3,4], digits=2); matSummer[idx, 3] <-  round(temp[3,5], digits=2)
    matFall[idx, 1] <- round(temp[4,2], digits=2); matFall[idx, 2] <- round(temp[4,4], digit=2); matFall[idx, 3] <-  round(temp[4,5], digits=2)
    matWet[idx, 1] <- round(temp[5,2], digits=2); matWet[idx, 2] <- round(temp[5,4], digits=2); matWet[idx, 3] <-  round(temp[5,5], digits=2)
    matDry[idx, 1] <- round(temp[6,2], digits=2); matDry[idx, 2] <- round(temp[6,4], digits=2); matDry[idx, 3] <-  round(temp[6,5], digits=2)
    matWinter[idx,4] <- matSpring[idx,4] <- matSummer[idx,4] <- matFall[idx,4] <- matWet[idx,4] <- matDry[idx,4] <- names[idx]
  }
  
  ## Seasons
  matDry$season <- "Dry"
  matWet$season <- "Wet"
  matSeasons <- bind_rows(matDry, matWet)
  matSeasons <- filter(matSeasons, !county %in% c("Bay", "Franklin", "Okaloosa"))
  p <- ggplot(matSeasons, aes(x=county, y=mean)) + 
    geom_pointrange(aes(ymin=lowPI, ymax=highPI), size=0.3) +
    labs( x ="County", y = "No. of outbreak generations") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=6)) +
    facet_wrap(~season, ncol = 2)
  pdf(paste0("figs/chains/OutbreakLength_season_", type[j], ".pdf"), height=5, width=9)
  print(p)
  dev.off()
  
  ## Year
  matWinter$season <- "Winter"
  matSpring$season <- "Spring"
  matSummer$season <- "Summer"
  matFall$season <- "Fall"
  matYear <- bind_rows(matWinter, matSpring, matSummer, matFall)
  matYear <- filter(matYear, !county %in% c("Bay", "Franklin", "Okaloosa"))
  p <- ggplot(matYear, aes(x=county, y=mean)) + 
    geom_pointrange(aes(ymin=lowPI, ymax=highPI), size=0.1, color = "#008080") +
    labs( x ="County", y = "No. of outbreak generations") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=12),
          axis.text.y = element_text(size=12),
          axis.title.x = element_text(size=19), axis.title.y = element_text(size=19),
          strip.text = element_text(size=20)) +
    facet_wrap(~season, ncol = 2)
  pdf(paste0("figs/chains/OutbreakLength_year_", type[j], ".pdf"), height=8, width=12)
  print(p)
  dev.off()
  
}

### Get outbreak size by county
for(j in 1:length(type)){
  chains <- all[[j]]
  
  matWinter <- matSpring <- matSummer <- matFall <- matWet <- matDry <- data.frame(matrix(NA, nrow=length(names), ncol=4))
  colnames(matWinter) <- colnames(matSpring) <- colnames(matSummer) <- colnames(matFall) <- colnames(matWet) <- colnames(matDry) <- c("mean", "highPI", "lowPI", "county")
  
  for (idx in 1:length(names)){
    temp <- chains[[idx]]
    temp[,2:7] <- temp[,2:7] %>% 
      mutate_if(is.character,as.numeric)
    matWinter[idx, 1] <- round(temp[1,3], digits=2); matWinter[idx, 2] <- round(temp[1,6], digits=2); matWinter[idx, 3] <-  round(temp[1,7], digits=2)
    matSpring[idx, 1] <- round(temp[2,3], digits=2); matSpring[idx, 2] <- round(temp[2,6], digits=2); matSpring[idx, 3] <-  round(temp[2,7], digits=2)
    matSummer[idx, 1] <- round(temp[3,3], digits=2); matSummer[idx, 2] <- round(temp[3,6], digits=2); matSummer[idx, 3] <-  round(temp[3,7], digits=2)
    matFall[idx, 1] <- round(temp[4,3], digits=2); matFall[idx, 2] <- round(temp[4,6], digit=2); matFall[idx, 3] <-  round(temp[4,7], digits=2)
    matWet[idx, 1] <- round(temp[5,3], digits=2); matWet[idx, 2] <- round(temp[5,6], digits=2); matWet[idx, 3] <-  round(temp[5,7], digits=2)
    matDry[idx, 1] <- round(temp[6,3], digits=2); matDry[idx, 2] <- round(temp[6,6], digits=2); matDry[idx, 3] <-  round(temp[6,7], digits=2)
    matWinter[idx,4] <- matSpring[idx,4] <- matSummer[idx,4] <- matFall[idx,4] <- matWet[idx,4] <- matDry[idx,4] <- names[idx]
  }
  
  ## Seasons
  matDry$season <- "Dry"
  matWet$season <- "Wet"
  matSeasons <- bind_rows(matDry, matWet)
  matSeasons <- filter(matSeasons, !county %in% c("Bay", "Franklin", "Okaloosa"))
  p <- ggplot(matSeasons, aes(x=county, y=mean)) + 
    geom_pointrange(aes(ymin=lowPI, ymax=highPI), size=0.3) +
    labs( x ="County", y = "Outbreak size") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=6)) +
    facet_wrap(~season, ncol = 2)
  pdf(paste0("figs/chains/OutbreakSize_season_", type[j], ".pdf"), height=5, width=9)
  print(p)
  dev.off()
  
  ## Year
  matWinter$season <- "Winter"
  matSpring$season <- "Spring"
  matSummer$season <- "Summer"
  matFall$season <- "Fall"
  matYear <- bind_rows(matWinter, matSpring, matSummer, matFall)
  matYear <- filter(matYear, !county %in% c("Bay", "Franklin", "Okaloosa"))
  p <- ggplot(matYear, aes(x=county, y=mean)) + 
    geom_pointrange(aes(ymin=lowPI, ymax=highPI), size=0.1, color = "#008080") +
    labs( x ="County", y = "Outbreak size") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=12),
          axis.text.y = element_text(size=12),
          axis.title.x = element_text(size=19), axis.title.y = element_text(size=19),
          strip.text = element_text(size=20)) +
    facet_wrap(~season, ncol = 2)
  pdf(paste0("figs/chains/OutbreakSize_year_", type[j], ".pdf"), height=8, width=12)
  print(p)
  dev.off()
  
}
