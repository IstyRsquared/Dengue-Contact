## Name: Isty Rysava
## Date: 4/6/2024
## Code: Simulate transmission chains for FL counties - alpha

rm(list=ls())

## Data
params <- read.csv("output/variables/FL_summaries.csv")
pops <- read.csv("data/county_pops.csv")

## Initialize
library(tidyverse)
source("R/simChains_fc.R")

## Prep names
names <- c("MiamiDade", "Martin", "Monroe", "Hillsborough", "Orange", "Escambia", "Broward", 
           "Palm Beach", "Duval", "Volusia", "Collier", "Alachua", "Bay", "Osceola", "Charlotte",
           "Brevard", "Citrus", "Clay", "Columbia", "Desoto",
           "Flagler", "Franklin", "Gulf", "Hendry", "Hernando", "Highlands", "Indian River", "Lake", "Lee", "Leon",
           "Manatee", "Marion", "Nassau", "Okaloosa", "Pasco", "Pinellas", "Polk", "Putnam", "Santa Rosa", "Sarasota",
           "Seminole", "St. Johns", "St. Lucie", "Sumter", "Suwannee", "Washington")
pops$county[which(pops$county=="Miami-Dade")] <- "MiamiDade"

## Initialize
I=1
outbreaks <- vector(mode = "list", length = length(names))
names(outbreaks) <- names

### Sim outbreaks 
set.seed(1666318)
for (idx in 1:length(names)){
  ## Get alpha & pop
  alpha <- params$meanA[which(params$county==names[idx])]
  pop <- pops$pop[pops$county==names[idx]]
  
  ## Get R0s 
  R=5;M=2 # medium
  R0.temp <- read.csv(paste0("data/test_EstR0perT_", names[idx], "_R", R, "M", M, ".csv"))
  # get annul mean
  R0.temp$jbw <- rep(1:26, nrow(R0.temp)/26)
  R0s <- R0.temp %>%
    mutate(biweek = as.factor(as.character(jbw))) %>%
    group_by(biweek) %>%
    summarize(meanR0 = mean(R0)) %>%
    mutate(biweek = as.numeric(as.character(biweek)))
  R0s <- R0s[order(R0s$biweek,decreasing = FALSE), ]
  
  # prep seasons
  winter <- c(23:26, 1:4)
  spring <- 5:10
  summer <- 11:16
  fall <- 17:22
  wet <- 10:22
  dry <- c(22:26, 1:11)
  
  R0seas <- c(mean(R0s$meanR0[winter]), mean(R0s$meanR0[spring]), mean(R0s$meanR0[summer]), 
              mean(R0s$meanR0[fall]), mean(R0s$meanR0[wet]), mean(R0s$meanR0[dry]))
  seasons <- c("Winter", "Spring", "Summer", "Fall", "Wet", "Dry")
  
  ## Sim case chains
  mat <- matrix(NA, ncol=7, nrow=length(seasons))
  for(i in 1:length(R0seas)){
    out <- replicate(1000, chain_alpha(gens = 7, I = 1, R0seas=R0seas[i], alpha=alpha, pop=pop))
    mat[i,2:3] <- rowMeans(out)
    mat[i,4:5] <- quantile(out[1,], c(0.975, 0.025))
    mat[i,6:7] <- quantile(out[2,], c(0.975, 0.025))
    rm(out)
  }
  
  mat[,1] <- seasons
  colnames(mat) <- c("season", "mean_length", "mean_size", "uQt_length", "lQt_length", "uQt_size", "lQt_size")
  outbreaks[[idx]] <- as_tibble(mat)
  rm(mat)
  
  print(idx)
}
save(outbreaks, file="./output/county_chains_alpha.Rdata")




