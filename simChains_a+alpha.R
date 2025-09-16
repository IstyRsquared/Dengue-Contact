## Name: Isty Rysava
## Date: 5/6/2024
## Code: Simulate transmission chains for FL counties - a + alpha

rm(list=ls())

## Data
pops <- read.csv("data/county_pops.csv")
params <- read.csv("output/variables/FL_summaries_a+alpha.csv")

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
rm=0.2
ec=0.000001
outbreaks <- vector(mode = "list", length = length(names))
names(outbreaks) <- names

### Sim outbreaks 
for (idx in 1:length(names)){
  
  ## Get alpha & pop
  a <- params$meanA[which(params$county==names[idx])]
  alpha <- params$meanAlpha[which(params$county==names[idx])]
  pop <- pops$pop[pops$county==names[idx]]
  
  ## Get R0s 
  R=5;M=2 # medium
  paras <- read.csv(paste0("data/R0perTparas_", names[idx], "_R", R, "M", M, ".csv"))
  # get annul mean
  paras$jbw <- rep(1:26, nrow(paras)/26)
  parsseas <- paras %>%
    mutate(biweek = as.factor(as.character(jbw))) %>%
    group_by(biweek) %>%
    summarize(meanb = mean(b), meanc = mean(c), meanmu = mean(mu), meanPDR = mean(PDR)) %>%
    mutate(biweek = as.numeric(as.character(biweek)))
  parsseas <- parsseas[order(parsseas$biweek,decreasing = FALSE), ]
  
  # prep seasons
  winter <- c(23:26, 1:4)
  spring <- 5:10
  summer <- 11:16
  fall <- 17:22
  wet <- 10:22
  dry <- c(22:26, 1:11)
  
  bseas <- c(mean(parsseas$meanb[winter]), mean(parsseas$meanb[spring]), mean(parsseas$meanb[summer]), 
             mean(parsseas$meanb[fall]), mean(parsseas$meanb[wet]), mean(parsseas$meanb[dry]))
  cseas <- c(mean(parsseas$meanc[winter]), mean(parsseas$meanc[spring]), mean(parsseas$meanc[summer]), 
             mean(parsseas$meanc[fall]), mean(parsseas$meanc[wet]), mean(parsseas$meanc[dry]))
  mumseas <- c(mean(parsseas$meanmu[winter]), mean(parsseas$meanmu[spring]), mean(parsseas$meanmu[summer]), 
               mean(parsseas$meanmu[fall]), mean(parsseas$meanmu[wet]), mean(parsseas$meanmu[dry]))
  PDRseas <- c(mean(parsseas$meanPDR[winter]), mean(parsseas$meanPDR[spring]), mean(parsseas$meanPDR[summer]), 
               mean(parsseas$meanPDR[fall]), mean(parsseas$meanPDR[wet]), mean(parsseas$meanPDR[dry]))
  seasons <- c("Winter", "Spring", "Summer", "Fall", "Wet", "Dry")
  
  ## Sim case chains
  mat <- matrix(NA, ncol=7, nrow=length(seasons))
  for(i in 1:length(bseas)){
    out <- replicate(1000, chain_a_alpha(gens = 7, I = I, a=a, bseas=bseas[i], cseas=cseas[i], 
                                        mumseas=mumseas[i], PDRseas=PDRseas[i], pop=pop, rm=rm, alpha=alpha))
    mat[i,2:3] <- rowMeans(out)
    mat[i,4:5] <- quantile(out[1,], c(0.975, 0.025))
    mat[i,6:7] <- quantile(out[2,], c(0.975, 0.025))
    rm(out)
  }
  
  mat[,1] <- seasons
  colnames(mat) <- c("season", "mean_length", "mean_size", "uQt_length", "lQt_length", "uQt_size", "lQt_size")
  outbreaks[[idx]] <- as_tibble(mat)
  rm(mat); rm(a); rm(alpha)

  print(idx)
}

save(outbreaks, file="./output/county_chains_a+alpha.Rdata")
