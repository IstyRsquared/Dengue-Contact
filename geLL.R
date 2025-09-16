## Name: Isty Rysava
## Date: 10/7/2024
## Code: Get LL values for Scenarios 1-4 across FL counties

rm(list=ls())

## Libraries
library(tidyverse)

## Initialize
NSEQ=5000 
NP= 500 
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

LLs <- matrix(NA, nrow=length(ranges), ncol=4)

## Loop and read
for(idx in 1:length(ranges)){
  range <- ranges[idx]
  y.dengue1 <- read_rds(paste0("output/param_search_fl_Poiss_", NSEQ, "param", NP, "npBiting2_", range, "_364.rds")) # a
  y.dengue2 <- read_rds(paste0("output/param_search_fl_MediumR0_", NSEQ, "param", NP, "npMalaria_", range, "_364.rds")) # alpha
  y.dengue3 <- read_rds(paste0("output/param_search_fl_Poiss_", NSEQ, "param", NP, "npBitingAlpha_", range, "_364.rds")) # a+alpha
  y.dengue4 <- read_rds(paste0("output/param_search_fl_NB_", NSEQ, "param", NP, "npBiting_", range, "_364.rds")) # a+k
  LLs[idx, 1] <- max(y.dengue1$logLik)
  LLs[idx, 2] <- max(y.dengue2$logLik)
  LLs[idx, 3] <- max(y.dengue3$logLik)
  LLs[idx, 4] <- max(y.dengue4$logLik)
}

LLs <- cbind(data.frame(names), data.frame(LLs))
colnames(LLs) <- c("County", "Sc1", "Sc2", "Sc3", "Sc4")
head(LLs); tail(LLs)

write.csv(LLs, "output/FL_county_LLs.csv", row.names=F)

plot(1:46, LLs$Sc1, type="l", col="red")
lines(1:46, LLs$Sc2, col="orange")
lines(1:46, LLs$Sc2, col="blue")
lines(1:46, LLs$Sc2, col="green")


