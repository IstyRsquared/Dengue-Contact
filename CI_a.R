## Name: Isty Rysava
## Date: 18/04/2024
## Code: Get confidence intervals for estimated a values for county

rm(list=ls())

## Libraries
library(tidyverse)

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
  
  ## Combine
  summaries.temp <- cbind(names[idx], CI_MDm)
  if(idx==1){
    summaries <- summaries.temp
  }else{
    summaries <- rbind(summaries, summaries.temp)
  }
}

summaries <- data.frame(summaries)
colnames(summaries)[1] <- "county"
head(summaries); tail(summaries)
# write.csv(summaries, "output/variables/FL_summaries_a.csv", row.names=F)

pdf("figs/a_CI_estimates.pdf", width=7, height=4)
plot(1:nrow(summaries), summaries$meanA, ylim = c(0, 1), pch = rep(1:17, 4), 
     col="darkred", ylab = "Estimated contact", xlab = " ", axes = FALSE, cex=.9, cex.lab=.8)
axis(side = 1, at = 1:nrow(summaries), labels=summaries$county, las=2, cex.axis = 0.6)
axis(side = 2, cex.axis = 0.6)
for(i in 1:nrow(summaries)){
  segments(i, summaries$lowerA[i], i, summaries$upperA[i], col="darkred")
}
dev.off()

# remove counties with no information
summaries <- filter(summaries, !county %in% c("Bay", "Franklin", "Okaloosa"))

pdf("figs/a_CI_estimates_rm.pdf", width=7, height=4)
plot(1:nrow(summaries), summaries$meanA, ylim = c(0, 0.6), pch = rep(1:17, 4), 
     col="darkred", ylab = "Estimated contact", xlab = " ", axes = FALSE, cex=.9, cex.lab=.8)
axis(side = 1, at = 1:nrow(summaries), labels=summaries$county, las=2, cex.axis = 0.6)
axis(side = 2, cex.axis = 0.6)
for(i in 1:nrow(summaries)){
  segments(i, summaries$lowerA[i], i, summaries$upperA[i], col="darkred")
}
dev.off()