## Name: Isty Rysava
## Date: 14/3/2024
## Code: Simulate outbreaks for FL counties

rm(list=ls())

## Data
params <- read.csv("output/variables/FL_summaries.csv")
pops <- read.csv("data/county_pops.csv")

## Initialize
library(tidyverse)

## Prep names
names <- c("MiamiDade", "Martin", "Monroe", "Hillsborough", "Orange", "Escambia", "Broward", 
           "Palm Beach", "Duval", "Volusia", "Collier", "Alachua", "Bay", "Osceola", "Charlotte",
           "Brevard", "Citrus", "Clay", "Columbia", "Desoto",
           "Flagler", "Franklin", "Gulf", "Hendry", "Hernando", "Highlands", "Indian River", "Lake", "Lee", "Leon",
           "Manatee", "Marion", "Nassau", "Okaloosa", "Pasco", "Pinellas", "Polk", "Putnam", "Santa Rosa", "Sarasota",
           "Seminole", "St. Johns", "St. Lucie", "Sumter", "Suwannee", "Washington")
pops$county[which(pops$county=="Miami-Dade")] <- "MiamiDade"

## Initialize
I=2
outbreaks <- vector(mode = "list", length = length(names))
names(outbreaks) <- names

### Sim outbreaks 
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
  
  ## Sim secondary cases 
  mat <- matrix(NA, nrow=1000, ncol=length(seasons))
  pdf(paste0("figs/outbreaks_alpha/", names[idx], ".pdf"), height=8, width=6)
  par(mfrow=c(3,2))
  for(i in 1:length(R0seas)){
    lambda <- R0seas[i] * (I^alpha)/pop # calculate FOI 
    mat[,i] <- rbinom(1000, pop, lambda) 
    barplot(table(mat[,i])/1000, main=paste0(seasons[i], " R0 = ", round(R0seas[i], digits=2), " alpha = ", round(alpha, digits=2)), 
            cex.lab=1, cex.main=1, cex.axis=0.8, cex.names=0.8, xlab="No. of secondary cases", ylab="Frequency", ylim=c(0,1),
            col=c("lightblue", "gold", "#008080", "orange", "darkblue", "darkred")[i], border=NA)
  }
  dev.off()
  colnames(mat) <- seasons
  outbreaks[[idx]] <- mat
  
  print(idx)
}

save(outbreaks, file="./output/county_outbreaks_alpha.Rdata")

### Get mean and 95% interval for each county
matWinter <- matSpring <- matSummer <- matFall <- matWet <- matDry <- data.frame(matrix(NA, nrow=length(names), ncol=4))
colnames(matWinter) <- colnames(matSpring) <- colnames(matSummer) <- colnames(matFall) <- colnames(matWet) <- colnames(matDry) <- c("mean", "highPI", "lowPI", "county")

for (idx in 1:length(names)){
  temp <- outbreaks[[idx]]
  matWinter[idx, 1] <- mean(temp[,1]); matWinter[idx, 2] <- quantile(temp[,1], 0.975); matWinter[idx, 3] <-  quantile(temp[,1], 0.025)
  matSpring[idx, 1] <- mean(temp[,2]); matSpring[idx, 2] <- quantile(temp[,2], 0.975); matSpring[idx, 3] <-  quantile(temp[,2], 0.025)
  matSummer[idx, 1] <- mean(temp[,3]); matSummer[idx, 2] <- quantile(temp[,3], 0.975); matSummer[idx, 3] <-  quantile(temp[,3], 0.025)
  matFall[idx, 1] <- mean(temp[,4]); matFall[idx, 2] <- quantile(temp[,4], 0.975); matFall[idx, 3] <-  quantile(temp[,4], 0.025)
  matWet[idx, 1] <- mean(temp[,5]); matWet[idx, 2] <- quantile(temp[,5], 0.975); matWet[idx, 3] <-  quantile(temp[,5], 0.025)
  matDry[idx, 1] <- mean(temp[,6]); matDry[idx, 2] <- quantile(temp[,6], 0.975); matDry[idx, 3] <-  quantile(temp[,6], 0.025)
  matWinter[idx,4] <- matSpring[idx,4] <- matSummer[idx,4] <- matFall[idx,4] <- matWet[idx,4] <- matDry[idx,4] <- names[idx]
}

## Seasons
matDry$season <- "Dry"
matWet$season <- "Wet"
matSeasons <- bind_rows(matDry, matWet)
pdf("figs/outbreaks_alpha/SecCases_season.pdf", height=5, width=9)
ggplot(matSeasons, aes(x=county, y=mean)) + 
  geom_pointrange(aes(ymin=lowPI, ymax=highPI), size=0.3) +
  labs( x ="County", y = "No. of secondary cases") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=6)) +
  facet_wrap(~season, ncol = 2)
dev.off()
  
## Year
matWinter$season <- "Winter"
matSpring$season <- "Spring"
matSummer$season <- "Summer"
matFall$season <- "Fall"
matYear <- bind_rows(matWinter, matSpring, matSummer, matFall)
pdf("figs/outbreaks_alpha/SecCases_year.pdf", height=8, width=10)
ggplot(matYear, aes(x=county, y=mean)) + 
  geom_pointrange(aes(ymin=lowPI, ymax=highPI), size=0.1) +
  labs( x ="County", y = "No. of secondary cases") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=6),
        axis.title.x = element_text(size=9), axis.title.y = element_text(size=9)) +
  facet_wrap(~season, ncol = 2)
dev.off()



