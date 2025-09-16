## Name: Isty Rysava
## Date: 22/04//2024
## Code: Simulate outbreaks for FL counties - with a instead of alpha

rm(list=ls())

## Data
pops <- read.csv("data/county_pops.csv")
params <- read.csv("output/variables/FL_summaries_a.csv")

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
I=1
rm=0.2
ec=0.000001
outbreaks <- vector(mode = "list", length = length(names))
names(outbreaks) <- names

### Sim outbreaks 
for (idx in 1:length(names)){
  
  ## Get alpha & pop
  a <- params$meanA[which(params$county==names[idx])]
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
  
  ## Sim secondary cases 
  mat <- matrix(NA, nrow=1000, ncol=length(seasons))
  pdf(paste0("figs/outbreaks_a/", names[idx], ".pdf"), height=8, width=6)
  par(mfrow=c(3,2))
  for(i in 1:length(bseas)){
    lambdaM = (a * bseas[i] * exp(-mumseas[i]/(PDRseas[i]+ec)))/(rm) * I
    M = rpois(1000, lambdaM)
    lambdaH = a * cseas[i] * 1/(mumseas[i]) * M * ((pop-1)/pop)
    mat[,i] = rpois(1000, lambdaH)
    barplot(table(mat[,i])/1000, main=paste0(seasons[i], " a = ", round(a, digits=2)),
            cex.lab=1, cex.main=1, cex.axis=0.8, cex.names=0.8, xlab="No. of secondary cases", ylab="Frequency",
            col=c("lightblue", "gold", "#008080", "orange", "darkblue", "darkred")[i], border=NA, ylim=c(0,1))
  }
  dev.off()
  
  colnames(mat) <- seasons
  outbreaks[[idx]] <- mat
  rm(a)
  print(idx)
}

save(outbreaks, file="./output/county_outbreaks_a.Rdata")

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
matSeasons <- filter(matSeasons, !county %in% c("Bay", "Franklin", "Okaloosa"))
pdf("figs/outbreaks_a/SecCases_season_a.pdf", height=5, width=9)
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
matYear <- filter(matYear, !county %in% c("Bay", "Franklin", "Okaloosa"))
pdf("figs/outbreaks_a/SecCases_year_a.pdf", height=8, width=10)
ggplot(matYear, aes(x=county, y=mean)) + 
  geom_pointrange(aes(ymin=lowPI, ymax=highPI), size=0.1) +
  labs( x ="County", y = "No. of secondary cases") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=6),
        axis.title.x = element_text(size=9), axis.title.y = element_text(size=9)) +
  facet_wrap(~season, ncol = 2)
dev.off()

### Percanatege of secondary cases from single intro MD
table(md$Winter)/10 # W: 98% 0, 1% 1, 0.3% 2, 0.1 3 
table(md$Spring)/10 # Sp: 93.7% 0, 4.2% 1, 1.6% 2, 0.4% 3, 0.1% 5 
table(md$Summer)/10 # Sm: 82.7% 0, 8.2% 1, 5.9% 2, 2.3% 3, 0.7% 4, 0.2% 5 
table(md$Fall)/10 # F: 84.6% 0, 7.8% 1, 4.7% 2, 2.3% 3, 0.3% 4, 0.2% 5, 0.1% 9

mean(c(98, 93.7, 82.7, 84.6)); sd(c(98, 93.7, 82.7, 84.6)) # 0 cases
mean(c(1, 4.2, 8.2, 7.8)); sd(c(1, 4.2, 8.2, 7.8)) # 1 case
mean(c(0.3, 1.6, 5.9, 4.7)); sd(c(0.3, 1.6, 5.9, 4.7)) # 2 cases
mean(c(0.1, 0.4, 2.3, 2.3)); sd(c(0.1, 0.4, 2.3, 2.3)) # 3 cases
mean(c(0, 0, 0.7, 0.3)); sd(c(0, 0, 0.7, 0.3)) # 4 cases
