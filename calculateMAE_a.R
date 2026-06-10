## Name: Isty Rysava
## Date: 16/04/2026
## Code: Calculate MAE & draw ts figs for all counties
## Miami-Dade medium R0: V->H and H->V transmission with species-specific lambda

rm(list=ls())
setwd("/Users/u1674940/Documents/Arboviruses_CDC/Florida-Transmission-new")

## Libraries
library(Metrics)
library(tidyverse)

## Loop through all counties and combine
surv <- "Poiss"
range <- "Biting2_MD_364"

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
range <- paste0("Biting2_", ranges, "_364")

## prep data base
ts_dat <- c()
mae_dat <- matrix(NA, nrow=length(ranges), ncol=2)
mae_dat[,1] <- ranges

for(idx in 1:length(ranges)){
  input <- read_rds(paste0("output/sims/allcases", range[idx], ".rds"))
  data <- input[[1]]
  sims <- input[[2]]
  
  ts_dat <- rbind(ts_dat, data.frame(data$cases, sims$cases, sims$upper, sims$lower))
  mae_dat[idx,2] <- mae(data$cases, sims$cases) # calculate MAE
}

## Save MAE output
mae_data <- data.frame(mae_dat)
colnames(mae_data) <- c("range", "MAE")
mae_data$county <- names
mae_data$roundMAE <- round(as.numeric(mae_data$MAE), digits=4)
head(mae_data); tail(mae_data)
# write.csv(mae_data, "output/MAE_county.csv", row.names=F)

ts_dat <- data.frame(ts_dat)
ts_dat$county <- rep(names, each=168)
colnames(ts_dat) <- c("obs_cases","sim_cases", "upper", "lower", "county")
head(ts_dat); tail(ts_dat)

## Draw timeseries figures
ts_dat$date <- rep(seq(as.Date("2009/1/1"), as.Date("2022/12/31"), "months"), length(ranges))

pdf(paste0("figs/biting/FL_", surv, "_param_AllCounties_5PERC.pdf"), width=9, height=8)
ggplot(ts_dat, aes(x = as.Date(date), y = sim_cases)) +
  geom_line(color = "#008080") +
  geom_ribbon(aes(ymax = upper, ymin = lower), alpha = 0.5, fill = "#008080") +
  geom_line(data=ts_dat, aes(x=as.Date(date), y=obs_cases), color = "deeppink2") +
  xlab("Time (monthly)") + ylab("Local dengue cases") +
  scale_x_date(date_breaks = "12 month", date_labels =  "%Y") +
  facet_wrap(~county, scales = "free_y", ncol=6) +
  theme_minimal(base_size = 6) +
  theme(axis.text.x = element_text(angle = 60, vjust = 0.5, hjust=1)) 
dev.off()

# plot.title = element_text(size=5,  face="bold"),
# axis.text = element_text(size = 1), axis.title = element_text(size = 9),
