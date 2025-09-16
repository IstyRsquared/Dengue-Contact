## Name: Isty Rysava
## Date: 09/10/2024
## Code: Fit explanatory variables to gamma regression analysis: interpretative

rm(list=ls())

### Libraries
library(tidyverse)
library(stringr)
library(SpatialEpi)
library(sf)
library(boot)
library(harrypotter)
library(StepReg)
library(car)

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
vars_all <- read.csv("output/variables/vars_all_a.csv")
travel_all <- read.csv("output/FL_biweekly_CountyTravel.csv")

## Remove Bay, Franklin, Okaloosa 
vars_rm <- filter(vars_all, !NAME %in% c("Bay", "Franklin", "Okaloosa"))
travel <- tibble(travel_all) %>% 
  group_by(county) %>% summarize(cases = sum(cases)) %>%
  filter(!county %in% c("Bay", "Franklin", "Okaloosa"))

### Run Gamma regressions
datreg <- vars_rm %>% pivot_wider(names_from = label, values_from = estimate) # make wide
travel$county[which(travel$county=="Desoto")] <- "DeSoto" # add intros
ind <- match(travel$county, datreg$NAME)
datreg$intro <- NA
datreg$intro[ind] <- travel$cases # add incursions
datreg$mosq <- as.numeric(datreg$mosq) # prep mosquitos

### FULL MODEL & WEIGHTED by incursions 
fullreg <- datreg
wghtmod <- glm(meanA ~ meanT+minT+meanP+totP+density+MedianFamilyInc+MedianHHInc+MedianNoRooms+PerCapitaInc+TotPop+Poccup+ 
                 Powned+PbornUS+PnatUS+PampUS+PnotUS+as.factor(mosq), data = datreg,family = Gamma(link = "logit"), 
               weights=datreg$intro, maxit=1000) 
summary(wghtmod) 
vif(wghtmod) # A lot co-linear!: minT, totP, PbornUS, PnatUS, (TotPop)
wghtmod2 <- glm(meanA ~ meanT+meanP+density+MedianFamilyInc+MedianHHInc+MedianNoRooms+PerCapitaInc+Poccup+ 
                 Powned+PampUS+PnotUS+as.factor(mosq), data = datreg,family = Gamma(link = "log"), 
               weights=datreg$intro, maxit=1000) 
summary(wghtmod2) 
vif(wghtmod2) 
modAICW <- MASS::stepAIC(wghtmod2, k = 2, steps=1000)
summary(modAICW)

fit <- data.frame(exp_coeff=round(exp(coef(modAICW)), digits=3), coeff=round(coef(modAICW), digits=3), 
                  effects=round(effects(modAICW)[1:13], digits=3), pval=round(coef(summary(modAICW))[,4], digits=3))
# A unit increase in the varx causes a exp(coef(modAICW))-1 relative change in the contact
(round(exp(coef(modAICW)), digits=3)-1)*100
effect=data.frame(a=c(round(exp(coef(modAICW)[1]), digits=3), 
               round(exp(coef(modAICW)[1]) * exp(coef(modAICW)[2:13]), digits=3)))
fit <- data.frame(cbind(fit, effect))
head(fit)
# write.csv(fit, "output/GammaModFitLog_effect.csv")

###########################################################################################################################
### Weighted regression output map: predict for missing counties
reducedmod=modAICW

vars_all <- read.csv("output/variables/vars_all_counties.csv")
vars_rm <- vars_all
travel <- tibble(travel_all) %>% group_by(county) %>% summarize(cases = sum(cases)) 

### Run Gamma regressions
datreg <- vars_rm %>% pivot_wider(names_from = label, values_from = estimate) # make wide
travel$county[which(travel$county=="Desoto")] <- "DeSoto" # add intros
ind <- match(travel$county, datreg$NAME)
datreg$intro <- 0
datreg$intro[ind] <- travel$cases
datreg$mosq <- as.factor(datreg$mosq)

preds <- cbind(datreg[,1], fitted=predict(reducedmod, datreg))
summary(exp(preds$fitted))
summary(inv.logit(preds$fitted))
breaks <- seq(0, 0.5, by=.05)
pal <- hp(n = length(breaks), house = "NewtScamander")

# match with the shapefile
setdiff(preds$NAME, fl$NAME_2) 
preds$NAME[which(preds$NAME=="St. Johns")] <- "Saint Johns"
preds$NAME[which(preds$NAME=="St. Lucie")] <- "Saint Lucie"
preds$NAME[which(preds$NAME=="DeSoto")] <- "Desoto"
fl$estA <- NA
meanA.inds=match(preds$NAME, fl$NAME_2)
fl$estA[meanA.inds] <- exp(preds$fitted)

pdf("figs/MeanEstimatedA_WGHTregressions_AllCountiesVIFlog.pdf", height=5, width=6)
plot(st_geometry(fl), col="white",border=NA)
plot(st_geometry(fl), add=T, border=NA, col=pal[findInterval(as.vector(as.character(fl$estA)),
                                                             breaks,all.inside=T)])
#plot(st_geometry(fl[which(fl$NAME_2 %in% c("Bay", "Franklin", "Okaloosa")),]), add=T, col="black")
plot(st_geometry(fl), add=T, border="grey70", lwd=.1)
legend("bottomleft", 
       legend=c(leglabs(breaks, under="<",over=">")),
       pch=c(rep(22,(length(breaks)+1))),
       lty=c(rep(NA, 12)),
       col = c(rep("black",(length(breaks)+1))),
       pt.bg=c(pal[-length(pal)]), pt.cex=2, lwd=.1,
       title=c("Mean heterogeneous contact (a): \n Gamma regression estimates"),
       cex=.8, bty="n")
dev.off()

