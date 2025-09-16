## Name: Isty Rysava
## Date: 12/1/2024
## Code: Estimate dengue transmission parameters using POMP & binomial surveillance draw
## Miami-Dade medium R0: alpha saturation (Scenario 2)

## Dengue POMP model
rm(list=ls())

## Prep for Windows
params <-
  list(prefix = "windows")
set.seed(1350254336)
options(pomp_cdir="./tmp")
surv <- "MediumR0"
range <- "Malaria_MD_364"

## Libraries
library(tidyverse)
library(pomp)
library(foreach) # loop for do parallel
library(doParallel) # backed for foreach
library(iterators)
library(doFuture)
library(pbapply)

## Data
if(surv=="MediumR0"){R=5;M=2}
travel_all <- read.csv("data/FL_biweekly_CountyTravel.csv")
local_all <- read.csv("data/FL_biweekly_CountyLocal.csv")
R0s <- read.csv(paste0("data/test_EstR0perT_MiamiDade_R", R, "M", M, ".csv"))
R0 <- R0s$R0

## Format monthly
travel <- tibble(travel_all) %>% filter(county=="Miami-Dade") %>% group_by(biweek) %>% summarize(cases = sum(cases))  
local <- tibble(local_all) %>% filter(county=="Miami-Dade") %>% group_by(biweek) %>% summarize(cases = sum(cases))

weeks <- seq(as.Date("2009/1/1"), as.Date("2022/12/31"), "weeks")
biweeks <- weeks[seq(1,length(weeks),2)][2:(nrow(travel)+1)]
months <- ((as.POSIXlt(strptime(biweeks, format="%Y-%m-%d"))$year-109)*12) + as.POSIXlt(strptime(biweeks, format="%Y-%m-%d"))$mon+1

R0 <- tibble(R0) %>% group_by(floor(seq(1, 335, length.out=dim(travel)[1]))) %>% summarize(R0 = mean(R0))
covar <- data.frame(times=seq(1, 168, by=.5), introsR=rowsum(travel, group=floor(seq(1, 335, length.out=dim(travel)[1])))$cases, R0=R0$R0)
data <- data.frame(times=unique(months), cases=rowsum(local, group=months)$cases)

## Process model
rproc <- Csnippet("
                  double intros;
                  double tol = 1.0e-25;

                  // deaths take precedent over all other rates
                  double deathsS = rbinom(S, 1-exp(-mu*dt));
                  double deathsI = rbinom(I, 1-exp(-mu*dt));
                  double deathsR = rbinom(R, 1-exp(-mu*dt));
                  double deathsT = deathsS + deathsI + deathsR;
                  
                  // births for all that are not dead
                  double births = rbinom(S + I + R - deathsT, 1-exp(-br*dt));

                  // incursions at each observation time
                  intros = rnbinom((introsR+tol), rhoT);


                  // transition from S to I  
                  // need to keep track of local and intros inf!!
                  double transmitters = I + intros;
                  double seasbet = R0;
                  double pop = S + transmitters + R;
                  double lambda = seasbet*R_pow(transmitters, alpha)/pop;
                  double infected = rbinom(S - deathsS, 1-exp(-lambda));
                  
                  // Rprintf(\"  %f \\t  %f \\t %f \\t %f \\n  \", t, -sf, seasbet, lambda); 

                  // transition from R to S: waning immunity
                  double waning = rbinom(R - deathsR, 1-exp(-gamma*dt));

                  // balance equations 
                  S += births - infected - deathsS + waning;
                  
                  R += infected - waning - deathsR; 

                  // infectious period bi-weekly (i.e. timestep of process model) so not additive
                  I = infected;
                  
                  C += infected;           // true infecteds
                  
                  N = S + I + R;
                  
                  FOI = lambda;

                  //Rprintf(\"  %f \\t  %f \\t %f \\t %f \\t %f \\n  \", t, S, R, I, C);                  

                  ")

## Initial conditions
initlz <- Csnippet("
                   N = 2500000+10;
                   S = 2500000;
                   I = 10;
                   R = 0;
                   C = 10;
                   FOI = 0;
                   ")

### Observation model density
dmeas <- Csnippet("
                  double tol = 1.0e-25;
                  lik = dbinom (cases, C, rhoL, 0) + tol ;
                  if (give_log) lik = log(lik);
                  ")

### Observation model process
rmeas <- Csnippet("
                  cases = rbinom (C, rhoL);
                  ")

### Make POMP object
dengue.MDmedium <- data %>%
  pomp(t0 = 1, 
       times = data$times, 
       rprocess = discrete_time(step.fun = Csnippet(rproc), delta.t = 1/2), # 1/2 (so dt = bi-weekly)
       dmeasure = dmeas,
       rmeasure = rmeas,
       rinit = initlz,
       partrans = parameter_trans(logit=c("alpha", "rhoL", "rhoT")),
       covar = covariate_table(covar, times="times"),
       accumvars = "C",       
       obsnames=c("cases"), 
       statenames = c("N", "S", "I", "R", "C", "FOI"),
       paramnames = c("alpha", "rhoL", "rhoT",  "gamma", "mu", "br"),
       cdir=".", cfile="dengue.fl"
  ) 

NSEQ=5000 
NP= 500 

### Prep guesses
guesses <- sobol_design(
  lower=c(alpha=0.1, 
          rhoL=0.01, 
          rhoT=0.01),
  upper=c(alpha=1, 
          rhoL=0.2, 
          rhoT=0.2),
  nseq=NSEQ) |> 
  mutate(gamma = 30.5/(6*30.5), mu = 1027/(12*100000), br = 11/(12*1000))

### Search function
search <- function(dengue, i){
  dengue %>%
    pfilter(
      params=guesses[i,],
      Np=NP,
      dmeasure=dmeas,
      statenames = c("N", "S", "I", "C", "R"),
      paramnames = c("alpha", "rhoL", "rhoT", "gamma", "mu", "br"), 
      partrans = parameter_trans(logit=c("alpha", "rhoL", "rhoT")
      )
    ) %>%
    logLik() %>%
    logmeanexp() %>%
    {
      tibble(
        "alpha" = guesses[i, "alpha"],
        "rhoL" = guesses[i, "rhoL"],
        "rhoT" = guesses[i, "rhoT"],
        "gamma" = guesses[i, "gamma"],
        "mu" = guesses[i, "mu"],
        "br" = guesses[i, "br"],
        "logLik" = .
      )
    }
}

### Estimate
y.dengue <- pblapply(1:nrow(guesses), function(x) search(dengue = dengue.MDmedium, i = x)) %>% bind_rows()
write_rds(y.dengue, paste0("output/param_search_fl_", surv, "_", NSEQ, "param", NP, "np", range, ".rds"))

### Check and plot (ts and pairs boxes)
y.dengue <- read_rds(paste0("output/param_search_fl_", surv, "_", NSEQ, "param", NP, "np", range, ".rds"))
summary(y.dengue$logLik) 

test <- y.dengue |>
  filter(logLik > -500)
test2 <- y.dengue[order(y.dengue$logLik, decreasing = TRUE),][1:250,] 

y.out2 <- tibble(data.frame(alpha=mean(test2$alpha), rhoL=mean(test2$rhoL), rhoT=mean(test2$rhoT),
                            gamma=mean(test2$gamma), mu=mean(test2$mu), br=mean(test2$br)))
y.dengue |>
  filter(logLik == max(logLik)) -> y.out

### Plot output: y.out2
dengue.MDmedium %>%
  simulate(params=y.out2,
           nsim = 1000,
           format = "data.frame"
  ) -> dengue_sim

dengue_sim %>%
  group_by(time) %>%
  summarize(upper = quantile(cases, 0.975),
            lower = quantile(cases, 0.025),
            cases = mean(cases)) -> sims_summed
sims_summed$date <- seq(as.Date("2009/1/1"), as.Date("2022/12/31"), "months")
data$date <- seq(as.Date("2009/1/1"), as.Date("2022/12/31"), "months")
library(scales)

pdf(paste0("figs/range/FL_", surv, "_", NSEQ, "param", NP, "np", range, "5PERC.pdf"), width=6, height=4)
ggplot(sims_summed, aes(x = as.Date(date), y = cases)) +
  geom_line(color = "#008080") +
  geom_ribbon(aes(ymax = upper, ymin = lower), alpha = 0.5, fill = "#008080") +
  geom_line(data=data, aes(x=as.Date(date), y=cases), color = "deeppink2") +
  ggtitle(paste("Log-likelihood =", round(y.out$logLik, digits=2))) +
  xlab("Time (monthly)") + ylab("Local dengue cases") +
  scale_x_date(date_breaks = "12 month", date_labels =  "%Y") +
  theme(plot.title = element_text(size=10,  face="bold"),
        axis.text = element_text(size = 6), axis.title = element_text(size = 8)) +
  theme_minimal() 
dev.off()


|