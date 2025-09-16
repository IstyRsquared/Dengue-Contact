## Name: Isty Rysava
## Date: 5/3/2024
## Code: Explore tidycensus to extract population variables for FL counties

rm(list=ls())

### Libraries
library(tidycensus)
library(tidyverse)
library(stringr)
library(sf)

### Identify variables of interest
key <-"07605ba121a86ee059266976f8f46e06b6aff604" # API key
vars <-load_variables(year=2009,
                      dataset = c("acs5"),
                      cache = FALSE) 
# write.csv(vars, "output/tidycensus_vars.csv", row.names=F)

vars_to_pull <-c("B01001_026", "B01001_002", # pop by Sex: female, male
                 "B01003_001", # Total population 
                 
                 "B05001_002",	# U.S. citizen, born in the United States
                 "B05001_003",	# U.S. citizen, born in Puerto Rico or U.S. Island Areas
                 "B05001_004",	# U.S. citizen, born abroad of American parent(s)
                 "B05001_005",	# U.S. citizen by naturalization
                 "B05001_006",	# Not a U.S. citizen
                 

                 "B11002_001","B11002_002", "B11002_012", # Living total, in family vs non family HHs
                 # "B17001_001", "B17001_002", "B17001_031", # Poverty level ?
                 "B19013_001", # Medium HH income in inflation adjusted dollars per year 
                 #"B19083_001", # GINI 
                 "B19113_001", # Median family income in inflation adjusted dollars per year 
                 "B19301_001", # Per capita income adjusted dollars per year 
                 
                 # "B19051_001", "B19051_002", # Income earnings in the last 12 months (total) and with earnings 
                 
                 "B25002_001", "B25002_002","B25002_003", # Occupancy status: vacant vs occupied
                 # "B23006_001", "B23006_002", "B23006_009", "B23006_016", "B23006_023", # education
                 
                 ### Additional vars added in Sept 2023
                 "B25008_001", "B25008_002", "B25008_003",  # Total population in occupied housing: owner and renter occupied
                 "B25010_001", # Average HH size of occupied housing units
                 "B25021_001", # Medium number of rooms 
                 "B25022_001", # Aggregate number of rooms
                 "B25014_001", # Occupants per room
                 "B25047_001", "B25047_002" # Plumbing (total) and complete
)

### Prep list to collect info from all years
years_to_pull <-c(2009:2022) 
acs_dat <-list()
vars_to_pull <- unlist(vars_to_pull) 

# pull
for(i in 1:length(years_to_pull)){
  acs_dat[[i]] <- get_acs(geography = "county", 
                         variables = vars_to_pull,
                         #key=key,
                         state = "Florida", 
                         year = years_to_pull[i]) %>%
    mutate(year=years_to_pull[i]) %>%
    left_join(., vars, c("variable"="name"))
}

acs_dat <- do.call(rbind.data.frame, acs_dat)
acs_dat$NAME <- gsub(" County, Florida", "", acs_dat$NAME)
head(acs_dat); tail(acs_dat)

## save
save(acs_dat, file="output/variables/acs_dat_2009_2022.Rdata")
write.csv(acs_dat, file="output/variables/acs_dat_2009_2022.csv", row.names=F)

### Get density
acs_density <-list()
for(i in 1:length(years_to_pull)){
  acs_density[[i]] <-get_acs(geography = "county", 
                             variables = "B01003_001", # Total pop
                             state = "Florida", 
                             geometry = TRUE,
                             keep_geo_vars = TRUE,
                             year = years_to_pull[i]) %>%
    mutate(year=years_to_pull[i])
  
}

land_2011 <- acs_density[[3]] %>%
  select(GEOID, ALAND) %>%
  st_drop_geometry()

acs_density[[1]] <- acs_density[[1]] %>%
  left_join(., land_2011, by=c("GEOID"))

acs_density[[2]] <- acs_density[[2]] %>%
  left_join(., land_2011, by=c("GEOID"))

for(i in 1:length(years_to_pull)){
  acs_density[[i]] <- acs_density[[i]] %>%
    mutate(density = 1000*(estimate/ALAND)) %>%
    dplyr::select(GEOID, density, year)
}

acs_density <-do.call(rbind.data.frame, acs_density) %>%
  st_drop_geometry()
head(acs_density); tail(acs_density)
save(acs_density, file="output/variables/pop_density_2009_2022.Rdata")
write.csv(acs_density, file="output/variables/pop_density_2009_2022.csv", row.names=F)

