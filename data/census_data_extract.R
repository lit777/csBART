install.packages("tidycensus")
library(tidyverse)
library(tidycensus)
library(data.table)
library(rio)

#census_api_key("YOUR API KEY GOES HERE")

variables <- load_variables(2010, "sf1", cache = TRUE)
# housing_unit="H001001"
# housing_unit_urban="H002002"
# population="P001001"
# population_white="P003002"
# population_black="P003003"
# population_male="P012002"
# population_under20="P014001"
# age_median="P013001"
variables_acs <- load_variables(2010, "acs5", cache = TRUE)
# population_under18="B09001_001"
# poverty_status="B14006_002"
# highschool_gradudate_male="B15002_011"
# highschool_gradudate_female="B15002_028"
# bachelors_male="B15002_015"
# bachelors_female="B15002_032"
# income_over200K="B19001_017"
# median_income="B19013_001"
# gini_index="B19083_001"
# median_housing_value="B25107_001"

confounders <- get_decennial(geography = "zip code tabulation area", 
              variables = c(housing_unit="H001001",housing_unit_urban="H002002",population="P001001",population_white="P003002",population_black="P003003",population_male="P012002",population_agelessthan20="P014001",age_median="P013001"), 
              year = 2010)
census <- as_tibble(confounders) %>% select(GEOID, variable,value) %>% spread(variable, value)
confounders_acs <- get_acs(geography = "zip code tabulation area", 
                variables = c(population_under18="B09001_001", poverty_status="B14006_002", highschool_gradudate_male="B15002_011", highschool_gradudate_female="B15002_028", bachelors_male="B15002_015",bachelors_female="B15002_032",income_over200K="B19001_017",median_income="B19013_001",gini_index="B19083_001",median_housing_value="B25107_001"), 
                geometry = TRUE)
census_acs <- as_tibble(confounders_acs) %>% select(GEOID, variable, estimate) %>% spread(variable, estimate)

CENSUS <- data.table(merge(census, census_acs, by="GEOID"))

crosswalk <- data.table(import("Zip_to_ZCTA_Crosswalk_JSI2014.xlsx"))
crosswalk <- subset(crosswalk, !is.na(ZCTA_crosswalk))
setkey(crosswalk, ZCTA_crosswalk)
setkey(CENSUS, GEOID)
subset_census_2010_with_zip <- CENSUS[crosswalk]

write.csv(subset_census_2010_with_zip,
          "subset_census_2010_with_zip.csv")
