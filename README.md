# csBART
BART confounder selection

To estimate the causal effect
- run MCMC_mar.R for the marginal model
- run MCMC_sep.R for the separate model

To extract the census data
- run census_data_extract.R in the 'data' subfolder

To extract the weather data
- run weather_data.R in 'data' subfolder. (But, the weather source files are not available here). 
- to run the main functions (e.g., MCMC_mar.R or MCMC_sep.R), post-processed data (weather2013 and weather2014) are provided in the 'data' subfolder and those are sufficient to run the functions.
