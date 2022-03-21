# csBART
BART confounder selection

To run the proposed models for each simulation scenario
- set the working directory to sim"X" where X denotes the scenario number
- run test_R_mar.R for the marginal model
- run sep_R_sep.R for the separate model

To run the proposed models for the application study (year 2013)
- set the working directory to "app2013"
- run MCMC_mar.R for the marginal model
- run MCMC_sep.R for the separate model

To extract the census data
- run census_data_extract.R in the 'data' subfolder

To extract the weather data
- run weather_data.R in 'data' subfolder. (But, the weather source files are not available here). 
- to run the main functions (e.g., test_R_mar.R and test_R_sep.R), post-processed data (data2013) is provided in the 'data' subfolder and those are sufficient to run the functions.
