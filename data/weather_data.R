library(readr)
library(tidyverse)
weather <- read_rds("weather/2013.rds")
weather_sub <- filter(weather, season %in% c("Summer", "Winter"))

direc <- c("W",  "SW", "NW", "S",  "C",  "N",  "NE", "SE", "E")

weather_summer_C <- filter(weather_sub, season == 'Summer', direction == 'C', radius == 0)
names(weather_summer_C)[8:19] <- c('temp_SC','stemp_SC', 'apcp_SC','cpcp_SC',     
                                           'tcdc_SC','dswrf_SC','hpbl_SC', 'rhum_SC',
                                          'vwnd_SC', 'uwnd_SC','wspd_SC','phi_SC')
weather_summer_C <- weather_summer_C[,c(1,8:19)] 

weather_winter_C <- filter(weather_sub, season == 'Winter', direction == 'C', radius == 0)
names(weather_winter_C)[8:19] <- c('temp_WC','stemp_WC', 'apcp_WC','cpcp_WC',     
                                       'tcdc_WC','dswrf_WC','hpbl_WC', 'rhum_WC',
                                       'vwnd_WC', 'uwnd_WC','wspd_WC','phi_WC')
weather_winter_C <- weather_winter_C[,c(1,8:19)] 

for(i in c(1,2,3,4,6,7,8,9)){
  eval(parse(text=(paste0("weather_summer_",direc[i],"_100 <- filter(weather_sub, season == 'Summer', direction == '",direc[i],"', radius == 100)"))))
  eval(parse(text=(paste0("names(weather_summer_",direc[i],"_100)[8:19] <- c('temp_S",i,"','stemp_S",i,"', 'apcp_S",i,"','cpcp_S",i,"',     
                                           'tcdc_S",i,"','dswrf_S",i,"','hpbl_S",i,"', 'rhum_S",i,"',
                                          'vwnd_S",i,"', 'uwnd_S",i,"','wspd_S",i,"','phi_S",i,"')"))))
  eval(parse(text=(paste0("weather_summer_",direc[i],"_100 <- weather_summer_",direc[i],"_100[,c(1,8:19)]"))))
}
for(i in c(1,2,3,4,6,7,8,9)){
  eval(parse(text=(paste0("weather_winter_",direc[i],"_100 <- filter(weather_sub, season == 'Winter', direction == '",direc[i],"', radius == 100)"))))
  eval(parse(text=(paste0("names(weather_winter_",direc[i],"_100)[8:19] <- c('temp_W",i,"','stemp_W",i,"', 'apcp_W",i,"','cpcp_W",i,"',     
                                           'tcdc_W",i,"','dswrf_W",i,"','hpbl_W",i,"', 'rhum_W",i,"',
                                          'vwnd_W",i,"', 'uwnd_W",i,"','wspd_W",i,"','phi_W",i,"')"))))
  eval(parse(text=(paste0("weather_winter_",direc[i],"_100 <- weather_winter_",direc[i],"_100[,c(1,8:19)]"))))
}
for(i in c(1,2,3,4,6,7,8,9)){
  eval(parse(text=(paste0("weather_summer_",direc[i],"_300 <- filter(weather_sub, season == 'Summer', direction == '",direc[i],"', radius == 300)"))))
  eval(parse(text=(paste0("names(weather_summer_",direc[i],"_300)[8:19] <- c('temp_300S",i,"','stemp_300S",i,"', 'apcp_300S",i,"','cpcp_300S",i,"',     
                                           'tcdc_300S",i,"','dswrf_300S",i,"','hpbl_300S",i,"', 'rhum_300S",i,"',
                                          'vwnd_300S",i,"', 'uwnd_300S",i,"','wspd_300S",i,"','phi_300S",i,"')"))))
  eval(parse(text=(paste0("weather_summer_",direc[i],"_300 <- weather_summer_",direc[i],"_300[,c(1,8:19)]"))))
}
for(i in c(1,2,3,4,6,7,8,9)){
  eval(parse(text=(paste0("weather_winter_",direc[i],"_300 <- filter(weather_sub, season == 'Winter', direction == '",direc[i],"', radius == 300)"))))
  eval(parse(text=(paste0("names(weather_winter_",direc[i],"_300)[8:19] <- c('temp_300W",i,"','stemp_300W",i,"', 'apcp_300W",i,"','cpcp_300W",i,"',     
                                           'tcdc_300W",i,"','dswrf_300W",i,"','hpbl_300W",i,"', 'rhum_300W",i,"',
                                          'vwnd_300W",i,"', 'uwnd_300W",i,"','wspd_300W",i,"','phi_300W",i,"')"))))
  eval(parse(text=(paste0("weather_winter_",direc[i],"_300 <- weather_winter_",direc[i],"_300[,c(1,8:19)]"))))
}

weather.list = list(weather_summer_C, weather_winter_C,
                    weather_summer_N_100, weather_winter_N_100,
                    weather_summer_S_100, weather_winter_S_100,
                    weather_summer_NW_100, weather_winter_NW_100,
                    weather_summer_SW_100, weather_winter_SW_100,
                    weather_summer_SE_100, weather_winter_SE_100,
                    weather_summer_NE_100, weather_winter_NE_100,
                    weather_summer_E_100, weather_winter_E_100,
                    weather_summer_W_100, weather_winter_W_100,
                    weather_summer_N_300, weather_winter_N_300,
                    weather_summer_S_300, weather_winter_S_300,
                    weather_summer_NW_300, weather_winter_NW_300,
                    weather_summer_SW_300, weather_winter_SW_300,
                    weather_summer_SE_300, weather_winter_SE_300,
                    weather_summer_NE_300, weather_winter_NE_300,
                    weather_summer_E_300, weather_winter_E_300,
                    weather_summer_W_300, weather_winter_W_300)
Weather <- Reduce(function(x,y) merge(x = x, y = y, by = "id"), weather.list)


save(Weather, file="weather2013.RData")
