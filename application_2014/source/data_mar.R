#n=300
#P=100

#------ Load required libraries
library(zipcodeR)
library(data.table)
library(MCMCpack)
library(rootSolve)
library(truncnorm)
library(msm)
library(ordinal)
library(RSpectra)
library(fields)
library(readr)
library(fst)
library(tidyverse)

#------ Load the main datasets
# PM2.5
#pm2013 <- read_rds("data/pm25_2013.rds")
#pm2013 <- as.data.table(pm2013[,c(2,6)])[,-3]
#setnames(pm2013, c("ZCTA5","zip_concentrations"),c("zip","pm2013"))
#pm2014 <- read_rds("data/pm25_2014.rds")
#pm2014 <- as.data.table(pm2014[,c(2,6)])[,-3]
#setnames(pm2014, c("ZCTA5","zip_concentrations"),c("zip","pm2014"))

# CENSUS 2010
#census <- read.csv("data/subset_census_2010_with_zip.csv")
#census$GEOID <- formatC(census$GEOID, width = 5, format = "d", flag = "0")
#setnames(census, "GEOID","zip")
#census <- census[,-c(1,36:39)]
#census <- census[!duplicated(census$zip), ]


# Weather data
#load("data/weather2014.RData")
#setnames(Weather, "id","zip")

# HyADS exposures
#hyads2013 <- read.fst("data/zips_exposures_total_2013.fst", as.data.table = T)
#hyads2013 <- hyads2013[,c(1,2)]
#setnames(hyads2013, c("ZIP","hyads"), c("zip","hyads2013"))

#hyads2014 <- read.fst("data/zips_exposures_total_2014.fst", as.data.table = T)
#hyads2014 <- hyads2014[,c(1,2)]
#setnames(hyads2014, c("ZIP","hyads"), c("zip","hyads2014"))

#Master <- Reduce(function(x,y) merge(x = x, y = y, by = "zip"), list(pm2014, hyads2014, census, Weather[,1:217]))
#Master <- na.omit(Master)

#------ Restriction on the main dataset
#Northeast = c("ME", "NH", "VT", "NY", "PA", "DE", "NJ", "MD", "DC", "VA", "MA", "CT", "RI")
#IndustrialMidwest = c("WV", "OH", "KY", "IN", "IL", "WI", "MI")
#Southeast = c("FL", "GA", "SC", "NC", "TN", "AL", "MS", "AR","LA")

#zipcode <- search_state(c(Northeast, IndustrialMidwest, Southeast))
#zipcode <- subset(zipcode, zipcode_type != "PO Box")
#zipcode <- zipcode[c("zipcode","state")]
#Master.red <- merge(Master, zipcode, by.x="zip", by.y="zipcode")
#Master.red <- subset(Master.red, !(state %in% c("ME", "FL", "AR", "LA", "MS","VT","NH","WI")))
#Master <- as.data.frame(dplyr::select(as_tibble(Master.red), -state))

#Master$A <- ifelse(Master$hyads2014 < 980249908, 1,0 )  # Maybe going back to Median

load("data/Master2014.RData")

### Data Preparation
Y_trt <- as.numeric(Master$A)
Y_out <- as.numeric(Master$pm2014)

Xpred <- as.matrix(Master[,4:252])
Xpred <-  Xpred[,!str_detect(colnames(Xpred), pattern="^denom")]  # discard "denom"
Xpred[,3] <- Xpred[,3] / Xpred[,2] * 100 # housing_unit_urban
Xpred <- Xpred[,-c(2,5,6,7,8,13,14,16,17,18,21,22)] # discard housing unit and count variable

Xpred <- Xpred[,-5] # discard population black (vs population white)
Xpred <- Xpred[,-12] # discard bachelors_male (vs bachelors_female)
Xpred <- Xpred[,-c(8, 13)] # discard income_over200K and median_housing_value (vs medain_income)
Xpred <- Xpred[,-c(9, 10)] # discard highschool graduate (vs bachelors)

Xpred <- Xpred[,!str_detect(colnames(Xpred), pattern="^stemp_")]   # discard stemp (vs temp)
Xpred <- Xpred[,!str_detect(colnames(Xpred), pattern="^cpcp_")]   # discard cpcp (vs apcp)
Xpred <- Xpred[,!str_detect(colnames(Xpred), pattern="^dswrf_")]  # dicard dswrf (vs tcdc)
Xpred <- Xpred[,!str_detect(colnames(Xpred), pattern="^vwnd_")]  # dicard both wind components (vs wspd, phi)
Xpred <- Xpred[,!str_detect(colnames(Xpred), pattern="^uwnd_")]  # dicard both wind components (vs wspd, phi)
Xpred <- Xpred[,!str_detect(colnames(Xpred), pattern="^phi_.\\d{1}")]  # discard wind angle for neighboring locations
Xpred <- Xpred[,!str_detect(colnames(Xpred), pattern="^wspd_.\\d{1}")]  # discard wind speed for neighboring locations

Xpred[,2] <- round(Xpred[,2],1)
Xpred[,3] <- round(Xpred[,3]/1000, 1)
Xpred[,4] <- round(Xpred[,4],1)
Xpred[,5] <- round(Xpred[,5],1)
Xpred[,6] <- round(Xpred[,6],1)
Xpred[,8] <- round(Xpred[,8]/1000, 1)
Xpred[,9] <- round(Xpred[,9],1)
Xpred[,10] <- round(Xpred[,10],1)

Xpred <- apply(Xpred, 2, function(x) round(x, 3))

# total 104 potential confounders
S <- sample(1:104, 104, replace=F)
Xpred <- Xpred[,S]

P <- dim(Xpred)[2]
Xcut <- lapply(1:dim(Xpred)[2], function(t) sort(unique(Xpred[,t]))) # e.g. unique values of


shift <- mean(Y_out)
Y_out <- Y_out - shift


Xpred1 <- cbind(Y_trt, Xpred)
Xcut1 <- lapply(1:dim(Xpred1)[2], function(t) sort(unique(Xpred1[,t]))) # e.g. unique values of predictors


Y1 <- Y_out
n1 <- length(Y1)
n <- n1 <- length(Y1)

Z <- rnorm(n, qnorm(mean(Y_trt)), 1)    # latent variable N(0,1)


#------ Initial Setup (priors, initial values and hyper-parameters)
p.grow <- 0.28            # Prob. of GROW
p.prune <- 0.28           # Prob. of PRUNE
p.change <- 0.44          # Prob. of CHANGE
sigma2_1 <- var(Y_out)    # Initial value of SD^2
nu <- 3                   # default setting (nu, q) = (3, 0.90) from Chipman et al. 2010
f <- function(lambda) invgamma::qinvgamma(0.90, nu/2, rate = lambda*nu/2, lower.tail = TRUE, log.p = FALSE) - sqrt(sigma2_1)
lambda_1 <- uniroot.all(f, c(0.1^5,10))
sigma_mu_1 <- max((min(Y_out)/(-2*sqrt(m)))^2, (max(Y_out)/(2*sqrt(m)))^2) # sigma_mu based on min/max of Y

sigma2 <- 1               # sigma for the exposure model
f <- function(lambda) invgamma::qinvgamma(0.90, nu/2, rate = lambda*nu/2, lower.tail = TRUE, log.p = FALSE) - sqrt(sigma2)
lambda <- uniroot.all(f, c(0.1^5,10))
sigma_mu <- max((min(Z)/(-2*sqrt(m)))^2, (max(Z)/(2*sqrt(m)))^2) # sigma_mu based on min/max of Y

dir.alpha  <- dir.alpha1 <- 0.3         # Hyper-parameter on selection probabilities
alpha <- 0.95             # alpha (1+depth)^{-beta} where depth=0,1,2,...
beta <- 2                 # default setting (alpha, beta) = (0.95, 2)


