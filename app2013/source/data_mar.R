### Load the master dataset
load("data/data2013.RData")

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

### Data Preparation
Y_trt <- as.numeric(Master$A)
Y_out <- as.numeric(Master$pm2013)

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



#####################################################################
############# Run main MCMC #########################################
#####################################################################

Effect <- NULL
Tree <- matrix(0,nrow=n, ncol=m)
Tree1 <- matrix(0,nrow=n, ncol=m)
Sigma2 <- Sigma2_1 <- Sigma2_0 <- NULL    
Sigma2[1] <- sigma2
Sigma2_1[1] <- sigma2_1

R <- Z  # Initial values of R
R1 <- Y_out  # Initial values of R1

post.dir.alpha <- rep(1, P)
post.dir.alpha1  <- rep(1, P+1) # Initial values for the selection probabilities

seq <- seq(n.iter/2, n.iter-1, by=10)  # thin = 10, burn-ins = n.iter/2
ind <- rep(0,P+1) # num. of inclusion for each confounder


#------ Place-holders for the posterior samples
Obs_list <- Obs1_list <- list()
for(i in 1:m){
  Obs_list[[i]] <-   Obs1_list[[i]] <- rep(NA, n)
}
Tree <- Tree1 <- matrix(0,nrow=n, ncol=m)
Tree11 <- matrix(0,nrow=n*2, ncol=m)
dt_list <- dt1_list <- list()
for(i in 1:m){
  dt_list[[i]] <- dt1_list[[i]] <- list(position=rep(1,1), parent=rep(NA,1), Terminal=rep(1,1), Split=rep(NA,1), Value=rep(NA,1), MU=rep(NA,1), begin=rep(1,1), end=rep(n,1))
}

prop.prob <- rdirichlet(1, rep(dir.alpha, P+1))
xpred.mult <- list()
xpred.ind <- sample(1:dim(Xpred)[1], n*2, replace=TRUE)

for(i in 2:(P+1)){
  xpred.mult[[i]] <- Xpred1[xpred.ind,i]
}
xpred.mult[[P+2]] <- 1:(2*n)


xpred_mult_mat = matrix(rep(0, (2*n*(P+2))), nrow=2*n)
for(i in 2:(P+2)){
  xpred_mult_mat[,i] = xpred.mult[[i]]
}
xpred_mult_mat[,(P+2)] = xpred_mult_mat[,(P+2)] - 1

Xpred.list <- Xpred1.list <- list()
for(i in 1:P){
  Xpred.list[[i]] <- Xpred[,i]
  Xpred1.list[[i]] <- Xpred1[,i]
}
Xpred1.list[[P+1]] <- Xpred1[,(P+1)]



