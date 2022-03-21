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
n <- length(Y_out)


Xpred1 <- Xpred[which(Y_trt==1),]  # potential confounders (A=1)
Xpred0 <- Xpred[which(Y_trt==0),]  # potential confounders (A=0)
Xcut1 <- lapply(1:dim(Xpred1)[2], function(t) sort(unique(Xpred1[,t]))) # e.g. unique values of predictors
Xcut0 <- lapply(1:dim(Xpred0)[2], function(t) sort(unique(Xpred0[,t]))) # e.g. unique values of predictors

Y1 <- Y_out[which(Y_trt==1)]
Y0 <- Y_out[which(Y_trt==0)]
n1 <- length(Y1)
n0 <- length(Y0)

Z <- rnorm(n, qnorm(mean(Y_trt)), 1)    # latent variable N(0,1)


### Initial Setup (priors, initial values and hyper-parameters)
p.grow <- 0.28            # Prob. of GROW
p.prune <- 0.28           # Prob. of PRUNE
p.change <- 0.44          # Prob. of CHANGE
#m <- 100                  # Num. of Trees: default setting 100
sigma2_1 <- var(Y1)       # Initial value of SD^2
sigma2_0 <- var(Y0) 

nu <- 3                   # default setting (nu, q) = (3, 0.90) from Chipman et al. 2010
f <- function(lambda) invgamma::qinvgamma(0.90, nu/2, rate = lambda*nu/2, lower.tail = TRUE, log.p = FALSE) - sqrt(sigma2_1)
lambda_1 <- uniroot.all(f, c(0.1^5,10))
f <- function(lambda) invgamma::qinvgamma(0.90, nu/2, rate = lambda*nu/2, lower.tail = TRUE, log.p = FALSE) - sqrt(sigma2_0)
lambda_0 <- uniroot.all(f, c(0.1^5,10))

sigma2 <- 1
f <- function(lambda) invgamma::qinvgamma(0.90, nu/2, rate = lambda*nu/2, lower.tail = TRUE, log.p = FALSE) - sqrt(sigma2)
lambda <- uniroot.all(f, c(0.1^5,10))

sigma_mu_1 <- max((min(Y1)/(-2*sqrt(m)))^2, (max(Y1)/(2*sqrt(m)))^2) # sigma_mu based on min/max of Y (A=1)
sigma_mu_0 <- max((min(Y0)/(-2*sqrt(m)))^2, (max(Y0)/(2*sqrt(m)))^2) # sigma_mu based on min/max of Y (A=0)
sigma_mu   <- max((min(Z) /(-2*sqrt(m)))^2, (max(Z) /(2*sqrt(m)))^2) # sigma_mu based on min/max of Y

dir.alpha <- dir.alpha0 <- dir.alpha1 <- 5         # Hyper-parameter on selection probabilities
alpha <- 0.95             # alpha (1+depth)^{-beta} where depth=0,1,2,...
beta <- 2                 # default setting (alpha, beta) = (0.95, 2)


#####################################################################
############# Run main MCMC #########################################
#####################################################################

Effect <- NULL
#n.iter <- 25000            # Num. of Iterations
Tree <- matrix(0,nrow=n, ncol=m)
Tree1 <- matrix(0,nrow=n1, ncol=m)
Tree0 <- matrix(0,nrow=n0, ncol=m)
Sigma2 <- Sigma2_1 <- Sigma2_0 <- NULL            # Variance parameter
Sigma2[1] <- sigma2
Sigma2_1[1] <- sigma2_1
Sigma2_0[1] <- sigma2_0

R  <- Z  # Initial values of R
R1 <- Y1
R0 <- Y0

post.dir.alpha <- post.dir.alpha1 <- post.dir.alpha0 <- rep(1,P) # Initial values for the selection probabilities

seq <- seq(n.iter/2, n.iter, by=10)  # thin = 10, burn-ins = n.iter/2
ind <- rep(0, P)

Obs_list <- Obs1_list <-  Obs0_list <- list()
for(i in 1:m){
  Obs_list[[i]] <- 1:n
  Obs1_list[[i]] <- 1:n1
  Obs0_list[[i]] <- 1:n0
}
Tree <- matrix(0,nrow=n, ncol=m)
Tree1 <- matrix(0,nrow=n1, ncol=m)
Tree0 <- matrix(0,nrow=n0, ncol=m)
Tree11 <- matrix(0,nrow=n*2, ncol=m)
Tree00 <- matrix(0,nrow=n*2, ncol=m)

dt_list <- dt1_list <- dt0_list <- list()
for(i in 1:m){
  dt_list[[i]] <- list( position=rep(1,1), parent=rep(NA,1), Terminal=rep(1,1), Split=rep(NA,1), Value=rep(NA,1), MU=rep(NA,1), begin=rep(1,1), end=rep(n,1))
  dt1_list[[i]] <- list(position=rep(1,1), parent=rep(NA,1), Terminal=rep(1,1), Split=rep(NA,1), Value=rep(NA,1), MU=rep(NA,1), begin=rep(1,1), end=rep(n1,1))
  dt0_list[[i]] <- list( position=rep(1,1), parent=rep(NA,1), Terminal=rep(1,1), Split=rep(NA,1), Value=rep(NA,1), MU=rep(NA,1), begin=rep(1,1), end=rep(n0,1))
}

prop.prob <- rdirichlet(1, rep(dir.alpha,P))
xpred.mult <- list()
xpred.ind <- sample(1:dim(Xpred)[1], n*2, replace=TRUE)
for(i in 1:P){
  xpred.mult[[i]] <- Xpred[xpred.ind,i]
}
xpred.mult[[P+1]] <- 1:(2*n)

xpred_mult_mat = matrix(rep(0, (2*n*(P+1))), nrow=2*n)
for(i in 1:(P+1)){
  xpred_mult_mat[,i] = xpred.mult[[i]]
}
xpred_mult_mat[,(P+1)] = xpred_mult_mat[,(P+1)] - 1

Xpred.list <- Xpred1.list <- Xpred0.list <- list()
for(i in 1:P){
  Xpred.list[[i]] <- Xpred[,i]
  Xpred1.list[[i]] <- Xpred1[,i]
  Xpred0.list[[i]] <- Xpred0[,i]
}
