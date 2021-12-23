#------ Load required libraries
library(zipcodeR)
library(data.table)
library(MCMCpack)
library(rootSolve)
library(RSpectra)
library(fst)
library(tidyverse)


#------ Import Files
source("source/PRUNE_mar.R")
source("source/CHANGE_mar.R")
source("source/Common.R")
source("source/Prediction_mar.R")
source("source/GROW_mar.R")


#------ Load the main datasets
# PM2.5
pm2013 <- read_rds("data/pm25_2013.rds")
pm2013 <- pm2013[,c(2,6)]
setnames(pm2013, c("ZCTA5","zip_concentrations"),c("zip","pm2013"))
pm2014 <- read_rds("data/pm25_2014.rds")
pm2014 <- pm2014[,c(2,6)]
setnames(pm2014, c("ZCTA5","zip_concentrations"),c("zip","pm2014"))

# CENSUS 2010
census <- read.csv("data/subset_census_2010_with_zip.csv") 
census$GEOID <- formatC(census$GEOID, width = 5, format = "d", flag = "0")
setnames(census, "GEOID","zip")
census <- census[,-c(1,21:24)]
census <- census[!duplicated(census$zip), ]

# Weather data (2013)
load("data/weather2013.RData")
# load("data/weather2014.RData")
setnames(Weather, "id","zip")

# HyADS exposures
hyads2013 <- read.fst("data/zips_exposures_total_2013.fst", as.data.table = T)
hyads2013 <- hyads2013[,c(1,2)]
setnames(hyads2013, c("ZIP","hyads"), c("zip","hyads2013"))

hyads2014 <- read.fst("data/zips_exposures_total_2014.fst", as.data.table = T)
hyads2014 <- hyads2014[,c(1,2)]
setnames(hyads2014, c("ZIP","hyads"), c("zip","hyads2014"))

Master <- Reduce(function(x,y) merge(x = x, y = y, by = "zip"), list(pm2013, hyads2013, census, Weather[,1:217]))
Master <- na.omit(Master)


#------ Restriction on the main dataset
Northeast = c("ME", "NH", "VT", "NY", "PA", "DE", "NJ", "MD", "DC", "VA", "MA", "CT", "RI")
IndustrialMidwest = c("WV", "OH", "KY", "IN", "IL", "WI", "MI")
Southeast = c("FL", "GA", "SC", "NC", "TN", "AL", "MS", "AR","LA")

zipcode <- search_state(c(Northeast, IndustrialMidwest, Southeast))
zipcode <- subset(zipcode, zipcode_type != "PO Box")
zipcode <- zipcode[c("zipcode","state")]
Master.red <- merge(Master, zipcode, by.x="zip", by.y="zipcode")
Master.red <- subset(Master.red, !(state %in% c("ME", "FL", "AR", "LA", "MS","VT","NH","WI")))
Master <- as.data.frame(select(Master.red, -state))[,-238]


#------ Binary exposure (dichotomized based on the mean)
Master$A <- ifelse(Master$hyads2013 < summary(Master$hyads2013)[4], 1,0 )  


#------ Data Preparation
Y_trt <- as.numeric(Master$A)
Y_out <- as.numeric(Master$pm2013)
n <- n1 <- length(Y_out) # num. of observations
Xpred <- as.matrix(Master[,4:237]) # potential confounders
P <- dim(Xpred)[2] # num. of potential confounders
Xcut <- lapply(1:dim(Xpred)[2], function(t) sort(unique(Xpred[,t]))) # e.g. unique values of potential confounders
shift <- mean(Y_out) # shifting the outcome variable
Y_out <- Y_out - shift

Xpred1 <- cbind(Y_trt, Xpred) # potential confounders + exposure
Xcut1 <- lapply(1:dim(Xpred1)[2], function(t) sort(unique(Xpred1[,t]))) # e.g. unique values of confounders + exposure

Z <- rnorm(n, qnorm(mean(Y_trt)), 1)    # latent variable N(0,1) for the exposure model

#------ Initial Setup (priors, initial values and hyper-parameters)
p.grow <- 0.28            # Prob. of GROW
p.prune <- 0.28           # Prob. of PRUNE
p.change <- 0.44          # Prob. of CHANGE
m <- 100                  # Num. of Trees: default setting 100
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
n.iter <- 25000            # Num. of Iterations
Tree <- matrix(0,nrow=n, ncol=m)
Tree1 <- matrix(0,nrow=n, ncol=m)
Sigma2 <- Sigma2_1 <- Sigma2_0 <- NULL    
Sigma2[1] <- sigma2
Sigma2_1[1] <- sigma2_1

R <- Z  # Initial values of R
R1 <- Y_out  # Initial values of R1

post.dir.alpha <- rep(1, P)
post.dir.alpha1  <- rep(1, P+1) # Initial values for the selection probabilities

seq <- seq(n.iter/2, n.iter, by=10)  # thin = 10, burn-ins = n.iter/2
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

Xpred.list <- Xpred1.list <- list()
for(i in 1:P){
  Xpred.list[[i]] <- Xpred[,i]
  Xpred1.list[[i]] <- Xpred1[,i]
}
Xpred1.list[[P+1]] <- Xpred1[,(P+1)]



#------ Run MCMC
for(j in 2:n.iter){

    #------ Exposure Model
    Ystar <- rnorm(n, mean=rowSums(Tree), 1)
    Z <- (Y_trt)*pmax(Ystar,0)+(1-Y_trt)*pmin(Ystar,0)

    for(t in 1:m){
        R <- Z - rowSums(Tree[,-t]) 

        #------ Find the depth of the tree (0 or 1 or 2)
        tree.length <- length(dt_list[[t]]$position)
        if(tree.length == 1){  # tree has no node yet
            grow.step <- GROW.first(sigma2=Sigma2[j-1], sigma_mu=sigma_mu, dt=dt_list[[t]], R=R, prop.prob=prop.prob, Obs=Obs_list[[t]], ind=2)
            dt_list[[t]] <- grow.step$dt
            Obs_list[[t]] <- grow.step$Obs
        } else {
            step <- sample(1:3, 1, prob=c(p.grow, p.prune, p.change))  # Pick a step
            if(step==3){  # CHANGE step
              change.step <- CHANGE(sigma2=Sigma2[j-1], sigma_mu=sigma_mu, dt=dt_list[[t]], R=R, prop.prob=prop.prob, Obs=Obs_list[[t]], ind=2)
              dt_list[[t]] <- change.step$dt
              Obs_list[[t]] <- change.step$Obs
            }else{
            if(step==2){   # PRUNE step
                prune.step <- PRUNE(sigma2=Sigma2[j-1], sigma_mu=sigma_mu, dt=dt_list[[t]], R=R, prop.prob=prop.prob, Obs=Obs_list[[t]], ind=2)
                dt_list[[t]] <- prune.step$dt
                Obs_list[[t]] <- prune.step$Obs
            }else{
            if(step==1){   # GROW step
              grow.step <- GROW(sigma2=Sigma2[j-1], sigma_mu=sigma_mu, dt=dt_list[[t]], R=R, prop.prob=prop.prob,  Obs=Obs_list[[t]], ind=2)
              dt_list[[t]] <- grow.step$dt
              Obs_list[[t]] <- grow.step$Obs
            }}}
        }
        Mean <- Mean.Parameter(Sigma2[j-1], sigma_mu, dt=dt_list[[t]],  Obs=Obs_list[[t]], R, ind=2)
        Tree[,t] <- Mean$T
        dt_list[[t]] <- Mean$dt
    }

    
    #------ Outcome Model
    for(t in 1:m){
      R1 <- Y_out - rowSums(Tree1[,-t])
      
      ### Find the depth of the tree (0 or 1 or 2)
      tree.length1 <- length(dt1_list[[t]]$position)
      if(tree.length1 == 1){  # tree has no node yet
        grow.step <- GROW.first(sigma2=Sigma2_1[j-1], sigma_mu=sigma_mu_1, dt=dt1_list[[t]], R=R1, prop.prob=prop.prob, Obs=Obs1_list[[t]], ind=1)
        dt1_list[[t]] <- grow.step$dt
        Obs1_list[[t]] <- grow.step$Obs
      } else {
        step <- sample(1:3, 1, prob=c(p.grow, p.prune, p.change))  # Pick a step
        if(step==3){  # CHANGE step
          change.step <- CHANGE(sigma2=Sigma2_1[j-1], sigma_mu=sigma_mu_1, dt=dt1_list[[t]], R=R1, prop.prob=prop.prob, Obs=Obs1_list[[t]], ind=1)
          dt1_list[[t]] <- change.step$dt
          Obs1_list[[t]] <- change.step$Obs
        }else{
        if(step==2){   # PRUNE step
          prune.step <- PRUNE(sigma2=Sigma2_1[j-1], sigma_mu=sigma_mu_1, dt=dt1_list[[t]], R=R1, prop.prob=prop.prob, Obs=Obs1_list[[t]], ind=1)
          dt1_list[[t]] <- prune.step$dt
          Obs1_list[[t]] <- prune.step$Obs
        }else{
        if(step==1){   # GROW step
          grow.step <- GROW(sigma2=Sigma2_1[j-1], sigma_mu=sigma_mu_1, dt=dt1_list[[t]], R=R1, prop.prob=prop.prob, Obs=Obs1_list[[t]], ind=1)
          dt1_list[[t]] <- grow.step$dt
          Obs1_list[[t]] <- grow.step$Obs
        }}}
      }
      Mean <- Mean.Parameter(Sigma2_1[j-1], sigma_mu_1, dt=dt1_list[[t]],  Obs=Obs1_list[[t]], R1, ind=1)
      Tree1[,t] <- Mean$T
      dt1_list[[t]] <- Mean$dt
    }


    # Sample variance parameters
    Sigma2[j] <- 1
    Sigma2_1[j] <- rinvgamma(1, nu/2+n/2, scale = nu*lambda_1/2 + sum((Y_out-rowSums(Tree1))^2)/2)

    # Num. of inclusion of each potential confounder
    DT <- unlist(lapply(dt_list, function(x) x$Split))
    add <- as.numeric(table(factor(DT[!is.na(DT)], levels=1:P)))
    DT1 <- unlist(lapply(dt1_list, function(x) x$Split))
    add1 <- as.numeric(table(factor(DT1[!is.na(DT1)], levels=1:(P+1))))

    # M.H. algorithm for the alpha parameter in the dirichlet distribution
    p.dir.alpha <- max(rnorm(1, dir.alpha, 0.1), 0.1^10)
    SumS <-  log(ifelse(prop.prob<0.1^300, 0.1^300, prop.prob))
    dir_lik.p <- sum(SumS*(rep(p.dir.alpha/(P+1), (P+1))-1)) + lgamma(sum(rep(p.dir.alpha/(P+1), (P+1))))-sum(lgamma(rep(p.dir.alpha/(P+1), (P+1)))) 
    dir_lik <- sum(SumS*(rep(dir.alpha/(P+1), (P+1))-1)) + lgamma(sum(rep(dir.alpha/(P+1), (P+1))))-sum(lgamma(rep(dir.alpha/(P+1), (P+1))))
    ratio <- dir_lik.p + log((p.dir.alpha/(p.dir.alpha+(P+1)))^(0.5-1)*abs(1/(p.dir.alpha+(P+1))-p.dir.alpha/(p.dir.alpha+(P+1))^2 ) ) + dnorm(dir.alpha, p.dir.alpha,0.1, log=TRUE) - dir_lik - log((dir.alpha/(dir.alpha+(P+1)))^(0.5-1)*abs(1/(dir.alpha+(P+1))-dir.alpha/(dir.alpha+(P+1))^2 ) ) - dnorm(p.dir.alpha, dir.alpha,0.1, log=TRUE)
    if (log(runif(1))<ratio) {
        dir.alpha <- p.dir.alpha
    }
    post.dir.alpha <- rep(dir.alpha/(P+1), (P+1))
   
    # M.H. algorithm for the inclusion probabilities
    add.max <- max(add)
    p.prop.prob <- rdirichlet(1, add1+post.dir.alpha+c(add.max,add))
    dir_lik.p <- sum(add)*log(1/(1-sum(p.prop.prob[1])))+(add.max+post.dir.alpha[1]-1)*log(pmax(p.prop.prob[1], 0.1^300))+sum((add1[-1]+add+post.dir.alpha[-1]-1)*log(pmax(p.prop.prob[-1],0.1^300)))
    dir_lik <- sum(add)*log(1/(1-sum(prop.prob[1])))+(add.max+post.dir.alpha[1]-1)*log(pmax(prop.prob[1], 0.1^300))+sum((add1[-1]+add+post.dir.alpha[-1]-1)*log(pmax(prop.prob[-1],0.1^300)))
    ratio <- dir_lik.p + sum((add1+post.dir.alpha+c(add.max,add)-1)*log(pmax(0.1^300,prop.prob))) - dir_lik - sum((add1+post.dir.alpha+c(add.max,add)-1)*log(pmax(0.1^300,p.prop.prob)))
    if (log(runif(1))<ratio){
      prop.prob <- p.prop.prob
    }

    # Sampling E(Y(1)-Y(0))
    if(j %in% seq){
      Tree11 <- sapply(1:m, function(x) Mean.Parameter_pred(dt1_list[[x]], 1)$T)
      Tree00 <- sapply(1:m, function(x) Mean.Parameter_pred(dt1_list[[x]], 0)$T)
      Effect[j] <- mean(rowSums(Tree11) - rowSums(Tree00))
      ind.temp<-ifelse(add1 > 0, 1, 0)
      ind <- rbind(ind, ind.temp) # indicators of whether confounders are included
    }
}


