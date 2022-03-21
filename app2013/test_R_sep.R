## Main 2013

library(Rcpp); library(MCMCpack); library(rootSolve)
rm(list=ls())
sourceCpp("src/MCMC_sep.cpp")

m=200; n.iter=100000

  # sep --------------
  source("source/data_sep.R")
  
  rcpp = MCMC_sep(Xpred, Y_trt, Y_out, p.grow, p.prune, p.change, m, nu, lambda, lambda_1, lambda_0, dir.alpha, dir.alpha1, dir.alpha0, alpha, beta, n.iter)
  mean(rcpp$Effect)
