## Main 2013

library(Rcpp); library(MCMCpack); library(rootSolve)
rm(list=ls())
sourceCpp("src/MCMC_mar.cpp")

m=200; n.iter=100000


  # mar ------------
  source("source/data_mar.R")
  # run rcpp
  rcpp = MCMC_mar(Xpred, Y_trt, Y_out, p.grow, p.prune, p.change, m, nu, lambda, lambda_1, dir.alpha, dir.alpha1, alpha, beta, n.iter)
  mean(rcpp$Effect)
