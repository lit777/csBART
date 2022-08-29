## SIM 5

library(Rcpp); library(MCMCpack); library(rootSolve)
rm(list=ls())
sourceCpp("src/MCMC_sep.cpp")

num_trial = 200
n=500; P=100; m=100; n.iter=20000

rcpp_sep = vector(mode = "numeric", length = num_trial)

for(test_case in 1:200) {
  cat("Testing ", test_case, "of", num_trial, "at", format(Sys.time(), "%H:%M:%S"), "\n")
  
  # sep --------------
  source("source/data_sep.R")
  
  rcpp = MCMC_sep(Xpred, Y_trt, Y_out, p.grow, p.prune, p.change, m, nu, lambda, lambda_1, lambda_0, dir.alpha, alpha, beta, n.iter)
  rcpp_sep[test_case] = mean(rcpp$Effect)
}

mean(rcpp_sep)
