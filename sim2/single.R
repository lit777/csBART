## SIM2


library(Rcpp); library(MCMCpack); library(rootSolve)
rm(list=ls())
sourceCpp("src/MCMC_mar.cpp")

num_trial = 200
n=300; P=100; m=100; n.iter=20000

rcpp_mar = vector(mode = "numeric", length = num_trial)

for(test_case in 1:200) {
  cat("Testing ", test_case, "of", num_trial, "at", format(Sys.time(), "%H:%M:%S"), "\n")
  # mar ------------
  source("source/data_mar.R")
  # run rcpp
  rcpp = MCMC_mar(Xpred, Y_trt, Y_out, p.grow, p.prune, p.change, m, nu, lambda, lambda_1, dir.alpha,  alpha, beta, n.iter)
  rcpp_mar[test_case] = mean(rcpp$Effect)
}

mean(rcpp_mar)

