## SIM 2

library(Rcpp); library(MCMCpack); library(rootSolve)
rm(list=ls())
sourceCpp("src/MCMC_sep.cpp")

num_trial = 200
n=300; P=100; m=100; n.iter=10000

rcpp_sep = vector(mode = "numeric", length = num_trial)
cover = rep(0, num_trial)

for(test_case in 1:200) {
  cat("Testing ", test_case, "of", num_trial, "at", format(Sys.time(), "%H:%M:%S"), "\n")
  
  # sep --------------
  source("source/data_sep.R")
  
  rcpp = MCMC_sep(Xpred, Y_trt, Y_out, p.grow, p.prune, p.change, m, nu, lambda, lambda_1, lambda_0, dir.alpha, dir.alpha1, dir.alpha0, alpha, beta, n.iter)
  rcpp_sep[test_case] = mean(rcpp$Effect)
  ci <- quantile(rcpp$Effect, c(0.025, 0.975))
  if(ci[1] < -1.3989 & ci[2] > -1.3989){cover[test_case] = 1}
}


true <- -1.3989
mean(rcpp_sep - true)
mean((rcpp_sep-true)^2)
mean(cover)
