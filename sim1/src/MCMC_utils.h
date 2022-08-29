// MCMC_utils.h
#pragma once

#include <Rcpp.h>

using namespace Rcpp;

IntegerVector init_seq(const int n_iter, const int thin, const int burn_in)
{
    IntegerVector seq;
    for (int i = 0; i < n_iter; i++)
    {
        int temp = i - n_iter / 2;
        if (i >= burn_in && temp % thin == 0)
        {
            seq.push_back(i);
        }
    }
    return seq;
}

void log_with_LB(NumericVector& res, const NumericVector& src)
{
    const double LB     = pow(0.1, 300);
    const double log_LB = -300 * log(10);
    for (int i = 0; i < src.length(); i++)
    {
        if (src(i) > LB)
            res(i) = log(src(i));
        else
            res(i) = log_LB;
    }
}

NumericVector rowSums_without(const NumericMatrix& src, const int idx)
{
    NumericVector res = rowSums(src);
    for (int i = 0; i < src.nrow(); i++)
    {
        res(i) -= src(i, idx);
    }
    return res;
}

void update_R(NumericVector& R, const NumericVector& Z, const NumericMatrix& Tree, const int t)
{
    NumericVector mu = rowSums(Tree);
    for (int i = 0; i < Tree.nrow(); i++)
    {
        R(i) = Z(i) - mu(i) + Tree(i, t);
    }
}

void update_Z(NumericVector& Z, const NumericVector& Y_trt, const NumericMatrix& Tree)
{
    NumericVector mu = rowSums(Tree);
    double Ystar;
    for (int i = 0; i < Tree.nrow(); i++)
    {
        Ystar = R::rnorm(mu(i), 1);
        Z(i)  = Y_trt(i) * std::max(Ystar, 0.0) + (1 - Y_trt(i)) * std::min(Ystar, 0.0);
    }
}
