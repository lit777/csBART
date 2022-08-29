// functions.h
#ifndef __FUNCTIONS_H__
#define __FUNCTIONS_H__

#include <Rcpp.h>

using namespace Rcpp;
using namespace std;
using namespace sugar;

// [[Rcpp::plugins(cpp17)]]
IntegerVector append(const IntegerVector& vector, const int a, const int b);
NumericVector append(const NumericVector& vector, const int a, const int b);
void discard_row(NumericMatrix& res, const NumericMatrix& src, const IntegerVector& idx);
void discard_row(NumericMatrix& res, const NumericMatrix& src, const LogicalVector& idx);
IntegerVector merge(const IntegerVector& x, const IntegerVector& y);
LogicalVector pred_enough_unique(const NumericMatrix& matrix, const IntegerVector& idx);
double sum_by_idx(const NumericVector& src, const IntegerVector& idx);
double sum_union(const NumericVector& src, const IntegerVector& idx1, const IntegerVector& idx2);
int which(const IntegerVector& vector, const int value);
int which(const NumericVector& vector, const double value);
  


IntegerVector append(const IntegerVector& vector, const int a, const int b) {
  // append push back two elements to vector
  int n = vector.length();
  IntegerVector output(n + 2);
  for (int i=0; i<n; i++) {
    output(i) = vector(i);
  }
  output(n)   = a; 
  output(n+1) = b;
  return output;
}
NumericVector append(const NumericVector& vector, const int a, const int b) {
  // append push back two elements to vector
  int n = vector.length();
  NumericVector output(n + 2);
  for (int i=0; i<n; i++) {
    output(i) = vector(i);
  }
  output(n)   = a; 
  output(n+1) = b;
  return output;
}

void discard_row(NumericMatrix& res, const NumericMatrix& src, const IntegerVector& idx) {
  // equivalent to X[-idx, :] in R
  if (res.nrow() != src.nrow() - idx.length()) {
    stop("Error : dimension does not match. \nres row is %i while idx length is %i (discard_row)", 
         res.nrow(), idx.length());
  }
  const int p = src.ncol();
  
  for (int j=0; j<p; j++) {
    for (int i=0; i<res.nrow(); i++) res(i, j) = src(idx(i), j);
  }
}
void discard_row(NumericMatrix& res, const NumericMatrix& src, const LogicalVector& idx) {
  // equivalent to X[-idx, :] in R
  // idx : keep true and remove false
  if (res.nrow() != sum(idx)) {
    stop("Error : dimension does not match. \nres row is %i while sum(idx) is %i (discard_row)", 
         res.nrow(), sum(idx));
  }
  const int p = src.ncol();
  IntegerVector idx_int = Rcpp::Range(0, idx.length()-1);
  idx_int = idx_int[idx];
  
  for (int j=0; j<p; j++) {
    for (int i=0; i<res.nrow(); i++) res(i, j) = src(idx_int(i), j);
  }
}

IntegerVector merge(const IntegerVector& x, const IntegerVector& y) {
  const int nx = x.length();
  const int ny = y.length();
  if(nx == 0) return y;
  if(ny == 0) return x;
  
  IntegerVector res(nx + ny);
  res[Rcpp::Range(0, nx-1)] = x;
  res[Rcpp::Range(nx, nx+ny-1)] = y;
  
  return res;
}

LogicalVector pred_enough_unique(const NumericMatrix& matrix, const IntegerVector& idx) {
  // find idx and their probs of covariates with enough unique obs
  const int P = matrix.ncol();
  LogicalVector flag = rep(false, P);
  for (int j=0; j<P; j++) { // col idx
    double x = matrix(idx(0), j);
    int i = 1;
    do {
      if (x != matrix(idx(i), j)) {
        flag(j) = true;
        break;
      } else {
        i++;
      }
    } while (i < idx.length());
  }
  return flag;
}

double sum_by_idx(const NumericVector& src, const IntegerVector& idx) {
  // equivalent to sum(x[idx]) in R
  if (src.length() <= max(idx)) {
    stop("Error : Incorrect idx, \nmax(idx) %i is greater than src.length() %i (sum_by_idx)", 
         max(idx), src.length());
  }
  NumericVector temp = src[idx];
  return sum(temp);
}

double sum_union(const NumericVector& src, const IntegerVector& idx1, const IntegerVector& idx2) {
  
  IntegerVector idx_union = union_(idx1, idx2);
  if (src.length() <= max(idx_union)) {
    stop("Error : Incorrect idx, \nmax(idx1, idx2) %i is greater than src.length() %i (sum_union)", 
         max(idx_union), src.length());
  }
  NumericVector temp = src[idx_union];
  return sum(temp);
}

int which(const IntegerVector& vector, const int value) {
  // find idx where vector[idx] = value
  for (int i=0; i<vector.length(); i++) {
    if (vector(i) == value) return i;
  }
  return 0;
} 
int which(const NumericVector& vector, const double value) {
  // find idx where vector[idx] = value
  for (int i=0; i<vector.length(); i++) {
    if (vector(i) == value) return i;
  }
  return 0;
} 

#endif