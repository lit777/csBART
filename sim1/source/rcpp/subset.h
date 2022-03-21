//subset.h
#ifndef __subset_h__
#define __subset_h__

#include <Rcpp.h>

using namespace Rcpp;
using namespace std;
using namespace sugar;

// this header includes functions related to subsetting

void subset_by_idx(IntegerVector& res, const IntegerMatrix& src, const IntegerVector& idx, const int& col) {
  if (res.length() != idx.length()) {
    stop("Error : dimension does not match. \n  res length is %i while src length is %i (subset_by_idx)", 
         res.length(), src.length());
  }
  for (int i=0; i<idx.length(); i++) res(i) = src(idx(i), col);
};

void subset_by_idx(NumericVector& res, const NumericMatrix& src, const IntegerVector& idx, const int& col) {
  if (res.length() != idx.length()) {
    stop("Error : dimension does not match. \n  res length is %i while src length is %i (subset_by_idx)", 
         res.length(), src.length());
  }
  for (int i=0; i<idx.length(); i++) res(i) = src(idx(i), col);
};

void subset_by_range(IntegerVector& res, const IntegerMatrix& src, const int start, const int end, const int col) {
  if (res.length() != end-start+1) {
    stop("Error : dimension does not match. \n  res length is %i while range is %i (subset_by_range)", 
         res.length(), end-start+1);
  }
  if (end < start) stop("Error : start is greater than end. \nstart is %i while end is %i (subset by range)", start, end);
  for (int i=0; i<end-start+1; i++) res(i) = src(start+i, col);
}

void subset_by_range(NumericVector& res, const NumericMatrix& src, const int start, const int end, const int col) {
  if (res.length() != end-start+1) {
    stop("Error : dimension does not match. \n  res length is %i while range is %i (subset_by_range)", 
         res.length(), end-start+1);
  }
  if (end<start) stop("Error : start is greater than end. \n  start is %i while end is %i (subset by range)", start, end);
  for (int i=0; i<end-start+1; i++) res(i) = src(start+i, col);
}

NumericMatrix subset_by_row(const NumericMatrix& src, const LogicalVector& idx) {
  if(sum(idx) < 0) Rcout << idx << endl;
  NumericMatrix res(sum(idx), src.ncol());
  if (res.nrow() != sum(idx)) {
    stop("Error : dimension does not match. \n  res length is %i while sum(idx) is %i (subset_by_row)", 
              res.length(), sum(idx));
  }
  IntegerVector idx_int = Rcpp::Range(0, idx.length()-1);
  idx_int = idx_int[idx];
  
  for (int i=0; i<idx_int.length(); i++) res(i, _) = src(idx_int(i), _);
  return res;
}
#endif