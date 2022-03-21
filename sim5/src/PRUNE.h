//PRUNE.h
#ifndef __PRUNE_H__
#define __PRUNE_H__

#include "decision_tree.h"
#include "functions.h"
#include "subset.h"

using namespace std;
using namespace Rcpp;
using namespace sugar;

// Fun. of PRUNE alteration
void DecisionTree::PRUNE(
    const NumericMatrix& xpred, 
    const NumericVector xcut[],
    double sigma2, double sigma_mu, 
    const NumericVector& R,
    IntegerMatrix& Obs_list, 
    
    double p_prune, double p_grow,
    double alpha, double beta, 
    const NumericVector& prop_prob) {
  
  // PRUNE()
  const int P = xpred.ncol();
  IntegerVector col_idx = Rcpp::Range(0, P-1);
  //singly.position <- as.numeric(names(which(table(dt$parent[which(dt$Terminal==1)])==2)))
  IntegerVector singly_position = this->singly_position();
  int singly_inode; 
  switch(singly_position.length()) {
  case 0: 
    return;
  case 1: 
    singly_inode = singly_position(0);
    break;
  default: 
    singly_inode = sample(singly_position, 1)(0);
  }
  
  int subset_ind = which(this->position, singly_inode);
  int prop_pred  = this->Split(subset_ind);
  int begin      = this->begin(subset_ind);
  int end        = this->end(subset_ind);
  IntegerVector Obs_begin_end (end - begin + 1);
  subset_by_range(Obs_begin_end, Obs_list, begin, end, this->id);  // Obs[begin:end]
  
  subset_ind  = which(this->position, 2*singly_inode);
  int begin_L = this->begin(subset_ind);
  int end_L   = this->end(subset_ind);
  IntegerVector temp_L (end_L - begin_L + 1);
  subset_by_range(temp_L, Obs_list, begin_L, end_L, this->id);
  
  subset_ind  = which(this->position, 2*singly_inode+1);
  int begin_R = this->begin(subset_ind);
  int end_R   = this->end(subset_ind);
  IntegerVector temp_R (end_R - begin_R + 1);
  subset_by_range(temp_R, Obs_list, begin_R, end_R, this->id);
  
  LogicalVector flag = pred_enough_unique(xpred, Obs_begin_end);
  IntegerVector enough_unique      = col_idx[flag];
  NumericVector enough_unique_prob = prop_prob[flag];
  
  double prop_pred_prob = prop_prob(prop_pred) / sum(enough_unique_prob);
  NumericVector xpred_prop_pred(Obs_begin_end.length());
  subset_by_idx(xpred_prop_pred, xpred, Obs_begin_end, prop_pred);
  int unique_len = unique(xpred_prop_pred).length();
  
  // transition ratio
  double TRANS = log(p_grow) - (log(this->terminal_nodes().length() - 1)) 
    + log(max(prop_pred_prob, 0.0)) - log(unique_len)
    - log(p_prune) + log(singly_position.length());
  
  // likelihood ratio
  int nlL = temp_L.length();
  int nlR = temp_R.length();
  
  double sum_R_temp_L = sum_by_idx(R, temp_L);
  double sum_R_temp_R = sum_by_idx(R, temp_R);
  double sum_R_union  = sum_union(R, temp_L, temp_R);
  
  double LH = 
    log(sqrt((sigma2+nlL*sigma_mu) * (sigma2+nlR*sigma_mu))
      / sqrt(sigma2*(sigma2+(nlL+nlR)*sigma_mu)))
    + (sigma_mu / (2*sigma2)
    *(- pow(sum_R_temp_L, 2) / (sigma2 + nlL*sigma_mu)
      - pow(sum_R_temp_R, 2) / (sigma2 + nlR*sigma_mu) 
      + pow(sum_R_union,  2) / (sigma2 + (nlR+nlL) * sigma_mu)));
  
  // structure ratio
  int d = 1;
  while(singly_inode >= pow(2,d)) d += 1;
  d -= 1;
  
  double STR = -log(alpha) - 2*log(1 - alpha/pow(2+d, beta))
    + log(pow(1+d, beta) - alpha)
    - log(max(prop_pred_prob, 0.0)) + log(unique_len);
  
  double r = TRANS + LH + STR;
  
  if (r > log(R::runif(0, 1))) {
    
    int subset_inds = which(this->parent, singly_inode);
    
    this->position.erase( subset_inds, 2+subset_inds);
    this->parent.erase(   subset_inds, 2+subset_inds);
    this->Terminal.erase( subset_inds, 2+subset_inds);
    this->Split.erase(    subset_inds, 2+subset_inds);
    this->Value.erase(    subset_inds, 2+subset_inds);
    this->MU.erase(       subset_inds, 2+subset_inds);
    this->begin.erase(    subset_inds, 2+subset_inds);
    this->end.erase(      subset_inds, 2+subset_inds);
    
    subset_ind = which(this->position, singly_inode);
    this->Split(subset_ind) = NA_INTEGER;
    this->Value(subset_ind) = NA_INTEGER;
    this->Terminal(subset_ind) = 1;
    
    Obs_begin_end.sort();
    for (int i=0; i<Obs_begin_end.length(); i++) Obs_list(begin + i, this->id) = Obs_begin_end(i);
  }
}
#endif