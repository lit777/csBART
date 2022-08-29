#include <Rcpp.h>
using namespace Rcpp;
using namespace std;

#include "subset.h"
#include "functions.h"

IntegerVector terminal_nodes(IntegerVector Terminal) {
  // method to find idx of terminal nodes
  IntegerVector tnode_idx = Rcpp::Range(0, Terminal.length()-1);
  tnode_idx = tnode_idx[Terminal == 1];
  return tnode_idx;
};

IntegerVector find_singly_position(IntegerVector Terminal, IntegerVector parent) {
  // return index of singly positions
  IntegerVector singly_position;
  map<int, int> table;
  for (int i=0; i<Terminal.length(); i++) {  // equivalent to table() in R
    if (Terminal(i)==1) {
      table[parent(i)]++;                     
    }
  }
  for(map<int, int>::iterator i=table.begin(); i!=table.end(); i++){
    if (i->second == 2) {                       // count value with 2
      singly_position.push_back(i->first);      // i->first : key, i->second : value
    }
  }
  return singly_position;
}

int terminal_parents(IntegerVector tnode_idx, IntegerVector parent) {
  // method to count internal nodes with two terminal child nodes
  int count = 0; 
  map<int, int> table;                            // use map to count parent  
  for (int i = 0; i < tnode_idx.length(); i++) {  // equivalent to table() in R
    table[parent[tnode_idx(i)]]++;                     
  }
  for (map<int, int>::iterator i=table.begin(); i!=table.end(); i++) {
    if (i->second == 2) count += 1;               // i->first : key, i->second : value
  }                                               // count value with 2
  return count;
};


// [[Rcpp::export]]
List change(
    const NumericMatrix& xpred, 
    const List xcut,
    double sigma2, double sigma_mu, 
    const NumericVector& R, 
    IntegerVector& Obs_list, 
    
    double p_prune, double p_grow,
    double alpha, double beta, 
    const NumericVector& prop_prob,
    List dt,
    int ind
) {
  IntegerVector temp;
  IntegerVector dt_position = dt["position"];
  IntegerVector dt_parent   = dt["parent"];
  IntegerVector dt_Terminal = dt["Terminal"];
  temp                      = dt["Split"]; 
  IntegerVector dt_Split = clone(temp) - 1;
  temp    = dt["Value"]; 
  IntegerVector dt_Value = clone(temp) - 1;
  NumericVector dt_MU       = dt["MU"];
  temp    = dt["begin"]; 
  IntegerVector dt_begin = clone(temp) - 1;
  temp      = dt["end"];   
  IntegerVector dt_end   = clone(temp) - 1;
  
  // CHANGE()
  const int P = xpred.ncol();
  IntegerVector col_idx = Rcpp::Range(0, P-1);
  
  // find nodes with two terminal child nodes (singly internal parent nodes)
  IntegerVector singly_position = find_singly_position(dt_Terminal, dt_parent);
  int singly_inode;     // pick a internal parent node
  switch(singly_position.length()) {
  case 0: 
    return List::create(
      Named("dt") = dt,
      Named("Obs") = Obs_list
  );
  case 1: 
    singly_inode = singly_position(0);
    break;
  default: 
    singly_inode = sample(singly_position, 1)(0);
  }
  
  
  int subset_ind = which(dt_position, singly_inode);
  int begin      = dt_begin(subset_ind);
  int end        = dt_end(subset_ind);
  IntegerVector Obs_ind(end-begin+1);
  //subset_by_range(Obs_ind, Obs_list, begin, end, this->id);
  Obs_ind = Obs_list[Rcpp::Range(begin, end)];
  
  LogicalVector flag = pred_enough_unique(xpred, Obs_ind);
  IntegerVector enough_unique      = col_idx[flag];
  NumericVector enough_unique_prob = prop_prob[flag];
  
  // pick a predictor -> already index
  int prop_pred = sample(enough_unique, 1, false, enough_unique_prob)(0); 
  
  NumericVector xpred_prop_pred (Obs_ind.length());
  subset_by_idx(xpred_prop_pred, xpred, Obs_ind, prop_pred);
  
  NumericVector unique_xpred = unique(xpred_prop_pred);
  int unique_len = unique_xpred.length();    // number of unique values
  if (unique_len == 1) {
    return List::create(
      Named("dt") = dt,
      Named("Obs") = Obs_list
    );
  }
  
  // double value;
  // {
  //   NumericVector temp = clone(unique(xpred_prop_pred));
  //   temp.sort();
  //   int idx = sample(unique_len-1, 1)(0);
  //   value = temp(idx);
  // }
  // int prop_rule = which(xcut[prop_pred], value);
  
  double min_xpred = min(xpred_prop_pred);
  double value;
  do {
    value = sample(unique_xpred, 1)(0);
  } while(value == min_xpred);
  int prop_rule = which(xcut[prop_pred], value);
  
  IntegerVector RL_star = Obs_ind[xpred_prop_pred <  value];
  IntegerVector RR_star = Obs_ind[xpred_prop_pred >= value];
  
  int begin_idx, end_idx;
  begin_idx = dt_begin[which(dt_position, 2*singly_inode)];
  end_idx   = dt_end[  which(dt_position, 2*singly_inode)];
  IntegerVector RL (end_idx - begin_idx + 1);
  //subset_by_range(RL, Obs_list, begin_idx, end_idx, this->id);
  RL = Obs_list[Rcpp::Range(begin_idx, end_idx)];
  
  begin_idx = dt_begin[which(dt_position, 2*singly_inode + 1)];
  end_idx   = dt_end[  which(dt_position, 2*singly_inode + 1)];
  IntegerVector RR (end_idx - begin_idx + 1);
  //subset_by_range(RR, Obs_list, begin_idx, end_idx, this->id);
  RR = Obs_list[Rcpp::Range(begin_idx, end_idx)];
  
  // Likelihood ratio
  int nlL = RL.length();
  int nlR = RR.length();
  int nlL_star = RL_star.length();
  int nlR_star = RR_star.length();
  double sum_R_RL_star = sum_by_idx(R, RL_star);
  double sum_R_RR_star = sum_by_idx(R, RR_star);
  double sum_R_RL = sum_by_idx(R, RL);
  double sum_R_RR = sum_by_idx(R, RR);
  
  
  double LH = log(sqrt((sigma2/sigma_mu + nlL)     *(sigma2/sigma_mu + nlR)))
    - log(sqrt((sigma2/sigma_mu + nlL_star)*(sigma2/sigma_mu + nlR_star)))
    + (0.5 / sigma2
    *(pow(sum_R_RL_star, 2) / (sigma2/sigma_mu + nlL_star)
    + pow(sum_R_RR_star, 2) / (sigma2/sigma_mu + nlR_star)
    - pow(sum_R_RL, 2)      / (sigma2/sigma_mu + nlL) 
    - pow(sum_R_RR, 2)      / (sigma2/sigma_mu + nlR) ));
    
  double r = LH;
  
  if(r > log(R::runif(0, 1))) {
    // New tree structure
    subset_ind = which(dt_position, singly_inode);
    dt_Terminal(subset_ind) = 0;
    dt_Split(subset_ind) = prop_pred;
    dt_Value(subset_ind) = prop_rule;
    
    subset_ind = which(dt_position, 2*singly_inode);
    dt_begin(subset_ind) = begin;
    dt_end(subset_ind)   = begin + nlL_star - 1;
    
    subset_ind = which(dt_position, 2*singly_inode + 1);
    dt_begin(subset_ind) = begin + nlL_star;
    dt_end(subset_ind)   = end;
    
    RL_star.sort(); RR_star.sort();
    for (int i=0; i<nlL_star; i++) Obs_list(i + begin) = RL_star(i);
    for (int i=0; i<nlR_star; i++) Obs_list(i + begin + nlL_star) = RR_star(i);
    
    List dt_new = List::create(
      Named("position") = dt_position,
      Named("parent")   = dt_parent,
      Named("Terminal") = dt_Terminal,
      Named("Split")    = dt_Split+1,
      Named("Value")    = dt_Value+1,
      Named("MU")       = dt_MU,
      Named("begin")    = dt_begin+1,
      Named("end")      = dt_end+1
    );
    
    return List::create(
      Named("dt") = dt_new,
      Named("Obs") = Obs_list
    );
  }
  return List::create(
    Named("dt") = dt,
    Named("Obs") = Obs_list
  );
}