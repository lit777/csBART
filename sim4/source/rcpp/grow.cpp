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
List grow(
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
  
  // GROW()
  const int P = xpred.ncol();
  IntegerVector col_idx = Rcpp::Range(0, P-1);
  
  // which(dt$Terminal==1)
  IntegerVector tnode_idx = terminal_nodes(dt_Terminal);
  int prop_tnode = sample(tnode_idx, 1)(0); // randomly choose terminal node
  
  // find begin and end for chosen terminal node
  int begin = dt_begin(prop_tnode);
  int end   = dt_end(prop_tnode);
  if (end - begin + 1 < 2) {
    // return the current tree if there is no covariate with enough unique values
    return List::create(
      Named("dt") = dt,
      Named("Obs") = Obs_list
    );
  }
  
  IntegerVector Obs_ind (end - begin + 1);
  if (end == begin) 
    Obs_ind = Obs_list(begin);
  else
    Obs_ind = Obs_list[Rcpp::Range(begin, end)];
  
  // find covariate with enough unique values
  LogicalVector flag = pred_enough_unique(xpred, Obs_ind);
  IntegerVector enough_unique      = col_idx[flag];
  NumericVector enough_unique_prob = prop_prob[flag];
  
  // pick a predictor -> already index
  int prop_pred = sample(enough_unique, 1, false, enough_unique_prob)(0); 
  NumericVector xpred_prop_pred(Obs_ind.length());
  subset_by_idx(xpred_prop_pred, xpred, Obs_ind, prop_pred);
  
  NumericVector unique_xpred = unique(xpred_prop_pred);
  
  int unique_len = unique_xpred.length();    // number of unique values
  if (unique_len == 1) {
    return List::create(
      Named("dt") = dt,
      Named("Obs") = Obs_list
    );
  }
  
  double min_xpred = min(xpred_prop_pred);
  double value;
  do {
    value = sample(unique_xpred, 1)(0);
  } while(value == min_xpred);
  // double value;
  // {
  //   NumericVector temp = clone(unique(xpred_prop_pred));
  //   temp.sort();
  //   int idx = sample(unique_len-1, 1)(0);
  //   value = temp(idx);
  // }
  
  // sample from  unique values of xpred_prop_pred and it should not be minimum value
  
  
  int prop_rule = which(xcut[prop_pred], value);         // index of xcut[prop_pred] with value
  
  IntegerVector R_L = Obs_ind[xpred_prop_pred <  value];
  IntegerVector R_R = Obs_ind[xpred_prop_pred >= value];
  int nlL = R_L.length(); 
  int nlR = R_R.length();
  
  
  // Temporary new tree structure
  IntegerVector dt_new_position = clone(dt_position);
  IntegerVector dt_new_parent   = clone(dt_parent);
  IntegerVector dt_new_Terminal = clone(dt_Terminal);
  IntegerVector dt_new_Split    = clone(dt_Split);
  IntegerVector dt_new_Value    = clone(dt_Value);
  NumericVector dt_new_MU       = clone(dt_MU);
  IntegerVector dt_new_begin    = clone(dt_begin);
  IntegerVector dt_new_end      = clone(dt_end);
  
  int new_position = dt_new_position(prop_tnode);
  dt_new_Terminal(prop_tnode) = 0;
  
  dt_new_position = append(dt_new_position, 2*new_position, 2*new_position + 1);
  dt_new_parent   = append(dt_new_parent,   new_position,   new_position);
  dt_new_Terminal = append(dt_new_Terminal, 1,              1);
  dt_new_Split    = append(dt_new_Split,    NA_INTEGER,     NA_INTEGER);
  dt_new_Value    = append(dt_new_Value,    NA_INTEGER,     NA_INTEGER);
  dt_new_MU       = append(dt_new_MU,       NA_REAL,        NA_REAL);
  dt_new_begin    = append(dt_new_begin,    begin,          begin+nlL);
  dt_new_end      = append(dt_new_end,      begin+nlL-1,    end);
  
  
  // Prune Step
  // Find internal nodes with two terminal child nodes
  // equivalent expression in R is  
  // count <- length(which(table(dt.new$parent[which(dt.new$Terminal==1)])==2)) 
  IntegerVector tnode_idx_new = terminal_nodes(dt_new_Terminal);
  int count = terminal_parents(tnode_idx_new, dt_new_parent);
  
  // Transition ratio
  double sum_prop_prob_enough_unique = sum(enough_unique_prob);
  double sum_R_L = sum_by_idx(R, R_L);
  double sum_R_R = sum_by_idx(R, R_R);
  double sum_R_union = sum_union(R, R_L, R_R);
  
  double TRANS = log(p_prune) + log(tnode_idx.length())
    - log(max(prop_prob(prop_pred)/sum_prop_prob_enough_unique, 0.0))
    + log(unique_len-1) - log(p_grow) - log(count);
    
  double LH = 
  log(sqrt( sigma2*(sigma2+(nlL+nlR)*sigma_mu))
        / sqrt((sigma2+nlL*sigma_mu)*(sigma2+nlR*sigma_mu)))
    + (sigma_mu / (2*sigma2) 
    * (pow(sum_R_L, 2)     / (sigma2 + sigma_mu*nlL)
    +  pow(sum_R_R, 2)     / (sigma2 + sigma_mu*nlR)
    -  pow(sum_R_union, 2) / (sigma2 + sigma_mu*(nlR+nlL))));
    
  // Structure ratio
  int d = 1;
  while (dt_position(prop_tnode) >= pow(2,d)) {
    d += 1;
  }
  d -= 1;
  double STR = log(alpha) + 2*log(1 - alpha/pow(2+d, beta))
    - log(pow(1+d, beta) - alpha)
    + log(max(prop_prob(prop_pred)/sum_prop_prob_enough_unique, 0.0))
    - log(unique_len-1);
    
  double r = TRANS + LH + STR;
  
  if (r > log(R::runif(0, 1))) {
    dt_new_Split(prop_tnode) = prop_pred;
    dt_new_Value(prop_tnode) = prop_rule;
    
    List dt_new = List::create(
      Named("position") = dt_new_position,
      Named("parent")   = dt_new_parent,
      Named("Terminal") = dt_new_Terminal,
      Named("Split")    = dt_new_Split+1,
      Named("Value")    = dt_new_Value+1,
      Named("MU")       = dt_new_MU,
      Named("begin")    = dt_new_begin+1,
      Named("end")      = dt_new_end+1
    );
    
    // update obs
    R_L.sort(); R_R.sort();
    for (int i=0; i<nlL; i++) Obs_list(i + begin) = R_L(i);
    for (int i=0; i<nlR; i++) Obs_list(i + begin + nlL) = R_R(i);
    
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