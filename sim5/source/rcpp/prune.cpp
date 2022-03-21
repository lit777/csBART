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
List prune(
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
  temp                   = dt["Split"]; 
  IntegerVector dt_Split = clone(temp) - 1;
  temp                   = dt["Value"]; 
  IntegerVector dt_Value = clone(temp) - 1;
  NumericVector dt_MU    = dt["MU"];
  temp                   = dt["begin"];
  IntegerVector dt_begin = clone(temp) - 1;
  temp                   = dt["end"];  
  IntegerVector dt_end   = clone(temp) - 1;

  // PRUNE()
  const int P = xpred.ncol();
  IntegerVector col_idx = Rcpp::Range(0, P-1);
  //singly.position <- as.numeric(names(which(table(dt$parent[which(dt$Terminal==1)])==2)))
  IntegerVector singly_position = find_singly_position(dt_Terminal, dt_parent);
  int singly_inode; 
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
  int prop_pred  = dt_Split(subset_ind);
  int begin      = dt_begin(subset_ind);
  int end        = dt_end(subset_ind);
  IntegerVector Obs_begin_end (end - begin + 1);
  //subset_by_range(Obs_begin_end, Obs_list, begin, end, this->id);  // Obs[begin:end]
  Obs_begin_end = Obs_list[Rcpp::Range(begin, end)];
  

  subset_ind  = which(dt_position, 2*singly_inode);
  int begin_L = dt_begin(subset_ind);
  int end_L   = dt_end(subset_ind);
  IntegerVector temp_L (end_L - begin_L + 1);
  //subset_by_range(temp_L, Obs_list, begin_L, end_L, this->id);
  temp_L = Obs_list[Rcpp::Range(begin_L, end_L)];
  
  subset_ind  = which(dt_position, 2*singly_inode+1);
  int begin_R = dt_begin(subset_ind);
  int end_R   = dt_end(subset_ind);
  IntegerVector temp_R (end_R - begin_R + 1);
  //subset_by_range(temp_R, Obs_list, begin_R, end_R, this->id);
  temp_R = Obs_list[Rcpp::Range(begin_R, end_R)];
  
  LogicalVector flag = pred_enough_unique(xpred, Obs_begin_end);
  IntegerVector enough_unique      = col_idx[flag];
  NumericVector enough_unique_prob = prop_prob[flag];
  
  
  double prop_pred_prob = prop_prob(prop_pred) / sum(enough_unique_prob);
  NumericVector xpred_prop_pred(Obs_begin_end.length());
  subset_by_idx(xpred_prop_pred, xpred, Obs_begin_end, prop_pred);
  int unique_len = unique(xpred_prop_pred).length();
  
  // transition ratio
  double TRANS = log(p_grow) - (log(sum(dt_Terminal)- 1)) 
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
    
    int subset_inds = which(dt_parent, singly_inode);
    
    dt_position.erase( subset_inds, 2+subset_inds);
    dt_parent.erase(   subset_inds, 2+subset_inds);
    dt_Terminal.erase( subset_inds, 2+subset_inds);
    dt_Split.erase(    subset_inds, 2+subset_inds);
    dt_Value.erase(    subset_inds, 2+subset_inds);
    dt_MU.erase(       subset_inds, 2+subset_inds);
    dt_begin.erase(    subset_inds, 2+subset_inds);
    dt_end.erase(      subset_inds, 2+subset_inds);
    
    subset_ind = which(dt_position, singly_inode);
    dt_Split(subset_ind) = NA_INTEGER;
    dt_Value(subset_ind) = NA_INTEGER;
    dt_Terminal(subset_ind) = 1;
    
    
    Obs_begin_end.sort();
    for (int i=0; i<Obs_begin_end.length(); i++) Obs_list(begin + i) = Obs_begin_end(i);
    
    
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
