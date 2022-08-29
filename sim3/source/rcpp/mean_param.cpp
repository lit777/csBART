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

// [[Rcpp::export]]
List Mean_Parameter(
    double sigma2,
    double sigma_mu, 
    const NumericVector& R, 
    const IntegerVector& Obs_list, 
    List dt,
    int n) {
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
  
  
  NumericVector T(n);
  IntegerVector tnode_idx = terminal_nodes(dt_Terminal);
  for (int i=0; i<tnode_idx.length(); i++) {
    int begin = dt_begin(tnode_idx(i));       // match index
    int end   = dt_end(tnode_idx(i));
    
    IntegerVector Obs_ind(end-begin+1);
    if (begin == end) 
      Obs_ind = Obs_list(begin);
    else 
      Obs_ind = Obs_list[Rcpp::Range(begin, end)];
      //subset_by_range(Obs_ind, Obs_list, begin, end, this->id); // Obs_ind = Obs_list[begin:end, ]
    
    
    //double Var  = 1 / (1/sigma_mu + Obs_ind.length()/sigma2);
    //double Mean = Var * (sum_by_idx(R, Obs_ind) / sigma2);
    
    double logsd = - 0.5 * log(1/sigma_mu + Obs_ind.length()/sigma2);
    double Mean = exp(2*logsd) * (sum_by_idx(R, Obs_ind) / sigma2);
    
    //this->MU(tnode_idx(i)) = R::rnorm(Mean, sqrt(Var));
    dt_MU(tnode_idx(i)) = R::rnorm(Mean, exp(logsd));
    T[Obs_ind] = dt_MU(tnode_idx(i));
  }
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
    Named("T") = T
  );
}