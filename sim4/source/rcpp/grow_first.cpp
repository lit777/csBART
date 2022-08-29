#include <Rcpp.h>
using namespace Rcpp;
using namespace std;

// [[Rcpp::export]]
List grow_first(
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
  // GROW FIRST
  const int n = xpred.nrow();
  const int P = xpred.ncol();
  const IntegerVector row_idx = Rcpp::Range(0, n-1);
  int prop_pred, prop_rule;
  double value;
  prop_pred = sample(P, 1, false, prop_prob)(0) - 1;             // pick a predictor and match index with cpp
  NumericVector xcut_prop_pred = xcut[prop_pred];
  if (ind == 1 && prop_pred == 1) {
    prop_rule = 2;
    value = 1;
  } else {
    prop_rule = sample(xcut_prop_pred.length()-1, 1)(0);   // sample U(2, length(xcut(prop_pred)))
    value = xcut_prop_pred(prop_rule);                  // value for separation 
  }
  
  IntegerVector R_L = row_idx[xpred(_, prop_pred) <  value];
  IntegerVector R_R = row_idx[xpred(_, prop_pred) >= value];
  
  NumericVector xpred_prop_pred = xpred(_, prop_pred);
  
  // Transition ratio (log scale)
  double TRANS = log(p_prune) - log(max(prop_prob(prop_pred), 0.0))
    + log(xcut_prop_pred.length()-1) - log(p_grow);
  
  // Likelihood ratio (log scale)
  int nlL = R_L.length();
  int nlR = R_R.length();
  NumericVector temp;
  temp = R[R_L];
  double sum_R_L = sum(temp); //sum_by_idx(R, R_L); 
  temp = R[R_R];
  double sum_R_R = sum(temp);//sum_by_idx(R, R_R);
  
  double LH = 
    log(sqrt(sigma2 * (sigma2+(nlL+nlR) * sigma_mu)) 
          / sqrt( (sigma2 + nlL*sigma_mu)*(sigma2 + nlR*sigma_mu)))
    + (sigma_mu / (2*sigma2)
    * (pow(sum_R_L,2) / (sigma2 + sigma_mu* nlL)
    +  pow(sum_R_R,2) / (sigma2 + sigma_mu* nlR)
    -  pow(sum(R), 2) / (sigma2 + sigma_mu*(nlR+nlL))));
    
    // Structure ratio (log scale)
    int d = 0;
    double STR = log(alpha) + 2*log((1-alpha / pow(2+d, beta)))
      - log(pow(1+d, beta) - alpha)
      + log(max(prop_prob(prop_pred), 0.0))
      - log(xcut_prop_pred.length()-1);
      
    double r = TRANS + LH + STR;
    
    if (r > log(R::runif(0, 1))) {
      // New tree structure
      IntegerVector position = dt["position"];
      IntegerVector parent   = dt["parent"];
      IntegerVector Terminal = dt["Terminal"];
      IntegerVector Split    = dt["Split"];
      IntegerVector Value    = dt["Value"];
      NumericVector MU       = dt["MU"];
      IntegerVector temp;
      temp = dt["begin"];
      IntegerVector begin    = clone(temp) -1;
      temp = dt["end"];
      IntegerVector end      = clone(temp) -1;
      
      Split    = prop_pred;
      Terminal = rep(0, 1);
      Value    = prop_rule;
      position.push_back(2);       position.push_back(3);
      parent.push_back(1);         parent.push_back(1);
      Terminal.push_back(1);       Terminal.push_back(1);
      Split.push_back(NA_INTEGER); Split.push_back(NA_INTEGER);
      Value.push_back(NA_INTEGER); Value.push_back(NA_INTEGER);
      MU.push_back(NA_REAL);       MU.push_back(NA_REAL);
      begin.push_back(0);          begin.push_back(nlL);
      end.push_back(nlL-1);        end.push_back(n-1);
      
      List dt_new = List::create(
        Named("position") = position,
        Named("parent") = parent,
        Named("Terminal") = Terminal,
        Named("Split") = Split +1,
        Named("Value") = Value +1,
        Named("MU") = MU,
        Named("begin") = begin +1,
        Named("end") = end+1
      );
      
      // update obs
      R_L.sort(); R_R.sort();
      for (int i=0; i<nlL; i++) Obs_list(i)       = R_L(i);
      for (int i=0; i<nlR; i++) Obs_list(i + nlL) = R_R(i);
      
      return List::create(
        Named("dt") = dt_new,
        Named("Obs") = Obs_list
      );
    }
  return List::create(
    Named("dt") = dt,
    Named("Obs") = Obs_list
  );
} // end of GROW_first()
