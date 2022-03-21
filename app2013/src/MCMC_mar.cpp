#include <Rcpp.h>

#include "decision_tree.h"
#include "MCMC_utils.h"

using namespace std;
using namespace R;
using namespace Rcpp;
using namespace sugar;

// [[Rcpp::export]]
List MCMC_mar(
    const NumericMatrix& Xpred, 
    const NumericVector& Y_trt, 
    NumericVector& Y_out,
    const double p_grow,    // Prob. of GROW
    const double p_prune,   // Prob. of PRUNE 
    const double p_change,  // Prob. of CHANGE
    const int m,            // Num. of Trees: default setting 100
    const int nu,
    double lambda, double lambda_1, 
    double dir_alpha, double dir_alpha1, 
    double alpha, double beta,
    const int n_iter) {
  
  // Data preparation
  const int P = Xpred.ncol();
  const int n = Y_out.length();
  //const int n1 = n;
  //Y_out = Y_out - mean(Y_out);
  
  NumericMatrix Xpred1 = cbind(Y_trt, Xpred);
  
  // Xcut <- lappy(1:dim(Xpred)[2], function(t) sort(unique(Xpred[,t])))
  NumericVector Xcut[P], Xcut1[P+1];        // e.g. unique value of potential confounders
  for (int j=0; j<P; j++) {     
    NumericVector temp = unique(Xpred(_, j));
    temp.sort();
    Xcut[j] = temp;
  }
  for (int j=0; j<P+1; j++) {
    NumericVector temp = unique(Xpred1(_, j)); 
    temp.sort();
    Xcut1[j] = temp;
  }
  
  NumericVector Z = Rcpp::rnorm(n, R::qnorm(mean(Y_trt), 0, 1, true, false), 1);
  
  // Initial Setup
  // Priors, initial values and hyper-parameters
  NumericVector prob = {p_grow, p_prune, p_change};
  
  // double Sigma2[n_iter];    Sigma2[0]   = 1;
  // double Sigma2_1[n_iter];  Sigma2_1[0] = var(Y_out);
  double sigma2 = 1;
  NumericVector sigma2_1 (n_iter); sigma2_1(0) = var(Y_out);
  
  // sigma_mu based on min/max of Y, Y (A=1) and Y (A=0)
  double sigma_mu   = max(pow(min(Z)     / (-2*sqrt(m)), 2), pow(max(Z)     / (2*sqrt(m)), 2));
  double sigma_mu_1 = max(pow(min(Y_out) / (-2*sqrt(m)), 2), pow(max(Y_out) / (2*sqrt(m)), 2));
  
  // Initial values of R
  NumericVector R  = clone(Z);
  NumericVector R1 = clone(Y_out);
  
  // Initial values for the selection probabilities
  NumericVector post_dir_alpha  = rep(1.0, P);
  NumericVector post_dir_alpha1 = rep(1.0, P+1);
  
  // thin = 10, burn-ins = n_iter/2
  IntegerVector seq = init_seq(n_iter, 10, n_iter/2);
  IntegerVector::iterator seq_pt = seq.begin();
  
  IntegerMatrix ind (1+seq.length(), P+1);
  ind(0, _) = rep(0, P+1);
  int ind_idx = 0;
  
  IntegerMatrix Obs_list  (n, m);           //changed list to matrix
  IntegerMatrix Obs1_list (n, m);
  
  NumericVector Effect (seq.length());
  NumericVector::iterator Effect_pt = Effect.begin();
  
  // Place-holder for the posterior samples
  NumericMatrix Tree   (n,   m);
  NumericMatrix Tree1  (n,   m);
  NumericMatrix Tree11 (2*n, m);
  NumericMatrix Tree00 (2*n, m);
  
  DecisionTree dt_list[m];                  // changed to array of trees
  DecisionTree dt1_list[m];
  for(int t=0; t<m; t++){
    dt_list[t]  = DecisionTree(n, t, true,  false);
    dt1_list[t] = DecisionTree(n, t, false, true);
  }
  
  // Obtaining namespace of MCMCpack package
  Environment MCMCpack = Environment::namespace_env("MCMCpack");
  
  // Picking up rinvgamma() and rdirichlet() function from MCMCpack package
  Function rinvgamma  = MCMCpack["rinvgamma"];
  Function rdirichlet = MCMCpack["rdirichlet"];
  NumericVector prop_prob = rdirichlet(1, rep(dir_alpha, P+1));
  
  NumericMatrix xpred_mult (2*n, P+2);            // changed list to matrix
  IntegerVector xpred_ind = sample(Xpred.nrow(), 2*n, true) - 1;
  for (int j=1; j<P + 1; j++) {                     // first col is rep(0, 2n)
    for(int i=0; i<2*n; i++) {
      xpred_mult(i,j) = Xpred1(xpred_ind(i), j);
    }
  }
  for (int i=0; i<2*n; i++) xpred_mult(i, P+1) = i;
  
  
  
  // Xpred_list skipped
  
  
  //=================================================
  //=================================================
  NumericVector dir_alpha_hist (n_iter); dir_alpha_hist(0) = dir_alpha;
  NumericVector Tree11_hist;
  NumericVector Tree00_hist;
  //=================================================
  //=================================================
  
  ////////////////////////////////////////
  //////////   Run main MCMC    //////////
  ////////////////////////////////////////
  
  // Run MCMC
  for (int iter=1; iter<n_iter; iter++){
    
    //if (iter % 100 == 0) Rcout << "Rcpp iter : " << iter << " of separation model" << endl;
    
    // Ystar = rnorm(rowSums(Tree), 1);
    // Z     = Y_trt * pmax(Ystar, 0) + (1-Y_trt) * pmin(Ystar,0);
    update_Z(Z, Y_trt, Tree);
    
    // ------ Exposure Model
    // Rcout << "Exposure Model " << endl;
    for (int t=0; t<m; t++) {
      // decision trees
      //R = Z - rowSums_without(Tree, t);
      update_R(R, Z, Tree, t);
      
      NumericVector prop_prob_exp = prop_prob[Rcpp::Range(1, P)];
      prop_prob_exp = prop_prob_exp / (sum(prop_prob) - prop_prob(0));
      // create new tree instance
      if (dt_list[t].length() == 1) {     // tree has no node yet
        // grow first step
        dt_list[t].GROW_first(
            Xpred, Xcut, sigma2, sigma_mu, R, Obs_list,
            p_prune, p_grow, alpha, beta, prop_prob_exp
        );
      } else {
        int step = sample(3, 1, false, prob)(0);
        switch (step){
        case 1:   // GROW step
          dt_list[t].GROW(
              Xpred, Xcut, sigma2,sigma_mu, R, Obs_list,
              p_prune, p_grow, alpha, beta, prop_prob_exp
          );
          break;
          
        case 2:   // PRUNE step
          dt_list[t].PRUNE(
              Xpred, Xcut, sigma2, sigma_mu, R, Obs_list,
              p_prune, p_grow, alpha, beta, prop_prob_exp
          );
          break;
          
        case 3:   // CHANGE step
          dt_list[t].CHANGE(
              Xpred, Xcut, sigma2,sigma_mu, R, Obs_list,
              p_prune, p_grow, alpha, beta, prop_prob_exp
          );
          break;
          
        default: {}
        } // end of switch
      } // end of tree instance
      dt_list[t].Mean_Parameter(Tree, sigma2, sigma_mu, R, Obs_list);
    } // end of Exposure Model
    
    // ------ Outcome Model
    for (int t=0; t<m; t++) {
      //R1 = Y_out - rowSums_without(Tree1, t);
      update_R(R1, Y_out, Tree1, t);
      if (dt1_list[t].length() == 1) {     // tree has no node yet
        // grow_first step
        dt1_list[t].GROW_first(
            Xpred1, Xcut1, sigma2_1(iter-1), sigma_mu_1, R1, Obs1_list,
            p_prune, p_grow, alpha, beta, prop_prob
        );
      } else {
        int step = sample(3, 1, false, prob)(0);
        switch (step){
        case 1:   // GROW step
          dt1_list[t].GROW(
              Xpred1, Xcut1, sigma2_1(iter-1), sigma_mu_1, R1, Obs1_list,
              p_prune, p_grow, alpha, beta, prop_prob
          );
          break;
          
        case 2:   // PRUNE step
          dt1_list[t].PRUNE(
              Xpred1, Xcut1, sigma2_1(iter-1), sigma_mu_1, R1, Obs1_list,
              p_prune, p_grow, alpha, beta, prop_prob
          );
          break;
          
        case 3:   // CHANGE step
          dt1_list[t].CHANGE(
              Xpred1, Xcut1, sigma2_1(iter-1), sigma_mu_1, R1, Obs1_list,
              p_prune, p_grow, alpha, beta, prop_prob
          );
          break;
          
        default: {}
        } // end of switch
      } // end of tree instance
      dt1_list[t].Mean_Parameter(Tree1, sigma2_1(iter-1), sigma_mu_1, R1, Obs1_list);
    } // end of Outcome Model (A=1)
    
    // Sample variance parameter
    {
      //sigma2 = 1;
      NumericVector sigma2_1_temp = rinvgamma(1, nu/2+n/2, nu*lambda_1/2 + sum(pow(Y_out-rowSums(Tree1),2))/2);
      sigma2_1(iter) = sigma2_1_temp(0);
    }
    
    // Num. of inclusion of each potential confounder
    NumericVector add(P), add1(P+1);
    for (int t=0; t<m; t++) {
      add  += dt_list[t].num_included(P);
      add1 += dt1_list[t].num_included(P+1);
    }
    
    // Rcout << "M.H. algorithm for alpha" << endl;
    // M.H. algorithm for the alpha parameter in the dirichlet distribution
    {
      double p_dir_alpha, dir_lik_p, dir_lik, ratio;
      
      p_dir_alpha = max(rnorm(dir_alpha, 0.1), pow(0.1, 10));
      
      NumericVector SumS(P+1);
      log_with_LB(SumS, prop_prob);
      
      dir_lik_p = 
        sum(SumS*(   rep(p_dir_alpha/(P+1), (P+1))-1))
        + lgamma(sum(rep(p_dir_alpha/(P+1), (P+1))))
        - sum(lgamma(rep(p_dir_alpha/(P+1), (P+1))));
        
      dir_lik = 
        sum(SumS*(   rep(dir_alpha/(P+1), (P+1))-1))
        + lgamma(sum(rep(dir_alpha/(P+1), (P+1))))
        - sum(lgamma(rep(dir_alpha/(P+1), (P+1))));
      
      ratio = dir_lik_p
        + log(pow(p_dir_alpha/(p_dir_alpha + (P+1)),0.5-1)
          * abs(1/(p_dir_alpha+(P+1))
          - p_dir_alpha/pow(p_dir_alpha+(P+1),2)))
        + dnorm(dir_alpha, p_dir_alpha, 0.1, true)
        - dir_lik
        - log(pow(dir_alpha/(dir_alpha+(P+1)), 0.5-1)
          * abs(1/(dir_alpha+(P+1))
          - dir_alpha/pow(dir_alpha+(P+1),2)))
        - dnorm(p_dir_alpha, dir_alpha, 0.1, true);
                    
      if (ratio > log(R::runif(0,1))) {
        dir_alpha = p_dir_alpha;
      }
      dir_alpha_hist(iter) = dir_alpha;
      post_dir_alpha = rep(dir_alpha/(P+1), P+1);
    }
    
    // M.H. algorithm for the inclusion probabilities
    {
      double dir_lik_p, dir_lik, ratio;
      
      NumericVector add_max(add.length() + 1);
      add_max(0) = add1(0)/sum(add1[Rcpp::Range(1, add1.length())])*sum(add);
      add_max[Rcpp::Range(1, add.length())] = add;
      
      NumericVector p_prop_prob = rdirichlet(1, add1 + post_dir_alpha + add_max);
      NumericVector log_p_prop_prob(P+1), log_prop_prob(P+1);
      log_with_LB(log_p_prop_prob, p_prop_prob);
      log_with_LB(log_prop_prob,   prop_prob);
      
      dir_lik_p = sum(add) * log(1/(1-p_prop_prob(0)))
        + (add_max(0) + post_dir_alpha(0) - 1.0) * log_p_prop_prob(0)
        + sum((add1[Rcpp::Range(1,P)] + add + post_dir_alpha[Rcpp::Range(1,P)] - 1.0)
        * log_p_prop_prob[Rcpp::Range(1,P)]);
        
      dir_lik = sum(add) * log(1/(1-prop_prob(0)))
        +(add_max(0) + post_dir_alpha(0) - 1.0) * log_prop_prob(0)
        + sum((add1[Rcpp::Range(1,P)] + add + post_dir_alpha[Rcpp::Range(1,P)] - 1.0)
        * log_prop_prob[Rcpp::Range(1,P)]);
          
      ratio = dir_lik_p
        + sum((add1 + post_dir_alpha + add_max - 1.0) * log_prop_prob)
        - dir_lik   
        - sum((add1 + post_dir_alpha + add_max - 1.0) * log_p_prop_prob);
            
      if (ratio > log(R::runif(0,1))) {
        prop_prob = clone(p_prop_prob);
      }
    }
    
    // Sampling E[Y(1)-Y(0)]
    if(iter == *seq_pt) {
      for (int i=0; i<m; i++) {
        // Rcout << "\n\n" << endl;
        // dt1_list[i].print();
        // Rcout << "\nPred with Tree11\n\n";
        dt1_list[i].Mean_Parameter_pred(Tree11, Xcut1, xpred_mult, n, 1);
        
        dt1_list[i].Mean_Parameter_pred(Tree00, Xcut1, xpred_mult, n, 0);
        // return List::create(
        //   Named("dt") = dt1_list[i].export_tree(),
        //   Named("Tree11") = Tree11,
        //   Named("Tree00") = Tree00,
        //   Named("xpred_mult") = xpred_mult
        // );
      }
      *Effect_pt = mean(rowSums(Tree11) - rowSums(Tree00));
      Tree11_hist.push_back(mean(rowSums(Tree11)));
      Tree00_hist.push_back(mean(rowSums(Tree00)));
      Effect_pt++;
      IntegerVector ind_temp = ifelse(add1 > 0.0, 1, 0);
      ind(ind_idx + 1, _) = ind_temp; // indicator of whether confounders are included
      seq_pt++; ind_idx++;
    }
    
    Rcpp::checkUserInterrupt(); // check for break in R
    
  } // end of MCMC iterations
  
  List L = List::create(
    Named("Effect") = Effect, 
    Named("ind") = ind,
    Named("prop_prob") = prop_prob,
    Named("sigma2_1") = sigma2_1,
    Named("dir_alpha_hist") = dir_alpha_hist,
    Named("Tree") = Tree,
    Named("Tree1") = Tree1,
    Named("Tree11") = Tree11,
    Named("Tree00") = Tree00,
    Named("Tree11_hist") = Tree11_hist,
    Named("Tree00_hist") = Tree00_hist
  );
  
  return L;
}


