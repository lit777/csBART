#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
double update_sigma(int n, int nu, int lambda_1, NumericVector Y_out, NumericMatrix Tree1) {
  // Obtaining namespace of MCMCpack package
  Environment MCMCpack = Environment::namespace_env("MCMCpack");
  
  // Picking up rinvgamma() and rdirichlet() function from MCMCpack package
  Function rinvgamma  = MCMCpack["rinvgamma"];
  
  double sigma2_1;
  NumericVector sigma2_1_temp = rinvgamma(1, nu/2+n/2, nu*lambda_1/2 + sum(pow(Y_out-rowSums(Tree1),2))/2);
  sigma2_1 = sigma2_1_temp(0);
  
  return sigma2_1;
}
