// decision_tree.h
#ifndef __DECISION_TREE_H__
#define __DECISION_TREE_H__

#include <Rcpp.h>
#include <Rcpp/Benchmark/Timer.h>

//[[Rcpp::plugins(cpp11)]]

using namespace Rcpp;
using namespace std;
using namespace sugar;

class DecisionTree {
  
private:
  IntegerVector position = rep(1,1);
  IntegerVector parent   = rep(NA_INTEGER, 1);
  IntegerVector Terminal = rep(1,1);
  IntegerVector Split    = rep(NA_INTEGER, 1);
  IntegerVector Value    = rep(NA_INTEGER, 1);
  NumericVector MU       = rep(NA_REAL, 1);
  IntegerVector begin    = rep(0, 1);
  IntegerVector end      = rep(0, 1);
  int n;
  int id;
  bool mar_exp;
  bool mar_out;
  
public:
  DecisionTree() {};     // default constructor
  DecisionTree(int n, int id, bool mar_exp=false, bool mar_out=false) {
    this->n = n;
    this->id = id;
    this->end = rep(n-1, 1);
    this->mar_exp = mar_exp;
    this->mar_out = mar_out;
  }
  void copy_from(const DecisionTree &dt) {
    this->position = clone(dt.position);
    this->parent   = clone(dt.parent);
    this->Terminal = clone(dt.Terminal);
    this->Split    = clone(dt.Split);
    this->Value    = clone(dt.Value);
    this->MU       = clone(dt.MU);
    this->begin    = clone(dt.begin);
    this->end      = clone(dt.end);
    this->n        = dt.n;
    this->id       = dt.id;
    this->mar_exp  = dt.mar_exp;
    this->mar_out  = dt.mar_out;
  }
  
  void copy_to(DecisionTree &dt) {
    dt.position = clone(this->position);
    dt.parent   = clone(this->parent);
    dt.Terminal = clone(this->Terminal);
    dt.Split    = clone(this->Split);
    dt.Value    = clone(this->Value);
    dt.MU       = clone(this->MU);
    dt.begin    = clone(this->begin);
    dt.end      = clone(this->end);
    dt.n        = this->n;
    dt.id       = this->id;
    dt.mar_exp  = this->mar_exp;
    dt.mar_out  = this->mar_out;
  }
  
  DecisionTree(int n, int id, bool mar_exp, bool mar_out, const List dt) {
    this->n = n;
    this->id = id;
    this->end = rep(n-1, 1);
    this->mar_exp = mar_exp;
    this->mar_out = mar_out;
    IntegerVector temp = dt[0];
    this->position = temp;
    temp = dt[1];
    this->parent   = temp;
    temp = dt[2];
    this->Terminal = temp;
    temp = dt[3];
    this->Split    = temp - 1;
    temp = dt[4];
    this->Value    = temp - 1;
    temp = dt[6];
    this->begin    = temp - 1;
    temp = dt[7];
    this->end      = temp - 1;
    NumericVector temp2 = dt[5];
    this->MU = temp2;
  }
  
  List export_tree() {
    List L = List::create(
      Named("position") = this->position,
      Named("parent")   = this->parent,
      Named("Terminal") = this->Terminal,
      Named("Split")    = this->Split + 1,
      Named("Value")    = this->Value + 1,
      Named("MU")       = this->MU,
      Named("begin")    = this->begin + 1,
      Named("end")      = this->end   + 1
    );
    
    return L;
  }
  
  inline int  length() {return this->position.length();};
  inline bool is_mar() {return this->mar_exp | this->mar_out;};
  inline bool is_mar_exposure() {return this->mar_exp;};
  inline bool is_mar_outcome()  {return this->mar_out;};
  
  NumericVector num_included(const int P);
  void get_hist(IntegerVector Split_hist[], IntegerVector Value_hist[], IntegerVector Side_hist[], bool verbose=false);
  IntegerVector singly_position();
  IntegerVector terminal_nodes();
  int terminal_parents(IntegerVector tnode_idx);
  void update_xpred_temp(
      NumericMatrix& xpred_temp, 
      const NumericVector xcut[],
      const IntegerVector Split_hist[], 
      const IntegerVector Value_hist[], 
      const IntegerVector Side_hist[], 
      const int i
  );
  void update_tree(
      NumericMatrix& Tree,
      const NumericVector xcut[],
      const NumericMatrix& xpred_mult, 
      const IntegerVector Split_hist[], 
      const IntegerVector Value_hist[], 
      const IntegerVector Side_hist[]
  );
  void print(){
    Rcout << "position : "; print(this->position); Rcout << endl;
    Rcout << "parent   : "; print(this->parent);   Rcout << endl;
    Rcout << "Terminal : "; print(this->Terminal); Rcout << endl;
    Rcout << "Split    : "; print(this->Split);    Rcout << endl;
    Rcout << "Value    : "; print(this->Value);    Rcout << endl;
    Rcout << "begin    : "; print(this->begin);    Rcout << endl;
    Rcout << "end      : "; print(this->end);      Rcout<< endl << endl;
  }
  void print(IntegerVector x) {
    for (int i=0; i<x.length(); i++) {
      if (x(i) == NA_INTEGER) Rcout << "NA ";
      else if (x(i) < 10) Rcout << " " << x(i) << " ";
      else Rcout << x(i) << " ";
    }
  }
  void print_MU() {
    for(int i=0; i<this->Terminal.length(); i++) {
      if(this->Terminal(i) == 0) Rcout << this->MU(i) << endl;
    }
  }
  
  void GROW_first(
      const NumericMatrix& xpred, 
      const NumericVector xcut[],
      double sigma2, double sigma_mu, 
      const NumericVector& R, 
      IntegerMatrix& Obs_list, 
      
      double p_prune, double p_grow,
      double alpha, double beta, 
      const NumericVector& prop_prob
  ); 
  
  void GROW(
      const NumericMatrix& xpred, 
      const NumericVector xcut[],
      double sigma2, double sigma_mu, 
      const NumericVector& R, 
      IntegerMatrix& Obs_list, 
      
      double p_prune, double p_grow,
      double alpha, double beta, 
      const NumericVector& prop_prob
  ); 
  
  void PRUNE(
      const NumericMatrix& xpred, 
      const NumericVector xcut[],
      double sigma2, double sigma_mu, 
      const NumericVector& R,
      IntegerMatrix& Obs_list, 
      
      double p_prune, double p_grow,
      double alpha, double beta, 
      const NumericVector& prop_prob
  );
  
  void CHANGE(
      const NumericMatrix& xpred, 
      const NumericVector xcut[],
      double sigma2, double sigma_mu, 
      const NumericVector& R,
      IntegerMatrix& Obs_list, 
      
      double p_prune, double p_grow,
      double alpha, double beta, 
      const NumericVector& prop_prob
  );
  
  void Mean_Parameter(
      NumericMatrix& Tree,
      double sigma2, double sigma_mu, 
      const NumericVector& R, 
      const IntegerMatrix& Obs_list
  );
  
  void Mean_Parameter_pred(
      NumericMatrix& Tree, 
      const NumericVector xcut[], 
      const NumericMatrix& xpred_mult, 
      const int n,
      const int ind=0
  );
};

#endif