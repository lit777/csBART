// Prediction.h
#ifndef __PREDICTION_H__
#define __PREDICTION_H__

#include "decision_tree.h"
#include "functions.h"
#include "subset.h"

using namespace Rcpp;
using namespace sugar;

// Prediction (binary exposure)
void DecisionTree::Mean_Parameter_pred(
    NumericMatrix& Tree,
    const NumericVector xcut[], 
    const NumericMatrix& xpred_mult,
    const int n,
    const int ind) {
  
  IntegerVector tnode_idx = this->terminal_nodes();
  const int Terminal_len = tnode_idx.length();
  if (Terminal_len > 1) {
    // placeholder for hist
    IntegerVector Split_hist[Terminal_len], Value_hist[Terminal_len], Side_hist[Terminal_len];
    this->get_hist(Split_hist, Value_hist, Side_hist);
    
    LogicalVector ind_list = rep(true, 2*n); // flag variable for matrix subsetting
    NumericVector T(2*n);
    
    NumericMatrix xpred_Mult;
    const int P = xpred_mult.ncol() - 1; // idx of last col
    if (this->is_mar()) { // mar model
      xpred_Mult = clone(xpred_mult);
      if (ind == 1) {
        xpred_Mult(_, 0) = rep(1, 2*n);
      }
    } else {              // sep model
      xpred_Mult = xpred_mult;
    }
    NumericMatrix xpred_temp = clone(xpred_Mult);
    
    for (int i=0; i<Terminal_len; i++) {
      this->update_xpred_temp(xpred_temp, xcut, Split_hist, Value_hist, Side_hist, i);
      NumericVector temp = xpred_temp(_, P);
      if (xpred_temp.nrow() == 0) {
        xpred_temp = clone(xpred_Mult);
      } else {
        ind_list[xpred_temp(_, P)] = false; // will remove rows with idx = xpred_temp(_, P)
        T[xpred_temp(_, P)] = this->MU(tnode_idx(i));
        xpred_temp = subset_by_row(xpred_Mult, ind_list);
      }
    }
    Tree(_, this->id) = clone(T);
  } else {
    Tree(_, this->id) = rep(this->MU(0), 2*n);
  }
}

void DecisionTree::get_hist(IntegerVector Split_hist[], IntegerVector Value_hist[], IntegerVector Side_hist[], bool verbose) {
  IntegerVector tnode_idx = this->terminal_nodes();
  const int Terminal_len = tnode_idx.length();
  
  for (int i=0; i<Terminal_len; i++) {
    int parent_node, current_node;
    IntegerVector Split_temp, Value_temp, Side_temp;
    
    parent_node  = this->parent(tnode_idx(i));
    current_node = this->position(tnode_idx(i));
    int idx = which(this->position, parent_node);
    Split_temp.push_front(this->Split(idx));
    Value_temp.push_front(this->Value(idx));
    Side_temp.push_front(current_node % 2);
    while (parent_node != 1) {
      idx = which(this->position, parent_node);
      current_node = this->position(idx);
      parent_node  = this->parent(idx);
      
      idx = which(this->position, parent_node);
      Split_temp.push_front(this->Split(idx));
      Value_temp.push_front(this->Value(idx));
      Side_temp.push_front(current_node % 2);
    }
    Split_hist[i] = clone(Split_temp);
    Value_hist[i] = clone(Value_temp);
    Side_hist[i]  = clone(Side_temp);
  } // end for
  if (verbose == true) {
    for (int i=0; i<Terminal_len; i++) {
      Rcout << "i          : " << i << endl;
      Rcout << "Split_hist : " << Split_hist[i] << endl;
      Rcout << "Value_hist : " << Value_hist[i] << endl;
      Rcout << "Side_hist  : " << Side_hist[i] << endl << endl;
    }
  }
}

void DecisionTree::update_xpred_temp(
    NumericMatrix& xpred_temp, 
    const NumericVector xcut[],
    const IntegerVector Split_hist[], 
    const IntegerVector Value_hist[], 
    const IntegerVector Side_hist[], 
    const int i) {
  int count = 0;
  while (count < Split_hist[i].length()) {
    if (xpred_temp.nrow() != 0) {
      LogicalVector sub_ind = rep(false, xpred_temp.nrow());
      int idx = Split_hist[i](count); 
      double value = xcut[idx](Value_hist[i](count));
      // Rcout << " i    : " << i << endl;
      // Rcout << "idx   : " << idx << endl;
      // Rcout << "value : " << value << endl;
      if (Side_hist[i][count] == 0) {
        for (int i=0; i<xpred_temp.nrow(); i++) {
          if (xpred_temp(i, idx) <  value) sub_ind(i) = true;
        }
      } else {
        for (int i=0; i<xpred_temp.nrow(); i++) {
          if (xpred_temp(i, idx) >= value) sub_ind(i) = true;
        }
      }
      xpred_temp = subset_by_row(xpred_temp, sub_ind);
    }  
    count++;
  } // end of while
}


#endif