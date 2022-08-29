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

void update_xpred_temp(
    NumericMatrix& xpred_temp, const NumericVector xcut[], 
    const IntegerVector Split_hist[], const IntegerVector Value_hist[],
    const IntegerVector Side_hist[], const int i
);

void get_hist(
    IntegerVector Split_hist[], IntegerVector Value_hist[], IntegerVector Side_hist[], bool verbose,
    IntegerVector dt_parent, IntegerVector dt_position, IntegerVector dt_Split, IntegerVector dt_Value, IntegerVector dt_Terminal
);

// [[Rcpp::export]]
NumericVector predict(
    const List xcut_list,
    const NumericMatrix& xpred_mult,
    const int n,
    const int ind,
    List dt,
    int ncol,
    bool is_mar=true
  ) {
  IntegerVector temp;
  IntegerVector dt_position = dt["position"];
  IntegerVector dt_parent   = dt["parent"];
  IntegerVector dt_Terminal = dt["Terminal"];
  temp                      = dt["Split"]; 
  IntegerVector dt_Split = clone(temp) - 1;
  temp                   = dt["Value"]; 
  IntegerVector dt_Value = clone(temp) - 1;
  NumericVector dt_MU    = dt["MU"];
  temp                   = dt["begin"]; 
  IntegerVector dt_begin = clone(temp) - 1;
  temp                   = dt["end"];   
  IntegerVector dt_end   = clone(temp) - 1;
  
  NumericVector xcut[ncol];
  for (int j=0; j<ncol; j++) {
    NumericVector xcut_temp = xcut_list[j];
    xcut[j] = clone(xcut_temp);
  }
  
  IntegerVector tnode_idx = terminal_nodes(dt_Terminal);
  
  // Rcout << "tnode_idx : " << tnode_idx << endl;
  
  const int Terminal_len = tnode_idx.length();
  if (Terminal_len > 1) {
    // placeholder for hist
    IntegerVector Split_hist[Terminal_len], Value_hist[Terminal_len], Side_hist[Terminal_len];
    get_hist(Split_hist, Value_hist, Side_hist, false, 
             dt_parent, dt_position, dt_Split, dt_Value, dt_Terminal);
    
    
    LogicalVector ind_list = rep(true, 2*n); // flag variable for matrix subsetting
    NumericVector T(2*n);
    
    NumericMatrix xpred_Mult;
    const int P = xpred_mult.ncol() - 1; // idx of last col
    if (is_mar) { // mar model
      xpred_Mult = clone(xpred_mult);
      if (ind == 1) {
        xpred_Mult(_, 0) = rep(1, 2*n);
      }
    } else {              // sep model
      xpred_Mult = xpred_mult;
    }
    NumericMatrix xpred_temp = clone(xpred_Mult);
    
    
    for (int i=0; i<Terminal_len; i++) {
      int count = 0;
      while (count < Split_hist[i].length()) {
        
        // Rcout << "while count = " << count << " < " << Split_hist[i].length() << endl;
        // Rcout << "i = " << i << endl;
        // Rcout << "xpred_temp_nrow : " << xpred_temp.nrow() << endl;
        
        if (xpred_temp.nrow() != 0) {
          
          int idx = Split_hist[i](count); 
          double value = xcut[idx](Value_hist[i](count));
          
          // Rcout << "idx : " << idx << endl;
          // Rcout << "value : " << value << endl;
          // Rcout << "side  : " << Side_hist[i][count] << endl;
          
          LogicalVector sub_ind = rep(false, xpred_temp.nrow());
          if (Side_hist[i][count] == 0) {
            for (int k=0; k<xpred_temp.nrow(); k++) {
              if (xpred_temp(k, idx) <  value) sub_ind(k) = true;
            }
          } else {
            for (int k=0; k<xpred_temp.nrow(); k++) {
              if (xpred_temp(k, idx) >= value) sub_ind(k) = true;
            }
          }
          
          // NumericVector temp123 = xpred_temp(_, idx);
          // Rcout << "\nsub_ind\n" << temp123 << endl<< endl;
          // Rcout << "remaining row : " << sum(sub_ind) << endl;
          
          xpred_temp = subset_by_row(xpred_temp, sub_ind);
        }  
        count++;
      } // end of while
      //Rcout << endl;
      if (xpred_temp.nrow() == 0) {
        xpred_temp = clone(xpred_Mult);
      } else {
        ind_list[xpred_temp(_, P)] = false; // will remove rows with idx = xpred_temp(_, P)
        
        // Rcout << "i : " << i << endl;
        // NumericVector temp = xpred_temp(_, P);
        // Rcout << "T idx : " << temp << endl;
        
        T[xpred_temp(_, P)] = dt_MU(tnode_idx(i));
        xpred_temp = subset_by_row(xpred_Mult, ind_list);
        //Rcout << endl;
      }
    }
    return T;
  } else {
    return rep(dt_MU(0), 2*n);
  }
}

void get_hist(
    IntegerVector Split_hist[], IntegerVector Value_hist[], IntegerVector Side_hist[], bool verbose,
    IntegerVector dt_parent, IntegerVector dt_position, IntegerVector dt_Split, 
    IntegerVector dt_Value, IntegerVector dt_Terminal
  ) {
  IntegerVector tnode_idx = terminal_nodes(dt_Terminal);
  const int Terminal_len = tnode_idx.length();
  
  for (int i=0; i<Terminal_len; i++) {
    int parent_node, current_node;
    IntegerVector Split_temp, Value_temp, Side_temp;
    
    parent_node  = dt_parent(tnode_idx(i));
    current_node = dt_position(tnode_idx(i));
    int idx = which(dt_position, parent_node);
    Split_temp.push_front(dt_Split(idx));
    Value_temp.push_front(dt_Value(idx));
    Side_temp.push_front(current_node % 2);
    while (parent_node != 1) {
      idx = which(dt_position, parent_node);
      current_node = dt_position(idx);
      parent_node  = dt_parent(idx);
      
      idx = which(dt_position, parent_node);
      Split_temp.push_front(dt_Split(idx));
      Value_temp.push_front(dt_Value(idx));
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

void update_xpred_temp(
    NumericMatrix& xpred_temp, 
    const NumericVector xcut[],
    const IntegerVector Split_hist[], 
    const IntegerVector Value_hist[], 
    const IntegerVector Side_hist[], 
    const int i) {
  
}