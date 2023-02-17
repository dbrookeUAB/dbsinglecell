#include <RcppArmadillo.h>
// [[Rcpp::depends("RcppArmadillo")]]
using namespace Rcpp;

using namespace arma;


// [[Rcpp::export]]
NumericVector test(arma::sp_mat x){;
  int size = x.n_cols;
  x.as_col();
  NumericVector tmp(size);
  for(int i=0; i<tmp.size(); i++){
    tmp[i] = x[0,i];
    }
  return tmp;
}


/***R
test(tmp2)
*/
