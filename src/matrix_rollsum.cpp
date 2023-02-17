#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
Rcpp::NumericMatrix matrix_rollsum2(SEXP x, int n, int margin) {
  Rcpp::NumericMatrix y(x);
  int NR = y.nrow(), NC = y.ncol();
  Rcpp::NumericMatrix result(NR,NC);
  std::fill( result.begin(), result.end(), NumericVector::get_na() ) ;

  double s=0.0;

  if(margin==1){
    for(int i = 0; i < NR; ++i){
      Rcpp::NumericMatrix::Row tmpvec = y(i,_);
      for(int j = 0; j < NC-n+1;++j){
        for(int q=j; q<j+n;q++){
          s+=tmpvec[q];
        }
        result(i,j+n-1) = s;
        s = 0.0;
      }}}

  if(margin==2){
    for(int i = 0; i < NC; ++i){
      Rcpp::NumericMatrix::Column tmpvec = y(_,i);
      for(int j = 0; j < NR-n+1;++j){
        for(int q=j; q<j+n;q++){
          s+=tmpvec[q];
        }
        result(j+n-1,i) = s;
        s = 0.0;
      }}}

  return result;
}
