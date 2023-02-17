#include <Rcpp.h>
#include <iostream>
#include <string>

using namespace Rcpp;

// [[Rcpp::export]]
std::string substr_r(std::string str, int n) {
  return str.substr(str.length()-n);
}

// [[Rcpp::export]]
std::string substr_l(std::string str, int n) {
  return str.substr(0, n-1);
}

  
  
/*** R
substr_r('apple', 2)
substr_l('apple', 2)
*/
