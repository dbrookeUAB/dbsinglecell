// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// matrix_rollsum2
Rcpp::NumericMatrix matrix_rollsum2(SEXP x, int n, int margin);
RcppExport SEXP _dbsinglecell_matrix_rollsum2(SEXP xSEXP, SEXP nSEXP, SEXP marginSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type x(xSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< int >::type margin(marginSEXP);
    rcpp_result_gen = Rcpp::wrap(matrix_rollsum2(x, n, margin));
    return rcpp_result_gen;
END_RCPP
}
// substr_r
std::string substr_r(std::string str, int n);
RcppExport SEXP _dbsinglecell_substr_r(SEXP strSEXP, SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type str(strSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    rcpp_result_gen = Rcpp::wrap(substr_r(str, n));
    return rcpp_result_gen;
END_RCPP
}
// substr_l
std::string substr_l(std::string str, int n);
RcppExport SEXP _dbsinglecell_substr_l(SEXP strSEXP, SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type str(strSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    rcpp_result_gen = Rcpp::wrap(substr_l(str, n));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_dbsinglecell_matrix_rollsum2", (DL_FUNC) &_dbsinglecell_matrix_rollsum2, 3},
    {"_dbsinglecell_substr_r", (DL_FUNC) &_dbsinglecell_substr_r, 2},
    {"_dbsinglecell_substr_l", (DL_FUNC) &_dbsinglecell_substr_l, 2},
    {NULL, NULL, 0}
};

RcppExport void R_init_dbsinglecell(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
