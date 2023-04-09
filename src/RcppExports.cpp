// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// cd_score
float cd_score(Eigen::MatrixXd X, Eigen::MatrixXd X0);
RcppExport SEXP _TSDFGS_cd_score(SEXP XSEXP, SEXP X0SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type X(XSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type X0(X0SEXP);
    rcpp_result_gen = Rcpp::wrap(cd_score(X, X0));
    return rcpp_result_gen;
END_RCPP
}
// pev_score
float pev_score(Eigen::MatrixXd X, Eigen::MatrixXd X0);
RcppExport SEXP _TSDFGS_pev_score(SEXP XSEXP, SEXP X0SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type X(XSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type X0(X0SEXP);
    rcpp_result_gen = Rcpp::wrap(pev_score(X, X0));
    return rcpp_result_gen;
END_RCPP
}
// r_score
float r_score(Eigen::MatrixXd X, Eigen::MatrixXd X0);
RcppExport SEXP _TSDFGS_r_score(SEXP XSEXP, SEXP X0SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type X(XSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type X0(X0SEXP);
    rcpp_result_gen = Rcpp::wrap(r_score(X, X0));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_TSDFGS_cd_score", (DL_FUNC) &_TSDFGS_cd_score, 2},
    {"_TSDFGS_pev_score", (DL_FUNC) &_TSDFGS_pev_score, 2},
    {"_TSDFGS_r_score", (DL_FUNC) &_TSDFGS_r_score, 2},
    {NULL, NULL, 0}
};

RcppExport void R_init_TSDFGS(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
