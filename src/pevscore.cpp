#define EIGEN_PERMANENTLY_DISABLE_STUPID_WARNINGS
#include <Rcpp.h>
#include <RcppEigen.h>
#include <Eigen/Core>
#include <Eigen/Dense>

// [[Rcpp::depends(RcppEigen)]]
using namespace Rcpp;
using Eigen::DiagonalMatrix;
using Eigen::Map;
using Eigen::MatrixXd;
using Eigen::VectorXd;


//' PEV score
//' 
//' This function calculate prediction error variance (PEV) score <doi:10.1186/s12711-015-0116-6> by given training set and test set. 
//' 
//' @param X A numeric matrix. The training set genotypic information matrix can be given as genotype matrix (coded as -1, 0, 1) or principle component matrix (row: sample; column: marker).
//' @param X0 A numeric mareix. The test set genotypic information matrix can be given as genotype matrix (coded as -1, 0, 1) or principle component matrix (row: sample; column: marker).
//' 
//' @return A floating-point number, PEV score.
//' 
//' @author Jen-Hsiang Ou
//' 
//' @import Rcpp
//' @export
//' @examples
//' data(geno)
//' \dontrun{pev_score(geno[1:50, ], geno[51:100])}
//' @rawNamespace useDynLib(TSDFGS); import(RcppEigen); importFrom(Rcpp, evalCpp)
//' 
// [[Rcpp::export]]
float pev_score(Eigen::MatrixXd X, Eigen::MatrixXd X0)
{
    int p = X.cols();
    int n = X.rows();
    int n0 = X0.rows();
    Eigen::MatrixXd x(n, p + 1);
    x << (Eigen::ArrayXd::Zero(n) + 1), X;
    Eigen::MatrixXd x0(n0, p + 1);
    x0 << (Eigen::ArrayXd::Zero(n0) + 1), X0;
    return (x0 * ((x.transpose() * x + (Eigen::MatrixXd::Identity(p + 1, p + 1) * (1.0 / p))).inverse()) * x0.transpose()).trace();
}