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


//' CD-score
//' 
//' This function calculate CD-score <doi:10.1186/1297-9686-28-4-359> by given training set and test set. 
//' 
//' @author Jen-Hsiang Ou
//' 
//' @param X A numeric matrix. The training set genotypic information matrix can be given as genotype matrix (coded as -1, 0, 1) or principle component matrix (row: sample; column: marker).
//' @param X0 A numeric mareix. The test set genotypic information matrix can be given as genotype matrix (coded as -1, 0, 1) or principle component matrix (row: sample; column: marker).
//' @return A floating-point number, CD score.
//' 
//' @import Rcpp
//' 
//' @export
//' 
//' @examples
//' data(geno)
//' \dontrun{cd_score(geno[1:50, ], geno[51:100])}
//' 
// [[Rcpp::export]]
float cd_score(Eigen::MatrixXd X, Eigen::MatrixXd X0)
{
    int p = X.cols();
    Eigen::MatrixXd cd = X0 * (X.transpose() * X + (Eigen::MatrixXd::Identity(p, p) * 1.0)).inverse() * X.transpose();
    return (((cd * cd.transpose()).diagonal().array()) / ((X0 * X0.transpose()).diagonal().array())).mean();
}
