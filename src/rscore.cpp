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


//' r-score
//' 
//' This function calculate r-score <doi:10.1007/s00122-019-03387-0> by given training set and test set. 
//' 
//' @author Jen-Hsiang Ou
//' 
//' @param X A numeric matrix. The training set genotypic information matrix can be given as genotype matrix (coded as -1, 0, 1) or principle component matrix (row: sample; column: marker).
//' @param X0 A numeric mareix. The test set genotypic information matrix can be given as genotype matrix (coded as -1, 0, 1) or principle component matrix (row: sample; column: marker).
//' @return A floating-point number, r-score.
//' 
//' @import Rcpp
//' 
//' @export
//' @examples
//' data(geno)
//' \dontrun{r_score(geno[1:50, ], geno[51:100])}
//' 
// [[Rcpp::export]]
float r_score(Eigen::MatrixXd X, Eigen::MatrixXd X0)
{
    int nr = X.rows();
    int nc = X.cols();
    int n0r = X0.rows();
    Eigen::MatrixXd A = X.transpose() * ((X * X.transpose() + MatrixXd::Identity(nr, nr) * (1.0 / nc)).inverse());
    Eigen::MatrixXd IJ = MatrixXd::Identity(n0r, n0r) - (MatrixXd::Identity(n0r, n0r) * (1.0 / n0r));
    float q1 = n0r - 1 + (IJ * X0).array().square().sum();
    Eigen::MatrixXd IJX0A = IJ * X0 * A;
    Eigen::MatrixXd IJX0AX = IJX0A * X;
    float q2 = IJX0A.array().square().sum() + IJX0AX.array().square().sum();
    float q12 = (X0.transpose() * IJX0AX).trace();
    return q12 / (sqrt(q1 * q2));
}
