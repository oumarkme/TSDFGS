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


//' r-score
//' 
//' This function calculate r-score <doi:10.1007/s00122-019-03387-0> by given training set and test set. 
//' 
//' @author Jen-Hsiang Ou
//' 
//' @param X A numeric matrix. The training set genotypic information matrix can be given as genotype matrix (coded as -1, 0, 1) or principle component matrix (row: sample; column: marker).
//' @param X0 A numeric mareix. The test set genotypic information matrix can be given as genotype matrix (coded as -1, 0, 1) or principle component matrix (row: sample; column: marker).
//' @param lambda A float number. The default value is set as 1.0 and you may change it with the knowledge of heritability.
//' @return A floating-point number, r-score.
//' 
//' @import Rcpp
//' 
//' @export
//' @examples
//' data(geno)
//' \dontrun{r_score(geno[1:50, ], geno[51:100])}
//' @rawNamespace useDynLib(TSDFGS); import(RcppEigen); importFrom(Rcpp, evalCpp)
//' 
// [[Rcpp::export]]
float r_score(Eigen::MatrixXd X, Eigen::MatrixXd X0, float lambda = 1.0)
{
    int nr = X.rows();
    int n0r = X0.rows();
    Eigen::MatrixXd A = X.transpose() * ((X * X.transpose() + MatrixXd::Identity(nr, nr) * lambda).inverse());
    Eigen::MatrixXd IJ = MatrixXd::Identity(n0r, n0r) - (MatrixXd::Ones(n0r, n0r) / n0r);
    float q1 = n0r - 1 + (IJ * X0).array().square().sum();
    Eigen::MatrixXd IJX0A = IJ * X0 * A;
    Eigen::MatrixXd IJX0AX = IJX0A * X;
    float q2 = IJX0A.array().square().sum() + IJX0AX.array().square().sum();
    float q12 = (X0.transpose() * IJX0AX).trace();
    return q12 / (sqrt(q1 * q2));
}
