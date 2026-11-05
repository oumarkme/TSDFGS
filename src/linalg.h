// Small dense linear algebra helpers on Eigen (RcppEigen).
#ifndef TSDFGS_LINALG_H_
#define TSDFGS_LINALG_H_

#include <RcppEigen.h>

#include <stdexcept>
#include <vector>

namespace tsdfgs {

using Eigen::Index;
using Eigen::Map;
using Eigen::MatrixXd;
using Mat = Eigen::MatrixXd;

// Copies the rows `idx` (0-based) of a matrix.
inline Mat GatherRows(const Map<const MatrixXd>& g, const std::vector<int>& idx) {
  Mat out(static_cast<Index>(idx.size()), g.cols());
  for (Index j = 0; j < g.cols(); ++j)
    for (Index i = 0; i < out.rows(); ++i) out(i, j) = g(idx[i], j);
  return out;
}

// A'A with a rank update on one triangle (half the work of a full product).
inline Mat CrossProd(const Mat& a) {
  Mat c = Mat::Zero(a.cols(), a.cols());
  c.selfadjointView<Eigen::Lower>().rankUpdate(a.transpose());
  return Mat(c.selfadjointView<Eigen::Lower>());
}

inline double SumProd(const Mat& a, const Mat& b) { return a.cwiseProduct(b).sum(); }

// Cholesky solver for the symmetric positive definite matrix G + lambda I.
class SpdSolver {
 public:
  SpdSolver(const Mat& g, double lambda) : llt_(g + lambda * Mat::Identity(g.rows(), g.cols())) {
    if (llt_.info() != Eigen::Success) {
      throw std::runtime_error(
          "the training set matrix is not positive definite; check the input for missing or "
          "non-finite values");
    }
  }
  // Returns (G + lambda I)^-1 B.
  Mat Solve(const Mat& b) const { return llt_.solve(b); }
  // Returns L^-1 B, where L L' = G + lambda I.
  Mat LSolve(const Mat& b) const { return llt_.matrixL().solve(b); }
  // log det(G + lambda I) = 2 sum(log diag(L)).
  double LogDet() const { return 2.0 * llt_.matrixL().toDenseMatrix().diagonal().array().log().sum(); }

 private:
  Eigen::LLT<Mat> llt_;
};

}  // namespace tsdfgs

#endif  // TSDFGS_LINALG_H_
