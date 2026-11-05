// r-score (Ou & Liao, 2019, doi:10.1007/s00122-019-03387-0), PEV, and CD criteria.
//
// The values equal those of v2.4.2. Matrices of size n0 x n0 are never formed:
//   p <= n: p x p matrices (primal form, B = (X'X + lambda I)^-1).
//   p >  n: n x n matrices (dual form, K = (XX' + lambda I)^-1).
// The r-score keeps the scalar c = 1 - 1/n0 of v2.4.2 where a centering matrix might be
// intended (see NEWS.md).
#include "scores.h"

#include <cmath>
#include <limits>

namespace tsdfgs {
namespace {

float AsFloat(double x) { return static_cast<float>(x); }

double RScore(const Mat& x, const Mat* x0, Index n0, const TestCache* cache, double lam, bool legacy) {
  const Index n = x.rows(), p = x.cols();
  const double lambda = lam > 0 ? lam : 1.0 / static_cast<double>(p);
  const double c = 1.0 - 1.0 / static_cast<double>(n0);
  double q1, q2, q12;
  if (p <= n) {
    const Mat gx = CrossProd(x);
    Mat g0_local;
    const Mat* g0 = (cache && cache->g0.size() > 0) ? &cache->g0 : nullptr;
    if (!g0) {
      g0_local = (cache && cache->total.size() > 0) ? Mat(cache->total - gx) : CrossProd(*x0);
      g0 = &g0_local;
    }
    const SpdSolver solver(gx, lambda);
    const Mat m = solver.Solve(gx);      // M = B X'X
    const Mat mb = solver.Solve(m);      // M B (B and X'X commute)
    const Mat g0m = *g0 * m;
    q1 = n0 - 1.0 + c * c * g0->trace();
    q2 = c * c * (SumProd(mb, *g0) + SumProd(g0m, m));
    q12 = c * SumProd(*g0, m);
  } else {
    const Mat gm = x * x.transpose();
    const Mat z = *x0 * x.transpose();
    const SpdSolver solver(gm, lambda);
    const Mat zk = solver.Solve(z.transpose()).transpose();  // Z K
    const Mat t = zk * gm;
    q1 = n0 - 1.0 + c * c * SumProd(*x0, *x0);
    q2 = c * c * (SumProd(zk, zk) + SumProd(t, zk));
    q12 = c * SumProd(zk, z);
  }
  if (legacy) {
    // v2.4.2 stored q1, q2, and q12 as float and returned a float.
    const float f1 = AsFloat(q1), f2 = AsFloat(q2), f12 = AsFloat(q12);
    const double s = std::sqrt(static_cast<double>(f1 * f2));
    return static_cast<double>(AsFloat(static_cast<double>(f12) / s));
  }
  return q12 / std::sqrt(q1 * q2);
}

Mat WithIntercept(const Mat& x) {
  Mat out(x.rows(), x.cols() + 1);
  out.col(0).setOnes();
  out.rightCols(x.cols()) = x;
  return out;
}

double PevScore(const Mat& xm, const Mat* x0m, const TestCache* cache, double lam, bool legacy) {
  const Index p = xm.cols();
  const double lambda = lam > 0 ? lam : 1.0 / static_cast<double>(p);
  const Mat x = WithIntercept(xm);
  double value;
  if (p + 1 <= x.rows()) {
    const Mat gx = CrossProd(x);
    const SpdSolver solver(gx, lambda);
    if (cache && cache->g0.size() > 0) {
      value = solver.Solve(cache->g0).trace();
    } else if (cache && cache->total.size() > 0) {
      value = solver.Solve(Mat(cache->total - gx)).trace();
    } else {
      value = solver.Solve(CrossProd(WithIntercept(*x0m))).trace();
    }
  } else {
    const Mat x0 = WithIntercept(*x0m);
    const Mat gm = x * x.transpose();
    const Mat z = x0 * x.transpose();
    const SpdSolver solver(gm, lambda);
    const Mat zk = solver.Solve(z.transpose()).transpose();
    value = (SumProd(x0, x0) - SumProd(zk, z)) / lambda;
  }
  return legacy ? static_cast<double>(AsFloat(value)) : value;
}

double CdScore(const Mat& x, const Mat& x0, double lam, bool legacy) {
  const double lambda = lam > 0 ? lam : 1.0;
  const Index n = x.rows(), p = x.cols();
  Eigen::ArrayXd num;
  if (p <= n) {
    const Mat gx = CrossProd(x);
    const SpdSolver solver(gx, lambda);
    const Mat q = solver.Solve(solver.Solve(gx));  // B X'X B
    num = ((x0 * q).cwiseProduct(x0)).rowwise().sum().array();
  } else {
    const Mat gm = x * x.transpose();
    const Mat z = x0 * x.transpose();
    const SpdSolver solver(gm, lambda);
    const Mat zk = solver.Solve(z.transpose()).transpose();
    num = zk.cwiseProduct(zk).rowwise().sum().array();
  }
  const Eigen::ArrayXd den = x0.cwiseProduct(x0).rowwise().sum().array();
  const double value = (num / den).mean();
  return legacy ? static_cast<double>(AsFloat(value)) : value;
}


}  // namespace

TestCache MakeTestCache(Criterion crit, const Mat& x0) {
  TestCache cache;
  if (crit == kRScore) cache.g0 = CrossProd(x0);
  if (crit == kPev || crit == kGOpt) cache.g0 = CrossProd(WithIntercept(x0));
  if (crit == kMspe) cache.g0 = CrossProd(x0);
  return cache;
}

TestCache MakeCandidateCache(Criterion crit, const Mat& candidates) {
  TestCache cache;
  if (crit == kRScore) cache.total = CrossProd(candidates);
  if (crit == kPev || crit == kGOpt) cache.total = CrossProd(WithIntercept(candidates));
  if (crit == kMspe) cache.total = CrossProd(candidates);
  return cache;
}

bool IsSetOnly(Criterion crit) {
  return crit == kAOpt || crit == kDOpt || crit == kEOpt || crit == kGvOverall ||
         crit == kGvAverage || crit == kAvgGrmSelf;
}

bool TargetIsCandidates(Criterion crit) {
  return crit == kMspe || crit == kAvgGrm || crit == kCdMean || crit == kCdMeanV2;
}

bool CanUseComplement(Criterion crit, Index n, Index p) {
  if (crit == kRScore) return p <= n;
  if (crit == kPev || crit == kGOpt) return p + 1 <= n;
  if (crit == kMspe) return p <= n;
  return false;
}

bool UseDual(Criterion crit, Index n, Index p) {
  if (crit == kPev) return n <= p;
  if (crit == kRScore || crit == kCd) return n < p;
  return false;
}

double ScoreGramDual(Criterion crit, const Mat& gram, const std::vector<int>& train,
                     const std::vector<int>& test, double lam, bool legacy) {
  const Index n = static_cast<Index>(train.size()), n0 = static_cast<Index>(test.size());
  Mat gm(n, n), z(n0, n);
  for (Index j = 0; j < n; ++j) {
    for (Index i = 0; i < n; ++i) gm(i, j) = gram(train[static_cast<size_t>(i)], train[static_cast<size_t>(j)]);
    for (Index i = 0; i < n0; ++i) z(i, j) = gram(test[static_cast<size_t>(i)], train[static_cast<size_t>(j)]);
  }
  double diag_sum = 0.0;
  Eigen::ArrayXd dg(n0);
  for (Index i = 0; i < n0; ++i) {
    dg(i) = gram(test[static_cast<size_t>(i)], test[static_cast<size_t>(i)]);
    diag_sum += dg(i);
  }
  const double lambda = lam;  // already resolved by the caller (default 1/p, or 1 for CD)
  if (crit == kPev) {
    gm.array() += 1.0;
    z.array() += 1.0;
    const SpdSolver solver(gm, lambda);
    const Mat zk = solver.Solve(z.transpose()).transpose();
    const double value = (diag_sum + static_cast<double>(n0) - SumProd(zk, z)) / lambda;
    return legacy ? static_cast<double>(AsFloat(value)) : value;
  }
  const SpdSolver solver(gm, lambda);
  const Mat zk = solver.Solve(z.transpose()).transpose();
  if (crit == kCd) {
    const Eigen::ArrayXd num = zk.cwiseProduct(zk).rowwise().sum().array();
    const double value = (num / dg).mean();
    return legacy ? static_cast<double>(AsFloat(value)) : value;
  }
  const double c = 1.0 - 1.0 / static_cast<double>(n0);
  const Mat t = zk * gm;
  const double q1 = n0 - 1.0 + c * c * diag_sum;
  const double q2 = c * c * (SumProd(zk, zk) + SumProd(t, zk));
  const double q12 = c * SumProd(zk, z);
  if (legacy) {
    const float f1 = AsFloat(q1), f2 = AsFloat(q2), f12 = AsFloat(q12);
    const double s = std::sqrt(static_cast<double>(f1 * f2));
    return static_cast<double>(AsFloat(static_cast<double>(f12) / s));
  }
  return q12 / std::sqrt(q1 * q2);
}

double Score(Criterion crit, const Mat& x, const Mat* x0, Index n0, const TestCache* cache,
             double lambda, bool legacy) {
  if (x0 && x.cols() != x0->cols()) throw std::runtime_error("X and X0 must have the same number of columns.");
  if (crit >= kAOpt) {
    if (x.rows() < 1 || x.cols() < 1) throw std::runtime_error("X must have at least one row and one column.");
    if (!IsSetOnly(crit) && !x0 && !(cache && (cache->total.size() > 0 || cache->g0.size() > 0))) {
      throw std::runtime_error("internal error: the test matrix is missing");
    }
    return ScoreV3(crit, x, x0, n0, cache, lambda);
  }
  if (x.rows() < 1 || n0 < 1 || x.cols() < 1) {
    throw std::runtime_error("X and X0 must have at least one row and one column.");
  }
  const bool complement = cache && cache->total.size() > 0 && cache->g0.size() == 0;
  if (!x0 && !(complement && CanUseComplement(crit, x.rows(), x.cols()))) {
    throw std::runtime_error("internal error: the test matrix is missing");
  }
  switch (crit) {
    case kRScore: return RScore(x, x0, n0, cache, lambda, legacy);
    case kPev: return PevScore(x, x0, cache, lambda, legacy);
    default: return CdScore(x, *x0, lambda, legacy);
  }
}

}  // namespace tsdfgs
