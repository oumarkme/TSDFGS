// Training set criteria added in v3. Each one is written from its published definition:
//   A-, D-, E-, G-optimality, PEVmax    Akdemir & Isidro-Sanchez (2019, doi:10.1038/s41598-018-38081-6);
//                                       the STPGA manual lists the same definitions
//   MSPE (ridge), GVoverall, GVaverage, CDmean(v2)
//                                       Chen, Sung & Liao (2024, doi:10.1007/s00122-024-04766-y)
//   CDmean (contrasts)                  Rincent et al. (2012, doi:10.1534/genetics.112.141473)
//   Avg_GRM, Avg_GRM_self, Avg_GRM_MinMax
//                                       Fernandez-Gonzalez et al. (2023, doi:10.1007/s00122-023-04265-6)
//   PAM (k-medoids cost)                Kaufman & Rousseeuw; Guo et al. (2019); Rio et al. (2021)
// No source code of other packages was read or copied.
//
// X holds the training rows. The target is the test set X0, or, without a test set, the
// candidates (criteria in TargetIsCandidates()) or the candidates left out of the training set.
#include <cmath>
#include <limits>

#include "scores.h"

namespace tsdfgs {
namespace {

double Lam(double lambda, double fallback) { return lambda > 0 ? lambda : fallback; }

Mat WithIntercept(const Mat& x) {
  Mat out(x.rows(), x.cols() + 1);
  out.col(0).setOnes();
  out.rightCols(x.cols()) = x;
  return out;
}

Mat StackRows(const Mat& a, const Mat& b) {
  Mat out(a.rows() + b.rows(), a.cols());
  out.topRows(a.rows()) = a;
  out.bottomRows(b.rows()) = b;
  return out;
}

// The target population: the test set, or all candidates if the criterion asks for them.
Mat Target(Criterion crit, const Mat& x, const Mat& x0, const TestCache* cache) {
  return (TargetIsCandidates(crit) && cache && cache->untargeted) ? StackRows(x, x0) : x0;
}

double Norm2(const Eigen::VectorXd& v) { return v.squaredNorm(); }

// ---- ridge information matrix: tr, log det, smallest eigenvalue of X'X + lambda I ----
double AOpt(const Mat& x, double lam) {
  const Index n = x.rows(), p = x.cols();
  const double lambda = Lam(lam, 1.0 / static_cast<double>(p));
  if (p <= n) return SpdSolver(CrossProd(x), lambda).Solve(Mat::Identity(p, p)).trace();
  // the p - n zero eigenvalues of X'X each add 1 / lambda
  return SpdSolver(x * x.transpose(), lambda).Solve(Mat::Identity(n, n)).trace() +
         static_cast<double>(p - n) / lambda;
}

double DOpt(const Mat& x, double lam) {
  const Index n = x.rows(), p = x.cols();
  const double lambda = Lam(lam, 1.0 / static_cast<double>(p));
  if (p <= n) return SpdSolver(CrossProd(x), lambda).LogDet();
  return SpdSolver(x * x.transpose(), lambda).LogDet() +
         static_cast<double>(p - n) * std::log(lambda);
}

// largest eigenvalue of (X'X + lambda I)^-1 = 1 / (smallest eigenvalue of X'X + lambda)
double EOpt(const Mat& x, double lam) {
  const Index n = x.rows(), p = x.cols();
  const double lambda = Lam(lam, 1.0 / static_cast<double>(p));
  if (p > n) return 1.0 / lambda;  // X'X is singular
  Eigen::SelfAdjointEigenSolver<Mat> es(CrossProd(x), Eigen::EigenvaluesOnly);
  return 1.0 / (es.eigenvalues().minCoeff() + lambda);
}

// ---- prediction error variance of the target, with an intercept column ----
// PEV matrix of the target rows: X0e (Xe'Xe + lambda I)^-1 X0e' (primal),
// or (X0e X0e' - Z (Xe Xe' + lambda I)^-1 Z') / lambda with Z = X0e Xe' (dual).
double PevMax(const Mat& xm, const Mat& x0m, double lam) {
  const Index p = xm.cols();
  const double lambda = Lam(lam, 1.0 / static_cast<double>(p));
  const Mat x = WithIntercept(xm), x0 = WithIntercept(x0m);
  if (p + 1 <= x.rows()) {
    const SpdSolver solver(CrossProd(x), lambda);
    const Mat q = solver.Solve(x0.transpose());
    return (x0.transpose().cwiseProduct(q)).colwise().sum().maxCoeff();
  }
  const Mat gm = x * x.transpose();
  const Mat z = x0 * x.transpose();
  const SpdSolver solver(gm, lambda);
  const Mat zk = solver.Solve(z.transpose()).transpose();
  const Eigen::VectorXd d = x0.rowwise().squaredNorm() - zk.cwiseProduct(z).rowwise().sum();
  return d.maxCoeff() / lambda;
}

double GOpt(const Mat& xm, const Mat* x0m, const TestCache* cache, double lam) {
  const Index p = xm.cols();
  const double lambda = Lam(lam, 1.0 / static_cast<double>(p));
  const Mat x = WithIntercept(xm);
  if (p + 1 <= x.rows()) {
    const Mat gx = CrossProd(x);
    Mat g0_local;
    const Mat* g0 = (cache && cache->g0.size() > 0) ? &cache->g0 : nullptr;
    if (!g0) {
      g0_local = (cache && cache->total.size() > 0) ? Mat(cache->total - gx)
                                                    : CrossProd(WithIntercept(*x0m));
      g0 = &g0_local;
    }
    // eigenvalues of X0 B X0' equal those of L^-1 G0 L^-T with L L' = X'X + lambda I
    const SpdSolver solver(gx, lambda);
    const Mat t1 = solver.LSolve(*g0);
    const Mat s = solver.LSolve(t1.transpose());
    Eigen::SelfAdjointEigenSolver<Mat> es(s, Eigen::EigenvaluesOnly);
    return es.eigenvalues().maxCoeff();
  }
  const Mat x0 = WithIntercept(*x0m);
  const Mat gm = x * x.transpose();
  const Mat z = x0 * x.transpose();
  const SpdSolver solver(gm, lambda);
  const Mat zk = solver.Solve(z.transpose()).transpose();
  const Mat pev = (x0 * x0.transpose() - zk * z.transpose()) / lambda;
  Eigen::SelfAdjointEigenSolver<Mat> es(pev, Eigen::EigenvaluesOnly);
  return es.eigenvalues().maxCoeff();
}

// ---- MSPE of ridge regression, approximation of Chen et al. (2024), Eq. 4 ----
// MSPE = 1 + tr(Xc A A' Xc')/nc + || Xc - Xc A Xt ||_F^2, A = Xt'(Xt Xt' + lambda I)^-1.
double Mspe(const Mat& x, const Mat* x0, Index n0, const TestCache* cache, double lam) {
  const Index n = x.rows(), p = x.cols();
  const double lambda = Lam(lam, 1.0);
  const bool untargeted = cache && cache->untargeted;
  const double nc = static_cast<double>(untargeted ? n + n0 : n0);
  if (p <= n) {
    const Mat gx = CrossProd(x);
    Mat gc;
    if (untargeted) {
      gc = (cache && cache->total.size() > 0) ? cache->total : Mat(gx + CrossProd(*x0));
    } else if (cache && cache->g0.size() > 0) {
      gc = cache->g0;
    } else {
      gc = CrossProd(*x0);
    }
    const SpdSolver solver(gx, lambda);
    const Mat s1 = solver.Solve(gx);                       // B X'X
    const Mat bgxb = solver.Solve(s1);                     // B X'X B
    const Mat b2 = solver.Solve(solver.Solve(Mat::Identity(p, p)));  // B B
    return 1.0 + SumProd(bgxb, gc) / nc + lambda * lambda * SumProd(b2, gc);
  }
  const Mat xc = untargeted ? StackRows(x, *x0) : *x0;
  const Mat gm = x * x.transpose();
  const SpdSolver solver(gm, lambda);
  const Mat zc = xc * x.transpose();
  const Mat zk = solver.Solve(zc.transpose()).transpose();  // Xc A
  const Mat resid = xc - zk * x;
  return 1.0 + zk.squaredNorm() / nc + resid.squaredNorm();
}

// ---- relationship matrix G = X X' / p ----
double GvOverall(const Mat& x) {
  const Index n = x.rows(), p = x.cols();
  Mat k = (x * x.transpose()) / static_cast<double>(p);
  Eigen::LLT<Mat> llt(k);
  if (llt.info() != Eigen::Success) return -std::numeric_limits<double>::infinity();
  return 2.0 * llt.matrixL().toDenseMatrix().diagonal().array().log().sum();
}

double ColSumDot(const Mat& a, const Mat& b) {
  const Eigen::VectorXd sa = a.colwise().sum().transpose();
  const Eigen::VectorXd sb = b.colwise().sum().transpose();
  return sa.dot(sb);
}

double AvgGrm(Criterion crit, const Mat& x, const Mat& x0, const TestCache* cache) {
  const double n = static_cast<double>(x.rows()), p = static_cast<double>(x.cols());
  const Eigen::VectorXd s = x.colwise().sum().transpose();
  const double self = Norm2(s) / (n * n * p);
  if (crit == kAvgGrmSelf) return self;
  const bool untargeted = cache && cache->untargeted;
  if (crit == kAvgGrm) {
    // target: the test set, or all candidates (training set included) without a test set
    if (untargeted) {
      const Eigen::VectorXd sc = s + Eigen::VectorXd(x0.colwise().sum().transpose());
      return s.dot(sc) / (n * (n + static_cast<double>(x0.rows())) * p);
    }
    return ColSumDot(x, x0) / (n * static_cast<double>(x0.rows()) * p);
  }
  // Avg_GRM_MinMax: target is the test set, or the candidates left out of the training set
  return ColSumDot(x, x0) / (n * static_cast<double>(x0.rows()) * p) - self;
}

// ---- CDmean (contrasts, Rincent et al. 2012) and CDmean(v2) (Chen et al. 2024) ----
// With K = XX'/p, M = I - J/n_t and W = (M K_t + lambda I)^-1 M, the covariance of the
// estimated and true genotypic values of the target rows r is K_rT W K_Tr.
// CDmean:    mean over targets i of  v_i' W v_i / (c_i' K c_i),  v_i = K_TP c_i,
//            c_i = e_i - 1/n_P (contrast with the target mean).
// CDmean(v2): mean over targets of  diag(K_PT W K_TP) / diag(K_PP).
double CdMean(Criterion crit, const Mat& x, const Mat& target, double lam) {
  const Index nt = x.rows();
  const double p = static_cast<double>(x.cols());
  const double lambda = Lam(lam, 1.0);
  const Mat k = (x * x.transpose()) / p;
  const Mat m = Mat::Identity(nt, nt) - Mat::Constant(nt, nt, 1.0 / static_cast<double>(nt));
  const Mat a = m * k + lambda * Mat::Identity(nt, nt);
  const Eigen::PartialPivLU<Mat> lu(a);
  Mat v;        // columns: v_i
  Eigen::ArrayXd den;
  if (crit == kCdMean) {
    Mat xc = target;
    xc.rowwise() -= target.colwise().mean();  // rows of the contrasts c_i' X_P
    v = (x * xc.transpose()) / p;
    den = (xc.rowwise().squaredNorm() / p).array();
  } else {
    v = (x * target.transpose()) / p;
    den = (target.rowwise().squaredNorm() / p).array();
  }
  const Mat wv = lu.solve(m * v);
  const Eigen::ArrayXd num = v.cwiseProduct(wv).colwise().sum().transpose().array();
  return (num / den).mean();
}

// ---- CDmax: largest per-individual CD of the ridge form used by cd_score (lambda = 1) ----
double CdMax(const Mat& x, const Mat& x0, double lam) {
  const double lambda = Lam(lam, 1.0);
  const Index n = x.rows(), p = x.cols();
  Eigen::ArrayXd num;
  if (p <= n) {
    const Mat gx = x.transpose() * x;
    const SpdSolver solver(gx, lambda);
    const Mat q = solver.Solve(solver.Solve(gx));
    num = ((x0 * q).cwiseProduct(x0)).rowwise().sum().array();
  } else {
    const Mat gm = x * x.transpose();
    const SpdSolver solver(gm, lambda);
    const Mat zk = solver.Solve((x0 * x.transpose()).transpose()).transpose();
    num = zk.cwiseProduct(zk).rowwise().sum().array();
  }
  const Eigen::ArrayXd den = x0.cwiseProduct(x0).rowwise().sum().array();
  return (num / den).maxCoeff();
}

// ---- k-medoids cost: sum over target rows of the distance to the nearest training row ----
double PamCost(const Mat& x, const Mat& target) {
  const Eigen::VectorXd rt = target.rowwise().squaredNorm();
  const Eigen::VectorXd rx = x.rowwise().squaredNorm();
  Mat d2 = (-2.0 * target * x.transpose()).eval();
  d2.colwise() += rt;
  d2.rowwise() += rx.transpose();
  return d2.cwiseMax(0.0).rowwise().minCoeff().cwiseSqrt().sum();
}

}  // namespace

double ScoreV3(Criterion crit, const Mat& x, const Mat* x0, Index n0, const TestCache* cache,
               double lambda) {
  switch (crit) {
    case kAOpt: return AOpt(x, lambda);
    case kDOpt: return DOpt(x, lambda);
    case kEOpt: return EOpt(x, lambda);
    case kGvAverage: return x.squaredNorm() / static_cast<double>(x.cols());
    case kGvOverall:
      if (x.cols() < x.rows()) {
        throw std::runtime_error(
            "GVoverall needs at least as many columns as training individuals (use all PCs or the markers)");
      }
      return GvOverall(x);
    case kAvgGrmSelf: return AvgGrm(crit, x, x, cache);
    case kGOpt: return GOpt(x, x0, cache, lambda);
    case kMspe: return Mspe(x, x0, n0, cache, lambda);
    default: break;
  }
  if (!x0) throw std::runtime_error("internal error: the test matrix is missing");
  switch (crit) {
    case kPevMax: return PevMax(x, *x0, lambda);
    case kAvgGrm:
    case kAvgGrmMinMax: return AvgGrm(crit, x, *x0, cache);
    case kCdMean:
    case kCdMeanV2: return CdMean(crit, x, Target(crit, x, *x0, cache), lambda);
    case kPam: return PamCost(x, *x0);
    case kCdMax: return CdMax(x, *x0, lambda);
    default: throw std::runtime_error("unknown criterion");
  }
}

}  // namespace tsdfgs
