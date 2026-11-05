// Exchange (hill-climbing) search for the training set, and score entry points for R.
//
// The moves and the order of random draws are those of the R code of v2.4.2, drawn from
// R's own generator (R_unif_index, as used by sample.int), so that set.seed() gives the
// same training set. Rcpp saves and restores R's random-number state around each call.
#include <R_ext/Random.h>
#include <Rcpp.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <memory>
#include <random>
#include <string>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "scores.h"

using namespace tsdfgs;

namespace {

int UnifIndex(int n) { return static_cast<int>(R_unif_index(static_cast<double>(n))); }

std::vector<int> ToZeroBased(const Rcpp::IntegerVector& v, const char* what, int n_all) {
  std::vector<int> out(static_cast<size_t>(v.size()));
  for (R_xlen_t i = 0; i < v.size(); ++i) {
    const int r = v[i];
    if (r == NA_INTEGER || r < 1 || r > n_all) {
      throw std::runtime_error(std::string(what) + " has a row number outside geno");
    }
    out[static_cast<size_t>(i)] = r - 1;
  }
  return out;
}

struct Problem {
  Map<const MatrixXd> geno;
  Criterion crit;
  double lambda;
  bool legacy;
  std::vector<int> cand;
  bool fixed_test = false;
  Mat x0_fixed;
  TestCache cache;
  bool use_complement = false;  // untargeted r-score/PEV: test cross-product = total - training
  std::vector<char> in_sol;     // by row of geno
  // Dual form (more columns than training rows): Gram matrix of all rows, built once.
  bool use_gram = false;
  std::shared_ptr<const Mat> gram;  // shared by the copies used in parallel searches
  std::vector<int> test_idx;    // fixed test rows (0-based), if any
  double gram_lambda = 0.0;

  // Uses a Gram matrix of all rows that the caller built once (e.g. with tcrossprod in R).
  void SetGram(const Mat& g) {
    if (crit > kCd || IsSetOnly(crit)) return;
    gram = std::make_shared<const Mat>(g);
    gram_lambda = lambda > 0 ? lambda : (crit == kCd ? 1.0 : 1.0 / static_cast<double>(geno.cols()));
    use_gram = true;
  }

  // Builds the Gram matrix when the dual form will be used for training sets of n rows.
  void EnableGram(int n_train) {
    if (crit > kCd || IsSetOnly(crit)) return;
    if (!UseDual(crit, n_train, geno.cols())) return;
    if (geno.rows() > 5000) return;  // the Gram matrix would take over 200 MB
    gram = std::make_shared<const Mat>(Mat(geno * geno.transpose()));
    gram_lambda = lambda > 0 ? lambda : (crit == kCd ? 1.0 : 1.0 / static_cast<double>(geno.cols()));
    use_gram = true;
  }

  Problem(const Rcpp::NumericMatrix& g, Criterion c, double lam, bool leg)
      : geno(&g[0], g.nrow(), g.ncol()), crit(c), lambda(lam), legacy(leg) {}

  // Score of the training set `sol`; `in_sol` flags the rows of `sol` (used when the test set is the
  // candidates left out). A score that is not finite stops the search: it happens when the test
  // set has fewer than two individuals, or when the data are degenerate.
  double Evaluate(const std::vector<int>& sol, const std::vector<char>& in_sol) const {
    const double v = EvaluateRaw(sol, in_sol);
    if (!std::isfinite(v)) {
      throw std::runtime_error(
          "the score is not finite; check that the test set (or the candidates left out of the "
          "training set) has at least 2 individuals and that geno has no constant rows or columns");
    }
    return v;
  }
  double Evaluate(const std::vector<int>& sol) const { return Evaluate(sol, in_sol); }

  double EvaluateRaw(const std::vector<int>& sol, const std::vector<char>& in_sol) const {
    if (use_gram && UseDual(crit, static_cast<Index>(sol.size()), geno.cols())) {
      if (fixed_test) return ScoreGramDual(crit, *gram, sol, test_idx, gram_lambda, legacy);
      std::vector<int> rest;
      rest.reserve(cand.size());
      for (int row : cand)
        if (!in_sol[static_cast<size_t>(row)]) rest.push_back(row);
      if (rest.empty()) throw std::runtime_error("no candidates are left outside the training set");
      return ScoreGramDual(crit, *gram, sol, rest, gram_lambda, legacy);
    }
    const Mat x = GatherRows(geno, sol);
    if (IsSetOnly(crit)) return Score(crit, x, nullptr, 0, nullptr, lambda, legacy);
    if (fixed_test) return Score(crit, x, &x0_fixed, x0_fixed.rows(), &cache, lambda, legacy);
    const Index n0 = static_cast<Index>(cand.size()) - x.rows();
    if (n0 < 1) throw std::runtime_error("no candidates are left outside the training set");
    if (use_complement && CanUseComplement(crit, x.rows(), x.cols())) {
      return Score(crit, x, nullptr, n0, &cache, lambda, legacy);
    }
    std::vector<int> rest;
    rest.reserve(cand.size());
    for (int row : cand)
      if (!in_sol[static_cast<size_t>(row)]) rest.push_back(row);
    if (rest.empty()) throw std::runtime_error("no candidates are left outside the training set");
    const Mat x0 = GatherRows(geno, rest);
    return Score(crit, x, &x0, x0.rows(), &untargeted_cache, lambda, legacy);
  }
  TestCache untargeted_cache = MakeUntargetedCache();
  static TestCache MakeUntargetedCache() {
    TestCache c;
    c.untargeted = true;
    return c;
  }
};

}  // namespace

// [[Rcpp::export]]
double score_matrices(int crit, Rcpp::NumericMatrix x, Rcpp::NumericMatrix x0, double lambda, bool legacy) {
  const Map<const MatrixXd> a(&x[0], x.nrow(), x.ncol());
  const Map<const MatrixXd> b(&x0[0], x0.nrow(), x0.ncol());
  const Mat xa(a), xb(b);
  return Score(static_cast<Criterion>(crit), xa, &xb, xb.rows(), nullptr, lambda, legacy);
}

// [[Rcpp::export]]
double score_rows(int crit, Rcpp::NumericMatrix geno, Rcpp::IntegerVector train,
                  Rcpp::IntegerVector test, double lambda, bool legacy) {
  const Map<const MatrixXd> g(&geno[0], geno.nrow(), geno.ncol());
  const Mat x = GatherRows(g, ToZeroBased(train, "train", geno.nrow()));
  const Mat x0 = GatherRows(g, ToZeroBased(test, "test", geno.nrow()));
  return Score(static_cast<Criterion>(crit), x, &x0, x0.rows(), nullptr, lambda, legacy);
}

// Scores of the columns of `sets` (each a training set, 1-based rows of geno) against one fixed
// test set. The test-set terms and, for more columns than rows, the Gram matrix are built once,
// and the sets can be scored by several threads. Same values as score_rows() up to rounding.
// [[Rcpp::export]]
Rcpp::NumericVector score_many(int crit, Rcpp::NumericMatrix geno, Rcpp::IntegerMatrix sets,
                               Rcpp::IntegerVector test, double lambda, bool legacy, int n_threads,
                               Rcpp::Nullable<Rcpp::NumericMatrix> gram) {
  Problem pr(geno, static_cast<Criterion>(crit), lambda, legacy);
  const int n_all = geno.nrow();
  const int k = sets.ncol();
  pr.fixed_test = true;
  pr.test_idx = ToZeroBased(test, "test", n_all);
  pr.x0_fixed = GatherRows(pr.geno, pr.test_idx);
  pr.cache = MakeTestCache(pr.crit, pr.x0_fixed);
  if (gram.isNotNull()) {
    const Rcpp::NumericMatrix gm(gram.get());
    pr.SetGram(Eigen::Map<const MatrixXd>(&gm[0], gm.nrow(), gm.ncol()));
  } else {
    pr.EnableGram(sets.nrow());
  }
  std::vector<std::vector<int>> idx(static_cast<size_t>(k));
  for (int j = 0; j < k; ++j) {
    idx[static_cast<size_t>(j)] = ToZeroBased(Rcpp::IntegerVector(sets(Rcpp::_, j)), "sol", n_all);
  }
  std::vector<double> out(static_cast<size_t>(k), NA_REAL);
  std::string error;
  int threads = n_threads < 1 ? 1 : n_threads;
#ifndef _OPENMP
  threads = 1;
#endif
#ifdef _OPENMP
#pragma omp parallel for num_threads(threads) schedule(static)
#endif
  for (int j = 0; j < k; ++j) {
    try {
      out[static_cast<size_t>(j)] = pr.Evaluate(idx[static_cast<size_t>(j)]);
    } catch (const std::exception& e) {
#ifdef _OPENMP
#pragma omp critical(tsdfgs_error_many)
#endif
      error = e.what();
    }
  }
  if (!error.empty()) throw std::runtime_error(error);
  return Rcpp::NumericVector(out.begin(), out.end());
}

// Hill climbing from the initial set `sol0`.
//   test given:    the criterion compares the training set with this fixed test set.
//   test NULL:     the test set is the candidates outside the training set.
//   group given:   integer group (sub-population) of each candidate, 1..n_groups; a move
//                  replaces a member by a candidate of the same, randomly chosen group.
// [[Rcpp::export]]
Rcpp::List exchange_search(Rcpp::NumericMatrix geno, Rcpp::IntegerVector cand,
                           Rcpp::IntegerVector sol0, Rcpp::Nullable<Rcpp::IntegerVector> test,
                           Rcpp::Nullable<Rcpp::IntegerVector> group, int n_groups, int crit,
                           bool maximize, double min_iter, double tol, double window,
                           double lambda, bool legacy, bool console) {
  if (!std::isfinite(min_iter) || min_iter < 1) throw std::runtime_error("min_iter must be a finite number of at least 1");
  Problem pr(geno, static_cast<Criterion>(crit), lambda, legacy);
  const int n_all = geno.nrow();
  pr.cand = ToZeroBased(cand, "cand", n_all);
  std::vector<int> sol = ToZeroBased(sol0, "sol", n_all);
  pr.in_sol.assign(static_cast<size_t>(n_all), 0);
  for (int row : sol) pr.in_sol[static_cast<size_t>(row)] = 1;
  if (IsSetOnly(pr.crit)) {
    // the training set alone defines the score; no test cache is needed
  } else if (test.isNotNull()) {
    pr.fixed_test = true;
    pr.test_idx = ToZeroBased(Rcpp::IntegerVector(test.get()), "test", n_all);
    pr.x0_fixed = GatherRows(pr.geno, pr.test_idx);
    pr.cache = MakeTestCache(pr.crit, pr.x0_fixed);
  } else {
    if (pr.crit != kCd) pr.cache = MakeCandidateCache(pr.crit, GatherRows(pr.geno, pr.cand));
    pr.cache.untargeted = true;
    pr.use_complement = true;
  }
  pr.EnableGram(static_cast<int>(sol.size()));
  const bool use_groups = group.isNotNull();
  std::vector<int> gpos, row_pos;
  if (use_groups) {
    Rcpp::IntegerVector g(group.get());
    if (static_cast<size_t>(g.size()) != pr.cand.size()) {
      throw std::runtime_error("group must have one entry per candidate");
    }
    gpos.assign(g.begin(), g.end());
    row_pos.assign(static_cast<size_t>(n_all), -1);
    for (size_t i = 0; i < pr.cand.size(); ++i) row_pos[static_cast<size_t>(pr.cand[i])] = static_cast<int>(i);
  }
  const int n = static_cast<int>(sol.size());

  std::vector<double> top, iter_score;
  double score = pr.Evaluate(sol);
  top.push_back(score);
  iter_score.push_back(score);
  std::vector<int> next, pool, slots;
  int iter = 1, status = 0;
  while (status == 0) {
    next = sol;
    if (!use_groups) {
      pool.clear();
      for (int row : pr.cand)
        if (!pr.in_sol[static_cast<size_t>(row)]) pool.push_back(row);
      if (pool.empty()) throw std::runtime_error("no candidates are left outside the training set");
      const int value = pool[static_cast<size_t>(UnifIndex(static_cast<int>(pool.size())))];
      next[static_cast<size_t>(UnifIndex(n))] = value;
    } else {
      const int g = UnifIndex(n_groups);
      pool.clear();
      slots.clear();
      for (size_t i = 0; i < pr.cand.size(); ++i)
        if (gpos[i] == g + 1 && !pr.in_sol[static_cast<size_t>(pr.cand[i])]) pool.push_back(pr.cand[i]);
      for (int k = 0; k < n; ++k)
        if (gpos[static_cast<size_t>(row_pos[static_cast<size_t>(sol[static_cast<size_t>(k)])])] == g + 1) slots.push_back(k);
      if (!pool.empty() && !slots.empty()) {
        const int value = pool[static_cast<size_t>(UnifIndex(static_cast<int>(pool.size())))];
        next[static_cast<size_t>(slots[static_cast<size_t>(UnifIndex(static_cast<int>(slots.size())))])] = value;
      }
    }
    // The complement of the evaluated set is the test set when there is no fixed test set.
    for (int row : sol) pr.in_sol[static_cast<size_t>(row)] = 0;
    for (int row : next) pr.in_sol[static_cast<size_t>(row)] = 1;
    const double new_score = pr.Evaluate(next);
    const bool better = maximize ? new_score > score : new_score < score;
    if (better) {
      sol = next;
      score = new_score;
    }
    for (int row : next) pr.in_sol[static_cast<size_t>(row)] = 0;
    for (int row : sol) pr.in_sol[static_cast<size_t>(row)] = 1;
    iter_score.push_back(new_score);
    top.push_back(score);
    if (console) Rcpp::Rcout << iter << "..";
    ++iter;
    if (iter > min_iter) {
      const int back = static_cast<int>(std::trunc(iter - window));
      if (back >= 1 && std::fabs(top[static_cast<size_t>(iter - 1)] - top[static_cast<size_t>(back - 1)]) < tol) status = 1;
    }
    if (status == 0 && iter > min_iter * 2) status = 2;
    if ((iter & 255) == 0) Rcpp::checkUserInterrupt();
  }

  Rcpp::IntegerVector sol_r(sol.size());
  for (size_t i = 0; i < sol.size(); ++i) sol_r[static_cast<R_xlen_t>(i)] = sol[i] + 1;
  return Rcpp::List::create(
      Rcpp::Named("sol") = sol_r,
      Rcpp::Named("top") = Rcpp::NumericVector(top.begin() + 1, top.end()),
      Rcpp::Named("iter") = Rcpp::NumericVector(iter_score.begin() + 1, iter_score.end()),
      Rcpp::Named("n.iter") = iter - 1,
      Rcpp::Named("converged") = status == 1);
}

// ---- Independent exchange searches in parallel (OpenMP) ----
// One search per element of `sols` (untargeted: the candidates left out of the training set are
// the test set; or against a fixed `test`). Each search has its own generator, seeded from
// `seeds`, so the result of search i does not depend on the number of threads or on the
// order in which the threads finish. These random numbers are not those of R's generator, so a
// search here does not return the training set of optTrain() with the same seed.
namespace {

inline int Draw(std::mt19937_64& rng, int n) {
  const double u = static_cast<double>(rng() >> 11) * (1.0 / 9007199254740992.0);
  const int k = static_cast<int>(u * n);
  return k < n ? k : n - 1;
}

// Hill climbing without sub-populations. Returns the best score; `best_sol` gets the set.
double ExchangeOnce(const Problem& pr, std::vector<int> sol, bool maximize, double min_iter, double tol,
                    double window, std::uint64_t seed, int* n_iter, std::vector<int>* best_sol,
                    bool* converged) {
  std::mt19937_64 rng(seed);
  const int n = static_cast<int>(sol.size());
  std::vector<char> in_sol(pr.in_sol.size(), 0);
  for (int row : sol) in_sol[static_cast<size_t>(row)] = 1;
  std::vector<double> top;
  double score = pr.Evaluate(sol, in_sol);
  top.push_back(score);
  std::vector<int> next, pool;
  int iter = 1, status = 0;
  while (status == 0) {
    next = sol;
    pool.clear();
    for (int row : pr.cand)
      if (!in_sol[static_cast<size_t>(row)]) pool.push_back(row);
    if (pool.empty()) throw std::runtime_error("no candidates are left outside the training set");
    const int value = pool[static_cast<size_t>(Draw(rng, static_cast<int>(pool.size())))];
    next[static_cast<size_t>(Draw(rng, n))] = value;
    for (int row : sol) in_sol[static_cast<size_t>(row)] = 0;
    for (int row : next) in_sol[static_cast<size_t>(row)] = 1;
    const double new_score = pr.Evaluate(next, in_sol);
    const bool better = maximize ? new_score > score : new_score < score;
    if (better) {
      sol = next;
      score = new_score;
    }
    for (int row : next) in_sol[static_cast<size_t>(row)] = 0;
    for (int row : sol) in_sol[static_cast<size_t>(row)] = 1;
    top.push_back(score);
    ++iter;
    if (iter > min_iter) {
      const int back = static_cast<int>(std::trunc(iter - window));
      if (back >= 1 && std::fabs(top[static_cast<size_t>(iter - 1)] - top[static_cast<size_t>(back - 1)]) < tol) status = 1;
    }
    if (status == 0 && iter > min_iter * 2) status = 2;
  }
  *n_iter = iter - 1;
  *converged = status == 1;
  *best_sol = sol;
  return score;
}

}  // namespace

// [[Rcpp::export]]
Rcpp::List exchange_search_multi(Rcpp::NumericMatrix geno, Rcpp::IntegerVector cand,
                                 Rcpp::List sols, Rcpp::Nullable<Rcpp::IntegerVector> test,
                                 int crit, bool maximize, Rcpp::NumericVector min_iter, double tol,
                                 Rcpp::NumericVector window, double lambda, bool legacy,
                                 Rcpp::NumericVector seeds, int n_threads) {
  const R_xlen_t n_rep = sols.size();
  if (min_iter.size() != n_rep || window.size() != n_rep || seeds.size() != n_rep) {
    throw std::runtime_error("min_iter, window, and seeds need one entry per search");
  }
  for (R_xlen_t i = 0; i < n_rep; ++i) {
    if (!std::isfinite(min_iter[i]) || min_iter[i] < 1) throw std::runtime_error("min_iter must be a finite number of at least 1");
  }
  Problem base(geno, static_cast<Criterion>(crit), lambda, legacy);
  const int n_all = geno.nrow();
  base.cand = ToZeroBased(cand, "cand", n_all);
  base.in_sol.assign(static_cast<size_t>(n_all), 0);
  if (IsSetOnly(base.crit)) {
  } else if (test.isNotNull()) {
    base.fixed_test = true;
    base.test_idx = ToZeroBased(Rcpp::IntegerVector(test.get()), "test", n_all);
    base.x0_fixed = GatherRows(base.geno, base.test_idx);
    base.cache = MakeTestCache(base.crit, base.x0_fixed);
  } else {
    if (base.crit != kCd) base.cache = MakeCandidateCache(base.crit, GatherRows(base.geno, base.cand));
    base.cache.untargeted = true;
    base.use_complement = true;
  }
  if (n_rep > 0) {
    int smallest = n_all;
    for (R_xlen_t i = 0; i < n_rep; ++i) smallest = std::min(smallest, static_cast<int>(Rcpp::IntegerVector(sols[i]).size()));
    base.EnableGram(smallest);
  }
  std::vector<std::vector<int>> start(static_cast<size_t>(n_rep));
  for (R_xlen_t i = 0; i < n_rep; ++i) {
    start[static_cast<size_t>(i)] = ToZeroBased(Rcpp::IntegerVector(sols[i]), "sol", n_all);
  }
  std::vector<double> best(static_cast<size_t>(n_rep), NA_REAL);
  std::vector<int> iters(static_cast<size_t>(n_rep), 0);
  std::vector<char> conv(static_cast<size_t>(n_rep), 0);
  std::vector<std::vector<int>> best_sol(static_cast<size_t>(n_rep));
  std::vector<std::uint64_t> sd(static_cast<size_t>(n_rep));
  for (R_xlen_t i = 0; i < n_rep; ++i) sd[static_cast<size_t>(i)] = static_cast<std::uint64_t>(seeds[i]);
  std::string error;
  int threads = n_threads < 1 ? 1 : n_threads;
#ifndef _OPENMP
  threads = 1;
#endif
  // No Rcpp or R API call inside the parallel loop; errors are collected and rethrown after it.
#ifdef _OPENMP
#pragma omp parallel for num_threads(threads) schedule(dynamic, 1)
#endif
  for (R_xlen_t i = 0; i < n_rep; ++i) {
    try {
      const size_t k = static_cast<size_t>(i);
      bool c = false;
      best[k] = ExchangeOnce(base, start[k], maximize, min_iter[i], tol, window[i], sd[k], &iters[k],
                             &best_sol[k], &c);
      conv[k] = c;
    } catch (const std::exception& e) {
#ifdef _OPENMP
#pragma omp critical(tsdfgs_error)
#endif
      error = e.what();
    }
  }
  if (!error.empty()) throw std::runtime_error(error);
  Rcpp::NumericVector score(n_rep);
  Rcpp::IntegerVector n_iter(n_rep);
  Rcpp::LogicalVector converged(n_rep);
  Rcpp::List sol_out(n_rep);
  for (R_xlen_t i = 0; i < n_rep; ++i) {
    const size_t k = static_cast<size_t>(i);
    score[i] = best[k];
    n_iter[i] = iters[k];
    converged[i] = conv[k] != 0;
    Rcpp::IntegerVector s(best_sol[k].size());
    for (size_t j = 0; j < best_sol[k].size(); ++j) s[static_cast<R_xlen_t>(j)] = best_sol[k][j] + 1;
    sol_out[i] = s;
  }
  return Rcpp::List::create(Rcpp::Named("score") = score, Rcpp::Named("n.iter") = n_iter,
                            Rcpp::Named("converged") = converged, Rcpp::Named("sol") = sol_out,
                            Rcpp::Named("threads") = threads);
}

// Number of threads Eigen may use inside its own products (1 keeps every call single-threaded
// unless exchange_search_multi() runs searches in parallel).
// [[Rcpp::export]]
void set_eigen_threads(int n) {
  Eigen::setNbThreads(n < 1 ? 1 : n);
}

// Number of threads that may be used by default: the processors available to this process,
// limited by OMP_NUM_THREADS and OMP_THREAD_LIMIT when they are set; 1 without OpenMP.
// [[Rcpp::export]]
int n_procs() {
#ifdef _OPENMP
  return std::max(1, std::min(omp_get_num_procs(), omp_get_max_threads()));
#else
  return 1;
#endif
}

// [[Rcpp::export]]
bool has_openmp() {
#ifdef _OPENMP
  return true;
#else
  return false;
#endif
}
