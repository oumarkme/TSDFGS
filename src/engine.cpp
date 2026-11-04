// Exchange (hill-climbing) search for the training set, and score entry points for R.
//
// The moves and the order of random draws are those of the R code of v2.4.2, drawn from
// R's own generator (R_unif_index, as used by sample.int), so that set.seed() gives the
// same training set. Rcpp saves and restores R's random-number state around each call.
#include <R_ext/Random.h>
#include <Rcpp.h>

#include <algorithm>
#include <cmath>
#include <string>
#include <vector>

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
  bool legacy;
  std::vector<int> cand;
  bool fixed_test = false;
  Mat x0_fixed;
  TestCache cache;
  bool use_complement = false;  // untargeted r-score/PEV: test cross-product = total - training
  std::vector<char> in_sol;     // by row of geno

  Problem(const Rcpp::NumericMatrix& g, Criterion c, bool leg)
      : geno(&g[0], g.nrow(), g.ncol()), crit(c), legacy(leg) {}

  double Evaluate(const std::vector<int>& sol) const {
    const Mat x = GatherRows(geno, sol);
    if (fixed_test) return Score(crit, x, &x0_fixed, x0_fixed.rows(), &cache, legacy);
    const Index n0 = static_cast<Index>(cand.size()) - x.rows();
    if (n0 < 1) throw std::runtime_error("no candidates are left outside the training set");
    if (use_complement && CanUseComplement(crit, x.rows(), x.cols())) {
      return Score(crit, x, nullptr, n0, &cache, legacy);
    }
    std::vector<int> rest;
    rest.reserve(cand.size());
    for (int row : cand)
      if (!in_sol[static_cast<size_t>(row)]) rest.push_back(row);
    if (rest.empty()) throw std::runtime_error("no candidates are left outside the training set");
    const Mat x0 = GatherRows(geno, rest);
    return Score(crit, x, &x0, x0.rows(), nullptr, legacy);
  }
};

}  // namespace

// [[Rcpp::export]]
double score_matrices(int crit, Rcpp::NumericMatrix x, Rcpp::NumericMatrix x0, bool legacy) {
  const Map<const MatrixXd> a(&x[0], x.nrow(), x.ncol());
  const Map<const MatrixXd> b(&x0[0], x0.nrow(), x0.ncol());
  const Mat xa(a), xb(b);
  return Score(static_cast<Criterion>(crit), xa, &xb, xb.rows(), nullptr, legacy);
}

// [[Rcpp::export]]
double score_rows(int crit, Rcpp::NumericMatrix geno, Rcpp::IntegerVector train,
                  Rcpp::IntegerVector test, bool legacy) {
  const Map<const MatrixXd> g(&geno[0], geno.nrow(), geno.ncol());
  const Mat x = GatherRows(g, ToZeroBased(train, "train", geno.nrow()));
  const Mat x0 = GatherRows(g, ToZeroBased(test, "test", geno.nrow()));
  return Score(static_cast<Criterion>(crit), x, &x0, x0.rows(), nullptr, legacy);
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
                           bool legacy, bool console) {
  Problem pr(geno, static_cast<Criterion>(crit), legacy);
  const int n_all = geno.nrow();
  pr.cand = ToZeroBased(cand, "cand", n_all);
  std::vector<int> sol = ToZeroBased(sol0, "sol", n_all);
  pr.in_sol.assign(static_cast<size_t>(n_all), 0);
  for (int row : sol) pr.in_sol[static_cast<size_t>(row)] = 1;
  if (test.isNotNull()) {
    pr.fixed_test = true;
    pr.x0_fixed = GatherRows(pr.geno, ToZeroBased(Rcpp::IntegerVector(test.get()), "test", n_all));
    pr.cache = MakeTestCache(pr.crit, pr.x0_fixed);
  } else if (pr.crit != kCd) {
    pr.cache = MakeCandidateCache(pr.crit, GatherRows(pr.geno, pr.cand));
    pr.use_complement = true;
  }
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
