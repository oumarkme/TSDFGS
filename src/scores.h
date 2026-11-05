#ifndef TSDFGS_SCORES_H_
#define TSDFGS_SCORES_H_

#include <vector>

#include "linalg.h"

namespace tsdfgs {

// 0-2: criteria of v2.4.2. 3-17: added in v3 (one page per criterion on the package website).
enum Criterion {
  kRScore = 0, kPev = 1, kCd = 2,
  kAOpt = 3,           // A-optimality (set only)
  kDOpt = 4,           // D-optimality (set only)
  kEOpt = 5,           // E-optimality (set only)
  kPevMax = 6,         // largest prediction error variance of the target
  kGOpt = 7,           // largest eigenvalue of the PEV matrix of the target (G-optimality)
  kMspe = 8,           // MSPE of ridge regression (Chen et al., 2024)
  kGvOverall = 9,      // log det of the relationship matrix of the training set (set only)
  kGvAverage = 10,     // trace of the relationship matrix of the training set (set only)
  kAvgGrm = 11,        // mean relationship of the training set to the target
  kAvgGrmSelf = 12,    // mean relationship within the training set (set only)
  kAvgGrmMinMax = 13,  // relationship to the target minus relationship within the training set
  kCdMean = 14,        // CDmean, GBLUP form with contrasts (Rincent et al., 2012)
  kCdMeanV2 = 15,      // CDmean(v2) (Chen et al., 2024)
  kPam = 16,           // sum of distances of the target to the nearest training individual
  kCdMax = 17          // largest CD (ridge, PCs) among the target individuals
};

// True if the criterion uses the training set only (the test set, if any, is ignored).
bool IsSetOnly(Criterion crit);

// Criteria that treat the whole candidate set as the target when no test set is given.
// (The criteria of v2.4.2 use the candidates left out of the training set.)
bool TargetIsCandidates(Criterion crit);

// Optional per-test-set quantities that do not change during a search with a fixed test set.
struct TestCache {
  Mat g0;     // X0'X0 (r-score), or [1 X0]'[1 X0] (PEV) of a fixed test set; empty if not used
  Mat total;  // the same sum over all candidates; the test set is then "candidates - training set"
  bool untargeted = false;  // no test set: x0 holds the candidates left out of the training set
};

// Cache for a fixed test matrix X0 (empty cache if nothing to reuse).
TestCache MakeTestCache(Criterion crit, const Mat& x0);

// Cache for the untargeted search: the cross-product of all candidates. The complement of the
// training set then has the cross-product total - X'X, so its rows need not be copied.
TestCache MakeCandidateCache(Criterion crit, const Mat& candidates);

// True if Score() can run with x0 = nullptr for this criterion and training set size.
bool CanUseComplement(Criterion crit, Index n, Index p);

// Scores of the training matrix X against the test matrix X0 (same columns); n0 is the number
// of test rows. x0 may be nullptr only when cache->total is set and CanUseComplement() is true.
// lambda <= 0 uses the default penalty of the criterion.
// legacy = true reproduces the single-precision values of v2.4.2.
double Score(Criterion crit, const Mat& x, const Mat* x0, Index n0, const TestCache* cache,
             double lambda, bool legacy);

// Dual form of the r-score, PEV, and CD (more columns than training rows) from the Gram matrix
// of all rows of geno: the training and test rows are index sets, so no column of geno is
// touched during a search. Same values as Score() up to rounding.
double ScoreGramDual(Criterion crit, const Mat& gram, const std::vector<int>& train,
                     const std::vector<int>& test, double lambda, bool legacy);

// True if the dual form applies to a training set of n rows and p columns.
bool UseDual(Criterion crit, Index n, Index p);

// Criteria added in v3 (src/criteria_v3.cpp). lambda <= 0 selects the default of the criterion.
double ScoreV3(Criterion crit, const Mat& x, const Mat* x0, Index n0, const TestCache* cache,
               double lambda);

}  // namespace tsdfgs

#endif  // TSDFGS_SCORES_H_
