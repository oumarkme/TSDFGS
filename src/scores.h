#ifndef TSDFGS_SCORES_H_
#define TSDFGS_SCORES_H_

#include "linalg.h"

namespace tsdfgs {

enum Criterion { kRScore = 0, kPev = 1, kCd = 2 };

// Optional per-test-set quantities that do not change during a search with a fixed test set.
struct TestCache {
  Mat g0;     // X0'X0 (r-score), or [1 X0]'[1 X0] (PEV) of a fixed test set; empty if not used
  Mat total;  // the same sum over all candidates; the test set is then "candidates - training set"
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
// legacy = true reproduces the single-precision values of v2.4.2.
double Score(Criterion crit, const Mat& x, const Mat* x0, Index n0, const TestCache* cache,
             bool legacy);

}  // namespace tsdfgs

#endif  // TSDFGS_SCORES_H_
