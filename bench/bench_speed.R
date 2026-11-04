# Speed of optTrain(): v2.4.2 versus the current branch, same seed, 1 thread.
# Usage (inside env/TSDFGS.sif): Rscript bench/bench_speed.R <lib_old> <lib_new> <out.tsv>
args <- commandArgs(TRUE)
libs <- c(old = args[1], new = args[2])
out <- args[3]
res <- list()
for (v in names(libs)) {
  suppressMessages(library(TSDFGS, lib.loc = libs[[v]]))
  data(geno); g <- as.matrix(geno)
  for (p in c(10, 50, 404)) {
    gg <- g[, 1:p]
    # v2.4.2 fails when min.iter < 5 * sqrt(candidates * n.train); stay just above it.
    mi <- function(nc, n) ceiling(5 * sqrt(nc * n)) + 20
    cases <- list(
      list("r-score, untargeted, n=100",  function() optTrain(gg, 1:404, 100, method = "rScore", min.iter = mi(404, 100), console = FALSE)),
      list("PEV, untargeted, n=100",      function() optTrain(gg, 1:404, 100, method = "PEV", min.iter = mi(404, 100), console = FALSE)),
      list("CD, targeted, n=50, test=300", function() optTrain(gg, 1:100, 50, test = 101:400, method = "CD", min.iter = mi(100, 50), console = FALSE)))
    for (cs in cases) {
      set.seed(1)
      t <- system.time(r <- cs[[2]]())[["elapsed"]]
      n.it <- length(r$TOPscore)
      res[[length(res) + 1]] <- data.frame(version = v, p = p, case = cs[[1]], iterations = n.it,
                                           seconds = t, ms_per_iter = 1000 * t / n.it, final = tail(r$TOPscore, 1))
    }
  }
  detach("package:TSDFGS", unload = TRUE, character.only = TRUE)
}
write.table(do.call(rbind, res), out, sep = "\t", quote = FALSE, row.names = FALSE)
