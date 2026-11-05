# Benchmark of the parallel simulation of SSDFGS(): serial, parallel::mclapply (the method of
# v2.4.2), and OpenMP (threads = 1, 2, 5, 10, 25, 50).
# Usage: Rscript bench/bench_openmp.R <library> <out.tsv> [n_iter] [n_pcs] [openmp_only]
args <- commandArgs(TRUE)
library(TSDFGS, lib.loc = args[1])
out <- args[2]
n_iter <- if (length(args) >= 3) as.integer(args[3]) else 10L
data(geno)
npc <- if (length(args) >= 4) as.integer(args[4]) else 50L
g <- as.matrix(geno)[, seq_len(npc)]
nt <- round(seq(nrow(g) / 16, nrow(g) / 2, length.out = 10))
tasks <- rep(nt, each = n_iter)
cat("rows:", nrow(g), "columns:", ncol(g), "searches:", length(tasks), "\n")
cat("has_openmp:", TSDFGS:::has_openmp(), "\n")

openmp_only <- length(args) >= 5 && args[5] == "openmp_only"
one <- function(k) max(optTrain(g, seq_len(nrow(g)), k, console = FALSE)$TOPscore)
res <- data.frame()
add <- function(method, threads, expr) {
  t <- system.time(r <- expr)[["elapsed"]]
  res <<- rbind(res, data.frame(method = method, threads = threads, seconds = t,
                                mean_r = mean(r), n = length(r)))
  cat(sprintf("%-8s threads=%2d  %8.1f s  mean r=%.4f\n", method, threads, t, mean(r)))
  write.table(res, out, sep = "\t", quote = FALSE, row.names = FALSE)
  invisible(r)
}
if (!openmp_only) {
  set.seed(1); add("serial", 1L, vapply(tasks, one, 0))
}
for (th in if (openmp_only) integer() else c(10L, 50L)) {
  set.seed(1); add("mclapply", th, unlist(parallel::mclapply(tasks, one, mc.cores = th)))
}
ref <- NULL
for (th in if (openmp_only) c(10L, 25L, 50L) else c(1L, 2L, 5L, 10L, 25L, 50L)) {
  set.seed(1)
  r <- add("openmp", th, TSDFGS:::.simulate_rscore(g, tasks, th))
  if (is.null(ref)) ref <- r else cat("  identical to threads=1:", identical(r, ref), "\n")
}
