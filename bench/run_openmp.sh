#!/bin/bash
# Run with Apptainer (env/TSDFGS.sif): installs the current source into a private library and benchmarks the OpenMP
# simulation (up to 50 threads). Usage: bench/run_openmp.sh <work_dir> [n_iter] [n_pcs] [openmp_only]
set -eu
d=$1
[ -d "$d" ] || { echo "d: directory not found: $d" >&2; exit 2; }
proj=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
export APPTAINER_QUIET=1
export BENCH_NITER=${2:-10} BENCH_NPC=${3:-50} BENCH_MODE=${4:-all}
apptainer exec --bind "$d:/scr" --bind "$proj:/proj:ro" --pwd /scr "$proj/env/TSDFGS.sif" bash -c '
  set -e; export TMPDIR=/scr/tmp OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1; mkdir -p $TMPDIR /scr/lib_new
  rm -rf /scr/new_src /scr/lib_new/TSDFGS; mkdir /scr/new_src
  (cd /proj && tar --exclude=env/*.sif --exclude=.git --exclude=docs --exclude=site --exclude=src/*.o --exclude=src/*.so -cf - .) | tar -xf - -C /scr/new_src
  R CMD INSTALL --no-docs --library=/scr/lib_new /scr/new_src > /scr/install_new.log 2>&1
  Rscript /proj/bench/bench_openmp.R /scr/lib_new /scr/bench_openmp.tsv $BENCH_NITER $BENCH_NPC $BENCH_MODE'
echo BENCH_DONE
