#!/bin/bash
# Run with Apptainer (env/TSDFGS.sif): installs v2.4.2 and the current source into a private library, then benchmarks.
# Usage: bench/run_bench.sh <work_dir>  (an empty, writable directory)
# The v2.4.2 source must be unpacked in <work_dir>/old_src (git archive a802b81).
set -eu
d=$1
[ -d "$d" ] || { echo "d: directory not found: $d" >&2; exit 2; }
proj=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
export APPTAINER_QUIET=1
apptainer exec --bind "$d:/scr" --bind "$proj:/proj" --pwd /scr "$proj/env/TSDFGS.sif" bash -c '
  set -e; export TMPDIR=/scr/tmp; mkdir -p $TMPDIR /scr/lib_old /scr/lib_new
  rm -f /scr/old_src/src/*.o /scr/old_src/src/*.so
  [ -d /scr/lib_old/TSDFGS ] || R CMD INSTALL --no-docs --library=/scr/lib_old /scr/old_src > /scr/install_old.log 2>&1
  rm -rf /scr/new_src /scr/lib_new/TSDFGS; mkdir /scr/new_src; (cd /proj && tar --exclude=env/*.sif --exclude=.git --exclude=docs --exclude=src/*.o --exclude=src/*.so -cf - .) | tar -xf - -C /scr/new_src
  R CMD INSTALL --no-docs --library=/scr/lib_new /scr/new_src > /scr/install_new.log 2>&1
  Rscript /proj/bench/bench_speed.R /scr/lib_old /scr/lib_new /scr/bench_speed.tsv'
echo BENCH_DONE
