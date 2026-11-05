#!/bin/bash
# Build the source package, run R CMD check --as-cran (no manual: the image has no LaTeX),
# inside env/TSDFGS.sif (built with env/build_sif.sh). The pkgdown site is built by GitHub Actions.
# Usage: env/check_pkg.sh <work_dir>  (an empty, writable directory, for example from mktemp -d)
# Results: <work_dir>/out/TSDFGS.Rcheck/00check.log
set -eu
d=$1
[ -d "$d" ] || { echo "d: directory not found: $d" >&2; exit 2; }
proj=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
export APPTAINER_QUIET=1
mkdir -p "$d/out" "$d/tmp"
apptainer exec --bind "$d:/scr" --bind "$proj:/proj:ro" --pwd /scr/out "$proj/env/TSDFGS.sif" bash -c '
  set -e; export TMPDIR=/scr/tmp
  rm -rf /scr/src_copy; mkdir /scr/src_copy
  (cd /proj && tar --exclude=env/*.sif --exclude=.git --exclude=src/*.o --exclude=src/*.so -cf - .) | tar -xf - -C /scr/src_copy
  R CMD build --no-manual /scr/src_copy
  _R_CHECK_CRAN_INCOMING_REMOTE_=false R CMD check --as-cran --no-manual TSDFGS_*.tar.gz || true'
echo CHECK_DONE
