#!/bin/bash
# Build the source package, run R CMD check --as-cran (no manual: the image has no LaTeX),
# inside env/TSDFGS.sif. Run on gs4. The pkgdown site is built on wsl (gs4 cannot reach the web).
# Usage: env/check_pkg.sh <scratch_dir>  (made with mktemp -d /raid/tmp/TSDFGS.XXXXXX)
# Results: <scratch_dir>/out/TSDFGS.Rcheck/00check.log
set -eu
export PATH=/home/user/miniforge-pypy3/envs/mark/bin:$PATH
d=$1
case "$d" in /raid/tmp/TSDFGS.*) ;; *) echo "bad scratch $d" >&2; exit 2;; esac
proj=/lab_nas_home_user/Mark/TSDFGS
export APPTAINER_QUIET=1
mkdir -p "$d/out" "$d/tmp"
apptainer exec --bind "$d:/scr" --bind "$proj:/proj:ro" --pwd /scr/out "$proj/env/TSDFGS.sif" bash -c '
  set -e; export TMPDIR=/scr/tmp
  rm -rf /scr/src_copy; mkdir /scr/src_copy
  (cd /proj && tar --exclude=env/*.sif --exclude=.git --exclude=src/*.o --exclude=src/*.so -cf - .) | tar -xf - -C /scr/src_copy
  R CMD build --no-manual /scr/src_copy
  _R_CHECK_CRAN_INCOMING_REMOTE_=false R CMD check --as-cran --no-manual TSDFGS_*.tar.gz || true'
echo CHECK_DONE
