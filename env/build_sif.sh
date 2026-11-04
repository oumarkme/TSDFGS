#!/bin/bash
# Build env/TSDFGS.sif on gs4. Usage: env/build_sif.sh <scratch_dir>
# <scratch_dir> must be a private directory made with mktemp -d /raid/tmp/TSDFGS.XXXXXX.
# The repository copies of the .def and .yml use CRLF; Apptainer needs LF, so the
# build runs on LF copies made in the scratch directory.
set -eu
export PATH=/home/user/miniforge-pypy3/envs/mark/bin:$PATH
scratch=$1
case "$scratch" in /raid/tmp/TSDFGS.*) ;; *) echo "bad scratch $scratch" >&2; exit 2;; esac
here=/lab_nas_home_user/Mark/TSDFGS/env
mkdir -p "$scratch/build" "$scratch/tmp" "$scratch/cache"
sed 's/\r$//' "$here/TSDFGS.def" > "$scratch/build/TSDFGS.def"
sed 's/\r$//' "$here/TSDFGS.yml" > "$scratch/build/TSDFGS.yml"
export APPTAINER_TMPDIR=$scratch/tmp APPTAINER_CACHEDIR=$scratch/cache TMPDIR=$scratch/tmp
cd "$scratch/build"
start=$(date +%s)
if [ -e "$here/TSDFGS.sif" ]; then mv "$here/TSDFGS.sif" "$here/TSDFGS_$(date +%Y%m%d).sif"; fi
apptainer build --fakeroot --mksquashfs-args '-comp zstd -Xcompression-level 15 -b 256K -processors 8' \
    "$scratch/TSDFGS.sif" TSDFGS.def
mv "$scratch/TSDFGS.sif" "$here/TSDFGS.sif"
echo "BUILD_OK secs=$(( $(date +%s)-start )) size=$(du -m "$here/TSDFGS.sif" | cut -f1)MB"
