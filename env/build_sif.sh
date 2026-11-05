#!/bin/bash
# Build env/TSDFGS.sif with Apptainer. Usage: env/build_sif.sh <work_dir>
# <work_dir> is an empty, writable directory for temporary files (for example from mktemp -d).
# The repository copies of the .def and .yml use CRLF; Apptainer needs LF, so the
# build runs on LF copies made in the work directory.
set -eu
scratch=$1
[ -d "$scratch" ] || { echo "scratch: directory not found: $scratch" >&2; exit 2; }
here=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
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
