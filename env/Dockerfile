# Image with the R, compilers, and R packages used to test TSDFGS, plus TSDFGS itself.
# Build from the project root:  docker build -f env/Dockerfile -t tsdfgs .
FROM ubuntu:24.04

LABEL org.opencontainers.image.title="TSDFGS" \
      org.opencontainers.image.description="R 4.4.3, Rcpp, RcppEigen, GCC 16, and the TSDFGS package" \
      org.opencontainers.image.source="https://github.com/oumarkme/TSDFGS" \
      org.opencontainers.image.licenses="GPL-3.0-or-later"

ENV DEBIAN_FRONTEND=noninteractive \
    MAMBA_ROOT_PREFIX=/opt/conda \
    PATH=/opt/conda/envs/TSDFGS/bin:/usr/local/bin:/usr/bin:/bin \
    LC_ALL=C.UTF-8

COPY environment.yml /opt/environment.yml

# Same steps as env/TSDFGS.def (the Apptainer recipe). The compilers stay because
# R CMD INSTALL needs them at run time.
RUN set -eux; \
    apt-get update; \
    apt-get install -y --no-install-recommends ca-certificates curl bzip2; \
    curl -fsSL https://github.com/mamba-org/micromamba-releases/releases/download/2.5.0-2/micromamba-linux-64 \
        -o /usr/local/bin/micromamba; \
    chmod +x /usr/local/bin/micromamba; \
    micromamba create -y -f /opt/environment.yml; \
    micromamba list -n TSDFGS > /opt/TSDFGS.packages.txt; \
    micromamba clean -afy; \
    rm -rf /opt/conda/pkgs; \
    find /opt/conda -type d -name '__pycache__' -prune -exec rm -rf {} +; \
    apt-get purge -y curl bzip2; \
    apt-get autoremove -y; \
    apt-get clean; \
    rm -rf /var/lib/apt/lists/*

# Install TSDFGS from the checked-out source, then drop the source.
COPY . /tmp/TSDFGS_src
RUN set -eux; \
    R CMD INSTALL --no-docs /tmp/TSDFGS_src; \
    Rscript -e 'library(TSDFGS); cat("TSDFGS", as.character(packageVersion("TSDFGS")), "\n")'; \
    rm -rf /tmp/TSDFGS_src

WORKDIR /work
CMD ["R"]
