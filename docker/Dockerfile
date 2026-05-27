# syntax=docker/dockerfile:1.6
# Multi-stage build: solve the conda env once, copy a slim runtime image.
# Built on mambaforge so pertpy / pydeseq2 / scanpy resolve from conda-forge + bioconda.

ARG MAMBA_VERSION=23.3.1-1
FROM condaforge/mambaforge:${MAMBA_VERSION} AS builder

ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    MAMBA_NO_BANNER=1

WORKDIR /opt/perturbflow

# Solve env first so layer is cached when only source changes.
COPY docker/environment.yml /tmp/environment.yml
RUN mamba env create -p /opt/conda/envs/perturbflow -f /tmp/environment.yml \
    && mamba clean -afy \
    && find /opt/conda/envs/perturbflow -follow -type f -name '*.a' -delete \
    && find /opt/conda/envs/perturbflow -follow -type f -name '*.pyc' -delete

# Install the package itself (editable not needed in image).
COPY pyproject.toml README.md LICENSE ./
COPY src ./src
RUN /opt/conda/envs/perturbflow/bin/pip install --no-deps --no-cache-dir .

# ---------- Runtime image ----------
FROM condaforge/mambaforge:${MAMBA_VERSION} AS runtime

ENV PATH=/opt/conda/envs/perturbflow/bin:$PATH \
    CONDA_DEFAULT_ENV=perturbflow \
    PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    MPLBACKEND=Agg \
    NUMBA_CACHE_DIR=/tmp/numba_cache \
    OPENBLAS_NUM_THREADS=1 \
    OMP_NUM_THREADS=1

# Pull the solved env over.
COPY --from=builder /opt/conda/envs/perturbflow /opt/conda/envs/perturbflow

# Pipeline files needed at runtime (Snakefile, configs).
WORKDIR /workspace
COPY workflow /workspace/workflow

# Non-root user (AWS Batch friendly).
RUN groupadd -r perturb && useradd -r -g perturb -m -d /home/perturb perturb \
    && chown -R perturb:perturb /workspace
USER perturb

ENTRYPOINT ["perturbflow"]
CMD ["--help"]
