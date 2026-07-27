# ghcr.io/dhslab/docker-casoffinder-bulge
#
# Cas-OFFinder (CPU/OpenCL) + the cas-offinder-bulge wrapper (DNA/RNA bulge enumeration,
# whose output columns the pipeline's combine_offtarget_results.py parses) + UCSC
# faToTwoBit (so PREP_CASOFFINDER_REF can build the .2bit from params.fasta in the same
# image). Satisfies `container` in prep_casoffinder_ref.nf and casoffinder.nf.
#
# CPU-only: cas-offinder runs on the 'C' (CPU) OpenCL device via pocl, so no GPU is needed
# on the SLURM nodes.
FROM ubuntu:22.04

ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update && apt-get install -y --no-install-recommends \
        ca-certificates wget unzip git g++ cmake make \
        opencl-headers ocl-icd-opencl-dev pocl-opencl-icd \
        python3 python3-pip python-is-python3 \
    && rm -rf /var/lib/apt/lists/*
# cas-offinder-bulge ships a `#!/usr/bin/python` shebang; python-is-python3 provides
# /usr/bin/python -> python3 so the wrapper runs (its body is Python-3 clean).

# --- Cas-OFFinder (3.0.0b3) built against pocl for CPU execution -----------------
ARG CASOFFINDER_VERSION=3.0.0b3
RUN git clone --depth 1 --branch ${CASOFFINDER_VERSION} \
        https://github.com/snugel/cas-offinder.git /opt/cas-offinder \
    && cmake -S /opt/cas-offinder -B /opt/cas-offinder/build \
    && cmake --build /opt/cas-offinder/build --target cas-offinder \
    && install -m755 /opt/cas-offinder/build/cas-offinder /usr/local/bin/cas-offinder

# --- cas-offinder-bulge wrapper (DNA/RNA bulges; emits the Bulge type/Size columns) ---
RUN git clone --depth 1 https://github.com/hyugel/cas-offinder-bulge.git /opt/cas-offinder-bulge \
    && install -m755 /opt/cas-offinder-bulge/cas-offinder-bulge /usr/local/bin/cas-offinder-bulge

# --- UCSC faToTwoBit (for PREP_CASOFFINDER_REF) --------------------------------------
RUN wget -q http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/faToTwoBit \
        -O /usr/local/bin/faToTwoBit \
    && chmod +x /usr/local/bin/faToTwoBit

CMD ["/bin/bash"]
