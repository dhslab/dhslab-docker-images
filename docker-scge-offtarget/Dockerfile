FROM python:3.11-slim

LABEL version="1.0"
LABEL description="SCGE off-target analysis: pinned Python ML/genomics stack (pysam, pandas, numpy, scipy, scikit-learn, joblib, matplotlib, pyranges, edlib, openpyxl) + shipped analysis scripts"

ARG DEBIAN_FRONTEND=noninteractive
ENV TZ=America/Chicago

# System layer: CA certs for pip over https, procps for pipeline runtimes, tzdata to set CDT
# (matches the other dhslab images), and build-essential so pyranges<1's Cython deps
# (ncls, sorted_nearest) and edlib can compile if no cp311 wheel is available.
RUN apt-get update -y && \
    apt-get install -y --no-install-recommends \
    ca-certificates \
    procps \
    tzdata \
    build-essential \
    && ln -sf /usr/share/zoneinfo/$TZ /etc/localtime \
    && echo $TZ > /etc/timezone \
    && dpkg-reconfigure --frontend noninteractive tzdata \
    && rm -rf /var/lib/apt/lists/*

# Pinned Python stack.
#  - scikit-learn==1.8.0: loads wgs_shape_model.pkl (pickled under 1.8.0) exactly and
#    forward-loads site14_site5_combined_model.pkl (1.6.1) with an InconsistentVersionWarning.
#    Do NOT drop below the newest model's pickle version or HistGradientBoosting load can raise.
#  - "pyranges<1": find_edited_reads.py uses the v0 .df API, which pyranges 1.x removed.
#  - scipy/edlib/pyranges are imported directly by the off-target scripts (not just transitively).
#  - openpyxl is not imported by any off-target script; kept for parity with docker-scge.
RUN python -m pip install --no-cache-dir --upgrade pip && \
    python -m pip install --no-cache-dir \
    pysam==0.24.0 \
    pandas==2.3.3 \
    numpy==2.3.5 \
    scipy \
    scikit-learn==1.8.0 \
    joblib==1.5.2 \
    matplotlib==3.10.8 \
    "pyranges<1" \
    edlib \
    openpyxl==3.1.5

# Fail the build early if any pin becomes unresolvable / uninstallable on this base.
RUN python -c "import pysam,pandas,numpy,scipy,sklearn,joblib,matplotlib,pyranges,edlib,openpyxl; print('scge-offtarget py deps OK')"

# Off-target analysis scripts (kept alongside the deps that run them).
COPY bin/ /opt/scge-offtarget/bin/
RUN chmod -R +x /opt/scge-offtarget/bin
ENV PATH=/opt/scge-offtarget/bin:$PATH

# Versioned model assets + a build-time load check. Fails the build here (not mid-run)
# if the pinned scikit-learn cannot unpickle a shipped model — the exact failure the
# 1.7.2 pin would have caused for wgs_shape_model.pkl (1.8.0). The forward-compatible
# load of site14_site5_combined_model.pkl (1.6.1) only warns and is expected to pass.
# This mirrors the runtime guard in bin/features.py::check_sklearn_version.
COPY models/ /opt/scge-offtarget/models/
RUN python -c "import joblib; \
[joblib.load(p) for p in ['/opt/scge-offtarget/models/wgs_shape_model.pkl', '/opt/scge-offtarget/models/site14_site5_combined_model.pkl']]; \
print('models load OK under pinned scikit-learn')"

WORKDIR /data
ENTRYPOINT []
