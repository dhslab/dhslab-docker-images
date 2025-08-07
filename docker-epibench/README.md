FROM ubuntu:focal-20221130

LABEL version="1.0"
LABEL description="Image for EpiBench a software tool designed for predicting DNA methylation levels using genomic sequence data and histone modification marks"

ENV DEBIAN_FRONTEND=noninteractive

# Install essential packages
RUN apt-get update -y && \
    apt-get install -y --no-install-recommends \
    build-essential \
    curl \
    git \
    wget \
    unzip \
    && apt-get clean && \
    rm -rf /var/lib/apt/lists/*

#set timezone to CDT
RUN ln -sf /usr/share/zoneinfo/America/Chicago /etc/localtime
RUN echo "America/Chicago" > /etc/timezone
RUN dpkg-reconfigure --frontend noninteractive tzdata

# Install conda (mamba version)
ENV CONDA_DIR /opt/conda
RUN cd /tmp && mkdir -p $CONDA_DIR && \
    wget "https://github.com/conda-forge/miniforge/releases/download/24.3.0-0/Mambaforge-24.3.0-0-Linux-x86_64.sh" && \
    bash Mambaforge-24.3.0-0-Linux-x86_64.sh -f -b -p $CONDA_DIR && \
    rm -f Mambaforge-24.3.0-0-Linux-x86_64.sh && \
    $CONDA_DIR/bin/mamba install -y --channel conda-forge --channel bioconda python=3.11 pandas openpyxl cython globus-cli SQLAlchemy biotite pyodbc boto3 && \
    $CONDA_DIR/bin/conda clean -y --all && \
    ln -s /opt/conda/bin/* /usr/local/bin/

RUN curl "https://awscli.amazonaws.com/awscli-exe-linux-x86_64.zip" -o "awscliv2.zip" && \
    unzip awscliv2.zip && ./aws/install

# Add conda to PATH
ENV PATH=$CONDA_DIR/bin:$PATH

# Install EpiBench dependencies and other useful tools
RUN mamba install -y -c pytorch -c bioconda -c conda-forge \
    'pytorch>=1.12.0' cpuonly \
    'captum>=0.5.0' \
    'numpy>=1.21.0' \
    'pandas>=1.3.0' \
    'scipy>=1.7.0' \
    'scikit-learn>=1.0.0' \
    'h5py>=3.1.0' \
    'matplotlib>=3.5.0' \
    'seaborn>=0.11.0' \
    'plotly>=5.3.0' \
    'biopython>=1.79' \
    'pybigwig>=0.3.18' \
    'pyfaidx>=0.6.4' \
    'shap>=0.40.0' \
    'optuna>=3.0.0' \
    'tensorboard>=2.10.0' \
    'tqdm>=4.62.0' \
    'pyyaml>=6.0' \
    'joblib>=1.1.0' \
    'pydantic>=1.9.0' \
    'jupyter>=1.0.0' \
    'statsmodels>=0.13.0' \
    'jinja2>=3.0.0' \
    'psutil>=5.9.0' \
    'jsonschema>=4.0.0' \
    'rich>=12.0.0' \
    'pytest>=7.0.0' \
    'tabulate>=0.8.9' \
    bedtools && \
    mamba clean -y --all

# ucsc utilities
RUN mkdir -p /tmp/ucsc && \
    cd /tmp/ucsc && \
    wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bedGraphToBigWig http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bedToBigBed http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bigBedToBed http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bigWigAverageOverBed http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bigWigToBedGraph http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/wigToBigWig && \
    chmod ugo+x * && \
    mv * /usr/local/bin/ && \
    rm -rf /tmp/ucsc
    
# Clone and install EpiBench
RUN git clone https://github.com/Bonney96/epibench.git /opt/epibench
WORKDIR /opt/epibench
RUN pip install .

# Set default work directory
WORKDIR /data

ENV PATH=/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:/opt/conda/bin 