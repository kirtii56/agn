# Dockerfile for UHDM AGN Analysis
# Reproducible environment for running all analysis code

FROM python:3.9-slim

LABEL maintainer="UHDM AGN Team"
LABEL description="Docker container for UHDM AGN paper analysis"

# Install system dependencies
RUN apt-get update && apt-get install -y \
    build-essential \
    gfortran \
    git \
    wget \
    texlive-full \
    texlive-latex-extra \
    && rm -rf /var/lib/apt/lists/*

# Set working directory
WORKDIR /workspace

# Copy requirements and install Python packages
COPY code/requirements.txt /workspace/requirements.txt
RUN pip install --no-cache-dir -r requirements.txt

# Copy all code
COPY code/ /workspace/code/
COPY manuscript/ /workspace/manuscript/
COPY figures/ /workspace/figures/
COPY data/ /workspace/data/
COPY workflow/ /workspace/workflow/

# Set environment variables
ENV PYTHONPATH=/workspace/code:$PYTHONPATH
ENV OMP_NUM_THREADS=4

# Default command
CMD ["/bin/bash"]

# To build:
#   docker build -t uhdm-agn -f workflow/Dockerfile .
#
# To run:
#   docker run -it -v $(pwd):/workspace uhdm-agn
#
# To run analysis:
#   docker run -it uhdm-agn bash workflow/run_analysis.sh
