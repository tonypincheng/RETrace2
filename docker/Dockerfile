FROM python:3.9-slim

# Install system dependencies
RUN apt-get update && apt-get install -y \
    git \
    make \
    g++ \
    wget \
    curl \
    bcftools \
    samtools \
    tabix \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    procps \
    default-jre \
    bowtie2 \
    python3-pyqt5 \
    xvfb \
    && apt-get clean && rm -rf /var/lib/apt/lists/*

# Install Picard
RUN wget https://github.com/broadinstitute/picard/releases/download/3.1.1/picard.jar -O /usr/local/bin/picard.jar && \
    echo '#!/bin/bash\njava -jar /usr/local/bin/picard.jar "$@"' > /usr/local/bin/picard && \
    chmod +x /usr/local/bin/picard

# Install Python dependencies from environment.yml (Python-related only)
RUN pip install --no-cache-dir \
    pysam==0.22.1 \
    matplotlib==3.9.4 \
    seaborn==0.13.2 \
    more-itertools==10.7.0 \
    scikit-bio==0.6.3 \
    ete3==3.1.3 \
    biopython==1.85 \
    pandas \
    numpy \
    tqdm \
    psutil \
    methylpy==1.4.7 \
    setuptools \
    PyQt5

# Install HipSTR
RUN git clone https://github.com/HipSTR-Tool/HipSTR.git && \
    cd HipSTR && \
    make && \
    mv HipSTR /usr/local/bin/ && \
    cd .. && rm -rf HipSTR

# Set working directory
WORKDIR /app

# Default command
CMD ["/bin/bash"] 