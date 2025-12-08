FROM ubuntu:22.04

RUN apt-get update && \
    apt-get install -y \
    curl gzip python3-pip \
    samtools bedtools bcftools tabix && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

WORKDIR /tools

RUN curl -L "https://cdn.oxfordnanoportal.com/software/analysis/dorado-1.2.0-linux-x64.tar.gz" -o dorado-1.2.0-linux-x64.tar.gz && \
    tar -xzf dorado-1.2.0-linux-x64.tar.gz && \
    rm dorado-1.2.0-linux-x64.tar.gz && \
    curl -L "https://github.com/nanoporetech/modkit/releases/download/v0.5.0/modkit_v0.5.0_u16_x86_64.tar.gz" -o modkit_v0.5.0_u16_x86_64.tar.gz && \
    tar -xzf modkit_v0.5.0_u16_x86_64.tar.gz && \
    rm modkit_v0.5.0_u16_x86_64.tar.gz && \
    curl -L "https://github.com/lh3/minimap2/releases/download/v2.30/minimap2-2.30_x64-linux.tar.bz2" -o minimap2-2.30_x64-linux.tar.bz2 && \
    tar -xvf minimap2-2.30_x64-linux.tar.bz2 && \
    rm minimap2-2.30_x64-linux.tar.bz2 && \
    pip3 install multiqc pod5

ENV PATH="/tools/dorado-1.2.0-linux-x64/bin:/tools/dist_modkit_v0.5.0_5120ef7:/tools/minimap2-2.30_x64-linux:$PATH"
