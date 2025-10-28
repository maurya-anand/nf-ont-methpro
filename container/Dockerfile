FROM mambaorg/micromamba:2.3.3-ubuntu22.04

USER root

RUN apt-get update && \
    apt-get install -y \
    curl && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

USER $MAMBA_USER

WORKDIR /tools

RUN curl "https://cdn.oxfordnanoportal.com/software/analysis/dorado-1.2.0-linux-x64.tar.gz" -o dorado-1.2.0-linux-x64.tar.gz && \
    tar -xzf dorado-1.2.0-linux-x64.tar.gz

ARG MAMBA_DOCKERFILE_ACTIVATE=1

ENV PATH="/tools/dorado-1.2.0-linux-x64/bin:$PATH"