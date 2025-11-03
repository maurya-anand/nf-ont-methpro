FROM mambaorg/micromamba:2.3.3-ubuntu22.04

USER root

RUN apt-get update && \
    apt-get install -y \
    curl gzip && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

USER $MAMBA_USER

COPY --chown=$MAMBA_USER:$MAMBA_USER ./container/environment.yml /tmp/

RUN micromamba install -y -n base -f /tmp/environment.yml && \
    micromamba clean --all --yes

WORKDIR /tools

RUN curl -L "https://cdn.oxfordnanoportal.com/software/analysis/dorado-1.2.0-linux-x64.tar.gz" -o dorado-1.2.0-linux-x64.tar.gz && \
    tar -xzf dorado-1.2.0-linux-x64.tar.gz && \
    rm dorado-1.2.0-linux-x64.tar.gz && \
    curl -L "https://github.com/nanoporetech/modkit/releases/download/v0.5.0/modkit_v0.5.0_u16_x86_64.tar.gz" -o modkit_v0.5.0_u16_x86_64.tar.gz && \
    tar -xzf modkit_v0.5.0_u16_x86_64.tar.gz && \
    rm modkit_v0.5.0_u16_x86_64.tar.gz

ARG MAMBA_DOCKERFILE_ACTIVATE=1

ENV PATH="/tools/dorado-1.2.0-linux-x64/bin:/tools/modkit_v0.5.0_u16_x86_64/bin:$PATH"