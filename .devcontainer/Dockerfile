FROM ubuntu:22.04

ENV PYTHONPATH=/usr/local/lib

ARG DEPS="\
        clang \
        cmake \
        g++ \
        gemmi-dev \
        git \
        libboost-filesystem-dev \
        libboost-program-options-dev \
        libboost-system-dev \
        libeigen3-dev \
        libfmt-dev \
        libfmt-dev \
        libnanoflann-dev \
        libomp-dev \
        nlohmann-json3-dev \
        python3-dev \
        python3-pybind11 \
        sudo \
        tao-pegtl-dev"

RUN apt-get update && \
        apt-get install -y ${DEPS} && \
        rm -rf /var/lib/apt/lists/* && \
        apt-get clean && \
        apt-get autoremove
RUN apt-get update && \
        apt-get install -y ${DEPS} && \
        rm -rf /var/lib/apt/lists/* && \
        apt-get clean && \
        apt-get autoremove

# non-root user
ARG UNAME=nonroot
ARG UNAME=nonroot
ARG UID=1000
ARG GID=1000
RUN groupadd -g ${GID} -o ${UNAME} && \
        useradd -m -u ${UID} -g ${GID} -o -s /bin/bash ${UNAME} && \
        echo "${UNAME}:${UNAME}" | chpasswd && adduser nonroot sudo

USER ${UNAME}
