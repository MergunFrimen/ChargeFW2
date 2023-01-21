FROM ubuntu:22.04

ENV PYTHONPATH=/usr/local/lib

RUN apt-get update

ARG DEPS="git cmake libnanoflann-dev \
        libomp-dev libeigen3-dev libfmt-dev \
        libboost-system-dev libboost-program-options-dev libboost-filesystem-dev \
        libfmt-dev nlohmann-json3-dev python3-pybind11 \
        python3-dev g++ gemmi-dev tao-pegtl-dev"

RUN apt-get install -y ${DEPS}

# non-root user
ARG UNAME=vscode
ARG UID=1000
ARG GID=1000
RUN groupadd -g ${GID} -o ${UNAME}
RUN useradd -m -u ${UID} -g ${GID} -o -s /bin/bash ${UNAME}
ENV HOME /home/${UNAME}
USER ${UNAME}

# WORKDIR ${HOME}

# COPY --chown=${UNAME}:${UNAME} . ./ChargeFW2

# RUN cd ChargeFW2 && mkdir build && cd build && \
#     cmake .. -DCMAKE_BUILD_TYPE=Debug -DCMAKE_INSTALL_PREFIX=. && \
#     make && make install

# CMD ["bash"]