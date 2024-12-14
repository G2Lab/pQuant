# Use Ubuntu 20.04 as the base image
FROM ubuntu:20.04

# Set non-interactive frontend
ENV DEBIAN_FRONTEND=noninteractive

# Install dependencies and add the PPA for GCC-11
RUN apt-get update && apt-get install -y software-properties-common \
    && add-apt-repository ppa:ubuntu-toolchain-r/test \
    && apt-get update \
    && apt-get install -y \
    build-essential \
    cmake \
    gcc-11 \
    g++-11 \
    llvm-9 \
    clang-9 \
    python3 \
    python3-pip \
    git \
    wget \
    curl \
    unzip \
    zlib1g-dev \
    && update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-11 100 \
    && update-alternatives --install /usr/bin/g++ g++ /usr/bin/g++-11 100 \
    && apt-get clean

# Install Snakemake
RUN pip3 install --upgrade pip && pip3 install snakemake==7.30.1

# Clone and install OpenFHE
RUN git clone https://github.com/openfheorg/openfhe-development.git /opt/openfhe-development \
    && cd /opt/openfhe-development \
    && git checkout v1.1.2 \
    && mkdir build && cd build \
    && cmake .. -DCMAKE_INSTALL_PREFIX=../../openfhe \
    && make -j$(nproc) && make install

# Clone the pQuant repository
RUN git clone https://github.com/g2lab/pQuant.git /opt/pQuant

# Build pQuant
RUN cd /opt/pQuant && mkdir build && cd build && cmake .. && make -j$(nproc)

# Set up the working directory
WORKDIR /opt/pQuant

# Entry point for the container (adjust based on usage)
CMD ["snakemake", "--help"]