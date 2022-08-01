ARG FUNCTION_DIR="/function"
ARG BINS_DIR="/function/bin"

# Set the base image to Ubuntu
FROM ubuntu:20.04

# Add timezone so that tzdata does not hang to ask for timezone
# https://grigorkh.medium.com/fix-tzdata-hangs-docker-image-build-cdb52cc3360d
ENV TZ=Europe/London
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone


# Create function work directories
ARG FUNCTION_DIR
ARG BINS_DIR
RUN mkdir -p ${FUNCTION_DIR} \
	&& mkdir -p ${BINS_DIR}

# Add basic utilities and compilers
RUN apt-get update \
	&& apt-get upgrade -y \
	&& apt-get install -y \
    apt-utils \
    zip \
    unzip \
    git \
    curl \
    gawk \
    wget \
    cmake \
    gcc \
    gfortran \
    apt-transport-https \
    ca-certificates \
    software-properties-common

# install libraries and dependencies needed for python/ R environments
RUN apt-get install -y \
  zlib1g-dev \
  libbz2-dev \
  libgsl0-dev \
  libncurses5-dev \
  libmysqlclient-dev \
  default-libmysqlclient-dev \
  autotools-dev \
  automake \
  redis-tools \
  libblas-dev \
  liblapack-dev \
  libcurl4-openssl-dev \
  libssl-dev \
  libxml2-dev \
  libpng-dev \
  libjpeg-dev \
  zlib1g-dev \
  libfontconfig1-dev \
  libpq-dev \
  libudunits2-dev \
  dirmngr \
  libglu1-mesa-dev \
  freeglut3-dev \
  mesa-common-dev \
	libffi-dev \
  libgdal-dev

############################################################
# Install and configure Julia
# version=current

RUN apt-get install -y julia

RUN julia -e 'using Pkg; Pkg.add("PyCall"); Pkg.add("Circuitscape")'

############################################################
# Install and configure R environment
# Version=latest available

# install r-base
RUN apt-get install -y \
	r-base \
	build-essential \
	&& apt-get update --fix-missing

# install R packages
RUN Rscript -e "install.packages('lme4', dependencies=TRUE)" \
    && Rscript -e "install.packages('MuMIn', dependencies=TRUE)" \
    && Rscript -e "install.packages('Matrix', dependencies=TRUE)"

############################################################
# Install and configure Python
# version=3.10.5

RUN wget "https://www.python.org/ftp/python/3.10.5/Python-3.10.5.tgz" \
  && tar -xvzf Python*.tgz \
	&& cd Python-* \
	&& ./configure --enable-shared \
		--enable-ipv6 \
		--with-optimizations \
		--with-lto \
		--prefix=/usr/local \
		LDFLAGS="-Wl,-rpath /usr/local/lib" \
	&& make \
	&& make install \
	&& export LD_LIBRARY_PATH=/usr/local/lib/

# link python to python3 binary
RUN ln -s /usr/bin/python3 /usr/bin/python

RUN pip3 install --no-cache --upgrade pip wheel setuptools \
	&& pip3 install --no-cache --upgrade \
		numpy \
		scipy \
		networkx \
		seaborn \
		matplotlib \
		pandas \
		sortedcontainers \
		julia \
		geopy \
		geopandas \
		shapely \
		scikit-learn \
		pysam \
		rpy2

RUN pip3 install "setuptools<58.0.0" && pip3 install deap


############################################################
# Finalize

# clean up
RUN apt-get autoremove && apt-get clean

# install resistnet from git
RUN git clone https://github.com/tkchafin/resistnet.git

# Update permissions for binaries and append location to PATH
RUN chmod -R ugoa+rwx ${BINS_DIR}
ENV PATH="${BINS_DIR}:${PATH}"
ENV HOME=/root
ENV LD_LIBRARY_PATH=/usr/local/lib/
