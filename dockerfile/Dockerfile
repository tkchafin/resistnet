ARG FUNCTION_DIR="/function"
ARG BINS_DIR="/function/bin"

# Set the base image to Ubuntu
FROM --platform=linux/amd64 ubuntu:22.04

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
RUN add-apt-repository -y universe \
	&& add-apt-repository -y multiverse \
	&& apt-get install -y \
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
		gdal-bin \
	  libgdal-dev

RUN export DEBIAN_FRONTEND=noninteractive \
	&& yes | apt install -y --force-yes libmkl-rt

############################################################
# Install and configure R environment
# Version=latest available

# install r-base
RUN apt-get install -y \
		r-base \
		build-essential \
	&& apt-get update --fix-missing \
	&& apt-get install -y \
		r-cran-lme4 \
		r-cran-matrix

# install R packages
RUN Rscript -e "install.packages('MuMIn', dependencies=TRUE)" \
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


############################################################
# Install and configure Julia
# version=1.6.7 (LTS)

# install julia=1.6.7 (LTS)
RUN wget https://julialang-s3.julialang.org/bin/linux/x64/1.6/julia-1.6.7-linux-x86_64.tar.gz \
	&& tar -xvzf julia-1.6.7-linux-x86_64.tar.gz \
	&& ln -s /julia-1.6.7/bin/julia /usr/bin/julia

# install and precompile packages
RUN /usr/bin/julia -e 'using Pkg; \
	ENV["PYTHON"] = "/usr/local/bin/python3"; \
	Pkg.add("PyCall"); \
	Pkg.add("Circuitscape"); \
	Pkg.add("Suppressor"); \
	using PyCall; using Circuitscape; using Suppressor'


############################################################
# Install python dependencies
#

RUN pip3 install --no-cache --upgrade pip wheel setuptools \
	&& pip3 install GDAL==$(gdal-config --version | awk -F'[.]' '{print $1"."$2}') \
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
RUN apt-get autoremove \
	&& apt-get clean \
	&& pip3 cache purge

# Update permissions for binaries and append location to PATH
RUN chmod -R ugoa+rwx ${BINS_DIR}
ENV PATH="${BINS_DIR}:${PATH}"
ENV HOME=/root
ENV LD_LIBRARY_PATH=/usr/local/lib/:/julia-1.6.7/lib:$LD_LIBRARY_PATH
