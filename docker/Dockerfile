FROM --platform=linux/amd64 condaforge/miniforge3

# autostreamtree
RUN conda install -y python=3.10 geopandas geopy
ENV PROJ_LIB=/opt/conda/share/proj
RUN conda install -c ecoevoinfo -c conda-forge -c bioconda -y resistnet
