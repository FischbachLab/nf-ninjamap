# Dockerfile - https://github.com/blekhmanlab/docker_dada2/blob/master/Dockerfile

# Base Image
FROM bioconductor/bioconductor_docker:RELEASE_3_12

LABEL authors="Xiandong Meng" \
      description="ninjaMap results aggregation" \
      software.name="ninjaMap R" \
      software.description="ninjaMap results aggregation" \
      software.website="https://github.com/FischbachLab/nf-ninjamap/" \
      container.category="tool"

USER root:root

# Install base utilities
RUN apt-get update && apt-get install -y build-essential wget && apt-get clean

#RUN apt-get update && \
#    apt-get install -y build-essentials  && \
#    apt-get install -y wget && \
#    apt-get clean 
    #rm -rf /var/lib/apt/lists/*

# install iniconda
RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda.sh && \
    /bin/bash ~/miniconda.sh -b -p /opt/conda

ENV PATH "/opt/conda/bin:$PATH"

# this "version" label reflects the installed version of DADA2:

#RUN R -e 'BiocManager::install("dada2", version="3.12", update=TRUE, ask=FALSE)'
#RUN R -e 'BiocManager::install("phyloseq", version="3.12", update=TRUE, ask=FALSE)'
#RUN R -e 'BiocManager::install("decontam", version="3.12", update=TRUE, ask=FALSE)'
RUN R -e 'install.packages("tidyverse")'
RUN R -e 'install.packages("foreach")'
RUN R -e 'install.packages("fs")'
RUN R -e 'install.packages("ggbeeswarm")'
RUN R -e 'install.packages("lubridate")'
RUN R -e 'install.packages("zip")'
RUN R -e 'install.packages("hrbrthemes")'
RUN R -e 'install.packages("waffle")'
RUN R -e 'install.packages("mdthemes")'
RUN R -e 'install.packages("scales")'
RUN R -e 'install.packages("purrr")'
RUN R -e 'install.packages("magrittr")'
RUN R -e 'install.packages("ggplot2")'
RUN R -e 'install.packages("remotes")'
RUN R -e 'remotes::install_github("hrbrmstr/waffle")'

# Update conda to latest version.
RUN conda update -n base -c defaults conda

# Add conda channels
RUN conda config --add channels defaults
RUN conda config --add channels bioconda
RUN conda config --add channels conda-forge

# Install software packages
#RUN conda install --freeze-installed -y numpy pandas python=3.7.6 \
RUN conda install --freeze-installed -c conda-forge -y awscli \
   && conda clean -afy

# Install time and memory usage tool (timem)
#RUN conda install -c jrmadsen -y timemory

#ADD fix_summary.py /usr/local/bin/fix_summary.py
#RUN chmod +x /usr/local/bin/fix_summary.py
#ADD 16s_wrapper.sh /usr/local/bin/16s_wrapper.sh
#RUN chmod +x /usr/local/bin/16s_wrapper.sh
#ADD 16S.R /usr/local/bin/16S.R
#RUN chmod +x /usr/local/bin/16S.R

RUN mkdir -p /work
WORKDIR /work

# Get Repo
COPY . .

RUN chmod -R +rx ./*.R

CMD ["pwd"]
