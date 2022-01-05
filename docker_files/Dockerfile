# Setup base image
FROM continuumio/miniconda3:latest
# 4.8.2

USER root:root
ENV PATH "/opt/conda/bin:$PATH"

# Install dependencies
RUN  apt-get update --allow-releaseinfo-change  && apt-get install -y libtbb2
RUN  conda install --freeze-installed -y nomkl numpy pandas \
      && conda install --freeze-installed -c conda-forge -y awscli biopython \
      && conda install --freeze-installed -c bioconda -y bowtie2 bedtools \
                              vcftools samtools sambamba pysam \
                              pysamstats pybedtools bbmap \
      && conda clean -afy

# && conda install -c jrmadsen -y timemory \
# Use connection.py from python version 3.8
COPY connection_py38.py /opt/conda/lib/python3.7/multiprocessing/connection.py

RUN mkdir -p /work
WORKDIR /work

# Get Repo
COPY . .

# Metadata
LABEL container.maintainer="Xian Meng" \
      container.base.image="continuumio/miniconda3:4.8.2" \
      software.name="ninjaMap" \
      software.description="Strain abundance pipeline" \
      software.website="https://github.com/FischbachLab/ninjaMap" \
      container.category="aligner"

RUN chmod -R +rx ./
