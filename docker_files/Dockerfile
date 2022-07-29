# Setup base image
FROM continuumio/miniconda3:latest

USER root:root
ENV PATH "/opt/conda/bin:$PATH"

# Install dependencies
RUN apt-get update --allow-releaseinfo-change  && apt-get install -y libtbb2

# Update conda to latest version.
RUN conda update -n base -c defaults conda

# channel order is important
RUN conda install --freeze-installed -c bioconda -c conda-forge -y nomkl numpy pandas \
       awscli biopython \
       bowtie2 bedtools vcftools samtools==1.11 sambamba pysam pysamstats pybedtools bbmap \
       && conda clean -afy

RUN conda install --freeze-installed -c conda-forge bc
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
