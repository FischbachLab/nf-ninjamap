#!/bin/bash -x

## Setup
# docker container run -it --rm -v ${HOME}/efs/docker/:/mnt continuumio/miniconda
# docker container run -it --rm -v ${HOME}/efs/docker/:/mnt sunitjain/ninjamap
# cd /mnt/NinjaMap/MultiSample/
# conda install -c conda-forge -y awscli vim htop
# conda install -c bioconda -y sourmash bbmap
# cd /mnt/NinjaMap/MultiSample/bbnorm

## Commands
# bbnorm.sh \
#     in=../fastq/W3/FMT2-W3-mouse11_S245_R1_001.fastq.gz \
#     in2=../fastq/W3/FMT2-W3-mouse11_S245_R2_001.fastq.gz \
#     out=highpass_FMT2-W3-mouse11.fastq \
#     outt=lowpass_FMT2-W3-mouse11.fastq \
#     passes=1 \
#     target=999999999 \
#     min=10 \
#     hist=hist_FMT2-W3-mouse11.txt &> bbnorm.log &