# Use an official Python image as a parent image.
# python:3.10-slim-bullseye is a great choice: it's smaller and based on Debian "Bullseye".
FROM python:3.10-slim-bullseye

# Avoid interactive prompts during installation
ENV DEBIAN_FRONTEND=noninteractive

# Install all required dependencies in a single layer for a smaller image
# Includes: git, ca-certificates (for downloading), bedtools, samtools, bc, and gnu-parallel
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
    ca-certificates \
    git \
    bedtools \
    samtools \
    bc \
    parallel && \
    rm -rf /var/lib/apt/lists/*

# Set a working directory
WORKDIR /app

# This makes the container executable by default, dropping into a shell
CMD ["/bin/bash"]