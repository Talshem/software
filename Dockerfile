# Use a base image with necessary build tools
FROM ubuntu:20.04

# Set environment variables
ENV DEBIAN_FRONTEND=noninteractive

# Install required packages and NUPACK
RUN apt-get update && \
    apt-get install -y build-essential gfortran wget && \
    wget https://www.nupack.org/download/nupack-4.0.tar.gz && \
    tar -xzf nupack-4.0.tar.gz && \
    cd nupack-4.0 && \
    ./configure && \
    make && \
    make install

# Copy your server code
COPY tool/server.py

# Set the working directory
WORKDIR /tool

# Set the command to run your server
CMD ["python3", "server.py"]
