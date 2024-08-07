# Use a base image with necessary build tools
FROM ubuntu:20.04

# Set environment variables
ENV DEBIAN_FRONTEND=noninteractive

# Install required packages
RUN apt-get update && \
    apt-get install -y build-essential gfortran wget python3 python3-pip

# Download and check the NUPACK file
RUN wget https://www.nupack.org/download/nupack-4.0.tar.gz -O /tmp/nupack-4.0.tar.gz && \
    file /tmp/nupack-4.0.tar.gz && \
    tar -xzf /tmp/nupack-4.0.tar.gz -C /tmp && \
    cd /tmp/nupack-4.0 && \
    ./configure && \
    make && \
    make install && \
    cd /tmp && \
    rm -rf nupack-4.0 nupack-4.0.tar.gz

# Set the working directory
WORKDIR /tool

# Copy server code and requirements.txt
COPY tool/server.py .
COPY requirements.txt .

# Install Python dependencies
RUN pip3 install --no-cache-dir -r requirements.txt

# Set the command to run your server
CMD ["python3", "server.py"]
