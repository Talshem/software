# Use an official Python runtime as a parent image
FROM python:3.9-slim

# Set the working directory in the container
WORKDIR /app

# Copy the requirements.txt file into the container at /app
COPY requirements.txt /app/

# Install any needed packages specified in requirements.txt
RUN pip install --no-cache-dir -r requirements.txt

# Install dependencies for NUPACK and Miniconda
RUN apt-get update && \
    apt-get install -y \
    wget \
    unzip \
    build-essential \
    libcurl4-openssl-dev \
    libssl-dev \
    libbz2-dev \
    liblzma-dev \
    zlib1g-dev \
    libncurses5-dev \
    libsqlite3-dev \
    libreadline-dev \
    libffi-dev

# Install Miniconda
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh && \
    /bin/bash Miniconda3-latest-Linux-x86_64.sh -b -p /opt/miniconda && \
    rm Miniconda3-latest-Linux-x86_64.sh && \
    /opt/miniconda/bin/conda update -n base -c defaults conda

# Set environment variables for Conda
ENV PATH /opt/miniconda/bin:$PATH

# Copy NUPACK zip file to the container
COPY nupack-VERSION.zip ./

# Unzip and install NUPACK
RUN unzip nupack-VERSION.zip && \
    cd nupack-VERSION && \
    pip install -U nupack -f ./package

# Copy credentials file to the container
COPY credentials.json /app/credentials.json

# Copy the rest of the application code to /app
COPY . /app

# Set environment variable for credentials
ENV GOOGLE_APPLICATION_CREDENTIALS="/app/credentials.json"

# Make port 8080 available to the world outside this container
EXPOSE 8080

# Run tool/server.py when the container launches
CMD ["python", "tool/server.py"]
