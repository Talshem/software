# Use an official Python runtime as a parent image
FROM python:3.9-slim

# Set the working directory in the container
WORKDIR /app

# Install dependencies for NUPACK, Miniconda, and Google Cloud SDK
RUN apt-get update && \
    apt-get install -y \ 
    wget \
    curl \
    unzip \
    gnupg \
    build-essential \
    libcurl4-openssl-dev \
    libssl-dev \
    libbz2-dev \
    liblzma-dev \
    zlib1g-dev \
    libncurses5-dev \
    libsqlite3-dev \
    libreadline-dev \
    libffi-dev \
    apt-transport-https \
    ca-certificates
    
COPY credentials.json /app/


RUN echo "deb [signed-by=/usr/share/keyrings/cloud.google.gpg] https://packages.cloud.google.com/apt cloud-sdk main" \
    | tee -a /etc/apt/sources.list.d/google-cloud-sdk.list
    
RUN curl https://packages.cloud.google.com/apt/doc/apt-key.gpg | apt-key --keyring /usr/share/keyrings/cloud.google.gpg add -

RUN apt-get update && apt-get install -y google-cloud-sdk

RUN gcloud auth activate-service-account --key-file=/app/credentials.json

RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh && \
    /bin/bash Miniconda3-latest-Linux-x86_64.sh -b -p /opt/miniconda && \
    rm Miniconda3-latest-Linux-x86_64.sh && \
    /opt/miniconda/bin/conda update -n base -c defaults conda

ENV PATH /opt/miniconda/bin:$PATH

RUN gsutil cp gs://spry-ivy-431810-v0.appspot.com/nupack-4.0.1.12.zip /app/nupack-4.0.1.12.zip && \
    cd /app && \
    unzip nupack-4.0.1.12.zip && \
    cd nupack-4.0.1.12 && \
    pip install -U nupack -f ./package

COPY . /app

COPY requirements.txt /app/

RUN pip install --no-cache-dir -r requirements.txt

EXPOSE 8080

CMD ["python", "tool/server.py"] 
