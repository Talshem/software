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

# Download the specific Google Cloud SDK public key and add it to the keyring
RUN curl -fsSL https://packages.cloud.google.com/apt/doc/apt-key.gpg | gpg --dearmor -o /usr/share/keyrings/cloud.google.gpg

# Add the Google Cloud SDK repository and associate it with the keyring
RUN echo "deb [signed-by=/usr/share/keyrings/cloud.google.gpg] https://packages.cloud.google.com/apt cloud-sdk main" | tee /etc/apt/sources.list.d/google-cloud-sdk.list > /dev/null

# Install Google Cloud SDK
RUN apt-get update && apt-get install -y google-cloud-sdk

# Install Miniconda
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh && \
    /bin/bash Miniconda3-latest-Linux-x86_64.sh -b -p /opt/miniconda && \
    rm Miniconda3-latest-Linux-x86_64.sh && \
    /opt/miniconda/bin/conda update -n base -c defaults conda

# Set environment variable for Google Application Credentials
ARG GOOGLE_APPLICATION_CREDENTIALS
COPY ${GOOGLE_APPLICATION_CREDENTIALS_PATH} /app/credentials.json

RUN gcloud auth activate-service-account --key-file=/app/credentials.json

# Download NUPACK zip file from the Google Cloud Storage bucket and install NUPACK
RUN gsutil cp gs://spry-ivy-431810-v0.appspot.com/nupack-4.0.1.7.zip /app/nupack-4.0.1.7.zip && \
    cd /app && \
    unzip nupack-4.0.1.7.zip && \
    cd nupack-4.0.1.7 && \
    pip install -U nupack -f ./package
# Copy the rest of the application code to /app
COPY . /app

# Set environment variables for Conda
ENV PATH /opt/miniconda/bin:$PATH

# Copy the requirements.txt file into the container at /app
COPY requirements.txt /app/

# Install any needed packages specified in requirements.txt
RUN pip install --no-cache-dir -r requirements.txt

# Make port 8080 available to the world outside this container
EXPOSE 8080

# Run tool/server.py when the container launches
CMD ["python", "tool/server.py"]
