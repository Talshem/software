# Use an official Python runtime as a parent image
FROM python:3.9-slim

# Set the working directory in the container
WORKDIR /app

# Copy the requirements.txt file into the container at /app
COPY requirements.txt /app/

# Install any needed packages specified in requirements.txt
RUN pip install --no-cache-dir -r requirements.txt

# Install dependencies for NUPACK
RUN apt-get update && \
    apt-get install -y build-essential cmake wget

# Download and install NUPACK
RUN wget https://github.com/Caltech-NUPACK/nupack/releases/download/4.0.0.27/nupack-4.0.0.27-linux64.tar.gz -O nupack.tar.gz && \
    ls -l && \
    tar -xvzf nupack.tar.gz && \
    ls -l && \
    mv nupack-4.0.0.27 /opt/nupack && \
    rm nupack.tar.gz

# Set NUPACK environment variables
ENV PATH="/opt/nupack/bin:$PATH"
ENV NUPACKHOME="/opt/nupack"

# Copy the rest of the application code to /app
COPY . /app

# Make port 8080 available to the world outside this container
EXPOSE 8080

# Run tool/server.py when the container launches
CMD ["python", "tool/server.py"]
