FROM ubuntu:20.04

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && \
    apt-get install -y build-essential gfortran wget && \
    wget https://www.nupack.org/download/nupack-4.0.tar.gz && \
    tar -xzf nupack-4.0.tar.gz && \
    cd nupack-4.0 && \
    ./configure && \
    make && \
    make install

COPY tool/server.py .
COPY requirements.txt .

RUN pip3 install --no-cache-dir -r requirements.txt

# Set the command to run your server
CMD ["python3", "server.py"]
