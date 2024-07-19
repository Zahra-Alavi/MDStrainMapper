# Start with a CUDA base image that supports CUDA 12.5.1
FROM nvidia/cuda:12.5.1-cudnn-devel-ubuntu20.04

# Set environment variables for CUDA and to prevent interactive prompts
ENV CUDA_VERSION=12.5.1
ENV DEBIAN_FRONTEND=noninteractive

# Install dependencies in a single RUN command to reduce the number of layers and ensure non-interactive mode
RUN apt-get update && apt-get install -y \
    software-properties-common \
    wget \
    && add-apt-repository ppa:deadsnakes/ppa \
    && apt-get update \
    && apt-get install -y \
    python3.11 \
    python3.11-venv \
    python3.11-dev \
    curl \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Set python3.11 as the default python
RUN update-alternatives --install /usr/bin/python3 python3 /usr/bin/python3.11 1

# Install the latest version of pip using get-pip.py
RUN curl https://bootstrap.pypa.io/get-pip.py -o get-pip.py && python3.11 get-pip.py && rm get-pip.py

# Set the working directory in the container
WORKDIR /app

# Set PYTHONPATH to include the /app directory
ENV PYTHONPATH=/app

# Copy the requirements.txt file into the container at /app
COPY requirements.txt /app/requirements.txt

# Install the dependencies
RUN pip install --no-cache-dir -r requirements.txt

# Copy the rest of the application code into the container at /app
COPY . /app

# Define an environment variable to prevent Python from writing .pyc files to disk
ENV PYTHONUNBUFFERED=1

# Define the command to run the application
ENTRYPOINT ["python3.11", "src/main.py"]
