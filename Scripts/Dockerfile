# Choose a base image
# Use a lightweight base image that matches your pipeline's requirements.
FROM ubuntu:22.04

# Set environment variables to prevent interactive prompts during installation
ENV DEBIAN_FRONTEND=noninteractive

# Install dependencies  
RUN apt-get update && apt-get install -y \
    build-essential wget unzip curl git cmake python3 python3-pip python3-venv \
    python3-all-dev clang autoconf automake make gcc g++ tabix gawk pandoc \
    r-base bzip2 zlib1g-dev libbz2-dev liblzma-dev libffi-dev libssl-dev \
    libncurses5-dev libdeflate-dev parallel \
    libcurl4-openssl-dev libxml2-dev \
    libharfbuzz-dev libboost-all-dev libfribidi-dev libfreetype6-dev libpng-dev \
    libtiff5-dev libjpeg-dev libcairo2-dev libfontconfig1-dev libx11-dev libpango1.0-dev \
    libxt-dev gfortran libreadline-dev && \
    ln -s /usr/bin/python3 /usr/bin/python && \
    apt-get clean && rm -rf /var/lib/apt/lists/*  

# Copy the required scripts into the Docker container
# Copy `setup.sh` and `ConstrunctValidation.sh` to `/app` in the container
WORKDIR /app
COPY . /app
 
# Make the scripts executable
RUN chmod +x /app/setup.sh /app/ConstructValidation.sh

# Create a virtual environment before running the setup script  
RUN python3 -m venv /opt/venv  

# Set the virtual environment as the default python environment  
ENV PATH="/opt/venv/bin:$PATH"  

# Run the setup.sh script to install required software and packages
RUN /app/setup.sh

# Set the main pipeline script (`ConstructValidation.sh`) as the default command
# This will run the pipeline when the container is started
#ENTRYPOINT ["/app/ConstructValidation.sh"]
ENTRYPOINT ["/bin/bash"] 

# If you want to pass additional arguments, you can configure CMD
# CMD ["arg1", "arg2"]
