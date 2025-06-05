bash
#!/bin/bash

# Activate the virtual environment  
source /opt/venv/bin/activate  

## Stop the script if any command fails
set -e

# Install dependencies
echo "=== Installing dependencies ==="
apt update && apt install -y \
  build-essential wget unzip curl git cmake python3 python3-pip python3-venv tabix gawk pandoc \
  r-base python3-all-dev clang autoconf automake make gcc g++ \
  bzip2 zlib1g-dev libbz2-dev liblzma-dev libffi-dev libssl-dev \
  libncurses5-dev libdeflate-dev parallel \
  # R-specific dependencies  
  libcurl4-openssl-dev libxml2-dev \
  libharfbuzz-dev libboost-all-dev libfribidi-dev libfreetype6-dev libpng-dev \
  libtiff5-dev libjpeg-dev libcairo2-dev libfontconfig1-dev libx11-dev libpango1.0-dev \
  libxt-dev gfortran libreadline-dev && \
  # Create symbolic link for Python  
  ln -s /usr/bin/python3 /usr/bin/python && \
  # Clean up APT  
  apt-get clean && rm -rf /var/lib/apt/lists/*


## Install R packages
echo "=== Installing R Base and R Packages ==="
chmod +x Install_R_Packages.R
Rscript Install_R_Packages.R

# Download minibar
# echo "=== Installing Minibar ==="
# wget https://raw.githubusercontent.com/calacademy-research/minibar/master/minibar.py
# mv minibar.py /usr/local/bin/minibar.py
# chmod +x /usr/local/bin/minibar.py
# pip3 install biopython

# Install Chopper (0.8.0)
echo "=== Installing Chopper (0.8.0) ==="
wget https://github.com/wdecoster/chopper/releases/download/v0.8.0/chopper-linux.zip
unzip chopper-linux.zip
chmod +x chopper
mv chopper /usr/local/bin
rm -f chopper-linux.zip

# Install minimap2 (2.24)
echo "=== Installing Minimap2 (2.24) ==="
wget https://github.com/lh3/minimap2/releases/download/v2.24/minimap2-2.24_x64-linux.tar.bz2
tar -jxvf minimap2-2.24_x64-linux.tar.bz2
mv minimap2-2.24_x64-linux/minimap2 /usr/local/bin
chmod +x /usr/local/bin/minimap2
rm -rf minimap2-2.24_x64-linux.tar.bz2 minimap2-2.24_x64-linux

# Install samtools (1.19.2)
echo "=== Installing samtools (1.19.2) ==="
wget https://github.com/samtools/samtools/releases/download/1.19.2/samtools-1.19.2.tar.bz2
tar -jxvf samtools-1.19.2.tar.bz2
cd samtools-1.19.2
make
make install
cd ..
rm -rf samtools-1.19.2 samtools-1.19.2.tar.bz2

# Install bcftools (1.19)
echo "=== Installing bcftools (1.19) ==="
wget https://github.com/samtools/bcftools/releases/download/1.19/bcftools-1.19.tar.bz2
tar -jxvf bcftools-1.19.tar.bz2
cd bcftools-1.19
make
make install
cd ..
rm -rf bcftools-1.19 bcftools-1.19.tar.bz2

# Install bedtools (2.30.0)
echo "=== Installing bedtools (2.30.0) ==="
wget https://github.com/arq5x/bedtools2/releases/download/v2.30.0/bedtools-2.30.0.tar.gz
tar -xzf bedtools-2.30.0.tar.gz
cd bedtools2
make
make install
cd ..
rm -rf bedtools2 bedtools-2.30.0.tar.gz

# Install Racon (1.5.0)
echo "=== Installing Racon (1.5.0) ==="
wget https://github.com/lbcb-sci/racon/archive/refs/tags/1.5.0.tar.gz
tar -xzf 1.5.0.tar.gz
cd racon-1.5.0
mkdir build && cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make
mv bin/racon /usr/local/bin/
cd ../..
rm -rf racon-1.5.0 1.5.0.tar.gz

# Install HTSlib  
echo "=== Installing HTSlib() ==="
wget https://github.com/samtools/htslib/releases/download/1.21/htslib-1.21.tar.bz2
tar -jxvf htslib-1.21.tar.bz2
cd htslib-1.21
make
make install
cd ..
rm -rf htslib-1.21 htslib-1.21.tar.bz2

# Install Medaka (2.0.1)
echo "=== Installing Medaka (2.0.1) ==="
# python3 -m venv medaka
# . ./medaka/bin/activate
# pip install --upgrade pip
# pip install medaka
pip3 install medaka==2.0.1 pyabpoa wheel
# wget https://github.com/nanoporetech/medaka/releases/download/v2.0.1/medaka-2.0.1.tar.gz
# tar -xzf medaka-2.0.1.tar.gz
# cd medaka-2.0.1
# make
# make install
# cd ..
# rm -rf medaka-2.0.1 medaka-2.0.1.tar.gz

# Install Seqtk (1.4)
echo "=== Installing Seqtk (1.4) ==="
wget https://github.com/lh3/seqtk/archive/refs/tags/v1.4.tar.gz
tar -xzf v1.4.tar.gz
cd seqtk-1.4
make
make install
cd ..
rm -rf seqtk-1.4 v1.4.tar.gz

# Install EMBOSS (6.6.0)  
echo "=== Installing EMBOSS (6.6.0) ==="
# wget -m 'ftp://emboss.open-bio.org/pub/EMBOSS/'
# cd emboss.open-bio.org/pub/EMBOSS
# tar -xzf EMBOSS-6.6.0.tar.gz
# cd EMBOSS-6.6.0
# ./configure
# make
# make install
# cd ../../..
wget https://github.com/kimrutherford/EMBOSS/archive/refs/tags/EMBOSS-6.6.0.tar.gz
tar -xzf EMBOSS-6.6.0.tar.gz
cd EMBOSS-EMBOSS-6.6.0/
./configure --prefix=/usr/local/emboss
make
make install
cd ..
rm -rf EMBOSS-EMBOSS-6.6.0 EMBOSS-6.6.0.tar.gz 
ln -s /usr/local/emboss/bin/* /usr/local/bin/

# Install MUSCLE3 (3.8.31)
echo "=== Installing MUSCLE3 (3.8.31) ==="
wget https://drive5.com/muscle/downloads3.8.31/muscle3.8.31_src.tar.gz
tar -xzf muscle3.8.31_src.tar.gz
make -C muscle3.8.31/src
cp muscle3.8.31/src/muscle /usr/local/bin/
rm -rf muscle3.8.31_src.tar.gz muscle3.8.31
