#!/bin/bash
set -e

# Global paths
INSTALL_DIR="/opt/conda"
MINIFORGE_SH="/tmp/Miniforge3.sh"

# Install Miniforge if not already installed
if [ ! -d "$INSTALL_DIR" ]; then
  echo "Installing Miniforge..."
  wget -O "$MINIFORGE_SH" "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh"
  bash "$MINIFORGE_SH" -b -p "$INSTALL_DIR"
else
  echo "Miniforge already installed at $INSTALL_DIR"
fi

# Set conda in PATH for all users
echo "Setting conda path for all users..."
cat <<EOF > /etc/profile.d/conda.sh
export PATH=$INSTALL_DIR/bin:\$PATH
source $INSTALL_DIR/etc/profile.d/conda.sh
EOF

chmod +x /etc/profile.d/conda.sh

# Use conda and mamba as root to create environment
source "$INSTALL_DIR/etc/profile.d/conda.sh"
source "$INSTALL_DIR/etc/profile.d/mamba.sh"

# Create the snakemake environment if it doesn't exist
if ! conda info --envs | grep -q "^wgsEnv\s"; then
  echo "Creating wgsEnv environment..."
  mamba create -y -n wgsEnv -c conda-forge -c bioconda \
    snakemake=7.19.1 sra-tools \
    samtools bedtools fastqc multiqc qualimap trimmomatic picard gatk4
else
  echo "wgsEnv already exists."
fi

# Fix permissions so all users can access the conda environments
echo "Adjusting permissions..."
chmod -R a+rx "$INSTALL_DIR"
chmod -R a+rX "$INSTALL_DIR/envs"
