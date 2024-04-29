#!/usr/bin/env bash

# Install FastQC
brew install fastqc

# Install multiQC
brew install multiqc

# Install Sickle
brew install sickle

# Install hisat2
CONDA_ SUBDIR=ox-64 conda install hisat2 -c bioconda

# Install samtools
brew install samtools

# Install subread to use featureCounts
wget "https://github.com/ShiLab-Bioinformatics/subread/releases/download/2.0.2/subread-2.0.2-macOS-x86_64.tar.gz"
export PATH="/Users/christophertarkaa/subread-2.0.2-macOS-x86_64/bin:$PATH"
source ~/.zshrc
