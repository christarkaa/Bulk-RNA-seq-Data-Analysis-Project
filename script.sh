#!/usr/bin/env bash

# Make a new directory for this project and change to the directory
mkdir bulk && cd '$_"

# Download datasets
fastq-dump --split-files SRR28420795

# Make a new directory for quality control
mkdir QC_Reports

# Quality control using FastQC
fastqc SRR28420795_1.fastq SRR28420795_2.fastq -o QC_Reports

# Summarising the QC results
multiqc .

# Trimming using Trimmomatic
