#!/usr/bin/env bash

# Make a new directory for this project
mkdir bulk

# Download datasets
fastq-dump --split-files SRR28420798

# Quality control using FastQC
fastqc 
