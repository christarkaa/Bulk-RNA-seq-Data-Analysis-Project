#!/usr/bin/env bash

# Make a new directory for this project and change to the directory
mkdir bulk && cd "$_"

# Download datasets
fastq-dump --split-files SRR28420795

# Make new directories for quality control and Mapping
mkdir QC_Reports
mkdir Mapping

# Quality control using FastQC
fastqc SRR28420795_1.fastq SRR28420795_2.fastq -o QC_Reports

# Summarizing the QC results
multiqc .

# Trimming using Sickle
sickle pe -f SRR28420795_1.fastq -r  SRR28420795_2.fastq -t sanger \
  -o trimmed.SRR28420795_1.fastq -p trimmed.SRR28420795_2.fastq \
  -s single.SRR28420795.fastq -q 20 -l 50

# Alignment with HISAT2
hisat2 -p 8 -x genome -1 paired1.fastq -2 paired2.fastq -S Mapping/SRR28420795.sam

# Convert to a BAM file
samtools view -@ 20 -S -b Mapping/SRR28420795.sam > Mapping/SRR28420795.bam

# Sort the BAM file
samtools sort -@ 32 -o Mapping/SRR28420795.sorted.bam Mapping/SRR28420795.bam

# Quantification of counts using featureCounts
featureCounts -T 8 -t exon -g gene_id -a annotations.gtf -o counts.SRR28420795.txt Mapping/SRR28420795.sorted.bam

