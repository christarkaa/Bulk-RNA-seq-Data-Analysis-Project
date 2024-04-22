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
java -jar ~/bin/trimmomatic/trimmomatic-0.39/trimmomatic-0.39.jar \
  PE -phred33 \
   SRR28420795_1.fastq SRR28420795_2.fastq \
  paired1.fastq unpaired1.fastq \
  paired2.fastq unpaired2.fastq \
  TRAILING:20 MINLEN:50

# Alignment with hisat2
hisat2 -x genome -1 paired1.fastq -2 paired2.fastq -S output.sam

# Convert to a bam file

# Sort the bam file

# Quatification of counts using featureCount 
featureCounts -T 8 -t exon -a annotations.gtf -o counts.txt output_sorted.bam

