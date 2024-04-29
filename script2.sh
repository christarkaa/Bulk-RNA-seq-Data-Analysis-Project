#!/usr/bin/env bash

# Define sample IDs
SAMPLES=("SRR28420795" "SRR28420796" "SRR28420797" "SRR28420798")

# Make a new directory for this project and change to the directory
mkdir bulk2 && cd "$_"

# Make common directories for quality control, Mapping and Quantification
mkdir QC_Reports Mapping Counts

# Loop over each sample
for SAMPLE in "${SAMPLES[@]}"; do
    # Download datasets
    fastq-dump --split-files "$SAMPLE"
    
    # Quality control using FastQC
    fastqc "${SAMPLE}_1.fastq" "${SAMPLE}_2.fastq" -o QC_Reports
    
    # Summarizing the QC results
    multiqc QC_Reports
    
    # Trimming using Sickle
    sickle pe -f "${SAMPLE}_1.fastq" -r "${SAMPLE}_2.fastq" -t sanger \
      -o "trimmed.${SAMPLE}_1.fastq" -p "trimmed.${SAMPLE}_2.fastq" \
      -s "single.${SAMPLE}.fastq" -q 20 -l 50
    
    # Alignment with HISAT2
    hisat2 -p 8 -x Genome/grch38/genome -1 "trimmed.${SAMPLE}_1.fastq" -2 "trimmed.${SAMPLE}_2.fastq" -S "Mapping/${SAMPLE}.sam"
    
    # Convert to a BAM file
    samtools view -@ 20 -S -b "Mapping/${SAMPLE}.sam" > "Mapping/${SAMPLE}.bam"
    
    # Sort the BAM file
    samtools sort -@ 32 -o "Mapping/${SAMPLE}.sorted.bam" "Mapping/${SAMPLE}.bam"
    
    # Quantification of counts using featureCounts
    featureCounts -T 8 -p --countReadPairs -B -t exon -g gene_id -a Genome/Homo_sapiens.GRCh38.111.gtf \
    -o "Counts/counts.${SAMPLE}.txt" "Mapping/${SAMPLE}.sorted.bam"
done

