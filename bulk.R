# Load the DESeq2 library
library(DESeq2)

# Set working directory
directory <- "~/Bioinformatics/bulk"
setwd(directory)

# Set the prefix for each output file name 
outputPrefix <- "bulk_DESeq2"

# List the sample files
sampleFiles <- c("SRR28420795.counts.txt",
                 "SRR28420796.counts.txt",
                 "SRR28420797.counts.txt",
                 "SRR28420798.counts.txt")

sampleNames <- c("Sample1", "Sample2", "Sample3", "Sample4")
sampleCondition <- c("Treated", "Control", "Treated", "Control") 
sampleTable <- data.frame(SampleName = sampleNames, fileName = sampleFiles, condition = sampleCondition)

treatments <- c("Treated", "Control")

ddsFeatureCounts <- DESeqDataSetFromFeatureCounts(sampleTable = sampleTable,
                                                 directory = directory,
                                                 design = ~ condition)

colData(ddsFeatureCounts)$condition <- factor(colData(ddsFeatureCounts
