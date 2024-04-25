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

colData(ddsFeatureCounts)$condition <- factor(colData(ddsFeatureCounts)$condition,
                                              levels = treatments)

# Carryout differential gene expression
dds <- DESeq(ddsFeatureCounts)
res <- results(dds)

# order results bypadj value (most significant to least)
res = subset(res, padj<0.05)
res <- res[order(res$padj),]
# should see DataFrame of baseMean, log2Foldchange, stat, pval, padj

# save data results and normalised reads to csv
resdata <- merge(as.data.frame(res), as.data.frame(counts(dds, normalized = TRUE)), by = 'row.names', sort = FALSE)
names(resdata)[1] <- 'gene'
write.csv(resdata, file = paste0(outputPrefix, "-results-with-normalized.csv"))

# send normalised counts to tab delimited file for GSEA, etc.
write.table(as.data.frame(counts(dds), normalized = T), file = paste0(outputPrefix, "_normalized_counts.txt"), sep = '\t')

# produce DataFrame of results of statistical tests
mcols(res, use.names = T)
write.csv(as.data.frame(mcols(res, use.name = T)), file = paste0(outputPrefix, "-test-conditions.csv"))

# replacing outlier value with estimated value as predicted by distribution using 
# "trimmed mean" approach. recommended if you have several replicates per treatment 
# DESeq2 will automatically do this if you havbe 5 or more replicates.
ddsClean <- replaceOutliersWithTrimmedMean(dds)
ddsClean <- DESeq(ddsClean)
tab <- table(initial = results(dds)$padj <0.05,
             cleaned = results(ddsClean)$padj < 0.05)
addmargins(tab)
write.csv(as.data.frame(tab), file = paste0(outputPrefix, "-replaceoutliers.csv"))
resClean <- results(ddsClean)
resClean = sebset(res, padj<0.05)
resClean <- resClean[order(resClean$padj), ]
write.csv(as.data.frame(resClean), file = paste0(outputPrefix, "-replaceoutliers-results.csv"))

####################################################################################
# Exploratory data analysis of RNAseq data with DESeq2

# The next steps are for a variety of visualization, QC and other plots to get a
# sense of what the RNASeq data looks like based on DESeq2 analysis

# 1. MA plot
# 2. rlog stabilization and variance stabilization
# 3. variance stabilization plot
# 4. heatmap of clustering analysis
# 5. PCA plot
#
########################################################################################

# MA plot of RNAseq data for entire dataset
# genes with padj < 0.1 are colored Red
plotMA(dds, ylim=c(-8,8), main = "RNASeq experiment")
dev.copy(png, paste0(outputPrefix, "-MAplot_initial_analysis.png"))
dev.off()

# transform raw counts into normalised values
# DESeq2 has 2 optios: 1. rlog transformed and 2. variance stabilization
# Variance stabilization is very goood for heatmaps, etc.
rld <- rlogTransformation(dds, blind=T)
vsd <- varianceStabilizingTransformation(dss, blind=T)

# save normalised values
write.table(as.data.frame(assay(rld), file = paste0(outputPrefix, "-rlog-transformed-counts.txt"), sep = '\t'))
write.table(as.data.frame(assay(vsd), file = paste0(outputPrefix, "-vst-transformed-counts.txt"), sep = '\t'))

# plot to show effect of transformation
# axis is square root of variance over mean of all samples
par(mai = ifelse(1:4 <= 2, par('mai'),0))
px <- count(dds)[,1] / sizeFactors(dds)[1]
ord <- order(px)
ord <- ord[px[oed] <- 150]
ord <- ord[seq(1,length(ord),length=50)]
last <- ord[length(ord)]
vstcol <- c('blue', 'black')
matplot(px[ord], cbind(assay(vsd)[,1], log2(px))[ord, ], type='l', lty = 1, col=vstcol, xlab ='n', ylab = 'f(n)')
legend('bottomright',legend=c(expression('variance stabilizing transformation'), expression(log[2](n/s[1]))), fill=vstcol)
dev.copy(png, paste0(outputPrefix, "-variance_stabilizing.png"))
dev.off()

# clustrering analysis
library("RColorBrewer")
library("ggplots")
distsRL <- dist(t(assay(rld)))
mat <- as.matrix(distsRL)
rownames(mat) <- colnames(mat) <- with(colData(dds), paste(condition, sampleNames,sep=" : "))
#or if you want to use condtions use:
#rownames(mat) <- colnames(mat) <- with(colData(dds),condition)
hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
dev.copy(png, paste0(outputPrefix, "-clustering.png"))
heatmap.2(mat, trace = "none", col = rev(hmcol), margin = c(13,13))

#Principal components plot shows additional but rough clustering of samples
library("genefilter")
library("ggplot2")
library("grDevices")

rv <- rowVars(assay(rld))
select <- order(rv, decreasing=T)[seq_len(min(500,length(rv)))]
pc <- prcomp(t(assay(vsd)[select,]))

# set condition
condition <- treatments
scores <- data.frame(pc$x, condition)
pcaplot <- ggplot(scores, aes(x = PC1, y = PC2, col = (factor(condition))) +
                  geom_point(size = 5) +
                  ggtitle("Principal Components") +
                  scale_colour_brewer(name = " ", palette = "Set1") +
                  theme(
                    Plot.title = element_text(face = 'bold'),
                    legend.position = c(.9,.2),
                    legend.key = element_rect(fill = 'NA'),
                    legend.text = element_text(size = 10, face = "bold"),
                    axis.text.y = element_text(colour = "Black"),
                    axis.text.x = element_text(colour = "Black"),
                    axis.title.x = element_text(face = "bold"),
                    axis.title.y = element_text(face = 'bold'),
                    panel.grid.major.x = element_blank(),
                    panel.grid.major.y = element_blank(),
                    panel.grid.minor.x = element_blank(),
                    panel.grid.minor.y = element_blank(),
                    panel.background = element_rect(color = 'black',fill = NA)
                    ))
ggsave(pcaplot, file = paste0(outputPrefix, "-ggplot2.pdf"))


# scatter plot of rlog transformations between sample conditions
# nice way to compare control and experimental samples
head(assay(rld))
# plot(log2(1+counts(dds,normalized=T)[,1:2]),col='black',pch=20,cex=0.3, main='Log2 transformed')
plot(assay(rld)[,1:3],col='#00000020',pch=20,cex=0.3, main = "rlog transformed")
plot(assay(rld)[,2:4],col='#00000020',pch=20,cex=0.3, main = "rlog transformed")

# heatmap of data
library("RcolorBrewer")
library("gplots")
# 1000 top expressed genes with heatmap.2
select <- order(rowMeans(counts(ddsClean, normalized=T)), decreasing=T)[1:1000]
my_palette <- colorRampPalette(c("blue", "white", "red"))(n=1000)
heatmap.2(assay(vsd)[select,], col=my_palette,
          scale="row", key=T, keysize=1, symkey=T,
          density.info="none", trace="none",
          cexCol=0.6, labRow=F,
          main="1000 Top Expressed Genes Heatmap")
dev.copy(png, paste0(outputPrefix, "-HEATMAP.png"))
dev.off()

##################################################################################
#Functional enrichment analysis
#

              
