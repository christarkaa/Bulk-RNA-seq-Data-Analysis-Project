getwd()
setwd("/Users/christophertarkaa/Desktop/bulk")

# Tutorial from https://gist.github.com/stephenturner/f60c1934405c127f09a6
# Load the count tables
counts_SRR28420795 <- read.table("counts.SRR28420795.txt",  header=TRUE, row.names = 1)
counts_SRR28420796 <- read.table("counts.SRR28420796.txt",  header=TRUE, row.names = 1)
counts_SRR28420797 <- read.table("counts.SRR28420797.txt",  header=TRUE, row.names = 1)
counts_SRR28420798 <- read.table("counts.SRR28420798.txt",  header=TRUE, row.names = 1)

# Combine count files into a single count file
countdata <- cbind(counts_SRR28420795, counts_SRR28420796, counts_SRR28420797, counts_SRR28420798)

# Remove first five columns (chr, start, end, strand, length) of every sample
countdata <- countdata[, seq(6, ncol(countdata), by = 6)]
View(countdata)

# Modify column names
colnames(countdata) <- gsub("^Mapping\\.|\\.sorted\\.bam$", "", colnames(countdata))
View(countdata)

# Convert to matrix
countdata <- as.matrix(countdata)
head(countdata)

# Assign condition
condition = c("Treated", "Control", "Treated", "Control")
 
# Analysis with DESeq2 ----------------------------------------------------
library(DESeq2)

# Create a coldata frame and instantiate the DESeqDataSet. See ?DESeqDataSetFromMatrix
(coldata <- data.frame(row.names=colnames(countdata), condition))
dds <- DESeqDataSetFromMatrix(countData=countdata, colData=coldata, design=~condition)
dds

# Run the DESeq pipeline
dds <- DESeq(dds)

# Plot dispersions
png("qc-dispersions.png", 1000, 1000, pointsize=20)
plotDispEsts(dds, main="Dispersion plot")
dev.off()

# Regularized log transformation for clustering/heatmaps, etc
rld <- rlogTransformation(dds)
head(assay(rld))
hist(assay(rld))

# Sample distance heatmap
sampleDists <- as.matrix(dist(t(assay(rld))))
library(gplots)
png("qc-heatmap-samples.png", w=1000, h=1000, pointsize=20)
heatmap.2(as.matrix(sampleDists), key=F, trace="none",
          col=colorpanel(100, "green", "yellow"),
          ColSideColors=mycols[condition], RowSideColors=mycols[condition],
          margin=c(15, 15), main="Sample Distance Matrix")
dev.off()

# Principal Component Analysis (PCA)
png("PCA.png", w=1000, h=1000)
DESeq2::plotPCA(rld, intgroup="condition")
dev.off()

# Let's customize it
library(ggplot2)
library(genefilter)
library(grDevices)

rv <- rowVars(assay(rld))
select <- order(rv, decreasing=T)[seq_len(min(500,length(rv)))]
pc <- prcomp(t(assay(vsd)[select,]))
scores <- data.frame(pc$x, condition)

png("PCA1.png", w=1000, h=1000)
(pcaplot <- ggplot(scores, aes(x = PC1, y = PC2, col = (factor(condition))))
  + geom_point(size = 5)
  + ggtitle("Principal Components")
  + scale_colour_brewer(name = " ", palette = "Set1") 
  + geom_text(aes(label = rownames(scores)), size = 4, vjust = -1.5)
  + theme(
    plot.title = element_text(face = 'bold'),
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
dev.off()

# Heatmap of data
select <- order(rowMeans(counts(dds,normalized=T)),decreasing=T)[1:1000]
png("Heatmap.png", w=1000, h=1000, pointsize=20)
my_palette <- colorRampPalette(c("blue",'white','red'))(n=1000)
heatmap.2(assay(vsd)[select,], col=my_palette,
          scale="row", key=T, keysize=1, symkey=T,
          density.info="none", trace="none",
          cexCol=0.6, labRow=F,
          main="1000 Top Expressed Genes Heatmap")
dev.off()

# Get differential expression results
res <- results(dds)
table(res$padj<0.05)
## Order by adjusted p-value
res <- res[order(res$padj), ]
## Merge with normalized count data
resdata <- merge(as.data.frame(res), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
names(resdata)[1] <- "Gene"
head(resdata)
## Write results
write.csv(resdata, file="diffexpr-results.csv")

## Examine plot of p-values
hist(res$pvalue, breaks=50, col="grey")

## MA plot
## Could do with built-in DESeq2 function:
## DESeq2::plotMA(dds, ylim=c(-1,1), cex=1)
## Let's customise
install.packages("calibrate")
library(calibrate)
maplot <- function (res, thresh=0.05, labelsig=TRUE, textcx=1, ...) {
  with(res, plot(baseMean, log2FoldChange, pch=20, cex=.5, log="x", ...))
  with(subset(res, padj<thresh), points(baseMean, log2FoldChange, col="red", pch=20, cex=1.5))
  if (labelsig) {
    require(calibrate)
    with(subset(res, padj<thresh), textxy(baseMean, log2FoldChange, labs=Gene, cex=textcx, col=2))
  }
}
png("diffexpr-maplot.png", 1500, 1000, pointsize=20)
maplot(resdata, main="MA Plot")
dev.off()

## Volcano plot with "significant" genes labeled
volcanoplot <- function (res, lfcthresh=2, sigthresh=0.05, main="Volcano Plot", legendpos="bottomright", labelsig=TRUE, textcx=1, ...) {
  with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main=main, ...))
  with(subset(res, padj<sigthresh ), points(log2FoldChange, -log10(pvalue), pch=20, col="red", ...))
  with(subset(res, abs(log2FoldChange)>lfcthresh), points(log2FoldChange, -log10(pvalue), pch=20, col="orange", ...))
  with(subset(res, padj<sigthresh & abs(log2FoldChange)>lfcthresh), points(log2FoldChange, -log10(pvalue), pch=20, col="green", ...))
  if (labelsig) {
    require(calibrate)
    with(subset(res, padj<sigthresh & abs(log2FoldChange)>lfcthresh), textxy(log2FoldChange, -log10(pvalue), labs=Gene, cex=textcx, ...))
  }
  legend(legendpos, xjust=1, yjust=1, legend=c(paste("FDR<",sigthresh,sep=""), paste("|LogFC|>",lfcthresh,sep=""), "both"), pch=20, col=c("red","orange","green"))
}
png("diffexpr-volcanoplot.png", 1200, 1000, pointsize=20)
volcanoplot(resdata, lfcthresh=1, sigthresh=0.05, textcx=.8, xlim=c(-2.3, 2))
dev.off()
