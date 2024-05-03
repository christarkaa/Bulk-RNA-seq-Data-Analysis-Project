# Set working directory
getwd()
setwd("/Users/christophertarkaa/Desktop/bulk")

# Install packages if not already installed
BiocManager::install("biomaRt")
aBiocManager::install("clusterProfiler")
BiocManager::install("enrichplot")
BiocManager::install("biomartr")
BiocManager::install("org.Hs.eg.db")
BiocManager::install("DOSE")

# Load packages
library(biomaRt)
library(clusterProfiler)
library(enrichplot)
library(biomartr)
library(tidyverse)
library(DOSE)
suppressPackageStartupMessages(library("org.Hs.eg.db"))

# Load data
resdata <- read.csv("diffexpr-results.csv")

# Retrieve information about Homo sapiens
genes <- biomartr::organismBM(organism = "Homo sapiens")
genes
View(genes) 

# Retrieve gene attributes
hsapiens_attributes = 
  biomartr::organismAttributes("Homo sapiens") %>% 
  filter(dataset == "hsapiens_gene_ensembl")
hsapiens_attributes
View(hsapiens_attributes)

# Retrieve filters
filters <- biomartr::getFilters(mart    = "ENSEMBL_MART_ENSEMBL", 
                          dataset = "hsapiens_gene_ensembl")
View(filters)

# Query Biomart
gene_set <- resdata$Gene
result_BM <- biomartr::biomart(
  genes = gene_set,
  mart = "ENSEMBL_MART_ENSEMBL",
  dataset = "hsapiens_gene_ensembl",
  attributes = c("hgnc_symbol", "entrezgene_id"),
  filters = "ensembl_gene_id"
)
head(result_BM)

# Convert result_BM to a data frame
result_df <- as.data.frame(result_BM)

# Rename the columns for easier access
colnames(result_df) <- c("Gene", "HGNC_symbol", "Entrez_gene_id")

# Perform functional enrichment analysis for GO terms
geneID <- result_df$Entrez_gene_id
ego <- enrichGO(
  gene = geneID,
  OrgDb = org.Hs.eg.db,  # Human genome annotation database
  keyType = "ENTREZID",  # Use Entrez gene IDs as keys
  ont = "BP",            # Biological Process ontology. Can either be BP, CC or MF
  pAdjustMethod = "BH",  #  
  pvalueCutoff = 0.05,   # P-value cutoff for significance
  qvalueCutoff = 0.05,   # Adjusted p-value (FDR) cutoff for significance
  readable      = TRUE)
ego

# Visualize the GO enrichment results
png("GO_BP.png", w=1000, h=1000, pointsize=20)
dotplot(ego, showCategory = 15) + ggtitle("Gene Enrichment Plot") + theme_bw() # Show top 15 enriched terms in a dotplot
dev.off()

# Perform functional enrichment analysis for KEGG pathways
kegg_enrich <- enrichKEGG(
  gene = geneID,
  organism = "hsa",       # Specify organism code for Homo sapiens
  pvalueCutoff = 0.05,    # P-value cutoff for significance
  qvalueCutoff = 0.05,     # Adjusted p-value (FDR) cutoff for significance
)
View(kegg_enrich)

# Visualize the KEGG pathway enrichment results
png("KEEG.png", w=1000, h=1000, pointsize=20)
dotplot(kegg_enrich, showCategory = 10, color = "qvalue", size = "Count")  # Show top 10 enriched pathways in a dotplot
dev.off()

