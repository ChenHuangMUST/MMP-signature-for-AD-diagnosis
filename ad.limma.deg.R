# Load packages
library(ReactomePA)
library(tidyverse)
library(data.table)
library(org.Hs.eg.db)
library(clusterProfiler)
library(biomaRt)
library(enrichplot)
library(limma)
library(ggplot2)

# Set download method for clusterProfiler package
R.utils::setOption("clusterProfiler.download.method", "auto")

# Load data files
group <- read.delim("D:/RProject/AD/NMF/2cluster.consensus.csv", sep = ",")
load("D:/RProject/AD/NMF/dat.nmf.ad.RData")
load(".RData")
go_kegg_file <- "go_kegg.Rdata"
load(go_kegg_file)

# Limma differential analysis
exp <- combat_edata[, group$X]
group <- t(cbind(t(groupq), t(groupw)))
group <- gsub("1", "group1", group, perl = TRUE)
group <- gsub("2", "group2", group, perl = TRUE)
design <- model.matrix(~0+factor(group))
colnames(design) <- levels(factor(group))
rownames(design) <- colnames(exp)

# Create contrast matrix
contrast.matrix <- makeContrasts(group1 - group2, levels = design)

# Perform differential analysis
fit <- lmFit(exp, design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
tempOutput <- topTable(fit2, coef = 1, n = Inf)
nrDEG <- na.omit(tempOutput)

# Generate heat map
choose_gene <- names(tail(sort(apply(exp, 1, mad)), 200))
choose_matrix <- exp[choose_gene, ]
choose_matrix <- t(scale(t(choose_matrix)))
pheatmap(choose_matrix)

# Generate volcano plot
DEG <- nrDEG
logFC_cutoff <- with(DEG, mean(abs(logFC)) + 2 * sd(abs(logFC)))
DEG$change <- as.factor(ifelse(DEG$P.Value < 0.05 & abs(DEG$logFC) > logFC_cutoff,
                                ifelse(DEG$logFC > logFC_cutoff, 'UP', 'DOWN'), 'NOT'))
ggplot(data = DEG, aes(x = logFC, y = -log10(P.Value), color = change)) +
  geom_point(alpha = 0.4, size = 1.75) +
  theme_set(theme_bw(base_size = 20)) +
  xlab("log2 fold change") + ylab("-log10 p-value") +
  ggtitle(paste0('Cutoff for logFC is ', round(logFC_cutoff, 3),
                 '\nThe number of up gene is ', nrow(DEG[DEG$change == 'UP',]),
                 '\nThe number of down gene is ', nrow(DEG[DEG$change == 'DOWN',]))) +
  theme(plot.title = element_text(size = 15, hjust = 0.5)) +
  scale_colour_manual(values = c('blue', 'black', 'red'))

# Perform pathway analysis on top 1000 genes
gene <- head(rownames(nrDEG), 1000)
gene.df <- bitr(gene, fromType = "SYMBOL", toType = c("ENSEMBL", "ENTREZID"), OrgDb = org.Hs.eg.db)

kk <- enrichKEGG(gene = gene.df$ENTREZID, organism = 'hsa', pvalueCutoff = 0.05)
kk2 <- gseKEGG(geneList = aaa[, 2], organism = 'hsa', nPerm = 1000, minGSSize = 120, pvalueCutoff = 0.05, verbose = FALSE)
gseaplot(kk2, geneSetID = "hsa04740", by = "all")

# Save data files
save(aaa, choose_matrix, contrast.matrix, DEG, design, exp, fit, fit2, g, gene_map, gene.df, genelist_input, group, file = go_kegg_file)

# Perform GSEA analysis
Go_gseresult <- gseGO(geneList, 'org.Hs.eg.db', keyType = "ENTREZID", ont = "all",
                      nPerm = 1000, minGSSize = 10, maxGSSize = 1000, pvalueCutoff = 1)
KEGG_gseresult <- gseKEGG(geneList, nPerm = 1000, minGSSize = 10, maxGSSize = 1000, pvalueCutoff = 1)
Go_Reactomeresult <- gsePathway(geneList, nPerm = 1000, minGSSize = 10, maxGSSize = 1000, pvalueCutoff=0.05)
