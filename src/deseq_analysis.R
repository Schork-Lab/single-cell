library("DESeq2")
library("DESeq")
countData = read.csv('../../data/htseq_counts.R-ready.csv', sep=',', row.names=1)
samples = sort(colnames(countData))
sample_types = c("Neuronal","Neuronal","Neuronal","Neuronal","Neuronal","Neuronal","Non-neuronal","Non-neuronal","Non-neuronal","Non-neuronal","Total RNA","Total RNA", "Total RNA", "Total RNA")
colData = data.frame(samples, sample_types, row.names=1)
dds <- DESeqDataSetFromMatrix(countData = countData,
colData = colData,
design = ~ sample_types)
dds <- DESeq(dds)
res <- results(dds)
resMLE <- results(dds, addMLE=TRUE, contrast=c("sample_types", "Neuronal","Non-neuronal"))
rld <- rlog(dds)
DESeq2::plotPCA(rld, intgroup=c("sample_types"))
distsRL <- dist(t(assay(rld)))
mat <- as.matrix(distsRL)
rownames(mat) <- colnames(mat) <- with(colData(dds),
paste(samples, sep=" : "))
hc <- hclust(distsRL)
heatmap.2(mat, Rowv=as.dendrogram(hc),
symm=TRUE, trace="none",
col = rev(hmcol), margin=c(13, 13))
