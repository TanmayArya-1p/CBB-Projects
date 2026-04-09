library(clusterProfiler)
library(org.Hs.eg.db)

genes <- read.table("gene_list.txt", header=FALSE)$V1

ego <- enrichGO(gene          = genes,
                OrgDb         = org.Hs.eg.db,
                keyType       = "SYMBOL",
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.05)

png("gene_enrichment_plot.png", width=10, height=8, units="in", res=300)
dotplot(ego, showCategory=15, title="BioFunctions of CRE containing Genes")
dev.off()
