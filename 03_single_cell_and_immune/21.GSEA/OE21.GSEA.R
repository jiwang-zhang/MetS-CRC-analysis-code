# ============================================================
# Step 7c: GSEA -- hub-gene high/low expression group KEGG enrichment
# Adapted from biowolf template; aligned with 22.GSVA convention:
#   - expression matrix = project merged batch-removed matrix (GPL570, tumor _tre only)
#   - gene list = hubGenes.List.txt (GPAT3/PLAU; 下游聚焦两交叉复现 hub)
#   - setwd = correct BaiduSyncdisk path
# Depends: limma, org.Hs.eg.db, clusterProfiler, enrichplot
# ============================================================

library(limma)
library(org.Hs.eg.db)
library(clusterProfiler)
library(enrichplot)

expFile  <- "../../第二部分DEG/25.diagnostic_model/merged_after_batch_removal.csv"
geneFile <- "hubGenes.List.txt"
gmtFile  <- "c2.cp.kegg.Hs.symbols.gmt"
setwd("D:/BaiduSyncdisk/daixie-colorectal/4.免疫浸润单细胞/21.GSEA")

rt <- read.csv(expFile, header = TRUE, check.names = FALSE, row.names = 1)
rt <- as.matrix(rt)
exp <- rt
treCols <- grepl("_tre$", colnames(exp))
exp <- exp[, treCols, drop = FALSE]
cat("Expression matrix (tumor only):", nrow(exp), "genes x", ncol(exp), "samples\n")

exp <- avereps(exp)
exp <- exp[rowMeans(exp) > 0, ]

geneRT <- read.table(geneFile, header = FALSE, sep = "\t", check.names = FALSE)
hubGenes <- as.vector(geneRT[, 1])
cat("Hub genes:", paste(hubGenes, collapse = ", "), "\n")

gmt <- read.gmt(gmtFile)

for (gene in hubGenes) {
  if (!gene %in% rownames(exp)) {
    cat("Skip", gene, ": not in expression matrix\n")
    next
  }
  geneExp <- exp[gene, ]
  dataL <- exp[, geneExp <  median(geneExp), drop = FALSE]
  dataH <- exp[, geneExp >= median(geneExp), drop = FALSE]
  meanL <- rowMeans(dataL)
  meanH <- rowMeans(dataH)
  meanL[meanL < 0.00001] <- 0.00001
  meanH[meanH < 0.00001] <- 0.00001
  logFC <- meanH - meanL
  logFC <- sort(logFC, decreasing = TRUE)

  kk <- GSEA(logFC, TERM2GENE = gmt, pvalueCutoff = 1)
  kkTab <- as.data.frame(kk)
  kkTab <- kkTab[kkTab$pvalue < 0.05, ]
  write.table(kkTab, file = paste0(gene, ".result.txt"), sep = "\t", quote = FALSE, row.names = FALSE)

  termNum <- 5
  kkUp <- kkTab[kkTab$NES > 0, ]
  if (nrow(kkUp) > 0) {
    if (nrow(kkUp) < termNum) termNum <- nrow(kkUp)
    showTerm <- row.names(kkUp)[1:termNum]
    gseaplot <- gseaplot2(kk, showTerm, base_size = 8, title = paste0("Enriched in high ", gene, " group"))
    pdf(file = paste0(gene, ".highExp.pdf"), width = 6.5, height = 5.5)
    print(gseaplot)
    dev.off()
  }

  termNum <- 5
  kkDown <- kkTab[kkTab$NES < 0, ]
  if (nrow(kkDown) > 0) {
    if (nrow(kkDown) < termNum) termNum <- nrow(kkDown)
    showTerm <- row.names(kkDown)[1:termNum]
    gseaplot <- gseaplot2(kk, showTerm, base_size = 8, title = paste0("Enriched in low ", gene, " group"))
    pdf(file = paste0(gene, ".lowExp.pdf"), width = 6.5, height = 5.5)
    print(gseaplot)
    dev.off()
  }
  cat(gene, "GSEA done\n")
}

cat("\n=== GSEA all done ===\n")
