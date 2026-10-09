# install.packages(c("scTenifoldNet","scTenifoldKnk","Seurat","ggplot2","dplyr","ggrepel"))

library(scTenifoldNet)
library(scTenifoldKnk)
library(Seurat)
library(ggplot2)
library(dplyr)
library(igraph)
library(ggrepel)
set.seed(12345)
options(bitmapType = "cairo")    # Adobe Illustrator 兼容 PDF

# 注意: R 4.6+ 已移除 memory.limit() (调用只返回 Inf 并警告), 物理 RAM 才是真正限制, 故不再设置。

hubGenes <- c("PLAU")   # 2026-08-24: GPAT3 OE 已完成(Aug6)，仅跑 PLAU 避免重复计算
baseDir   <- "D:/BaiduSyncdisk/daixie-colorectal/4.免疫浸润单细胞/30.scTenifoldOE"
seuratFile <- "D:/BaiduSyncdisk/daixie-colorectal/4.免疫浸润单细胞/12.Seurat/GSE132465/Seurat.Rdata"
referFile <- file.path(baseDir, "refer.scTenifoldOE.R")

# 自定义 scTenifoldOE (含 scQC + scTenifoldKnk:::strictDirection/dRegulation 后处理) 在主进程定义即可。
# 2026-08-04 决策: 改回【串行 for 循环】(不再用 doParallel 并行) —— 与 28(KO) 完全对齐。
#   原因: 每个 scTenifoldOE 内部 pcNet+tensorDecomposition 会建 ~1GB 稠密张量; 2 个并行 worker 同时
#   跑会把 16GB 笔记本内存撑爆 (报错 cannot allocate vector of size 951.0 Mb)。28 串行 5000 基因已验证可跑,
#   30 的 nc_nNet=5/nc_nCells=300 比 28 重跑的 10/1000 还轻, 串行必然放得下。
source(referFile)

# ===== 主进程只 load 一次 542MB Rdata, 预先算好每个基因的输入矩阵(~27MB/基因), 避免循环内重复加载。
load(seuratFile)
pbmc <- subset(pbmc, subset = Type == "Disease")
countMat <- GetAssayData(pbmc, layer = "counts")
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 5000)
hvgs <- VariableFeatures(pbmc)
geneData <- setNames(lapply(hubGenes, function(g){
  as.data.frame(countMat[unique(c(g, hvgs)), ])
}), hubGenes)
rm(pbmc, countMat); gc()

results <- list()

for (target_gene in hubGenes) {
  cat("\n========================================\n")
  cat("Target gene (OE):", target_gene, "\n")
  cat("========================================\n")

  geneDir <- file.path(baseDir, target_gene)
  dir.create(geneDir, showWarnings = FALSE, recursive = TRUE)
  setwd(geneDir)

  set.seed(12345 + which(hubGenes == target_gene))   # 每基因独立种子, 保证可复现

  data <- geneData[[target_gene]]                    # 用预提取的小矩阵, 不再 load 542MB

  result <- scTenifoldOE(countMatrix = data, gOE = target_gene,
                         oe_factor = 10, nc_nNet = 5, nc_nCells = 300,
                         nc_nComp = 3, nCores = 1)

  df <- result$diffRegulation
  df <- df[df$gene != target_gene, ]
  outTab <- df[df$p.adj < 0.05, ]
  write.table(outTab, file = "sigDiff.txt", sep = "\t", quote = FALSE, row.names = FALSE)
  cat("  Significant perturbed genes:", nrow(outTab), "\n")

  # 火山图
  df$log_p.adj <- -log10(pmax(df$p.adj, 1e-300))   # 防止 p.adj=0 时 -log10 出 Inf
  df$significant <- ifelse(df$p.adj < 0.05, "Significant", "Not significant")
  # 仅标注效应最强的 Top 30 (按 |Z| 排序)，其余显著基因保留红点但不标名
  label_genes <- df[df$p.adj < 0.05, , drop = FALSE]
  label_genes <- label_genes[order(-abs(label_genes$Z)), , drop = FALSE][1:30, , drop = FALSE]
  y_upper <- quantile(df$log_p.adj, 0.999, na.rm = TRUE)
  p2 <- ggplot(df, aes(x = Z, y = log_p.adj, color = significant)) +
    geom_point(alpha = 0.7, size = 1) +
    scale_color_manual(values = c("Significant" = "red", "Not significant" = "gray50")) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
    geom_text_repel(data = label_genes, aes(label = gene), size = 3, max.overlaps = 50) +
    labs(title = paste0(target_gene, " OE"), x = "Z-score", y = "-log10(p.adj)") +
    theme_classic() + coord_cartesian(ylim = c(0, y_upper)) + theme(legend.position = "none")
  pdf(file = "vol.pdf", width = 6, height = 5); print(p2); dev.off()

  # 柱状图
  top_genes <- head(df[order(-df$FC), ], 20)
  p1 <- ggplot(top_genes, aes(x = reorder(gene, FC), y = FC)) +
    geom_bar(stat = 'identity', fill = '#E27070') + coord_flip() +
    labs(title = paste0(target_gene, " OE: Top 20"), x = "Gene", y = "FC") +
    theme_minimal() + theme(plot.title = element_text(hjust = 0.5))
  pdf(file = "barplot.pdf", width = 6, height = 5); print(p1); dev.off()

  # 饼图
  sig_count <- table(df$significant)
  sig_df <- as.data.frame(sig_count)
  colnames(sig_df) <- c("category", "count")
  sig_df$percentage <- paste0(round(sig_df$count / sum(sig_df$count) * 100, 1), "%")
  p3 <- ggplot(sig_df, aes(x = "", y = count, fill = category)) +
    geom_bar(stat = "identity", width = 1) +
    geom_text(aes(label = percentage), position = position_stack(vjust = 0.5), size = 4) +
    coord_polar("y", start = 0) +
    scale_fill_manual(values = c("Significant" = "red", "Not significant" = "lightgray")) +
    labs(title = paste0(target_gene, " OE: Sig Genes"), fill = "") +
    theme_minimal() + theme(axis.title.x = element_blank(), axis.title.y = element_blank(),
                             axis.text = element_blank(), panel.grid = element_blank(),
                             legend.position = "right", plot.title = element_text(hjust = 0.5))
  pdf(file = "pie.pdf", width = 6, height = 5); print(p3); dev.off()

  results[[target_gene]] <- data.frame(gene = target_gene, n_sig = nrow(outTab), stringsAsFactors = FALSE)
  rm(result, df, outTab, label_genes, top_genes, sig_df, sig_count,
     p1, p2, p3, y_upper, data); gc()   # 每个基因跑完清场, 防止内存累积
  cat("  Output saved to:", geneDir, "\n")
}

cat("\n=== RESULTS ===\n")
for (r in results) cat(r$gene, ": ", r$n_sig, " significant genes\n")
cat("\n>>> ALL DONE: 2 genes (GPAT3/PLAU) virtual OE completed.\n")
