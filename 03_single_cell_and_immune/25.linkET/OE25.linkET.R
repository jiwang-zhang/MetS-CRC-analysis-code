# ============================================================
# Step 7c: 7 hub 基因 × 免疫细胞 Spearman 相关性
# 输入: 合并表达矩阵 + CIBERSORT结果 + hubGenes.List.txt
# 输出: linkET.pdf (相关性图) / hub_immune_cor.csv (相关表)
#
# 依赖: limma, dplyr, tidyverse, ggplot2, linkET
#        install: pak::pak("Hy4m/linkET")
# ============================================================

library(limma)
library(dplyr)
library(tidyverse)
library(ggplot2)

expFile   <- "../../第二部分DEG/25.diagnostic_model/merged_after_batch_removal.csv"  # 表达数据
geneFile  <- "hubGenes.List.txt"          # 7 hub 基因列表
immFile   <- "../23.CIBERSORT/CIBERSORT-Results.txt"    # CIBERSORT 结果
setwd("D:/BaiduSyncdisk/daixie-colorectal/4.免疫浸润单细胞/25.linkET")

# ============================================================
# 1. 读取表达矩阵
# ============================================================
rt <- read.table(expFile, header = TRUE, sep = ",", check.names = FALSE)
rt <- as.matrix(rt)
rownames(rt) <- rt[, 1]
exp <- rt[, 2:ncol(rt)]
dn  <- list(rownames(exp), colnames(exp))
exp <- matrix(as.numeric(exp), nrow = nrow(exp), dimnames = dn)
exp <- avereps(exp)

# 读取基因列表, 提取表达量
geneRT <- read.table(geneFile, header = FALSE, sep = "\t", check.names = FALSE)
hubGenes <- as.vector(geneRT[, 1])
cat("7 hub 基因:", paste(hubGenes, collapse = ", "), "\n")

# 仅肿瘤样本 (列名含 _tre)
treCols <- grepl("_tre$", colnames(exp))
data <- exp[hubGenes, treCols, drop = FALSE]
data <- t(data)
cat("肿瘤样本数:", nrow(data), "\n")

# ============================================================
# 2. 读取 CIBERSORT 结果
# ============================================================
immune <- read.table(immFile, header = TRUE, sep = "\t", check.names = FALSE, row.names = 1)
# 对齐样本
sameSample <- intersect(rownames(data), rownames(immune))
data   <- data[sameSample, , drop = FALSE]
immune <- immune[sameSample, , drop = FALSE]
immune <- immune[, apply(immune, 2, sd) > 0, drop = FALSE]
cat("对齐后:", length(sameSample), "样本,", ncol(immune), "免疫细胞\n")

# ============================================================
# 3. Spearman 相关
# ============================================================
geneCor <- data.frame()
for (cell in colnames(immune)) {
  if (sd(immune[, cell]) == 0) next
  for (gene in colnames(data)) {
    x <- as.numeric(immune[, cell])
    y <- as.numeric(data[, gene])
    ct <- cor.test(x, y, method = "spearman")
    geneCor <- rbind(geneCor, cbind(spec = gene, env = cell,
                                    r = ct$estimate, p = ct$p.value))
  }
}
geneCor$r <- as.numeric(geneCor$r)
geneCor$p <- as.numeric(geneCor$p)
geneCor$pd <- ifelse(geneCor$p < 0.05,
                     ifelse(geneCor$r > 0, "Positive", "Negative"), "Not")
geneCor$absR <- abs(geneCor$r)
geneCor <- geneCor %>% mutate(rd = cut(absR, breaks = c(-Inf, 0.2, 0.4, 0.6, Inf),
                                       labels = c("< 0.2", "0.2 - 0.4", "0.4 - 0.6", ">= 0.6")))

write.csv(geneCor, "hub_immune_correlation.csv", row.names = FALSE)
cat("相关性表已保存: hub_immune_correlation.csv\n")

# 显著关联摘要
sigRows <- geneCor[geneCor$pd != "Not", ]
if (nrow(sigRows) > 0) {
  cat("显著关联 (p<0.05):\n")
  print(sigRows[order(sigRows$absR, decreasing = TRUE),
                c("spec", "env", "r", "p", "pd")], row.names = FALSE)
}

# ============================================================
# 4. linkET 图 (需要 linkET 包)
# ============================================================
hasLinkET <- requireNamespace("linkET", quietly = TRUE)
if (hasLinkET) {
  library(linkET)
  qcorPlot <- qcorrplot(correlate(immune, method = "spearman"),
                        type = "lower", diag = FALSE) +
    geom_tile() +
    geom_couple(aes(colour = pd, size = rd),
                data = geneCor,
                curvature = nice_curvature()) +
    scale_fill_gradientn(colours = rev(RColorBrewer::brewer.pal(11, "RdBu"))) +
    scale_size_manual(values = c(1, 2, 3, 4)) +
    scale_colour_manual(values = c("#6699FFFF", "#CCCCCC99", "#E6550DFF")) +
    guides(size = guide_legend(title = "abs(Cor)",
                               override.aes = list(colour = "grey35"), order = 2),
           colour = guide_legend(title = "pvalue",
                                 override.aes = list(size = 3, linewidth = 2), order = 1),
           fill = guide_colorbar(title = "Cell-cell cor", order = 3))
  
  pdf("linkET.pdf", width = 9, height = 7)
  print(qcorPlot)
  dev.off()
  cat("linkET.pdf 已保存\n")
} else {
  # 降级: pheatmap
  library(pheatmap)
  cat("linkET 未安装, 改用 pheatmap\n")
  corMat <- matrix(NA, nrow = ncol(data), ncol = ncol(immune),
                   dimnames = list(colnames(data), colnames(immune)))
  for (g in colnames(data)) {
    for (c in colnames(immune)) {
      ok <- complete.cases(data[, g], immune[, c])
      if (sum(ok) > 5) corMat[g, c] <- cor(data[ok, g], immune[ok, c], method = "spearman")
    }
  }
  pdf("hub_immune_heatmap.pdf", width = max(8, ncol(corMat) * 0.5), height = 6)
  pheatmap(corMat, main = "7 Hub Genes vs Immune Cells (Spearman)",
           fontsize = 8, na_col = "grey90",
           cellwidth = 12, cellheight = 18)
  dev.off()
  cat("hub_immune_heatmap.pdf 已保存\n")
}