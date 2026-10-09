# ============================================================
# 扩展免疫浸润分析：3 hub → 7 hub
# 补齐 PPARD / PLAU / LPCAT4 / TNFRSF21 的 GSEA + ssGSEA barplot + linkET
#
# 已有（不重跑）：
#   - CIBERSORT 结果 (23.CIBERSORT/CIBERSORT-Results.txt)
#   - ssGSEA 评分矩阵 (22.GSVA/gsva_score.csv) — 样本级，与基因无关
#   - GPAT3/PLAU 的 GSEA 结果和 barplot (下游聚焦两交叉复现 hub; 其余 7 hub 见各脚本)
#
# 本次新增：
#   1) ssGSEA barplot: 4 个新基因 (复用 gsva_score.csv)
#   2) GSEA: 4 个新基因 (PPARD/PLAU/LPCAT4/TNFRSF21)
#   3) linkET: 7 hub × 22 免疫细胞 Spearman 相关 (全重跑)
# ============================================================

library(reshape2)
library(ggpubr)
library(limma)
# GSEABase/GSVA 不需要 — gsva_score.csv 已存在，只读 CSV 画 barplot
library(org.Hs.eg.db)
library(clusterProfiler)
library(enrichplot)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(pheatmap)

# ============================================================
# 公共：读取表达矩阵
# ============================================================
expFile <- "D:/BaiduSyncdisk/daixie-colorectal/第二部分DEG/25.diagnostic_model/merged_after_batch_removal.csv"

cat("===== 读取表达矩阵 =====\n")
rt <- read.table(expFile, header = TRUE, sep = ",", check.names = FALSE)
rt <- as.matrix(rt)
rownames(rt) <- rt[, 1]
exp <- rt[, 2:ncol(rt)]
dn <- list(rownames(exp), colnames(exp))
exp <- matrix(as.numeric(exp), nrow = nrow(exp), dimnames = dn)
exp <- avereps(exp)

# 仅肿瘤样本
treCols <- grepl("_tre$", colnames(exp))
data.tre <- exp[, treCols, drop = FALSE]
cat("肿瘤样本:", ncol(data.tre), "| 基因:", nrow(data.tre), "\n")

# 7 hub 基因
hubGenes <- c("GPAT3", "SGK1", "TMCC3", "PPARD", "PLAU", "LPCAT4", "TNFRSF21")
newGenes <- c("PPARD", "PLAU", "LPCAT4", "TNFRSF21")  # 本次新增

# 检查 7 hub 是否都在矩阵中
for (g in hubGenes) {
  cat(g, "in matrix:", g %in% rownames(data.tre), "\n")
}

# ============================================================
# Part 1: ssGSEA barplot — 4 个新基因 (复用已有 gsva_score.csv)
# ============================================================
cat("\n===== Part 1: ssGSEA barplot (新基因) =====\n")

setwd("D:/BaiduSyncdisk/daixie-colorectal/4.免疫浸润单细胞/22.GSVA")

# 读取已有的 ssGSEA 评分矩阵
gsvaScore <- read.csv("gsva_score.csv", header = TRUE, row.names = 1, check.names = FALSE)
gsvaScore <- as.matrix(gsvaScore)
cat("ssGSEA 评分矩阵:", nrow(gsvaScore), "通路 x", ncol(gsvaScore), "样本\n")

normalize <- function(x) { return((x - min(x)) / (max(x) - min(x))) }
# 已归一化，但确保
if (max(gsvaScore) > 1.01) {
  gsvaScore <- t(apply(gsvaScore, 1, normalize))
  colnames(gsvaScore) <- colnames(read.csv("gsva_score.csv", header = TRUE, row.names = 1, check.names = FALSE))
}
gsvaScore <- gsvaScore[apply(gsvaScore, 1, sd, na.rm = TRUE) > 0.01, ]

for (gene in newGenes) {
  if (!gene %in% rownames(data.tre)) {
    cat("跳过", gene, ": 不在表达矩阵中\n")
    next
  }
  geneExp <- data.tre[gene, ]
  med <- median(geneExp, na.rm = TRUE)
  lowName  <- colnames(data.tre)[geneExp < med]
  highName <- colnames(data.tre)[geneExp >= med]

  # 确保样本对齐
  lowName  <- intersect(lowName, colnames(gsvaScore))
  highName <- intersect(highName, colnames(gsvaScore))

  lowScore  <- gsvaScore[, lowName,  drop = FALSE]
  highScore <- gsvaScore[, highName, drop = FALSE]
  scoreData <- cbind(lowScore, highScore)
  conNum    <- ncol(lowScore)
  treatNum  <- ncol(highScore)
  Type <- c(rep("Low", conNum), rep("High", treatNum))

  outTab <- data.frame()
  for (i in rownames(scoreData)) {
    test <- t.test(scoreData[i, ] ~ Type)
    t <- test$statistic
    pvalue <- test$p.value
    if (test$estimate[2] > test$estimate[1]) t <- abs(t) else t <- -abs(t)
    Sig <- ifelse(pvalue > 0.05, "Not", ifelse(t > 0, "Up", "Down"))
    outTab <- rbind(outTab, data.frame(Pathway = i, t = t, pvalue = pvalue, Sig = Sig))
  }

  notSigTab <- outTab[outTab$Sig == "Not", ]
  notSigTab <- notSigTab[order(notSigTab$t), ]
  sigTab <- outTab[outTab$Sig != "Not", ]
  sigTab <- sigTab[order(sigTab$t), ]

  if (nrow(sigTab) > 20) {
    outTab <- rbind(sigTab[c(1:10, ((nrow(sigTab) - 9):nrow(sigTab))), ],
                    notSigTab[c(1:5, ((nrow(notSigTab) - 4):nrow(notSigTab))), ])
  } else {
    outTab <- rbind(sigTab,
                    notSigTab[c(1:5, ((nrow(notSigTab) - 4):nrow(notSigTab))), ])
  }

  outTab$t <- as.numeric(outTab$t)
  outTab$Sig <- factor(outTab$Sig, levels = c("Down", "Not", "Up"))

  pdf(file = paste0(gene, ".barplot.pdf"), width = 10, height = 7)
  gg1 <- ggbarplot(outTab, x = "Pathway", y = "t", fill = "Sig", color = "white",
                   palette = c("royalblue1", "gray80", "tomato"),
                   sort.val = "asc", sort.by.groups = TRUE,
                   rotate = TRUE, title = gene, legend.title = "Group",
                   legend = "right", xlab = "", ylab = "t value of GSVA score",
                   x.text.angle = 60)
  print(gg1)
  dev.off()
  cat(gene, "ssGSEA barplot 已保存\n")
}

cat("ssGSEA barplot 全部完成\n")

# ============================================================
# Part 2: GSEA — 4 个新基因
# ============================================================
cat("\n===== Part 2: GSEA (新基因) =====\n")

setwd("D:/BaiduSyncdisk/daixie-colorectal/4.免疫浸润单细胞/21.GSEA")

gmtFile <- "c2.cp.kegg.Hs.symbols.gmt"
gmt <- read.gmt(gmtFile)

# 用肿瘤样本的表达矩阵
exp.gsea <- data.tre
exp.gsea <- exp.gsea[rowMeans(exp.gsea, na.rm = TRUE) > 0, ]
cat("GSEA 表达矩阵:", nrow(exp.gsea), "基因 x", ncol(exp.gsea), "样本\n")

for (gene in newGenes) {
  if (!gene %in% rownames(exp.gsea)) {
    cat("跳过", gene, ": 不在表达矩阵中\n")
    next
  }
  geneExp <- exp.gsea[gene, ]
  dataL <- exp.gsea[, geneExp <  median(geneExp), drop = FALSE]
  dataH <- exp.gsea[, geneExp >= median(geneExp), drop = FALSE]
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

  # highExp 图
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

  # lowExp 图
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
  cat(gene, "GSEA 完成 | 显著通路:", nrow(kkTab), "\n")
}

cat("GSEA 全部完成\n")

# ============================================================
# Part 3: linkET — 7 hub × 22 免疫细胞 Spearman 相关 (全重跑)
# ============================================================
cat("\n===== Part 3: linkET (7 hub) =====\n")

setwd("D:/BaiduSyncdisk/daixie-colorectal/4.免疫浸润单细胞/25.linkET")

# 读取 CIBERSORT 结果
immFile <- "../23.CIBERSORT/CIBERSORT-Results.txt"
immune <- read.table(immFile, header = TRUE, sep = "\t", check.names = FALSE, row.names = 1)

# 提取 7 hub 在肿瘤样本中的表达
hubExp <- data.tre[hubGenes, , drop = FALSE]
hubExp <- t(hubExp)

# 对齐样本
sameSample <- intersect(rownames(hubExp), rownames(immune))
hubExp <- hubExp[sameSample, , drop = FALSE]
immune <- immune[sameSample, , drop = FALSE]
immune <- immune[, apply(immune, 2, sd, na.rm = TRUE) > 0, drop = FALSE]
cat("对齐后:", length(sameSample), "样本,", ncol(immune), "免疫细胞\n")

# Spearman 相关
geneCor <- data.frame()
for (cell in colnames(immune)) {
  if (sd(immune[, cell], na.rm = TRUE) == 0) next
  for (gene in colnames(hubExp)) {
    x <- as.numeric(immune[, cell])
    y <- as.numeric(hubExp[, gene])
    ok <- complete.cases(x, y)
    if (sum(ok) < 5) next
    ct <- cor.test(x[ok], y[ok], method = "spearman")
    geneCor <- rbind(geneCor, data.frame(
      spec = gene, env = cell,
      r = ct$estimate, p = ct$p.value
    ))
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
cat("7 hub × 免疫细胞相关性表已保存: hub_immune_correlation.csv\n")

# 显著关联摘要
sigRows <- geneCor[geneCor$pd != "Not", ]
cat("显著关联 (p<0.05):", nrow(sigRows), "/", nrow(geneCor), "\n")
if (nrow(sigRows) > 0) {
  print(sigRows[order(sigRows$absR, decreasing = TRUE),
                c("spec", "env", "r", "p", "pd")], row.names = FALSE)
}

# 每个基因的显著关联数
cat("\n各基因显著关联数:\n")
for (g in hubGenes) {
  n <- sum(geneCor$spec == g & geneCor$pd != "Not")
  cat(g, ":", n, "/", sum(geneCor$spec == g), "\n")
}

# linkET 图 (需要 linkET 包)
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
  cat("linkET 未安装, 改用 pheatmap\n")
  corMat <- matrix(NA, nrow = length(hubGenes), ncol = ncol(immune),
                   dimnames = list(hubGenes, colnames(immune)))
  for (g in hubGenes) {
    for (c in colnames(immune)) {
      ok <- complete.cases(hubExp[, g], immune[, c])
      if (sum(ok) > 5) corMat[g, c] <- cor(hubExp[ok, g], immune[ok, c], method = "spearman")
    }
  }
  pdf("hub_immune_heatmap.pdf", width = max(8, ncol(corMat) * 0.5), height = 6)
  pheatmap(corMat, main = "7 Hub Genes vs Immune Cells (Spearman)",
           fontsize = 8, na_col = "grey90",
           cellwidth = 12, cellheight = 18)
  dev.off()
  cat("hub_immune_heatmap.pdf 已保存\n")
}

cat("\n===== 全部完成 =====\n")
