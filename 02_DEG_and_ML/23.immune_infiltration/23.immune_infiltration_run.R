# ============================================================
# 23.immune_infiltration — 免疫浸润分析（管线第7步）
# ------------------------------------------------------------
# 输入 : 合并去批次 CRC 表达矩阵
#         D:/BaiduSyncdisk/daixie-colorectal/第二部分DEG/08.两数据集合并并去除批次/
#           BatchEffect_output_20260720_142148/merged_after_batch_removal.csv
#         结构: 首列 geneSymbol, 其余 748 列为样本 (log2 标准化+去批次)
#         样本名后缀: _tre=肿瘤, _con=正常 (687 肿瘤 / 61 正常)
# 方法 : ESTIMATE + ssGSEA(GSVA) + MCPcounter  (+ 可选 CIBERSORT)
# 目的 : 看 7 个 ML 锁定的 hub 与免疫细胞丰度/基质-免疫评分的关联
# 运行 : RStudio 直接 Source（中文路径 OK，勿在命令行跑）
# 注意 : 输入文件就在 BaiduSyncdisk 同步盘副本里(已验证 21755 基因 x 748 样本,
#        7 hub 全在), 用同步盘绝对路径即可; 不要连到已失效的 D:/daixie-colorectal。
# ============================================================

## ---- 0. 路径与参数 ----
input_mat <- "D:/BaiduSyncdisk/daixie-colorectal/第二部分DEG/08.两数据集合并并去除批次/BatchEffect_output_20260720_142148/merged_after_batch_removal.csv"
outdir    <- "D:/BaiduSyncdisk/daixie-colorectal/第二部分DEG/23.immune_infiltration"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
setwd(outdir)

## 7 个 ML 锁定的 hub（LASSO / RF / SVM-RFE 三交叉）
hub7 <- c("GPAT3","SGK1","TMCC3","PPARD","PLAU","LPCAT4","TNFRSF21")

## ---- 0.5 安装依赖（首次跑取消注释即可；之后可注释掉）----
# if (!require("BiocManager", quietly=TRUE)) install.packages("BiocManager")
# BiocManager::install(c("GSVA","estimate"), ask=FALSE)
# install.packages(c("MCPcounter","pheatmap","ggplot2","reshape2"))
# MCPcounter 还依赖 Bioconductor 的 Biobase: BiocManager::install("Biobase")

## ---- 1. 读取表达矩阵 ----
suppressMessages({ library(utils) })
raw <- read.csv(input_mat, row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)
expr_mat <- as.matrix(raw)                      # 基因(行) x 样本(列)
cat("表达矩阵维度:", nrow(expr_mat), "基因 x", ncol(expr_mat), "样本\n")

## 表型: 从样本名后缀判断 肿瘤/正常
group <- ifelse(grepl("_tre$", colnames(expr_mat)), "Tumor",
         ifelse(grepl("_con$", colnames(expr_mat)), "Normal", "Other"))
cat("样本分组: Tumor=", sum(group=="Tumor"),
    " Normal=", sum(group=="Normal"),
    " Other=", sum(group=="Other"), "\n")

## 校验 7 hub 都在矩阵里
miss <- setdiff(hub7, rownames(expr_mat))
if (length(miss) > 0) warning("以下 hub 不在矩阵中: ", paste(miss, collapse=", "))
stopifnot(all(hub7 %in% rownames(expr_mat)))

## ---- 2. ESTIMATE (基质/免疫评分) ----
est_out <- NULL
if (require("estimate", quietly = TRUE)) {
  cat("运行 ESTIMATE...\n")
  est <- estimate::estimateScore(expr_mat, platform = "affymetrix")
  if (nrow(est) == 3) est <- t(est)            # 规范化成 样本 x 3
  est_out <- as.data.frame(est)
  colnames(est_out) <- c("StromalScore","ImmuneScore","ESTIMATEScore")
  cat("  ESTIMATE 完成, 维度:", nrow(est_out), "x", ncol(est_out), "\n")
} else {
  warning("estimate 包未安装 -> 跳过 ESTIMATE（可在 0.5 段安装后重跑）")
}

## ---- 3. ssGSEA (GSVA) — 免疫细胞浸润分辨率 ----
## 注: 以下是简化版免疫标记集(每类几个经典标记), 出版级请替换为
##     MSigDB C7 (immunologic signatures) 或 CIBERSORT LM22 完整基因集。
immune_gs <- list(
  Macrophage = c("CD68","CD163","MRC1","IL10","TGFB1","MSR1","CD80","CD86"),
  Neutrophil = c("S100A8","S100A9","CXCR2","FCGR3B","CEACAM1","ELANE"),
  CD8_T      = c("CD8A","CD8B","GZMB","GZMK","PRF1","CCR7","TBX21"),
  CD4_T      = c("CD4","IL7R","CD40LG","LTB","TRBC1"),
  Treg       = c("FOXP3","IL2RA","CTLA4","TIGIT","IKZF2"),
  B_cell     = c("CD19","MS4A1","CD79A","IGHG1","IGHD","CD22"),
  NK         = c("NCAM1","NCR1","KLRD1","KLRC1","XCL1"),
  Dendritic  = c("CD1A","CD1C","FCER1A","CLEC9A","IRF7","BDCA2"),
  Mast       = c("TPSAB1","TPSB2","KIT","MS4A2","CPA3"),
  Monocyte   = c("CD14","LYZ","FCN1","S100A12","ITGAM")
)
## 只保留矩阵里实际存在的基因, 丢弃空集
immune_gs <- lapply(immune_gs, function(g) g[g %in% rownames(expr_mat)])
immune_gs <- immune_gs[sapply(immune_gs, length) > 0]
cat("ssGSEA 使用细胞类型:", paste(names(immune_gs), collapse=", "), "\n")

ssg_out <- NULL
if (require("GSVA", quietly = TRUE)) {
  cat("运行 ssGSEA (GSVA)...\n")
  ssg <- gsva(expr_mat, gene.sets = immune_gs, method = "ssgsea",
              kcdf = "Gaussian", parallel.sz = 1, verbose = FALSE)
  if (nrow(ssg) == length(immune_gs)) ssg <- t(ssg)   # 规范化成 样本 x 细胞类型
  ssg_out <- as.data.frame(ssg)
  cat("  ssGSEA 完成, 维度:", nrow(ssg_out), "x", ncol(ssg_out), "\n")
} else {
  warning("GSVA 未安装 -> 跳过 ssGSEA")
}

## ---- 4. MCPcounter (10 类免疫细胞/基质群体) ----
mcp_out <- NULL
if (require("MCPcounter", quietly = TRUE)) {
  cat("运行 MCPcounter...\n")
  mcp <- MCPcounter::estimateMCF(expr_mat)
  if (nrow(mcp) != ncol(expr_mat)) mcp <- t(mcp)     # estimateMCF 返回 群体x样本 -> 转置成 样本x群体
  mcp_out <- as.data.frame(mcp)
  cat("  MCPcounter 完成, 维度:", nrow(mcp_out), "x", ncol(mcp_out), "\n")
} else {
  warning("MCPcounter 未安装 -> 跳过 (CRAN 安装: install.packages('MCPcounter'))")
}

## ---- 5. 合并免疫评分 + hub 表达 ----
hub_expr <- t(expr_mat[hub7, ])                     # 样本 x 7 hub
colnames(hub_expr) <- paste0("hub_", hub7)
immune_all <- do.call(cbind, Filter(Negate(is.null), list(est_out, ssg_out, mcp_out)))
dat <- data.frame(Sample = rownames(hub_expr), Group = group,
                  hub_expr, immune_all, stringsAsFactors = FALSE)
rownames(dat) <- dat$Sample
write.csv(dat, "immune_infiltration_scores.csv", row.names = FALSE)
cat("已写出 immune_infiltration_scores.csv, 维度:", nrow(dat), "x", ncol(dat), "\n")

## ---- 6. 关联分析 ----
## 6a. 肿瘤组内: 7 hub 表达 vs 各免疫指标的 Spearman 相关
tum <- dat[dat$Group == "Tumor", ]
hub_cols  <- grep("^hub_", colnames(tum), value = TRUE)
imm_cols  <- setdiff(colnames(tum), c("Sample","Group", hub_cols))
cor_mat <- matrix(NA, nrow = length(hub_cols), ncol = length(imm_cols),
                  dimnames = list(hub_cols, imm_cols))
p_mat   <- cor_mat
for (i in hub_cols) for (j in imm_cols) {
  x <- tum[[i]]; y <- tum[[j]]
  ok <- complete.cases(x, y)
  if (sum(ok) > 5) {
    r <- cor.test(x[ok], y[ok], method = "spearman")
    cor_mat[i, j] <- r$estimate
    p_mat[i, j]   <- r$p.value
  }
}
cor_df <- as.data.frame(cor_mat)
p_df   <- as.data.frame(p_mat)
write.csv(cor_df, "hub_immune_spearman_corr.csv", row.names = TRUE)
write.csv(p_df,   "hub_immune_spearman_p.csv",    row.names = TRUE)
cat("已写出 hub_immune_spearman_corr.csv / _p.csv\n")

## 6b. 肿瘤 vs 正常: hub 表达与免疫评分的差异 (Wilcoxon)
tw <- data.frame(
  Feature = c(hub_cols, imm_cols),
  Mean_Tumor   = NA, Mean_Normal = NA,
  log2FC       = NA, P_value      = NA, stringsAsFactors = FALSE
)
for (k in seq_along(tw$Feature)) {
  f <- tw$Feature[k]
  a <- dat[[f]][dat$Group == "Tumor"]
  b <- dat[[f]][dat$Group == "Normal"]
  if (length(b) > 1) {
    tw$Mean_Tumor[k]  <- mean(a, na.rm = TRUE)
    tw$Mean_Normal[k] <- mean(b, na.rm = TRUE)
    tw$log2FC[k]      <- log2(mean(a, na.rm = TRUE) / mean(b, na.rm = TRUE))
    tw$P_value[k]     <- wilcox.test(a, b)$p.value
  }
}
tw$FDR <- p.adjust(tw$P_value, method = "fdr")
write.csv(tw, "hub_and_immune_tumor_vs_normal.csv", row.names = FALSE)
cat("已写出 hub_and_immune_tumor_vs_normal.csv\n")
cat("肿瘤vs正常 显著差异(FDR<0.05)的特征数:", sum(tw$FDR < 0.05, na.rm = TRUE), "\n")

## ---- 7. 可视化: hub × 免疫指标 相关热图 ----
if (require("pheatmap", quietly = TRUE)) {
  pdf("hub_immune_corr_heatmap.pdf", width = 10, height = 7)
  pheatmap::pheatmap(cor_mat, main = "Spearman corr: 7 hub vs immune metrics (Tumor)",
                     cluster_rows = TRUE, cluster_cols = TRUE,
                     display_numbers = TRUE, number_format = "%.2f",
                     fontsize = 9, na_col = "grey90")
  dev.off()
  cat("已出图 hub_immune_corr_heatmap.pdf\n")
} else {
  warning("pheatmap 未安装 -> 跳过热图(相关矩阵已在 CSV 中)")
}

cat("\n=== 完成 ===\n")
cat("产物: immune_infiltration_scores.csv / hub_immune_spearman_corr.csv(_p)\n")
cat("      hub_and_immune_tumor_vs_normal.csv / hub_immune_corr_heatmap.pdf\n")
cat("把这几个 CSV + 控制台末几行贴给我, 我帮你判读 7 hub 与哪些免疫细胞显著关联。\n")
