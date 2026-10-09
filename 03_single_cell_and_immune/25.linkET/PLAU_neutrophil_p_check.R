# =============================================================
# PLAU–Neutrophils Spearman p 值算法复核（可选，跑完贴回 AI 核对）
# 目的：验证 linkET 输出的 p = 1.63e-24 是否可复现
# 原理：OE25.linkET.R 只用肿瘤样本（_tre 列，n = 687），
#       cor.test(method="spearman") 在有 ties 时走 t 近似
# 预期：r ≈ 0.376, p ≈ 1.2~1.6e-24（与正文一致）
# =============================================================

library(limma)

setwd("D:/BaiduSyncdisk/daixie-colorectal/4.免疫浸润单细胞/25.linkET")

# 1. 表达矩阵（与 OE25.linkET.R 同源）
rt <- read.table("../../第二部分DEG/25.diagnostic_model/merged_after_batch_removal.csv",
                 header = TRUE, sep = ",", check.names = FALSE)
rt <- as.matrix(rt); rownames(rt) <- rt[, 1]
exp <- rt[, 2:ncol(rt)]
dn  <- list(rownames(exp), colnames(exp))
exp <- matrix(as.numeric(exp), nrow = nrow(exp), dimnames = dn)
exp <- avereps(exp)

# 2. 只取肿瘤样本
data <- t(exp["PLAU", grepl("_tre$", colnames(exp)), drop = FALSE])

# 3. CIBERSORT 结果并对齐
immune <- read.table("../23.CIBERSORT/CIBERSORT-Results.txt",
                     header = TRUE, sep = "\t", check.names = FALSE, row.names = 1)
sameSample <- intersect(rownames(data), rownames(immune))
x <- as.numeric(immune[sameSample, "Neutrophils"])
y <- as.numeric(data[sameSample, "PLAU"])
cat("对齐后肿瘤样本数:", length(sameSample), "\n")

# 4. 三种算法对比
ct1 <- cor.test(x, y, method = "spearman", exact = FALSE)  # t 近似（linkET 同款）
cat(sprintf("cor.test exact=FALSE : r = %.4f, p = %.3e\n", ct1$estimate, ct1$p.value))

ct2 <- suppressWarnings(cor.test(x, y, method = "spearman"))  # 默认（有 ties 自动近似）
cat(sprintf("cor.test 默认        : p = %.3e\n", ct2$p.value))

if (requireNamespace("psych", quietly = TRUE)) {
  pt <- psych::corr.test(data.frame(PLAU = y, Neutrophils = x), method = "spearman")
  cat(sprintf("psych::corr.test     : r = %.4f, p = %.3e\n",
              pt$r["PLAU","Neutrophils"], pt$p["PLAU","Neutrophils"]))
}

cat("\n正文报告值: r = 0.376, p = 1.63e-24（BH = 2.51e-22 = 原始p × 154）\n")
