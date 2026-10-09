# GSE145412 差异分析 —— limma (MetS vs Healthy, 敏感性)
# 依赖: Bioconductor limma
# 数据: 19.databases/GSE145412/GSE145412_rnaseq.txt.gz
#       + 19.databases/GSE145412/GSE145412_series_matrix.txt.gz (解析 sample 分组)
# 注意: 该矩阵是归一化表达量(类TPM), 非原始计数, 故用 limma 而非 DESeq2
# 用法: 在 19.databases/ 目录下运行 Rscript deg_gse145412_limma.R

library(limma)

f  <- "GSE145412/GSE145412_rnaseq.txt.gz"
df <- read.delim(f, check.names = FALSE)

# --- 从 series_matrix 的 !Sample_title 解析 sample_X -> (健康状态, BMI) ---
sm <- read.delim("GSE145412/GSE145412_series_matrix.txt.gz", sep = "\t",
                 header = FALSE, comment.char = "", stringsAsFactors = FALSE)
title_row <- sm[sm$V1 == "!Sample_title", ]
titles <- as.character(title_row[1, -1])

sample_cols <- grep("^sample_", colnames(df), value = TRUE)
health <- bmi <- rep(NA, length(sample_cols))
for (i in seq_along(sample_cols)) {
  s   <- sub("sample_", "", sample_cols[i])
  ttx <- titles[grep(paste0("sample_", s, "\\b"), titles)]
  if (length(ttx) > 0) {
    tok <- trimws(strsplit(gsub('"', '', ttx), ",")[[1]])
    health[i] <- if ("D" %in% tok) "MetS" else if ("H" %in% tok) "Healthy" else NA
    bmi[i]    <- if ("O" %in% tok) "Obese" else if ("L" %in% tok) "Lean"   else NA
  }
}
cat("各组样本数:\n"); print(table(health, bmi))

expr  <- as.matrix(df[, sample_cols])
expr  <- log2(expr + 1)                              # 归一化数据 -> log2

# 敏感性主对比: MetS(疾病 D) vs Healthy(健康 H), 跨 BMI 合并
group <- factor(ifelse(health == "MetS", "MetS", "Healthy"),
                levels = c("Healthy", "MetS"))
design <- model.matrix(~ group)
fit <- lmFit(expr, design)
fit <- eBayes(fit)
res <- topTable(fit, coef = 2, number = Inf, adjust = "BH")
res$Symbol <- df$External_gene_name
res$geneID <- df$geneID

sig <- res[res$adj.P.Val < 0.05 & abs(res$log2FC) > 0.5, ]
message(sprintf("MetS=%d, Healthy=%d, 基因=%d", sum(group=="MetS"),
                sum(group=="Healthy"), nrow(expr)))
message(sprintf("DEG(|log2FC|>0.5 & adj.P.Val<0.05): %d (上调 %d, 下调 %d)",
                nrow(sig), sum(sig$log2FC>0), sum(sig$log2FC<0)))

write.csv(res, "GSE145412/DEG_MetS_vs_Healthy_R.csv", row.names = TRUE)
write.csv(sig, "GSE145412/DEG_MetS_vs_Healthy_R.sig.csv", row.names = TRUE)
