# ============================================================
# 从 HALLMARK_TNFA_SIGNALING_VIA_NFKB 创始集 GMT 中
# 提取出现 ≥2 次的高频基因  (R 版)
# 直接 Ctrl+A → Ctrl+Enter 运行
# ============================================================
setwd("D:/daixie-colorectal/第二部分DEG/10.获取研究的基因集")

gmt_file <- "HALLMARK_TNFA_SIGNALING_VIA_NFKB_FOUNDERS.v2026.1.Hs.gmt"
min_freq <- 2   # 可改：2/3/5/10

# 读 GMT，统计每个基因的创始集出现次数
lines   <- readLines(gmt_file)
cat("创始集总数:", length(lines), "\n")

all_genes <- unlist(lapply(lines, function(line) {
  parts <- strsplit(line, "\t")[[1]]
  unique(parts[-c(1,2)])  # 跳过基因集名和URL
}))
gene_freq <- sort(table(all_genes), decreasing = TRUE)

cat("去重基因总数:", length(gene_freq), "\n")

# 筛选 ≥ min_freq
selected <- gene_freq[gene_freq >= min_freq]
cat("出现 ≥", min_freq, "次的基因:", length(selected), "\n")

# 保存 TXT（每行一个基因）
writeLines(names(selected), paste0("TNF_NFKB_freq", min_freq, ".txt"))
cat("已保存: TNF_NFKB_freq", min_freq, ".txt\n", sep = "")

# 保存带频次的 CSV
df <- data.frame(Gene = names(selected), FounderSetCount = as.integer(selected))
write.csv(df, paste0("TNF_NFKB_freq", min_freq, "_with_counts.csv"), row.names = FALSE)
cat("已保存: TNF_NFKB_freq", min_freq, "_with_counts.csv\n", sep = "")

# 频次分布
cat("\n频次分布 (Top 10):\n")
freq_dist <- table(as.integer(gene_freq))
top <- head(sort(as.integer(names(freq_dist)), decreasing = TRUE), 10)
for (f in top) {
  cat(sprintf("  出现 %d 次: %d 个基因\n", f, freq_dist[as.character(f)]))
}
