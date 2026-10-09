# ============================================================
# plot_volcano_standalone.R
#   独立火山图重绘：仅读取已算出的 sigDiff.txt，不重跑 scTenifoldKnk
#   用法：修改下方 target_gene，直接 source / 在 RStudio 点 Source 即可
#   说明：sigDiff.txt 只含显著差异基因(p.adj<0.05)，故底图仅画显著红点
#         （不含非显著灰云，比原始 vol.pdf 更干净；如需灰云背景需重算并 saveRDS(df)）
# ============================================================

# ---------- 参数 ----------
target_gene <- "GPAT3"   # 改成 "PLAU" 即可复用（需该基因目录下已有 sigDiff.txt）

# ---------- 路径 ----------
base_dir <- "D:/BaiduSyncdisk/daixie-colorectal/4.免疫浸润单细胞/28.scTenifoldKnk"
out_dir  <- file.path(base_dir, target_gene)
setwd(out_dir)

# ---------- 读数据 ----------
sig_file <- "sigDiff.txt"
if (!file.exists(sig_file)) {
  stop(paste0(">>> 未找到 ", out_dir, "/", sig_file, " —— 请先确认 ", target_gene, " 的 scTenifoldKnk 已跑完"))
}
df <- read.delim(sig_file, stringsAsFactors = FALSE)

# 列检查（兼容列名差异）
required <- c("gene", "Z", "p.adj")
if (!all(required %in% names(df))) {
  stop(paste0(">>> sigDiff.txt 缺少必要列，现有列: ", paste(names(df), collapse = ", ")))
}

# ---------- 构造火山图数据 ----------
df$log_p.adj  <- -log10(df$p.adj)
df$significant <- "Significant"   # sigDiff.txt 内全部为显著基因

# 仅标注效应最强的 Top 30（按 |Z| 排序），其余显著基因保留红点但不标名
label_genes <- df[order(-abs(df$Z)), , drop = FALSE][1:30, , drop = FALSE]

# ---------- 绘图 ----------
library(ggplot2)
library(ggrepel)

p <- ggplot(df, aes(x = Z, y = log_p.adj)) +
  geom_point(color = "red", alpha = 0.7, size = 1.5) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
  geom_text_repel(data = label_genes, aes(label = gene), size = 3, max.overlaps = 30) +
  labs(title = paste0(target_gene, " KO (virtual)"),
       x = "Z-score",
       y = expression(-log[10](p.adj))) +
  theme_classic()

# 输出（覆盖原 vol.pdf，dpi 300 出版级）
ggsave("vol.pdf", p, width = 6, height = 5, dpi = 300)
cat(sprintf(">>> %s 火山图已重绘并保存: %s/vol.pdf （显著基因数 = %d, 标注 Top 30）\n",
            target_gene, out_dir, nrow(df)))
