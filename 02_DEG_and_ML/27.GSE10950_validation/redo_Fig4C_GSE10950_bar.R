# ============================================================
# Fig4C — GSE10950 independent hub differential-expression validation
# R (ggplot2) re-implementation of the original Python bar chart.
#
# 输入 : GSE10950_hub_validation.csv  (Python 已算好: 7 hub 的
#         log2FC / Welch FDR / 方向一致性)
# 输出 : Fig4C_GSE10950_validation.tif  (投稿 TIFF, 300 dpi, LZW)
#        Fig4C_GSE10950_validation.pdf  (cairo, Illustrator 兼容)
#
# 依赖 : ggplot2
# 运行 : 在 RStudio 直接 Source 本文件即可 (Windows 路径已用正斜杠)
# ============================================================
setwd("D:/BaiduSyncdisk/daixie-colorectal/第二部分DEG/27.GSE10950_validation")

library(ggplot2)

## 1. 读取 Python 写出的验证结果 -------------------------------------------
df <- read.csv("GSE10950_hub_validation.csv", stringsAsFactors = FALSE)

## 原 CSV 表头含 "<"，read.csv 会规范成 "significant_FDR.0.05"，
## 这里用正则稳健抓取显著性列，避免列名写死。
sig.col <- grep("^significant_FDR", names(df), value = TRUE)[1]
df$sig        <- toupper(trimws(as.character(df[[sig.col]]))) == "TRUE"
df$consistent <- toupper(trimws(as.character(df$direction_consistent_with_main))) == "TRUE"

## 2. 分类上色 (与 Python 版配色一致) -------------------------------------
##  显著且方向一致 -> 绿 ; 显著但方向反转 -> 红 ; 不显著 -> 灰
df$status <- ifelse(df$sig & df$consistent, "Significant & concordant",
             ifelse(df$sig,                  "Significant but discordant",
                                              "Not significant"))
df$status <- factor(df$status,
                    levels = c("Significant & concordant",
                               "Significant but discordant",
                               "Not significant"))
## 保持 ML 锁定顺序 (与发现集/主文一致)
df$Gene  <- factor(df$Gene, levels = df$Gene)

## 柱顶标签: log2FC 值, 显著加 "*"
df$lab   <- sprintf("%+.2f%s", df$log2FC_tumor_vs_normal,
                     ifelse(df$sig, "*", ""))
df$vjust <- ifelse(df$log2FC_tumor_vs_normal >= 0, -0.5, 1.5)

## 3. 绘图 ----------------------------------------------------------------
p <- ggplot(df, aes(x = Gene, y = log2FC_tumor_vs_normal, fill = status)) +
  geom_bar(stat = "identity", color = "black", linewidth = 0.6, width = 0.7) +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.8) +
  geom_text(aes(label = lab, vjust = vjust), size = 3.2) +
  scale_fill_manual(
    values = c("Significant & concordant"   = "#2ca02c",
               "Significant but discordant" = "#d62728",
               "Not significant"            = "#7f7f7f")) +
  labs(x = NULL,
       y = expression(log[2]~FC~(tumor~vs~normal)),
       fill = NULL,
       title = "GSE10950 independent validation:\n7 ML-identified hub genes in CRC") +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "top",
    legend.direction = "horizontal",
    plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
    axis.text.x = element_text(angle = 30, hjust = 1, size = 10),
    axis.title.y = element_text(size = 11)
  )

## 4. 导出 ----------------------------------------------------------------
## TIFF (投稿用, 300 dpi, LZW 压缩)
tiff("Fig4C_GSE10950_validation.tif", width = 8, height = 5,
     units = "in", res = 300, type = "cairo", compression = "lzw")
print(p); dev.off()

## PDF (cairo, 供 Illustrator 组版)
cairo_pdf("Fig4C_GSE10950_validation.pdf", width = 8, height = 5)
print(p); dev.off()

cat("Fig4C written -> Fig4C_GSE10950_validation.tif / .pdf\n")
cat("显著且方向一致:", paste(df$Gene[df$sig & df$consistent], collapse = ", "), "\n")
