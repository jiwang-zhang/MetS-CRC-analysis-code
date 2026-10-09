# ============================================================
# TCGA 多基因 Cox 回归森���图 (v3 — base R)
# 输入: GEPIA2 data.csv
# 输出: forest_plot.pdf + TCGA_multicox_summary.csv
# ============================================================

setwd("D:/BaiduSyncdisk/daixie-colorectal/GEPIA2")

# ============================================================
# 1. 读取 & 映射
# ============================================================
raw <- read.csv("data.csv", stringsAsFactors = FALSE)

ensg_map <- c("ENSG00000138678" = "GPAT3",
              "ENSG00000176454" = "LPCAT4",
              "ENSG00000122861" = "PLAU",
              "ENSG00000112033" = "PPARD",
              "ENSG00000118515" = "SGK1",
              "ENSG00000057704" = "TMCC3",
              "ENSG00000146072" = "TNFRSF21")
raw$Gene <- ifelse(raw$qname %in% names(ensg_map),
                   ensg_map[raw$qname], raw$qname)

# 按 HR 从大到小 (高风险在上)
raw <- raw[order(raw$HR, decreasing = TRUE), ]

cat("=== TCGA 多基因 Cox 回归结果 ===\n")
print(raw[, c("Gene", "HR", "CI95.lower", "CI95.upper", "p.value")])

# ============================================================
# 2. 格式化
# ============================================================
n <- nrow(raw)
genes <- raw$Gene
hr    <- raw$HR
lo    <- raw$CI95.lower
hi    <- raw$CI95.upper
pv    <- raw$p.value

fmt_p <- function(p) {
  if (p < 0.001) return("< 0.001")
  if (p < 0.01)  return(sprintf("%.3f", p))
  return(sprintf("%.2f", p))
}
p_str  <- sapply(pv, fmt_p)
sig    <- ifelse(pv < 0.001, "***",
          ifelse(pv < 0.01, "**",
          ifelse(pv < 0.05, "*", " ")))

hr_ci  <- sprintf("%.2f (%.2f - %.2f)", hr, lo, hi)
p_col  <- paste0(p_str, sig)

# 颜色
point_col <- ifelse(pv < 0.05 & hr < 1, "#2166AC",
             ifelse(pv < 0.05 & hr > 1, "#B2182B",
                    "#AAAAAA"))

# ============================================================
# 3. 画图 (base R)
# ============================================================
cairo_pdf("forest_plot.pdf", width = 11, height = 5.5)

# 布局: 左侧文字 + 右侧森林图
layout(matrix(c(1, 1, 2), nrow = 1), widths = c(3.5, 4, 3.5))
par(mar = c(0, 0, 0, 0))

# ---- 左面板: 文字 (Gene / HR / P value) ----
plot.new()
plot.window(xlim = c(0, 1), ylim = c(0.5, n + 1.5))
# 表头
text(0.02, n + 1.0, "Gene",      pos = 4, cex = 1.6, font = 2)
text(0.40, n + 1.0, "HR (95% CI)", pos = 4, cex = 1.6, font = 2)
text(0.78, n + 1.0, "P value",    pos = 4, cex = 1.6, font = 2)
# 横线
segments(0, n + 0.7, 1, n + 0.7, col = "grey60")
# 每行
for (i in 1:n) {
  y <- n - i + 1
  text(0.02, y, genes[i], pos = 4, cex = 1.5, font = 2,
       col = point_col[i])
  text(0.40, y, hr_ci[i], pos = 4, cex = 1.3)
  text(0.78, y, p_col[i], pos = 4, cex = 1.3,
       col = ifelse(pv[i] < 0.05, "#B2182B", "#666666"))
}

# ---- 右面板: 森林图 ----
par(mar = c(4, 2, 0, 2))
plot.new()
xr <- c(0.1, 3.0)
plot.window(xlim = xr, ylim = c(0.5, n + 0.5))
# 参考线
abline(v = 1, lty = 2, col = "grey60", lwd = 1.5)
# X 轴
axis(1, at = c(0.25, 0.5, 1, 1.5, 2, 2.5, 3),
     labels = c("0.25", "0.5", "1.0", "1.5", "2.0", "2.5", "3.0"),
     cex.axis = 1.3)
title(xlab = "Hazard Ratio", cex.lab = 1.5)
# CI 线和点
for (i in 1:n) {
  y <- n - i + 1
  # CI 线
  segments(lo[i], y, hi[i], y, col = point_col[i], lwd = 3)
  # 点 (HR)
  points(hr[i], y, pch = 19, cex = 2.2, col = point_col[i])
  points(hr[i], y, pch = 1,  cex = 2.2, col = "black", lwd = 0.5)
}
# 图例
legend("topright",
       legend = c("Protective (P<0.05)", "Risk (P<0.05)", "Non-significant"),
       pch = 19, col = c("#2166AC", "#B2182B", "#AAAAAA"),
       pt.cex = 2.0, cex = 1.2, bty = "n")

dev.off()
cat("\nforest_plot.pdf saved\n")

# ============================================================
# 4. CSV 汇总
# ============================================================
out <- raw[, c("Gene", "coef", "se.coef.", "HR", "CI95.lower", "CI95.upper", "z", "p.value")]
write.csv(out, "TCGA_multicox_summary.csv", row.names = FALSE)
cat("TCGA_multicox_summary.csv saved\n")
