# ============================================================
# Step 7b: CIBERSORT 结果可视化
# 输入: ../23.CIBERSORT/CIBERSORT-Results.txt
# 输出: barplot.pdf (堆叠柱状图) + immune.diff.pdf (组间箱线图)
#
# 依赖: reshape2, ggpubr
# ============================================================

library(reshape2)
library(ggpubr)

inputFile <- "../23.CIBERSORT/CIBERSORT-Results.txt"   # CIBERSORT 结果
diseaseName <- "CRC"                                     # 疾病名称
setwd("D:/BaiduSyncdisk/daixie-colorectal/4.免疫浸润单细胞/24.barplot")

# 读取免疫细胞浸润文件
rt <- read.table(inputFile, header = TRUE, sep = "\t", check.names = FALSE, row.names = 1)

# 分组: 列名含 _tre 为肿瘤, 否则为正常 (列名=样本ID)
con  <- !grepl("_tre$", rownames(rt))
treat <- grepl("_tre$", rownames(rt))
conData  <- rt[con, , drop = FALSE]
treatData <- rt[treat, , drop = FALSE]
conNum   <- nrow(conData)
treatNum <- nrow(treatData)
cat("有效样本: con =", conNum, ", tre =", treatNum, "\n")

data <- t(rbind(conData, treatData))

# ============================================================
# 堆叠柱状图
# ============================================================
pdf(file = "barplot.pdf", width = 13, height = 7)
col <- rainbow(nrow(data), s = 0.7, v = 0.7)
par(las = 1, mar = c(8, 5, 4, 16), mgp = c(3, 0.1, 0), cex.axis = 1.5)
a1 <- barplot(data, col = col, xaxt = "n", yaxt = "n",
              ylab = "Relative Percent", cex.lab = 1.5)
a2 <- axis(2, tick = FALSE, labels = FALSE)
axis(2, a2, paste0(a2 * 100, "%"))
par(srt = 0, xpd = TRUE)
# 蓝色 = 正常, 红色 = 肿瘤
if (conNum > 0) {
  rect(xleft = a1[1] - 0.5, ybottom = -0.01,
       xright = a1[conNum] + 0.5, ytop = -0.05, col = "#6699FFFF")
  text(a1[conNum] / 2, -0.025, "Normal", cex = 1.8)
}
if (treatNum > 0) {
  rect(xleft = a1[conNum + 1] - 0.5, ybottom = -0.01,
       xright = a1[length(a1)] + 0.5, ytop = -0.05, col = "#E6550DFF")
  text((a1[length(a1)] + a1[conNum + 1]) / 2, -0.025, diseaseName, cex = 1.8)
}
legend(par('usr')[2] * 0.98, par('usr')[4], legend = rownames(data),
       col = col, pch = 15, bty = "n", cex = 1)
dev.off()
cat("barplot.pdf 已保存\n")

# ============================================================
# 组间箱线图
# ============================================================
plotDf <- as.data.frame(t(data))
plotDf$Type <- c(rep("Normal", conNum), rep(diseaseName, treatNum))
plotDf$Type <- factor(plotDf$Type, levels = c("Normal", diseaseName))
plotLong <- reshape2::melt(plotDf, id.vars = "Type",
                           variable.name = "Immune", value.name = "Fraction")

p_box <- ggboxplot(plotLong, x = "Immune", y = "Fraction", fill = "Type",
                   xlab = "", ylab = "Fraction",
                   legend.title = "Type", width = 0.8,
                   palette = c("#6699FFFF", "#E6550DFF")) +
  rotate_x_text(50) +
  stat_compare_means(aes(group = Type),
                     symnum.args = list(cutpoints = c(0, 0.001, 0.01, 0.05, 1),
                                        symbols = c("***", "**", "*", "")),
                     label = "p.signif")

pdf(file = "immune.diff.pdf", width = 8, height = 6)
print(p_box)
dev.off()
cat("immune.diff.pdf 已保存\n")