# =====================================================================
# 10.GO_run.R  ——  适配版 GO 富集分析（200 co-DEGs）
# 适配自课程脚本 AI17.GO.R，仅改路径与输入文件，逻辑保持不变。
# 工作目录已指向百度同步盘副本（跨电脑前请确认本机路径一致，否则改第 3 节）。
# 运行方式：在 R / RStudio 里 source 本文件，或打开后 Ctrl+Alt+R 全跑。
# 依赖包：clusterProfiler, org.Hs.eg.db, enrichplot, ggplot2,
#         circlize, RColorBrewer, dplyr, ggpubr, ComplexHeatmap
# （首次使用若缺包，取消下面 install 区块注释先装）
# =====================================================================

# ---- 0. 安装依赖（仅首次，按需取消注释）-------------------------------
# install.packages(c("colorspace","stringi","ggplot2","circlize",
#                    "RColorBrewer","ggpubr","dplyr"))
# if (!requireNamespace("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
# BiocManager::install(c("org.Hs.eg.db","DOSE","clusterProfiler",
#                        "enrichplot","ComplexHeatmap","ReactomePA"))

# ---- 1. 引用包 -------------------------------------------------------
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)
library(circlize)
library(RColorBrewer)
library(dplyr)
library("ggpubr")
library(ComplexHeatmap)

# ---- 2. 参数 ---------------------------------------------------------
pvalueFilter  = 0.05     # p 值过滤条件
p.adjustFilter= 1        # 矫正后 p 值过滤条件（=1 表示不过滤，靠 pvalue 控）

colorSel = "p.adjust"
if (p.adjustFilter > 0.05) {
  colorSel = "pvalue"
}

# ---- 3. 工作目录（★ 同步盘副本路径；换电脑后若盘符/同步名不同改这里）----
setwd("D:/BaiduSyncdisk/daixie-colorectal/第二部分DEG/10.GO")

# ---- 4. 读入 200 个 co-DEGs 的基因名 ---------------------------------
# 原课程脚本读 interGenes.txt（单列、无表头、tab 分隔）；
# 本项目用 co_DEGs_CRC_MetS.csv，取 Gene 列即可。
rt   = read.csv("co_DEGs_CRC_MetS.csv", header = TRUE, stringsAsFactors = FALSE)
genes = unique(as.vector(rt$Gene))
cat(sprintf("读入基因数（去重前）：%d，去重后：%d\n", nrow(rt), length(genes)))

# 基因符号 -> Entrez ID
entrezIDs = mget(genes, org.Hs.egSYMBOL2EG, ifnotfound = NA)
entrezIDs = as.character(entrezIDs)
gene = entrezIDs[entrezIDs != "NA"]      # 去除映射不到 Entrez 的基因
cat(sprintf("成功映射到 Entrez 的基因数：%d\n", length(gene)))

# ---- 5. GO 富集分析（BP/CC/MF 全做）----------------------------------
kk  = enrichGO(gene = gene, OrgDb = org.Hs.eg.db,
               pvalueCutoff = 1, qvalueCutoff = 1,
               ont = "all", readable = TRUE)
GO  = as.data.frame(kk)
GO  = GO[(GO$pvalue < pvalueFilter & GO$p.adjust < p.adjustFilter), ]
# 输出显著富集结果
write.table(GO, file = "GO.txt", sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("显著 GO 条目数（p<%.3f）：%d\n", pvalueFilter, nrow(GO)))

# 柱状图
pdf(file = "barplot.pdf", width = 11, height = 7)
bar = barplot(kk, drop = TRUE, showCategory = 10, label_format = 100,
              split = "ONTOLOGY", color = colorSel) + facet_grid(ONTOLOGY ~ ., scale = 'free')
print(bar)
dev.off()

# 气泡图
pdf(file = "bubble.pdf", width = 11, height = 7)
bub = dotplot(kk, showCategory = 10, orderBy = "GeneRatio",
              label_format = 100, split = "ONTOLOGY", color = colorSel) +
  facet_grid(ONTOLOGY ~ ., scale = 'free')
print(bub)
dev.off()

# ---- 6. GO 圈图（空结果保护）-----------------------------------------
ontology.col = c("#3498db", "#2ecc71", "#e67e22")
data    = GO[order(GO$p.adjust), ]
datasig = data[data$pvalue < 0.05, , drop = FALSE]

if (nrow(datasig) >= 1) {
  BP = head(datasig[datasig$ONTOLOGY == "BP", , drop = FALSE], 6)
  CC = head(datasig[datasig$ONTOLOGY == "CC", , drop = FALSE], 6)
  MF = head(datasig[datasig$ONTOLOGY == "MF", , drop = FALSE], 6)
  data = rbind(BP, CC, MF)
  main.col = ontology.col[as.numeric(as.factor(data$ONTOLOGY))]

  # 整理圈图数据
  BgGene = as.numeric(sapply(strsplit(data$BgRatio,  "/"), '[', 1))
  Gene   = as.numeric(sapply(strsplit(data$GeneRatio, '/'), '[', 1))
  ratio  = Gene / BgGene
  logpvalue   = -log(data$pvalue, 10)
  logpvalue.col = brewer.pal(n = 8, name = "Reds")
  f = colorRamp2(breaks = c(0, 2, 4, 6, 8, 10, 15, 20), colors = logpvalue.col)
  BgGene.col = f(logpvalue)
  df   = data.frame(GO = data$ID, start = 1, end = max(BgGene))
  rownames(df) = df$GO
  bed2 = data.frame(GO = data$ID, start = 1, end = BgGene, BgGene = BgGene, BgGene.col = BgGene.col)
  bed3 = data.frame(GO = data$ID, start = 1, end = Gene,   BgGene = Gene)
  bed4 = data.frame(GO = data$ID, start = 1, end = max(BgGene), ratio = ratio, col = main.col)
  bed4$ratio = bed4$ratio / max(bed4$ratio) * 9.5

  # 绘制圈图主体
  pdf(file = "GO.circlize.pdf", width = 10, height = 10)
  par(omi = c(0.1, 0.1, 0.1, 1.5))
  circos.par(track.margin = c(0.01, 0.01))
  circos.genomicInitialize(df, plotType = "none")
  circos.trackPlotRegion(ylim = c(0, 1), panel.fun = function(x, y) {
    sector.index = get.cell.meta.data("sector.index")
    xlim = get.cell.meta.data("xlim")
    ylim = get.cell.meta.data("ylim")
    circos.text(mean(xlim), mean(ylim), sector.index, cex = 0.8,
                facing = "bending.inside", niceFacing = TRUE)
  }, track.height = 0.08, bg.border = NA, bg.col = main.col)

  for (si in get.all.sector.index()) {
    circos.axis(h = "top", labels.cex = 0.6, sector.index = si, track.index = 1,
                major.at = seq(0, max(BgGene), by = 100), labels.facing = "clockwise")
  }
  f = colorRamp2(breaks = c(-1, 0, 1), colors = c("green", "black", "red"))
  circos.genomicTrack(bed2, ylim = c(0, 1), track.height = 0.1, bg.border = "white",
                      panel.fun = function(region, value, ...) {
                        i = getI(...)
                        circos.genomicRect(region, value, ytop = 0, ybottom = 1,
                                           col = value[, 2], border = NA, ...)
                        circos.genomicText(region, value, y = 0.4, labels = value[, 1],
                                           adj = 0, cex = 0.8, ...)
                      })
  circos.genomicTrack(bed3, ylim = c(0, 1), track.height = 0.1, bg.border = "white",
                      panel.fun = function(region, value, ...) {
                        i = getI(...)
                        circos.genomicRect(region, value, ytop = 0, ybottom = 1,
                                           col = '#BA55D3', border = NA, ...)
                        circos.genomicText(region, value, y = 0.4, labels = value[, 1],
                                           cex = 0.9, adj = 0, ...)
                      })
  circos.genomicTrack(bed4, ylim = c(0, 10), track.height = 0.35, bg.border = "white",
                      bg.col = "grey90",
                      panel.fun = function(region, value, ...) {
                        cell.xlim = get.cell.meta.data("cell.xlim")
                        cell.ylim = get.cell.meta.data("cell.ylim")
                        for (j in 1:9) {
                          y = cell.ylim[1] + (cell.ylim[2] - cell.ylim[1]) / 10 * j
                          circos.lines(cell.xlim, c(y, y), col = "#FFFFFF", lwd = 0.3)
                        }
                        circos.genomicRect(region, value, ytop = 0, ybottom = value[, 1],
                                           col = value[, 2], border = NA, ...)
                      })
  circos.clear()
  # 中间图例
  middle.legend = Legend(
    labels = c('Number of Genes', 'Number of Select', 'Rich Factor(0-1)'),
    type = "points", pch = c(15, 15, 17),
    legend_gp = gpar(col = c('pink', '#BA55D3', ontology.col[1])),
    title = "", nrow = 3, size = unit(3, "mm"))
  circle_size = unit(1, "snpc")
  draw(middle.legend, x = circle_size * 0.42)
  # ONTOLOGY 图例
  main.legend = Legend(
    labels = c("Biological Process", "Cellular Component", "Molecular Function"),
    type = "points", pch = 15, legend_gp = gpar(col = ontology.col),
    title_position = "topcenter", title = "ONTOLOGY", nrow = 3,
    size = unit(3, "mm"), grid_height = unit(5, "mm"), grid_width = unit(5, "mm"))
  # -log10(P) 图例
  logp.legend = Legend(
    labels = c('(0,2]', '(2,4]', '(4,6]', '(6,8]', '(8,10]', '(10,15]', '(15,20]', '>=20'),
    type = "points", pch = 16, legend_gp = gpar(col = logpvalue.col),
    title = "-log10(Pvalue)", title_position = "topcenter",
    grid_height = unit(5, "mm"), grid_width = unit(5, "mm"), size = unit(3, "mm"))
  lgd = packLegend(main.legend, logp.legend)
  circle_size = unit(1, "snpc")
  draw(lgd, x = circle_size * 0.85, y = circle_size * 0.55, just = "left")
  dev.off()
  cat("圈图已输出：GO.circlize.pdf\n")
} else {
  warning("无显著 GO 条目（p<0.05），跳过圈图绘制。")
}

cat("GO 富集完成。\n")
