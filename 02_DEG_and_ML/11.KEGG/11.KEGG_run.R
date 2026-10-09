###############################################################
# 11.KEGG 富集分析 (200 co-DEGs)
# 适配自课程版 AI18.KEGG.R
# 输入: co_DEGs_CRC_MetS.csv (首列 Gene, 200 个 CRC∩MetS co-DEGs)
# 输出: KEGG.txt + barplot.pdf + bubble.pdf
###############################################################

## 1. 引用包 (首次运行前请先装好, 安装命令见文件末尾注释) --------------
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)
library(circlize)
library(RColorBrewer)
library(dplyr)
library(ggpubr)

## 2. 过滤参数 ---------------------------------------------------------
pvalueFilter   = 0.05     # p 值过滤条件
p.adjustFilter = 1        # 矫正后 p 值过滤条件 (1=不按 adjP 卡, 仅看 pvalue)

# 图形着色列: adjP 阈值宽松时按 pvalue 上色, 否则按 p.adjust
colorSel = "p.adjust"
if(p.adjustFilter > 0.05){
  colorSel = "pvalue"
}

## 3. 设置工作目录 (换电脑/换盘符只改这一行) --------------------------
setwd("D:/BaiduSyncdisk/daixie-colorectal/第二部分DEG/11.KEGG")

## 4. 读入 co-DEGs 并转 Entrez ID ------------------------------------
rt = read.csv("co_DEGs_CRC_MetS.csv", header=T, check.names=F)   # 首列 Gene
genes = unique(as.vector(rt[,1]))
entrezIDs = mget(genes, org.Hs.egSYMBOL2EG, ifnotfound=NA)
entrezIDs = as.character(entrezIDs)
rt = data.frame(genes, entrezID=entrezIDs)
gene = entrezIDs[entrezIDs != "NA"]          # 去除无法映射到 Entrez 的基因

## 5. KEGG 富集分析 --------------------------------------------------
kk <- enrichKEGG(gene=gene, organism="hsa", pvalueCutoff=1, qvalueCutoff=1)
KEGG = as.data.frame(kk)
# 把结果里的 Entrez ID 还原成基因 Symbol, 便于阅读
KEGG$geneID = as.character(sapply(KEGG$geneID, function(x)
  paste(rt$genes[match(strsplit(x,"/")[[1]], as.character(rt$entrezID))], collapse="/")))
KEGG = na.omit(KEGG)
KEGG = KEGG[(KEGG$pvalue < pvalueFilter & KEGG$p.adjust < p.adjustFilter),]

# 保存显著富集结果
write.table(KEGG, file="KEGG.txt", sep="\t", quote=F, row.names=F)

## 6. 空结果保护 ----------------------------------------------------
if(nrow(KEGG) == 0){
  cat("警告: 没有满足 pvalue<", pvalueFilter, " & p.adjust<", p.adjustFilter,
      " 的显著 KEGG 通路, 跳过绘图。\n")
} else {
  # 展示通路数目 (最多前 30)
  showNum = 30
  if(nrow(KEGG) < showNum){
    showNum = nrow(KEGG)
  }

  # 柱状图
  pdf(file="barplot.pdf", width=8, height=7)
  print(barplot(kk, drop=TRUE, showCategory=showNum, label_format=130, color=colorSel))
  dev.off()

  # 气泡图
  pdf(file="bubble.pdf", width=8, height=7)
  print(dotplot(kk, showCategory=showNum, orderBy="GeneRatio", label_format=130, color=colorSel))
  dev.off()

  cat("KEGG 完成: 共", nrow(KEGG), "条显著通路, 已输出 KEGG.txt / barplot.pdf / bubble.pdf\n")
}

###############################################################
# 首次运行前的包安装 (已装好可忽略):
# install.packages(c("ggplot2","circlize","RColorBrewer","ggpubr","dplyr"))
# if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager")
# BiocManager::install(c("clusterProfiler","org.Hs.eg.db","enrichplot"))
###############################################################
