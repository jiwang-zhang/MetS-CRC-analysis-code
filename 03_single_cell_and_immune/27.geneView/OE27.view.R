#install.packages("clustree")
#install.packages("Seurat")
#install.packages("harmony")

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("GSVA")
#BiocManager::install("GSEABase")
#BiocManager::install("limma")
#BiocManager::install("SingleR")
#BiocManager::install("celldex")
#BiocManager::install("monocle")

#install.packages("devtools")
#devtools::install_github('immunogenomics/presto')

#install.packages("assertthat")
#install.packages("SCpubr")



#引用包 (DotPlot 只需 Seurat；其余为教程原版，按需保留)
library(Seurat)
library(dplyr)
library(magrittr)
library(ggplot2)
options(bitmapType = "cairo")    # Adobe Illustrator 兼容 PDF

#设置工作目录 (改为你本地实际路径)
workDir="D:/BaiduSyncdisk/daixie-colorectal/4.免疫浸润单细胞/27.geneView"
setwd(workDir)

#导入单细胞的对象
# ⚠️ 2026-07-30 修正: 原 load("Seurat.Rdata") 读的是本目录 7/27 旧副本(cellType 被污染)。
#    改为绝对路径指向 12.Seurat 修复版(6 类细胞齐全, 0 NA)。本目录旧 Seurat.Rdata 请勿再用。
load("D:/BaiduSyncdisk/daixie-colorectal/4.免疫浸润单细胞/12.Seurat/GSE132465/Seurat.Rdata")

# 修复 bug: barcodes.tsv.gz 的 \r 导致细胞名与注释不匹配
# colnames 含 \r 后缀, 匹配时大部分细胞丢失
colnames(pbmc) <- gsub("\r", "", colnames(pbmc))
cat("Cells after fixing barcode names:", ncol(pbmc), "\n")

# 重新读取注释并匹配 (解决"Only Epithelial cells"问题)
annot_file <- "../12.Seurat/GSE132465/GSE132465_GEO_processed_CRC_10X_cell_annotation.txt.gz"
if (file.exists(annot_file)) {
  annot <- data.table::fread(annot_file, sep = "\t", header = TRUE)
  annot_matched <- annot[match(colnames(pbmc), annot[["Index"]]), ]
  # 检查匹配率
  match_rate <- sum(!is.na(annot_matched[["Cell_type"]])) / ncol(pbmc)
  cat(sprintf("Annotation match rate: %.1f%% (%d/%d)\n", 
              match_rate * 100, sum(!is.na(annot_matched[["Cell_type"]])), ncol(pbmc)))
  pbmc$cellType <- annot_matched[["Cell_type"]]
  pbmc$sample   <- annot_matched[["Sample"]]
  pbmc$Type     <- ifelse(annot_matched[["Class"]] == "Tumor", "Disease", "Normal")
  
  # 恢复 cellType 为 active identity
  Idents(pbmc) <- pbmc$cellType
  cat("Cell types:", paste(levels(Idents(pbmc)), collapse=", "), "\n")
}

#读取交集基因的文件
geneRT=read.table("interGenes.List.txt", header=F, sep="\t", check.names=F)
interGenes=as.vector(geneRT[,1])
interGenes=interGenes[interGenes %in% rownames(pbmc)]

#绘制气泡图
pdf(file="Bubble.pdf", width=7.5, height=5)
DotPlot(object = pbmc, features = interGenes) + theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10))
dev.off()


######生信自学网: https://www.biowolf.cn/
######课程链接1: https://shop119322454.taobao.com
######课程链接2: https://ke.biowolf.cn
######课程链接3: https://ke.biowolf.cn/mobile
######光俊老师邮箱: seqbio@foxmail.com
######光俊老师微信