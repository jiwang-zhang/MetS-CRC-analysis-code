# install.packages("remotes")
# remotes::install_github('cailab-tamu/scTenifoldKnk')
# install.packages("Seurat")
# install.packages("ggplot2")
# install.packages("dplyr")
# install.packages("ggrepel")


# 加载包
library(scTenifoldKnk)
library(Seurat)
library(ggplot2)
library(dplyr)
library(igraph)
library(ggrepel)
set.seed(12345)
options(bitmapType = "cairo")    # Adobe Illustrator 兼容 PDF

# GPAT3(已跑)+PLAU 虚拟敲除 (下游聚焦两交叉复现 hub)
hubGenes <- c("PLAU")  # 2026-08-24: GPAT3 KO 已稳定(166基因); 下游聚焦 GPAT3+PLAU, 仅补跑 PLAU KO
baseDir   <- "D:/BaiduSyncdisk/daixie-colorectal/4.免疫浸润单细胞/28.scTenifoldKnk"
seuratFile <- "D:/BaiduSyncdisk/daixie-colorectal/4.免疫浸润单细胞/12.Seurat/GSE132465/Seurat.Rdata"

for (target_gene in hubGenes) {
  cat("\n========================================\n")
  cat("Target gene:", target_gene, "\n")
  cat("========================================\n")

  # 每个基因独立子文件夹
  geneDir <- file.path(baseDir, target_gene)
  dir.create(geneDir, showWarnings = FALSE, recursive = TRUE)
  setwd(geneDir)

  # 导入Seurat的对象
  load(seuratFile)
  # 只保留疾病组的样品
  pbmc <- subset(pbmc, subset = Type == "Disease")
  cat("Disease cells:", ncol(pbmc), "\n")

  #提取表达矩阵
  countMat = GetAssayData(pbmc, layer = "counts")
  #提取高可变基因
  pbmc <- FindVariableFeatures(object=pbmc, selection.method="vst", nfeatures=5000)
  hvgs <- VariableFeatures(pbmc)
  data=as.data.frame(countMat[unique(c(target_gene,hvgs)),])
  rm(pbmc, countMat);gc()

  #执行虚拟敲除
  result <- scTenifoldKnk(countMatrix = data, 
                          gKO = target_gene,          #需要敲除的基因
                          qc_mtThreshold = 0.1,       #线粒体含量的阈值
                          qc_minLSize = 500,          #文库阈值(细胞测到的基因总数)
                          nc_nNet = 5,                #2026-08-03 调大: 5->10 子网络(更多扰动捕捉)
                          nc_nCells = 300,            #2026-08-03 调大: 300->1000 每网抽样细胞(降抽样噪声)
                          nCores = 1                   # Windows 多核无效
                          )

  #输出差异分析的结果
  df=result$diffRegulation
  df <- df[df$gene != target_gene, ]
  outTab=df[df$p.adj<0.05,]
  write.table(outTab, file="sigDiff.txt", sep="\t", quote=F, row.names=F)
  cat("Significant perturbed genes:", nrow(outTab), "\n")


  ###########################柱状图###########################
  #绘制柱状图
  top_genes <- head(df[order(-df$FC), ], 20)
  p1=ggplot(top_genes, aes(x=reorder(gene, FC), y=FC)) +
    geom_bar(stat='identity', fill='#5A9BD4') +
    coord_flip() + 
    labs(title=paste0(target_gene, " KO: Top 20 Differentially Regulated Genes"), x="Gene", y="FC") +
    theme_minimal() + theme(plot.title = element_text(hjust = 0.5))
  #输出图形
  pdf(file="barplot.pdf", width=6, height=5)
  print(p1)
  dev.off()


  ###########################火山图###########################
  #准备火山图的数据
  df$log_p.adj <- -log10(df$p.adj)
  df$significant <- ifelse(df$p.adj < 0.05, "Significant", "Not significant")
  # 仅标注效应最强的 Top 30（按 |Z| 排序），其余显著基因保留红点但不标名
  label_genes <- df[df$p.adj < 0.05, , drop = FALSE]
  label_genes <- label_genes[order(-abs(label_genes$Z)), , drop = FALSE][1:30, , drop = FALSE]

  # 绘制火山图
  y_upper <- quantile(df$log_p.adj, 0.999, na.rm = TRUE)
  p2=ggplot(df, aes(x=Z, y=log_p.adj, color=significant)) +
    geom_point(alpha=0.7, size=1) +
    scale_color_manual(values = c("Significant" = "red", "Not significant" = "gray50")) +
    geom_hline(yintercept=-log10(0.05), linetype="dashed", color="red") +
    geom_text_repel(data=label_genes, aes(label=gene), size=3, max.overlaps=30) +
    labs(title=paste0(target_gene, " KO"), x="Z-score", y="-log10(p.adj)") +
    theme_classic() + 
    coord_cartesian(ylim = c(0, y_upper)) + theme(legend.position = "none")
  pdf(file="vol.pdf", width=6, height=5)
  print(p2)
  dev.off()


  ###########################饼图###########################
  sig_count <- table(df$significant)
  sig_df <- as.data.frame(sig_count)
  colnames(sig_df) <- c("category", "count")
  sig_df$percentage <- paste0(round(sig_df$count / sum(sig_df$count) * 100, 1), "%")
  p3 <- ggplot(sig_df, aes(x="", y=count, fill=category)) +
    geom_bar(stat="identity", width=1) +
    geom_text(aes(label=percentage), position=position_stack(vjust=0.5), size=4) + 
    coord_polar("y", start=0) + 
    scale_fill_manual(values=c("Significant"="red", "Not significant"="lightgray")) + 
    labs(title=paste0(target_gene, " KO: Significant Genes"), fill="") +
    theme_minimal() + 
    theme(axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text = element_blank(),
      panel.grid = element_blank(), legend.position = "right", plot.title = element_text(hjust = 0.5))
  pdf(file="pie.pdf", width=6, height=5)
  print(p3)
  dev.off()

  cat("Output saved to:", geneDir, "\n")
}

cat("\n>>> ALL DONE: 2 genes (GPAT3/PLAU) virtual KO completed.\n")

######生信自学网: https://www.biowolf.cn/
######课程链接1: https://shop119322454.taobao.com
######课程链接2: https://ke.biowolf.cn
######课程链接3: https://ke.biowolf.cn/mobile
######光俊老师邮箱: seqbio@foxmail.com
######光俊老师微信: eduBio
