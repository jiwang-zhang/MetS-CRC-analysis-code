#install.packages("dplyr")
#install.packages("ggplot2")
#install.packages("ggsci")
#install.packages('devtools')
#devtools::install_github("sajuukLyu/ggunchull", type = "source")
#devtools::install_github('junjunlab/scRNAtoolVis')


#引用包
library(dplyr)
library(ggplot2)
library(Seurat)
library(ggunchull)
library(scRNAtoolVis)

logFCfilter=1           #logFC的过滤条件
adjPvalFilter=0.05     #矫正后的pvalue的过滤条件

#设置工作目录
workDir="D:/BaiduSyncdisk/daixie-colorectal/4.免疫浸润单细胞/14.scRNAtoolVis"
setwd(workDir)

#导入单细胞的数据
# ⚠️ 2026-07-30 修正: 改为绝对路径指向 12.Seurat 修复版(本目录无 Rdata 副本)。
load("D:/BaiduSyncdisk/daixie-colorectal/4.免疫浸润单细胞/12.Seurat/GSE132465/Seurat.Rdata")

# 2026-08-04 修复: 计算各细胞类型差异 marker 基因, 供 markerVolcano / jjVolcano 使用。
# 原脚本 allGroupMarkers 未定义, 直接调用会报 object 'allGroupMarkers' not found。
# 顺带打印 meta.data 列名, 用于核对 sample.name 列确为小写 "sample"。
cat("meta.data columns:", paste(colnames(pbmc@meta.data), collapse = ", "), "\n")
allGroupMarkers <- FindAllMarkers(pbmc, group.by = "cellType", only.pos = FALSE,
                                  min.pct = 0.1, logfc.threshold = 0.25)

#绘制细胞注释的图形
pdf(file="01.cellType.pdf", width=8, height=6)
clusterCornerAxes(object = pbmc,
                  reduction = 'tsne',      #降维的方法 (修复版 Rdata 由 geo12 的 RunTSNE 生成, 含 tsne 槽, 原样保留)
                  clusterCol = "cellType",      #根据细胞类型定义细胞的颜色
                  noSplit = T, cellLabel = T,    #是否展示细胞的标签
			    addCircle = T, cicDelta = 0.2, cicAlpha=0.2, cicLineSize=1.5,cicLineColor="grey50",
				themebg = 'bwCorner',nbin=200,
                  cellLabelSize = 5,
                  show.legend = F)
dev.off()

#绘制细胞比例的柱状图
pdf(file="02.barplot.pdf", width=6, height=6)
cellRatioPlot(
  object = pbmc,
  sample.name = "sample",
  celltype.name = "cellType",
  col.width = 0.7
)
dev.off()

#绘制marker基因的火山图
pdf(file="03.markerVolcano.pdf", width=10, height=6)
markerVolcano(markers = allGroupMarkers, topn = 5, labelCol = ggsci::pal_npg()(9), log2FC=logFCfilter)
dev.off()

#绘制火山图1
pdf(file="04.jjVolcano1.pdf", width=12, height=6)
jjVolcano(diffData = allGroupMarkers, legend.position=c(0.08,0.95), log2FC.cutoff=logFCfilter, adjustP.cutoff=adjPvalFilter)+labs(x = "Cell type")
dev.off()

#绘制转置的火山图2
pdf(file="04.jjVolcano2.pdf", width=9, height=7)
jjVolcano(diffData = allGroupMarkers, legend.position=c(0.05,0.95), flip = T, log2FC.cutoff=logFCfilter, adjustP.cutoff=adjPvalFilter)+labs(x = "Cell type")
dev.off()

#绘制环形的火山图3
pdf(file="04.jjVolcano3.pdf", width=8, height=7)
jjVolcano(diffData = allGroupMarkers, legend.position=c(0.1,0.85), polar = T, size=3,fontface='italic', log2FC.cutoff=logFCfilter, adjustP.cutoff=adjPvalFilter)+labs(x = "Cell type")
dev.off()


######生信自学网: https://www.biowolf.cn/
######课程链接1: https://shop119322454.taobao.com
######课程链接2: https://ke.biowolf.cn
######课程链接3: https://ke.biowolf.cn/mobile
######光俊老师邮箱: seqbio@foxmail.com
######光俊老师微信: eduBio


