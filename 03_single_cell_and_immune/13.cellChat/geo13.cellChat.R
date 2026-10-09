#if (!require("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")
#BiocManager::install("ComplexHeatmap")

#install.packages("dplyr")
#install.packages("NMF")
#install.packages("patchwork")
#install.packages("igraph")
#install.packages("ggplot2")
#install.packages("ggalluvial")
#install.packages("circlize")
#install.packages("svglite")

#install.packages("devtools")
library(devtools)
devtools::install_github("sqjin/CellChat")


#引用包
library(limma)
library(NMF)
library(ggplot2)
library(ggalluvial)
library(svglite)
library(CellChat)

#设置工作目录
setwd("D:/BaiduSyncdisk/daixie-colorectal/4.免疫浸润单细胞/13.cellChat")

#准备细胞通讯的数据
# ?? 2026-07-30 修正: 改为绝对路径指向 12.Seurat 修复版(本目录无 Rdata 副本)。
load("D:/BaiduSyncdisk/daixie-colorectal/4.免疫浸润单细胞/12.Seurat/GSE132465/Seurat.Rdata")
pbmc <- subset(pbmc, subset = Type == "Disease")
expMatrix=as.matrix(pbmc@assays$RNA$data)
meta=as.data.frame(pbmc$cellType)
colnames(meta)[1]="labels"
row.names(meta)=colnames(pbmc)   # 2026-08-04 修复: 直接用细胞barcode(=expMatrix列名)对齐, 避免 seurat_clusters 缺失时报错/错位

#创建cellchat对象
cellchat <- createCellChat(object = expMatrix, meta = meta, group.by = "labels")
cellchat <- setIdent(cellchat, ident.use="labels")
groupSize <- as.numeric(table(cellchat@idents))      #每个细胞类型的细胞数目

#导入配体和受体数据库
CellChatDB <- CellChatDB.human       #如果物种是小鼠改成CellChatDB.mouse
#CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling")
cellchat@DB <- CellChatDB

#查看配体和受体对的类型
pdf(file="COMM01.DatabaseCategory.pdf", width=7, height=5)
showDatabaseCategory(CellChatDB)
dev.off()

#对表达数据进行预处理
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)      #识别每种细胞中过表达的基因
cellchat <- identifyOverExpressedInteractions(cellchat)      #识别基因之间的相互作用
cellchat <- projectData(cellchat, PPI.human)  

#计算细胞通讯的概率
cellchat <- computeCommunProb(cellchat)
#过滤掉小于10个细胞的细胞通讯
cellchat <- filterCommunication(cellchat, min.cells = 10)
#输出细胞之间的通讯关系
df.net=subsetCommunication(cellchat)
write.table(file="COMM02.Comm.network.xls", df.net, sep="\t", row.names=F, quote=F)

#在信号通路的水平进一步推测细胞间的通讯, 推断通路水平的互作网络
cellchat <- computeCommunProbPathway(cellchat)

#对计算结果汇总整合，展示整体细胞通讯状态
cellchat <- aggregateNet(cellchat)
#输出细胞通讯的图形(相互作用数目的图形)
pdf(file="COMM03.cellNetworkCount.pdf", width=7, height=6)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
dev.off()
#输出细胞通讯的图形(相互作用强度的图形)
pdf(file="COMM04.cellNetworkWeight.pdf", width=7, height=6)
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction strength")
dev.off()

#分细胞类型展示(把单个细胞类型提取出来,观察这个细胞与其他细胞的通讯)
pdf(file="COMM05.singleCell.pdf", width=9, height=7.5)
weight_mat <- cellchat@net$weight
par(mfrow = c(2,3), mgp=c(0,0,0), xpd=TRUE)
for (cel in unique(cellchat@idents)){
  cir_mat <- matrix(0, nrow = nrow(weight_mat), ncol = ncol(weight_mat), dimnames = dimnames(weight_mat))
  cir_mat[cel, ] <- weight_mat[cel, ]
  netVisual_circle( cir_mat, vertex.weight= groupSize, weight.scale= T,edge.weight.max = max(weight_mat), vertex.label.cex=0.8,title.name=cel)
}
dev.off()

#对于每个细胞, 得到单独的文件
for (cel in unique(cellchat@idents)){
  cir_mat <- matrix(0, nrow = nrow(weight_mat), ncol = ncol(weight_mat), dimnames = dimnames(weight_mat))
  cir_mat[cel, ] <- weight_mat[cel, ]
  pdf(file=paste0("COMM05.", cel, ".pdf"), width=6.5, height=5.5)
  netVisual_circle( cir_mat, vertex.weight= groupSize, weight.scale= T,edge.weight.max = max(weight_mat), vertex.label.cex=0.8,title.name=cel)
  dev.off()
}

#绘制细胞通讯的气泡图
pdf(file="COMM06.bubble.pdf", width=7, height=5)
netVisual_bubble(cellchat, remove.isolate = FALSE, angle.x=45, font.size = 6, thresh = 0.001)
dev.off()


######生信自学网: https://www.biowolf.cn/
######课程链接1: https://shop119322454.taobao.com
######课程链接2: https://ke.biowolf.cn
######课程链接3: https://ke.biowolf.cn/mobile
######光俊老师邮箱: seqbio@foxmail.com
######光俊老师微信
rm(list = ls())   # 删除全局环境所有对象
gc()
# 27：重出 Bubble.pdf（之前的图疑似同步盘丢了，重跑生成）
source("D:/BaiduSyncdisk/daixie-colorectal/4.免疫浸润单细胞/27.geneView/OE27.view.R")
rm(list = ls())   # 删除全局环境所有对象
gc()
# 31：补齐 SGK1/TMCC3 最小富集（<1 门槛已改好，读 30 的 sigDiff.txt）
source("D:/BaiduSyncdisk/daixie-colorectal/4.免疫浸润单细胞/31.Enrich/OE31.Enrich.R")
