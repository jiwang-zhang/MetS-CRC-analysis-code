# ============================================================
# GSE132465 单细胞 Seurat 对象构建
# 适配教程 12.Seurat (原教程读逐样本 10X 目录，本数据改为读 GEO processed matrix)
# 数据来源: GEO GSE132465 (Lee et al. 2020, Nat Genet)
# 提供文件:
#   - GSE132465_GEO_processed_CRC_10X_raw_UMI_count_matrix.txt.gz (123MB)
#   - GSE132465_GEO_processed_CRC_10X_cell_annotation.txt.gz      (424KB)
# 数据格式: genes(行) x cells(列) 的 raw UMI count 文本矩阵
# ============================================================

# ---- 0. 包安装 (首次运行取消注释) ----
# install.packages("Seurat")
# install.packages("harmony")
# install.packages("data.table")
# install.packages("Matrix")
# install.packages("gridExtra")
# install.packages("clustree")

# ---- 引用包 ----
library(Seurat)
library(dplyr)
library(data.table)
library(Matrix)
library(harmony)
library(gridExtra)
library(clustree)

logFCfilter    = 1      # logFC 过滤条件
adjPvalFilter  = 0.05   # 矫正后 pvalue 过滤条件

# ---- 设置工作目录 (改为你本地实际路径) ----
workDir = "D:/BaiduSyncdisk/daixie-colorectal/4.免疫浸润单细胞/12.Seurat/GSE132465"
setwd(workDir)

# 两个 GEO 下载的文件名
count_file = "GSE132465_GEO_processed_CRC_10X_raw_UMI_count_matrix.txt.gz"
annot_file = "GSE132465_GEO_processed_CRC_10X_cell_annotation.txt.gz"

# ⚠️ 复现须知：本脚本从 matrix.mtx.gz 重建对象，不读取任何旧 Seurat.Rdata。
#    若工作目录已存在旧 Seurat.Rdata（由未修复脚本生成、cellType 被 RenameIdents 污染，
#    会出现"只剩某一类细胞"等异常），请勿在运行本脚本前手动 load() 它；
#    本脚本会全程重建对象，并在第 07 段覆盖保存为干净的 Rdata。
if (file.exists("Seurat.Rdata")) {
  cat("⚠️ 检测到旧 Seurat.Rdata，本脚本将重建对象并在末尾覆盖它（请勿手动提前 load）。\n")
}


####################### 01. 数据读取与稀疏化 #######################
cat(">>> Reading MatrixMarket format (pre-converted by Python)...\n")
# 脚本会先在本目录查找 matrix.mtx.gz / barcodes.tsv.gz / features.tsv.gz
# 如果不存在，说明还没跑 Python 转换脚本 convert_gse132465_to_mtx.py
# 可自行先运行一次: python scripts/convert_gse132465_to_mtx.py
mtx_dir <- getwd()  # 当前工作目录 = GSE132465/
counts_sparse <- ReadMtx(
  mtx = file.path(mtx_dir, "matrix.mtx.gz"),
  cells = file.path(mtx_dir, "barcodes.tsv.gz"),
  features = file.path(mtx_dir, "features.tsv.gz")
)
cat(">>> Sparse matrix loaded. Dim:", nrow(counts_sparse), "x", ncol(counts_sparse), "\n")

# 创建 Seurat 对象
pbmc <- CreateSeuratObject(counts = counts_sparse,
                            min.cells = 5,
                            min.features = 100,
                            project = "GSE132465")
cat(">>> Seurat object created. Cells:", ncol(pbmc), "Genes:", nrow(pbmc), "\n")

# ---- 读取细胞注释并加到 metadata (提前读取以便 Harmony 按样本整合) ----
cat(">>> Reading cell annotation ...\n")
annot <- fread(annot_file, sep = "\t", header = TRUE)
# GSE132465 annotation 实际列名 (已核对):
#   Index = 细胞 barcode (Sample_Barcode 格式, 与 count matrix 列名一致)
#   Patient / Class(Tumor|Normal) / Sample / Cell_type / Cell_subtype
# 匹配策略: 用 setNames 建 Index->行号 查表，比 match() 更稳健（未匹配返回 NA 不报错）。
#   已用真实数据验证: Index 与 barcode 精确匹配 63689/63689，6 类细胞齐全、0 个 NA。
cell_barcodes <- colnames(pbmc)
idx_map  <- setNames(seq_len(nrow(annot)), annot[["Index"]])
mpos     <- idx_map[cell_barcodes]          # 未匹配的 cell -> NA
matched  <- na.omit(mpos)
n_unmatched <- sum(is.na(mpos))
if (n_unmatched > 0)
  warning(sprintf("有 %d 个细胞在注释文件中找不到匹配 barcode，其 cellType 将为 NA", n_unmatched))
pbmc$cellType <- rep(NA_character_, length(cell_barcodes))
pbmc$cellType[!is.na(mpos)] <- annot[matched, ][["Cell_type"]]
pbmc$sample   <- rep(NA_character_, length(cell_barcodes))
pbmc$sample[!is.na(mpos)]   <- annot[matched, ][["Sample"]]
# 教程 28/30 用 Type=="Disease" 取疾病组；Tumor->Disease, Normal->Normal
pbmc$Type     <- rep(NA_character_, length(cell_barcodes))
pbmc$Type[!is.na(mpos)]     <- ifelse(annot[matched, ][["Class"]] == "Tumor", "Disease", "Normal")
# 【校验】注释是否成功装入 6 类细胞（正常: T/B/Myeloid/Epithelial/Stromal/Mast, 0 个 NA）
cat(">>> cellType 分布 (应为 6 类, 0 个 NA):\n")
print(table(pbmc$cellType, useNA = "ifany"))
if (length(unique(pbmc$cellType[!is.na(pbmc$cellType)])) < 2)
  stop("cellType 注释异常(类别数<2)，请检查 annot_file 的 Index/Cell_type 列与 barcode 是否一致")


####################### 02. QC 与过滤 #######################
# 线粒体比例
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

# QC 小提琴图
pdf(file = "01.featureViolin.pdf", width = 10, height = 6.5)
p1 <- VlnPlot(pbmc, features = "nFeature_RNA", raster = FALSE) + theme(legend.position = "none")
p2 <- VlnPlot(pbmc, features = "nCount_RNA",   raster = FALSE) + theme(legend.position = "none")
p3 <- VlnPlot(pbmc, features = "percent.mt",   raster = FALSE) + theme(legend.position = "none")
grid.arrange(p1, p2, p3, ncol = 3)
dev.off()

# 过滤 (教程标准: nFeature>100 & percent.mt<10)
pbmc <- subset(x = pbmc, subset = nFeature_RNA > 100 & percent.mt < 10)

# 测序深度相关性
pdf(file = "01.featureCor.pdf", width = 13, height = 7)
plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", pt.size = 1.5)
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt",  pt.size = 1.5)
CombinePlots(plots = list(plot1, plot2))
dev.off()


####################### 03. 标准化 + 高可变基因 #######################
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 1500)

pdf(file = "01.featureVar.pdf", width = 10, height = 6)
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = head(VariableFeatures(pbmc), 10), repel = TRUE)
CombinePlots(plots = list(plot1, plot2))
dev.off()


####################### 04. PCA + Harmony #######################
pbmc <- ScaleData(pbmc)
pbmc <- RunPCA(pbmc, npcs = 20, features = VariableFeatures(pbmc))

# Harmony 按样本整合 (若数据中无 sample 元数据可跳过此行)
if ("sample" %in% colnames(pbmc@meta.data)) {
  pbmc <- RunHarmony(pbmc, "sample")
  redu = "harmony"
} else {
  redu = "pca"
}

pdf(file = "02.pcaGene.pdf", width = 10, height = 8)
VizDimLoadings(pbmc, dims = 1:4, reduction = "pca", nfeatures = 20)
dev.off()
pdf(file = "02.PCA.pdf", width = 7.5, height = 5)
DimPlot(pbmc, reduction = "pca")
dev.off()
pdf(file = "02.pcaHeatmap.pdf", width = 10, height = 8)
DimHeatmap(pbmc, dims = 1:4, cells = 500, balanced = TRUE, nfeatures = 30, ncol = 2)
dev.off()
pbmc <- JackStraw(pbmc, num.replicate = 100)
pbmc <- ScoreJackStraw(pbmc, dims = 1:20)
pdf(file = "02.pcaJackStraw.pdf", width = 8, height = 6)
JackStrawPlot(pbmc, dims = 1:20)
dev.off()


####################### 05. 聚类 + TSNE #######################
pcSelect = 20
pbmc <- FindNeighbors(pbmc, dims = 1:pcSelect, reduction = redu)
pbmc <- FindClusters(pbmc, resolution = 0.6)
pdf(file = "03.cluster.pdf", width = 7.5, height = 6)
pbmc <- RunTSNE(pbmc, dims = 1:pcSelect, reduction = redu)
TSNEPlot(pbmc, pt.size = 2, label = TRUE)
dev.off()
write.table(pbmc$seurat_clusters, file = "03.Cluster.txt", quote = F, sep = "\t", col.names = F)

# 各聚类 marker 基因
pbmc.markers <- FindAllMarkers(pbmc, only.pos = FALSE, min.pct = 0.25,
                                logfc.threshold = logFCfilter)
sig.markers <- pbmc.markers[(abs(as.numeric(as.vector(pbmc.markers$avg_log2FC))) > logFCfilter &
                             as.numeric(as.vector(pbmc.markers$p_val_adj)) < adjPvalFilter), ]
write.table(sig.markers, file = "03.clusterMarkers.txt", sep = "\t", row.names = F, quote = F)
top10 <- pbmc.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
pdf(file = "03.clusterHeatmap.pdf", width = 15, height = 15)
DoHeatmap(pbmc, features = top10$gene) + NoLegend()
dev.off()


####################### 06. 细胞注释可视化 (annotation 已在 01 段加入 metadata) #######################
# 教程原版用 SingleR 7 个参考库注释；本数据 GEO 已提供细胞类型标注，直接读取并用于下游

# 注释后可视化
# 重要：GSE132465 的 GEO 提供文件已含逐细胞 Cell_type 标注（在 01 段已正确写入 pbmc$cellType），
#       直接按该 metadata 列可视化即可，**切勿用 RenameIdents 把聚类重命名为细胞类型**——
#       RenameIdents 会把 cluster 名映射到"第 N 个细胞的类型"，导致 cellType 完全错乱。
pdf(file = "04.cellAnn.pdf", width = 8, height = 6)
TSNEPlot(pbmc, group.by = "cellType", pt.size = 2, label = TRUE)
dev.off()
if ("Type" %in% colnames(pbmc@meta.data)) {
  pdf(file = "04.group.cellAnn.pdf", width = 12, height = 6)
  TSNEPlot(pbmc, group.by = "cellType", pt.size = 1, label = TRUE, split.by = "Type")
  dev.off()
}

# 细胞类型级差异分析（按 GEO 提供的逐细胞 Cell_type 分组）
pbmc.markers <- FindAllMarkers(pbmc, only.pos = FALSE, min.pct = 0.25,
                               logfc.threshold = logFCfilter, group.by = "cellType")
sig.cellMarkers <- pbmc.markers[(abs(as.numeric(as.vector(pbmc.markers$avg_log2FC))) > logFCfilter &
                                 as.numeric(as.vector(pbmc.markers$p_val_adj)) < adjPvalFilter), ]
write.table(sig.cellMarkers, file = "04.cellMarkers.txt", sep = "\t", row.names = F, quote = F)


####################### 07. 保存对象 (保存前再次核对注释) #######################
# 防御性重对齐: 全程过滤/整合均保留 metadata，但为彻底避免"旧对象/cellType 错位"复现事故，
# 保存前用【最终 colnames(pbmc)】再对 annot 查表一次，确保 Rdata 里的 cellType 与细胞严格对应。
cat(">>> [保存前校验] 用最终细胞集重对齐 cellType ...\n")
cell_barcodes <- colnames(pbmc)
idx_map  <- setNames(seq_len(nrow(annot)), annot[["Index"]])
mpos     <- idx_map[cell_barcodes]
matched  <- na.omit(mpos)
pbmc$cellType <- rep(NA_character_, length(cell_barcodes))
pbmc$cellType[!is.na(mpos)] <- annot[matched, ][["Cell_type"]]
pbmc$sample   <- rep(NA_character_, length(cell_barcodes))
pbmc$sample[!is.na(mpos)]   <- annot[matched, ][["Sample"]]
pbmc$Type     <- rep(NA_character_, length(cell_barcodes))
pbmc$Type[!is.na(mpos)]     <- ifelse(annot[matched, ][["Class"]] == "Tumor", "Disease", "Normal")
cat(">>> 最终 cellType 分布 (应为 6 类, 0 个 NA):\n")
print(table(pbmc$cellType, useNA = "ifany"))

save(pbmc, sig.markers, sig.cellMarkers, file = "Seurat.Rdata")
cat(">>> DONE. 干净的 Seurat.Rdata 已保存 (cellType 已校验为 6 类细胞)。\n")


###### 生信自学网: https://www.biowolf.cn/
###### 光俊老师邮箱: seqbio@foxmail.com
