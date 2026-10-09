# 仅重出 04.group.cellAnn.pdf（Disease vs Normal 分面 tSNE）
# 无需重跑整条 Seurat 流程：Seurat.Rdata 已含 cellType / Type / tsne 降维
setwd("D:/BaiduSyncdisk/daixie-colorectal/4.免疫浸润单细胞/12.Seurat/GSE132465")
load("Seurat.Rdata")   # 542MB，含 pbmc（cellType + Type + tsne）

# 原脚本 line 192 的生成方式（Seurat 旧版 TSNEPlot）
pdf(file = "04.group.cellAnn.pdf", width = 12, height = 6)
TSNEPlot(pbmc, group.by = "cellType", pt.size = 1, label = TRUE, split.by = "Type")
dev.off()

# 若上面 TSNEPlot 仍报错/为空，改用等价 DimPlot（新版 Seurat 推荐写法）
# pdf(file = "04.group.cellAnn.pdf", width = 12, height = 6)
# DimPlot(pbmc, reduction = "tsne", group.by = "cellType", pt.size = 1,
#         label = TRUE, split.by = "Type")
# dev.off()

cat(">>> 04.group.cellAnn.pdf 已重出；请确认页数应 > 0\n")
