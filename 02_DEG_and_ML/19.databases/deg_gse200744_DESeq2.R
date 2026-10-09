# GSE200744 差异分析 —— DESeq2 (MetS vs Normal)
# 依赖: Bioconductor DESeq2
# 数据: 19.databases/GSE200744/GSE200744_RNAseq_count_in_Mets.txt.gz
# 用法: 在 19.databases/ 目录下运行 Rscript deg_gse200744_DESeq2.R

library(DESeq2)

f <- "GSE200744/GSE200744_RNAseq_count_in_Mets.txt.gz"
df <- read.delim(f, check.names = FALSE)

# 分组直接来自计数矩阵的列名: mets_* / premets_* / normal_*
sample_cols <- grep("^(mets|premets|normal)_", colnames(df), value = TRUE)
grp <- ifelse(grepl("^mets_",   sample_cols), "MetS",
       ifelse(grepl("^normal_", sample_cols), "Normal", "PreMetS"))
keep <- sample_cols[grp %in% c("MetS", "Normal")]

counts <- df[, keep]
counts <- counts[rowSums(counts) >= 10, ]          # 低表达过滤
counts <- round(as.matrix(counts))                  # 确保整数计数

meta <- data.frame(
  condition = factor(ifelse(grepl("^mets_", keep), "MetS", "Normal"),
                     levels = c("Normal", "MetS")),    # 参照=Normal
  row.names = keep
)

dds <- DESeqDataSetFromMatrix(countData = counts, colData = meta, design = ~ condition)
dds <- DESeq(dds)
res <- results(dds, name = "condition_MetS_vs_Normal")   # log2FC = MetS - Normal
res <- as.data.frame(res)
res$Symbol   <- df[rownames(counts), "Gene Symbol"]
res$Ensembl  <- df[rownames(counts), "Ensembl Gene ID"]
res <- res[order(res$padj, na.last = TRUE), ]

sig <- res[!is.na(res$padj) & res$padj < 0.05 & abs(res$log2FoldChange) > 0.5, ]
message(sprintf("MetS=%d, Normal=%d, 基因=%d", sum(grp[grp%in%c("MetS","Normal")]=="MetS"),
                sum(grp[grp%in%c("MetS","Normal")]=="Normal"), nrow(counts)))
message(sprintf("DEG(|log2FC|>0.5 & padj<0.05): %d (上调 %d, 下调 %d)",
                nrow(sig), sum(sig$log2FoldChange>0), sum(sig$log2FoldChange<0)))

write.csv(res, "GSE200744/DEG_MetS_vs_Normal_R.csv", row.names = TRUE)
write.csv(sig, "GSE200744/DEG_MetS_vs_Normal_R.sig.csv", row.names = TRUE)
