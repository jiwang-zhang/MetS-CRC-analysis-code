# ============================================================
# 单细胞补充图：对标模板 Fig 9B / 9C / 9D（GPRC5B, IJBM 2026）
# 我们自己的编号 = Fig 5B / 5C / 5D
#   5B: 细胞比例 logFC 发散条图        (对标模板 Fig 9B)
#   5C: hub 基因细胞类型表达 分组柱状图 (对标模板 Fig 9C, 已改为柱状图)
#   5D: hub 基因 x 细胞类型 气泡图      (对标模板 Fig 9D, bubble summary)
# ============================================================

library(Seurat)
library(dplyr)
library(ggplot2)
library(tidyr)   # pivot_wider / pivot_longer

workDir <- "D:/BaiduSyncdisk/daixie-colorectal/4.免疫浸润单细胞"
setwd(workDir)

load("12.Seurat/GSE132465/Seurat.Rdata")
cat(">>> pbmc loaded:", ncol(pbmc), "cells,", nrow(pbmc), "genes\n")

# ---- 诊断（数据核对，避免“图空白但不知为何”）----
cat(">>> Type 水平:", paste(unique(pbmc$Type), collapse = " | "), "\n")
cat(">>> cellType 唯一类别数:", length(unique(pbmc$cellType)), "\n")
print(table(pbmc$cellType, useNA = "ifany"))

hub_genes <- c("GPAT3", "PLAU")   # 两交叉复现 hub（GSE10950 ∩ TCGA/GEPIA2）

# ============================================================
# Fig 5B: 细胞比例 logFC 发散条图（对标模板 Fig 9B）
#   x = log2(比例倍数变化 Disease/Normal)，右扩=Disease 富集，左缩=减少
# ============================================================

freq <- pbmc@meta.data %>%
  group_by(Type, cellType) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(Type) %>%
  mutate(freq = n / sum(n)) %>%
  ungroup()

# 关键修复：freq 里残留的 n(计数) 列必须剔除，否则 pivot_wider 会把 n 也当
# 成分组键 -> cellType 不再唯一 -> Disease/Normal 两列错位成 NA -> logFC 全 NA
# -> 条宽全 0 -> 空白 PDF（这正是此前 Fig 9B 空白的真正根因，非列名不匹配）。
freq_wide <- freq %>% select(-n) %>% pivot_wider(names_from = Type, values_from = freq)

# 动态识别 disease / normal 列：绝不硬编码 "Disease"/"Normal"，
# 防止 pivot 列名与预期不符 -> logFC 全 NA -> 条宽全 0 -> 空白 PDF。
type_cols <- setdiff(names(freq_wide), "cellType")
norm_col  <- type_cols[grepl("norm", tolower(type_cols))]
dis_col   <- setdiff(type_cols, norm_col)
if (length(dis_col) != 1 || length(norm_col) != 1) {
  stop("Type 列识别失败，实际列名: ", paste(type_cols, collapse = ", "))
}
cat(">>> 比例宽表前几行 (dis_col =", dis_col, ", norm_col =", norm_col, "):\n")
print(head(freq_wide))

eps <- 1e-9   # 防止任一分组计数为 0 时 logFC 出 Inf
freq_wide$logFC <- log2((freq_wide[[dis_col]] + eps) / (freq_wide[[norm_col]] + eps))

# y 轴按 logFC 升序；先转字符再用「排序后的唯一值」构造因子水平，
# 避免 cellType 因子带多余/重复水平触发 “factor level duplicated”。
freq_wide$cellType <- as.character(freq_wide$cellType)
lvl <- freq_wide$cellType[order(freq_wide$logFC)]
lvl <- lvl[!duplicated(lvl)]
freq_wide$cellType <- factor(freq_wide$cellType, levels = lvl)

p5b <- ggplot(freq_wide, aes(x = logFC, y = cellType, fill = logFC > 0)) +
  geom_col(width = 0.7) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
  scale_fill_manual(values = c("TRUE" = "#E41A1C", "FALSE" = "#377EB8"),
                    labels = c("FALSE" = "Normal enriched",
                               "TRUE"  = "Disease enriched"),
                    name = "") +
  labs(x = expression(log[2]~"(fold change) in cell proportion (Disease / Normal)"),
       y = NULL,
       title = "Cell-type composition changes: CRC vs Normal") +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank())

ggsave("29.fig5B_cellprop_diverging.pdf", p5b, width = 8, height = 5, dpi = 300)
cat(">>> Fig 5B 保存\n")


# ============================================================
# Fig 5C: hub 基因细胞类型特异性表达 —— 分组柱状图（对标模板 Fig 9C, 改为柱状图）
#   每个基因一个分面；x = cellType, fill = Type(Disease/Normal), y = 平均归一化表达
# ============================================================

# 分别取 Disease / Normal 子集，各自按 cellType 求平均表达，避免依赖
# AverageExpression 内部列名拼接规则，最稳。
subD <- subset(pbmc, subset = Type == dis_col)
subN <- subset(pbmc, subset = Type == norm_col)

avgD <- as.data.frame(t(AverageExpression(subD, features = hub_genes,
                                          group.by = "cellType", slot = "data")$RNA))
avgN <- as.data.frame(t(AverageExpression(subN, features = hub_genes,
                                          group.by = "cellType", slot = "data")$RNA))
avgD$cellType <- rownames(avgD); avgD$Type <- dis_col
avgN$cellType <- rownames(avgN); avgN$Type <- norm_col

avg_long <- bind_rows(avgD, avgN) %>%
  pivot_longer(cols = all_of(hub_genes), names_to = "gene", values_to = "avg_expr")

p5c <- ggplot(avg_long, aes(x = cellType, y = avg_expr, fill = Type)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.75) +
  facet_wrap(~ gene, scales = "free_y") +
  labs(x = NULL, y = "Average normalized expression",
       title = "Hub gene expression by cell type (CRC vs Normal)") +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom",
        panel.grid.minor = element_blank())

ggsave("29.fig5C_hub_expression_bar.pdf", p5c, width = 9, height = 5, dpi = 300)
cat(">>> Fig 5C 保存\n")

# 注：原 Fig 5D（DotPlot 气泡图，含 pct.exp 维度）已按用户决定删除，
# 单细胞主图只保留 Fig 5B(比例发散条图) + Fig 5C(柱状图) 两张，
# 省一个 Sci Rep display item；pct.exp 维度如需可挪 Supplementary。

cat(">>> DONE\n")
