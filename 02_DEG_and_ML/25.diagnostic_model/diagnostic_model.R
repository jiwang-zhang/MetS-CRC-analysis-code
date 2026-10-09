# ============================================================
# Step 8a: 诊断模型 — 7 Hub 基因 ROC-AUC
# ============================================================
# 训练集: CRC合并矩阵 (GSE21510+GSE32323+GSE39582, 748样本)
# 验证集: GSE10950 (外部验证, GPL6104平台)
# 方法: 逻辑回归 → ROC-AUC (单基因 + 7基因联合)
#
# 需要安装的包:
#   install.packages(c("pROC", "ggplot2", "reshape2"))
#   limma 已有
#
# 使用前: 将 merged_after_batch_removal.csv 拷贝到本目录
#          GSE10950.csv 也拷贝过来(或修改下方路径)
# ============================================================

# --- 引用包 ---
library(limma)     # avereps
library(pROC)      # ROC曲线
library(ggplot2)
library(reshape2)

# --- 参数设置 ---
setwd("D:\\BaiduSyncdisk\\daixie-colorectal\\第二部分DEG\\25.diagnostic_model")

expFile_train <- "merged_after_batch_removal.csv"
# 如果不拷贝GSE10950, 改成原始路径:
# expFile_valid <- "D:/BaiduSyncdisk/daixie-colorectal/第二部分DEG/7.databases - 副本/GSE10950/Final_Results/GSE10950.csv"
expFile_valid <- "GSE10950.csv"

hubGenes <- c("GPAT3", "SGK1", "TMCC3", "PPARD", "PLAU", "LPCAT4", "TNFRSF21")

# ============================================================
# 读取表达矩阵的通用函数 (与ML脚本格式一致)
# ============================================================
read_expr <- function(file) {
  rt <- read.table(file, header = T, sep = ",", check.names = F)
  rt <- as.matrix(rt)
  rownames(rt) <- rt[, 1]
  exp <- rt[, 2:ncol(rt)]
  dn <- list(rownames(exp), colnames(exp))
  data <- matrix(as.numeric(exp), nrow = nrow(exp), dimnames = dn)
  data <- avereps(data)
  return(data)
}

# ============================================================
# 1. 读取训练集
# ============================================================
cat("=== 读取训练集 ===\n")
train_mat <- read_expr(expFile_train)
cat("训练集维度:", nrow(train_mat), "基因 x", ncol(train_mat), "样本\n")

hub_train <- hubGenes[hubGenes %in% rownames(train_mat)]
hub_missing <- hubGenes[!hubGenes %in% rownames(train_mat)]
if (length(hub_missing) > 0) {
  cat("  警告: 训练集缺失基因:", paste(hub_missing, collapse = ", "), "\n")
}
cat("  找到 hub 基因:", paste(hub_train, collapse = ", "), "\n")

# 提取7 hub表达量, 转置为 样本×基因
expr_train <- t(train_mat[hub_train, ])
group_train <- sub(".*_", "", rownames(expr_train))  # "tre" 或 "con"
cat("  分组统计:\n")
print(table(group_train))

# ============================================================
# 2. 读取验证集
# ============================================================
cat("\n=== 读取验证集 ===\n")
valid_mat <- read_expr(expFile_valid)
cat("验证集维度:", nrow(valid_mat), "基因 x", ncol(valid_mat), "样本\n")

hub_valid <- hub_train[hub_train %in% rownames(valid_mat)]
if (length(hub_valid) < length(hub_train)) {
  cat("  警告: 验证集缺失基因:", 
      paste(setdiff(hub_train, hub_valid), collapse = ", "), "\n")
}

expr_valid <- t(valid_mat[hub_valid, ])
group_valid <- sub(".*_", "", rownames(expr_valid))
cat("  分组统计:\n")
print(table(group_valid))

# ============================================================
# 3. 训练集 — 单基因 ROC
# ============================================================
cat("\n=== 训练集单基因 ROC ===\n")
auc_single <- numeric(length(hub_train))
names(auc_single) <- hub_train
ci_single <- matrix(NA, nrow = length(hub_train), ncol = 2)
rownames(ci_single) <- hub_train
colnames(ci_single) <- c("lower", "upper")

roc_list <- list()
for (i in seq_along(hub_train)) {
  gene <- hub_train[i]
  roc_obj <- roc(group_train, expr_train[, gene],
                 levels = c("con", "tre"), direction = "auto", quiet = TRUE)
  roc_list[[gene]] <- roc_obj
  auc_single[i] <- as.numeric(auc(roc_obj))
  ci_val <- ci.auc(roc_obj, conf.level = 0.95)
  ci_single[i, ] <- as.numeric(ci_val)[c(1, 3)]
  cat(sprintf("  %s: AUC = %.3f (95%% CI: %.3f-%.3f)\n",
              gene, auc_single[i], ci_single[i, 1], ci_single[i, 2]))
}

# --- 绘制单基因ROC图 (2x4布局, 7基因 + 1联合) ---
pdf("ROC_train_single.pdf", width = 12, height = 9)
par(mfrow = c(2, 4))
colors <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3",
            "#FF7F00", "#A65628", "#999999")
for (i in seq_along(hub_train)) {
  gene <- hub_train[i]
  plot(roc_list[[gene]], col = colors[i], lwd = 2,
       main = paste0(gene, "\nAUC = ", round(auc_single[i], 3)),
       legacy.axes = TRUE, xlab = "1 - Specificity", ylab = "Sensitivity")
}
# 第8个面板留空或放联合ROC (后面画)
plot.new()
dev.off()

# ============================================================
# 4. 训练集 — 7基因联合逻辑回归 ROC
# ============================================================
cat("\n=== 训练集联合 ROC ===\n")
train_df <- data.frame(group = factor(group_train, levels = c("con", "tre")),
                       expr_train)
model <- glm(group ~ ., data = train_df, family = binomial)
pred_train <- predict(model, type = "response")

roc_train <- roc(group_train, pred_train,
                 levels = c("con", "tre"), direction = "<", quiet = TRUE)
auc_train <- as.numeric(auc(roc_train))
ci_train <- as.numeric(ci.auc(roc_train, conf.level = 0.95))
cat(sprintf("  联合 AUC = %.3f (95%% CI: %.3f-%.3f)\n",
            auc_train, ci_train[1], ci_train[3]))

# ============================================================
# 5. 验证集 — 应用训练模型, 计算 ROC
# ============================================================
cat("\n=== 验证集联合 ROC ===\n")
# 用训练模型的系数预测验证集
valid_df <- as.data.frame(expr_valid)
colnames(valid_df) <- hub_valid
pred_valid <- predict(model, newdata = valid_df, type = "response")

roc_valid <- roc(group_valid, pred_valid,
                 levels = c("con", "tre"), direction = "<", quiet = TRUE)
auc_valid <- as.numeric(auc(roc_valid))
ci_valid <- as.numeric(ci.auc(roc_valid, conf.level = 0.95))
cat(sprintf("  验证集 AUC = %.3f (95%% CI: %.3f-%.3f)\n",
            auc_valid, ci_valid[1], ci_valid[3]))

# ============================================================
# 6. 绘制联合 ROC 对比图 (训练 vs 验证)
# ============================================================
pdf("ROC_combined.pdf", width = 6, height = 6)
plot(roc_train, col = "#E41A1C", lwd = 2.5, legacy.axes = TRUE,
     main = "7-Gene Diagnostic Model",
     xlab = "1 - Specificity", ylab = "Sensitivity")
plot(roc_valid, col = "#377EB8", lwd = 2.5, add = TRUE)
abline(a = 0, b = 1, lty = 2, col = "gray70")
legend("bottomright",
       c(paste0("Training set (AUC = ", round(auc_train, 3), ")"),
         paste0("Validation set (AUC = ", round(auc_valid, 3), ")")),
       col = c("#E41A1C", "#377EB8"), lwd = 2.5, bty = "n")
dev.off()

# ============================================================
# 7. 箱线图 — 7 hub 基因 肿瘤 vs 正常 (训练集)
# ============================================================
plot_df <- data.frame(Group = factor(group_train, levels = c("con", "tre"),
                                     labels = c("Normal", "Tumor")),
                      expr_train)
plot_long <- reshape2::melt(plot_df, id.vars = "Group",
                  variable.name = "Gene", value.name = "Expression")

p_box <- ggplot(plot_long, aes(x = Gene, y = Expression, fill = Group)) +
  geom_boxplot(outlier.size = 0.3, width = 0.6) +
  scale_fill_manual(values = c("#4DAF4A", "#E41A1C")) +
  theme_bw() +
  labs(x = "", y = "Expression Level",
       title = "7 Hub Gene Expression: Tumor vs Normal") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 11),
        plot.title = element_text(hjust = 0.5, size = 13))
ggsave("boxplot_7hub_train.pdf", p_box, width = 8, height = 5)

# ============================================================
# 8. 输出 AUC 汇总表
# ============================================================
auc_table <- data.frame(
  Gene = c(hub_train, "7-gene combined", "7-gene combined"),
  Dataset = c(rep("Training (GEO merged)", length(hub_train)),
              "Training (GEO merged)", "Validation (GSE10950)"),
  AUC = c(round(auc_single, 4), round(auc_train, 4), round(auc_valid, 4)),
  AUC_lower = c(round(ci_single[, 1], 4), round(ci_train[1], 4), round(ci_valid[1], 4)),
  AUC_upper = c(round(ci_single[, 2], 4), round(ci_train[3], 4), round(ci_valid[3], 4))
)
write.csv(auc_table, "AUC_summary.csv", row.names = FALSE)

cat("\n=== AUC 汇总 ===\n")
print(auc_table, row.names = FALSE)

cat("\n=== 诊断模型完成 ===\n")
cat("输出文件:\n")
cat("  ROC_train_single.pdf   — 训练集单基因ROC (7个面板)\n")
cat("  ROC_combined.pdf       — 联合ROC (训练 vs 验证)\n")
cat("  boxplot_7hub_train.pdf — 7基因表达箱线图\n")
cat("  AUC_summary.csv         — AUC汇总表\n")
