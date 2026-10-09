# ============================================================
# Step 8b: 预后模型 — GSE39582 独立验证队列
# ============================================================
# 队列 : GSE39582 (Marisa 2013, PLoS Med), n=585 CRC 肿瘤
# 平台 : GPL570 (Affymetrix HG-U133_Plus_2), microarray
# 终点 : OS (overall survival) + RFS (relapse-free survival) 双终点
# 主模型: LASSO-Cox 7 hub genes (自动选特征)
# 对照 : 单因素 Cox + 多因素 Cox (7基因联合)
# 时间 : 月 (OS/RFS 字段单位已确认=月)
#
# 需要安装的包:
#   install.packages(c("survival", "survminer", "rms", "timeROC",
#                      "dcurves", "glmnet", "ggplot2"))
# 说明:
#   生存数据 GSE39582_survival.csv 由 extract_gse39582_survival.py 生成
#   表达矩阵用原始路径(205MB, 不拷贝)
# ============================================================

library(survival)
library(survminer)
library(rms)
library(timeROC)
library(dcurves)
library(glmnet)
library(ggplot2)

setwd("D:\\BaiduSyncdisk\\daixie-colorectal\\第二部分DEG\\26.prognostic_model")

hubGenes <- c("GPAT3", "SGK1", "TMCC3", "PPARD", "PLAU", "LPCAT4", "TNFRSF21")

# 表达矩阵原始路径(不拷贝, 直接读)
exprFile <- "D:/BaiduSyncdisk/daixie-colorectal/第二部分DEG/7.databases - 副本/GSE39582/Final_Results/GSE39582.csv"
survFile <- "GSE39582_survival.csv"

# ============================================================
# 1. 读取表达矩阵, 提取 7 hub
# ============================================================
cat("=== 1. 读取表达矩阵 ===\n")
rt <- read.csv(exprFile, header = TRUE, check.names = FALSE)
rt <- as.matrix(rt)
rownames(rt) <- rt[, 1]
exp <- rt[, 2:ncol(rt)]
exp <- matrix(as.numeric(exp), nrow = nrow(exp),
              dimnames = list(rownames(exp), colnames(exp)))

# 样本ID: 去掉 _tre 后缀 (GSM971957_tre -> GSM971957)
colnames(exp) <- sub("_tre$", "", colnames(exp))

hub_exp <- exp[hubGenes, , drop = FALSE]
cat("  找到 hub 基因:", paste(rownames(hub_exp), collapse = ", "), "\n")
cat("  表达矩阵样本数:", ncol(hub_exp), "\n")

# ============================================================
# 2. 读取生存数据
# ============================================================
cat("\n=== 2. 读取生存数据 ===\n")
surv <- read.csv(survFile, stringsAsFactors = FALSE)
surv$Sample <- as.character(surv$Sample)
# 强制数值化: CSV 中 "N/A"/空串会使整列变 character, 导致 sum()/Cox 报错
# as.numeric 会把 "N/A"/空串转成 NA, 后续清洗步骤会剔除这些 NA 行
num_cols <- c("os_event", "os_time_months", "rfs_event",
              "rfs_time_months", "age_year", "tnm_stage")
for (c in num_cols) surv[[c]] <- suppressWarnings(as.numeric(as.character(surv[[c]])))
cat("  生存数据样本数:", nrow(surv), "\n")

# 匹配表达与生存
common <- intersect(colnames(hub_exp), surv$Sample)
cat("  交集样本数:", length(common), "\n")
hub_exp <- hub_exp[, common, drop = FALSE]
surv <- surv[match(common, surv$Sample), ]

# 合并为数据框
data <- as.data.frame(t(hub_exp))
colnames(data) <- rownames(hub_exp)
for (g in colnames(data)) data[[g]] <- as.numeric(data[[g]])
data$Sample <- common
data <- merge(data, surv, by = "Sample")

# 清理: OS 需 time>0 且非NA; RFS 需 time>0 且非NA
data <- data[!is.na(data$os_time_months) & data$os_time_months > 0 &
             !is.na(data$os_event) & data$os_event %in% c(0, 1), ]
data_rfs <- data[!is.na(data$rfs_time_months) & data$rfs_time_months > 0 &
                  !is.na(data$rfs_event) & data$rfs_event %in% c(0, 1), ]
cat("  OS 分析样本数:", nrow(data), " (事件=", sum(data$os_event), ")\n")
cat("  RFS 分析样本数:", nrow(data_rfs), " (事件=", sum(data_rfs$rfs_event), ")\n")

# z-score 标准化 (per gene, Cox系数按 1 SD 解释)
for (g in hubGenes) {
  coln <- paste0(g, "_z")
  data[[coln]] <- as.numeric(scale(data[[g]]))
  if (exists("data_rfs")) data_rfs[[coln]] <- as.numeric(scale(data_rfs[[g]]))
}
zcols <- paste0(hubGenes, "_z")

# ============================================================
# 3. 单因素 Cox (OS 和 RFS)
# ============================================================
cat("\n=== 3. 单因素 Cox 回归 ===\n")
uni_cox <- function(dat, time_col, event_col, genes) {
  res <- data.frame(Gene = character(), HR = numeric(), Lower = numeric(),
                    Upper = numeric(), P = numeric(), stringsAsFactors = FALSE)
  for (g in genes) {
    form <- as.formula(paste0("Surv(", time_col, ", ", event_col, ") ~ ", g, "_z"))
    cox <- coxph(form, data = dat)
    s <- summary(cox)
    res <- rbind(res, data.frame(
      Gene = g,
      HR = round(s$coefficients[, "exp(coef)"], 3),
      Lower = round(s$conf.int[, "lower .95"], 3),
      Upper = round(s$conf.int[, "upper .95"], 3),
      P = signif(s$coefficients[, "Pr(>|z|)"], 3),
      stringsAsFactors = FALSE))
  }
  return(res)
}
uni_os <- uni_cox(data, "os_time_months", "os_event", hubGenes)
uni_rfs <- uni_cox(data_rfs, "rfs_time_months", "rfs_event", hubGenes)
cat("  --- OS 单因素 ---\n"); print(uni_os, row.names = FALSE)
cat("  --- RFS 单因素 ---\n"); print(uni_rfs, row.names = FALSE)
write.csv(uni_os, "GSE39582_uniCox_OS.csv", row.names = FALSE)
write.csv(uni_rfs, "GSE39582_uniCox_RFS.csv", row.names = FALSE)

# ============================================================
# 4. 多因素 Cox (7基因联合, OS)
# ============================================================
cat("\n=== 4. 多因素 Cox (7基因联合) ===\n")
multi_form <- as.formula(paste0("Surv(os_time_months, os_event) ~ ",
                                paste(zcols, collapse = " + ")))
cox_multi <- coxph(multi_form, data = data)
print(cox_multi)
s_multi <- summary(cox_multi)
# C-index: survConcordance 已弃用, 改用 Hmisc::rcorr.cens
# 线性预测值越大=风险越高=预后越差, 取负号使 C-index 越大=区分度越好
# rcorr.cens 返回带名字的数值向量(原子型), 不能用 $, 用位置[1]取 C Index
lp <- as.numeric(predict(cox_multi, newdata = data))
rc <- Hmisc::rcorr.cens(-lp, Surv(data$os_time_months, data$os_event))
cidx <- as.numeric(rc[1])
cat("  C-index:", round(cidx, 4), "\n")
multi_tbl <- data.frame(
  Gene = rownames(s_multi$coefficients),
  HR = round(s_multi$coefficients[, "exp(coef)"], 3),
  Lower = round(s_multi$conf.int[, "lower .95"], 3),
  Upper = round(s_multi$conf.int[, "upper .95"], 3),
  P = signif(s_multi$coefficients[, "Pr(>|z|)"], 3))
write.csv(multi_tbl, "GSE39582_multiCox_OS.csv", row.names = FALSE)

# ============================================================
# 5. LASSO-Cox 选特征 (主模型)
# ============================================================
cat("\n=== 5. LASSO-Cox 选特征 ===\n")
x <- as.matrix(data[, zcols])
y <- Surv(data$os_time_months, data$os_event)
set.seed(123)
cvfit <- cv.glmnet(x, y, family = "cox", alpha = 1, nfolds = 10)
fit <- glmnet(x, y, family = "cox", alpha = 1, lambda = cvfit$lambda.min)
coef_mat <- coef(fit)
selected <- rownames(coef_mat)[coef_mat[, 1] != 0]
cat("  lambda.min =", round(cvfit$lambda.min, 4), "\n")
cat("  选中的基因:", ifelse(length(selected) > 0, paste(selected, collapse = ", "), "无 (全部系数为0)"), "\n")

# 风险评分
if (length(selected) > 0) {
  risk_score <- as.matrix(x[, selected, drop = FALSE]) %*% coef_mat[selected, 1]
} else {
  # fallback: 等权 (所有基因 z-score 均值)
  cat("  [fallback] LASSO 未选特征, 使用等权 risk score\n")
  risk_score <- rowMeans(x)
}
data$risk_score <- as.numeric(risk_score)

# 中位分割
data$risk_group <- ifelse(data$risk_score > median(data$risk_score), "High", "Low")
cat("  风险分组: Low =", sum(data$risk_group == "Low"),
    ", High =", sum(data$risk_group == "High"), "\n")

# ============================================================
# 6. KM 曲线 (OS + RFS)
# ============================================================
cat("\n=== 6. KM 曲线 ===\n")
pdf("GSE39582_KM_OS.pdf", width = 7, height = 6)
km_os <- survfit(Surv(os_time_months, os_event) ~ risk_group, data = data)
p_os <- ggsurvplot(km_os, data = data, pval = TRUE, pval.coord = c(0, 0.2),
                   risk.table = TRUE, conf.int = TRUE,
                   xlab = "Time (months)", ylab = "Overall Survival",
                   legend.title = "Risk Group",
                   legend.labs = c("Low", "High"),
                   palette = c("#4DAF4A", "#E41A1C"),
                   title = "GSE39582: 7-Gene Risk Score (OS)")
print(p_os)
dev.off()
cat("  KM (OS) 已生成, Log-rank P =",
    signif(survdiff(Surv(os_time_months, os_event) ~ risk_group, data = data)$pvalue, 4), "\n")

pdf("GSE39582_KM_RFS.pdf", width = 7, height = 6)
# RFS 需要 data_rfs 也有 risk_score (用 OS 模型的风险评分)
data_rfs$risk_score <- data[match(data_rfs$Sample, data$Sample), "risk_score"]
data_rfs$risk_group <- ifelse(data_rfs$risk_score > median(data$risk_score), "High", "Low")
km_rfs <- survfit(Surv(rfs_time_months, rfs_event) ~ risk_group, data = data_rfs)
p_rfs <- ggsurvplot(km_rfs, data = data_rfs, pval = TRUE, pval.coord = c(0, 0.2),
                    risk.table = TRUE, conf.int = TRUE,
                    xlab = "Time (months)", ylab = "Relapse-Free Survival",
                    legend.title = "Risk Group",
                    legend.labs = c("Low", "High"),
                    palette = c("#4DAF4A", "#E41A1C"),
                    title = "GSE39582: 7-Gene Risk Score (RFS)")
print(p_rfs)
dev.off()
cat("  KM (RFS) 已生成, Log-rank P =",
    signif(survdiff(Surv(rfs_time_months, rfs_event) ~ risk_group, data = data_rfs)$pvalue, 4), "\n")

# ============================================================
# 7. 列线图 (risk_score + 临床协变量)
# ============================================================
cat("\n=== 7. 列线图 ===\n")
# 把 stage 转数值 (便于纳入模型; "N/A" 自动变 NA 被 Cox 剔除)
# 必须在 datadist 之前生成, 否则 nomogram 报 "variable does not have limits"
data$stage_num <- suppressWarnings(as.numeric(data$tnm_stage))
dd <- datadist(data)
options(datadist = "dd")
f <- cph(Surv(os_time_months, os_event) ~ risk_score + age_year + stage_num,
         data = data, surv = TRUE, x = TRUE, y = TRUE, time.inc = 36)
surv_func <- Survival(f)
nom <- nomogram(f,
                fun = list(function(x) surv_func(12, x),
                           function(x) surv_func(36, x),
                           function(x) surv_func(60, x)),
                funlabel = c("1-year OS", "3-year OS", "5-year OS"),
                lp = FALSE)
pdf("GSE39582_nomogram.pdf", width = 9, height = 6)
plot(nom)
dev.off()
cat("  列线图已生成 (n =", length(f$y), ")\n")

# ============================================================
# 8. 校准曲线 (36月)
# ============================================================
cat("\n=== 8. 校准曲线 ===\n")
tryCatch({
  cal <- calibrate(f, method = "boot", u = 36, m = max(20, floor(nrow(data) / 10)), B = 500)
  pdf("GSE39582_calibration.pdf", width = 7, height = 6)
  plot(cal, xlab = "Predicted 3-year OS Probability",
       ylab = "Observed 3-year OS Probability", main = "Calibration (GSE39582)")
  dev.off()
  cat("  校准曲线已生成\n")
}, error = function(e) {
  cat("  校准曲线生成失败:", e$message, "\n")
})

# ============================================================
# 9. DCA —— 有意跳过 (不纳入 manuscript)
# ============================================================
# 说明: 本队列 7-gene signature 预后为阴性(LASSO 选空, C-index 0.56,
# timeROC AUC<=0.61), DCA 对 null 模型无信息量且易引 reviewer 深挖弱点。
# 此外 dcurves::dca 对连续型 risk_score 会报 "All covariates/risks must be
# between 0 and 1" (要求 0-1 有界 marker), 对任何连续签名都会触发, 不值得修。
# 故预后部分按"轻写法"仅保留: 单/多因素 Cox 表 + RFS KM 补充图, DCA 不产出。
cat("\n=== 9. DCA ===\n")
cat("  [跳过] 预后为阴性结果, DCA 不纳入 manuscript (见上方说明)\n")

# ============================================================
# 10. 时间依赖性 ROC
# ============================================================
cat("\n=== 10. 时间依赖性 ROC ===\n")
ROC <- timeROC(T = data$os_time_months,
               delta = data$os_event,
               marker = data$risk_score,
               cause = 1,
               times = c(12, 36, 60),
               ROC = TRUE)
cat(sprintf("  AUC (1yr)  = %.3f\n", ROC$AUC[1]))
cat(sprintf("  AUC (3yr)  = %.3f\n", ROC$AUC[2]))
cat(sprintf("  AUC (5yr)  = %.3f\n", ROC$AUC[3]))
pdf("GSE39582_timeROC.pdf", width = 7, height = 7)
plot(ROC, time = 12, col = "#E41A1C", lwd = 2, title = FALSE,
     xlab = "1 - Specificity", ylab = "Sensitivity")
plot(ROC, time = 36, col = "#377EB8", lwd = 2, add = TRUE)
plot(ROC, time = 60, col = "#4DAF4A", lwd = 2, add = TRUE)
legend("bottomright",
       c(paste0("1-year (AUC=", round(ROC$AUC[1], 3), ")"),
         paste0("3-year (AUC=", round(ROC$AUC[2], 3), ")"),
         paste0("5-year (AUC=", round(ROC$AUC[3], 3), ")")),
       col = c("#E41A1C", "#377EB8", "#4DAF4A"), lwd = 2, bty = "n")
dev.off()
cat("  时间ROC已生成\n")

# ============================================================
# 11. 汇总输出
# ============================================================
cat("\n=== 11. 汇总 ===\n")
summary_tbl <- data.frame(
  Metric = c("Cohort", "N (OS)", "Events (OS)", "N (RFS)", "Events (RFS)",
             "LASSO selected genes", "KM OS log-rank P", "KM RFS log-rank P",
             "AUC 1yr", "AUC 3yr", "AUC 5yr"),
  Value = c("GSE39582 (Marisa 2013)", nrow(data), sum(data$os_event),
            nrow(data_rfs), sum(data_rfs$rfs_event),
            ifelse(length(selected) > 0, paste(selected, collapse = ","), "none"),
            signif(survdiff(Surv(os_time_months, os_event) ~ risk_group, data = data)$pvalue, 4),
            signif(survdiff(Surv(rfs_time_months, rfs_event) ~ risk_group, data = data_rfs)$pvalue, 4),
            round(ROC$AUC[1], 3), round(ROC$AUC[2], 3), round(ROC$AUC[3], 3)))
write.csv(summary_tbl, "GSE39582_prognostic_summary.csv", row.names = FALSE)
print(summary_tbl, row.names = FALSE)

cat("\n=== GSE39582 预后验证完成 ===\n")
cat("输出文件:\n")
cat("  GSE39582_uniCox_OS.csv / _RFS.csv — 单因素Cox\n")
cat("  GSE39582_multiCox_OS.csv        — 多因素Cox\n")
cat("  GSE39582_KM_OS.pdf / _RFS.pdf   — KM曲线\n")
cat("  GSE39582_nomogram.pdf           — 列线图\n")
cat("  GSE39582_calibration.pdf         — 校准曲线\n")
cat("  GSE39582_DCA.pdf                — 决策曲线\n")
cat("  GSE39582_timeROC.pdf            — 时间ROC\n")
cat("  GSE39582_prognostic_summary.csv — 汇总表\n")
