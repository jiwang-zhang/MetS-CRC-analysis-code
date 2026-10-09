# =============================================================
# Forest_MetS_CRC_bidirectional.R
# 用途：MetS<->FinnGen CRC 双向 MR 森林图（对齐 forest.pdf 风格：
#       底部 x 轴"标尺" + 组间/表头"横线" + OR=1 红虚线）
# 正向：3.MR/17.MR/MetS__FinnGen_CRC.table.MRresult.csv   (MetS->CRC, 185 SNP, 5 方法)
# 反向：第一部分nhanes+MR/17.MR/reverse/FinnGen_CRC__MetS.table.MRresult.csv (CRC->MetS, 7 SNP, 3 方法)
# 依赖：forestploter / grid / readr   install.packages(c("forestploter","grid","readr"))
# 输出：3.MR/19.forest/Forest_MetS_CRC_bidirectional.pdf / .png
# =============================================================

suppressPackageStartupMessages({
  library(grid); library(readr); library(forestploter)
})

# ---- 路径（绝对，避免硬编码） ----
fwd_path <- "D:/BaiduSyncdisk/daixie-colorectal/3.MR/17.MR/MetS__FinnGen_CRC.table.MRresult.csv"
rev_path <- "D:/BaiduSyncdisk/daixie-colorectal/第一部分nhanes+MR/17.MR/reverse/FinnGen_CRC__MetS.table.MRresult.csv"
out_dir  <- "D:/BaiduSyncdisk/daixie-colorectal/3.MR/19.forest"
if (!file.exists(fwd_path)) stop("正向结果缺失: ", fwd_path)
if (!file.exists(rev_path)) stop("反向结果缺失: ", rev_path)

# ---- 读 + 合并 ----
df_fwd <- read.csv(fwd_path, check.names = FALSE)
df_rev <- read.csv(rev_path, check.names = FALSE)
data   <- rbind(df_fwd, df_rev)

# ---- 只保留 5 主方法（反向天然只有前 3 种） ----
keep <- c("Inverse variance weighted","MR Egger","Weighted median",
          "Simple mode","Weighted mode")
data <- data[data$method %in% keep, ]
data$method <- factor(data$method, levels = keep)

# ---- 排序：先正向(MetS->CRC)，后反向(CRC->MetS)；组内按方法顺序 ----
data$dir <- ifelse(data$exposure == "MetS", 1, 2)   # 正向 exposure=MetS
data <- data[order(data$dir, data$method), ]
rownames(data) <- NULL

# ---- 显示列整理 ----
data$" "          <- paste(rep(" ", 5), collapse = " ")
data$"OR(95% CI)" <- sprintf("%.2f (%.2f to %.2f)",
                             data$or, data$or_lci95, data$or_uci95)
data$pval <- ifelse(data$pval < 0.001, "<0.001", sprintf("%.3f", data$pval))
data$nsnp <- as.integer(data$nsnp)

# 组内重复 exposure/outcome 置空（仅每组首行显示）
data2 <- data[, c("exposure","outcome")]
data[duplicated(data2), "exposure"] <- ""
data[duplicated(data2), "outcome"]  <- ""

# ---- 主题（对齐 forest.pdf） ----
tm <- forest_theme(base_size = 12,
                   ci_pch = 15, ci_lty = 1, ci_lwd = 1.5, ci_col = "black",
                   ci_fill = "#339900FF", ci_Theight = 0.2,
                   refline_gp = gpar(lty = "dashed", lwd = 1, col = "red"),
                   xaxis_gp   = gpar(fontsize = 12))

# ---- 出图（x 轴标尺：0.5–1.5 + 刻度；ci 画在第 6 列空白列） ----
plot <- forestploter::forest(
  data[, c("exposure","outcome","nsnp","method","pval"," ","OR(95% CI)")],
  est    = data$or,
  lower  = data$or_lci95,
  upper  = data$or_uci95,
  ci_column = 6,
  ref_line  = 1,
  xlim      = c(0.5, 1.5),
  ticks_at  = c(0.5, 1, 1.5),
  theme     = tm
)

# ---- 5 方法 5 色 box ----
boxcolor <- c("Inverse variance weighted" = "#E64B35",
              "MR Egger"                = "#4DBBD5",
              "Weighted median"         = "#00A087",
              "Simple mode"             = "#3C5488",
              "Weighted mode"           = "#F39B7F")
for (i in 1:nrow(data)) {
  plot <- edit_plot(plot, col = 6, row = i, which = "ci",
                    gp = gpar(fill = boxcolor[as.character(data$method[i])],
                              col  = "black"))
}

# ---- P<0.05 加粗 ----
pos_bold <- which(as.numeric(gsub("<", "", data$pval)) < 0.05)
for (i in pos_bold) {
  plot <- edit_plot(plot, col = 5, row = i, which = "text",
                    gp = gpar(fontface = "bold"))
}

# ---- 横线：共 4 根 —— 表头顶线 + 表头底线 + 正反分隔(数据行5底) + 图底粗线(数据行10底) ----
# ⚠ 不能用 table(exposure,outcome) 算 lineVec：空字符串水平会污染分组，
#   实测 lineVec 会变成 c(8,8,8,8,8,9,9,10,10) 导致 row 错位（线落到 0.214 上方）。
#   改用 data$dir（排序后正向在前、反向在后）直接取行数，对 header 行数稳健。
n_fwd <- sum(data$dir == 1)        # 正向行数（MetS->CRC）
lineVec <- c(n_fwd, nrow(data))    # -> c(5, 10)
# ① 表头顶线（header 上沿，粗黑）
plot <- add_border(plot, part = "header", row = 1, where = "top",
                   gp = gpar(lwd = 2))
# ①b 表头底线（header 下沿，粗黑）—— 表头上下框线，分隔表头与数据区
plot <- add_border(plot, part = "header", row = 1, where = "bottom",
                   gp = gpar(lwd = 2))
# ② 区分正反 MR（body 数据行5 底 = 正向末行下方 = 反向 IVW 0.073 上方）
plot <- add_border(plot, part = "body", row = lineVec[1], where = "bottom",
                   gp = gpar(lwd = 1))
# ③ 图底横线（body 数据行10 底 = 反向末行下方 = 标尺之下，粗黑）
plot <- add_border(plot, part = "body", row = lineVec[2], where = "bottom",
                   gp = gpar(lwd = 2))

# ---- 字号 + 居中 ----
plot <- edit_plot(plot, col = 1:ncol(data), row = 1:nrow(data),
                  which = "text", gp = gpar(fontsize = 12))
plot <- edit_plot(plot, col = 1:ncol(data), which = "text",
                  hjust = unit(0.5, "npc"), x = unit(0.5, "npc"))
plot <- edit_plot(plot, col = 1:ncol(data), which = "text",
                  hjust = unit(0.5, "npc"), x = unit(0.5, "npc"), part = "header")

# ---- 输出 ----
out_pdf <- file.path(out_dir, "Forest_MetS_CRC_bidirectional.pdf")
out_png <- file.path(out_dir, "Forest_MetS_CRC_bidirectional.png")
pdf(out_pdf, width = 10, height = 4); print(plot); dev.off()
png(out_png, width = 10, height = 4, units = "in", res = 300)
print(plot); dev.off()

cat("DONE ->", out_pdf, "and", out_png, "\n\n")
for (i in 1:nrow(data)) {
  cat(sprintf("%-20s %-18s n=%-4d %-22s p=%s\n",
              ifelse(data$exposure[i]=="", "  ", data$exposure[i]),
              ifelse(data$outcome[i]=="", "  ", data$outcome[i]),
              data$nsnp[i], data$method[i], data$pval[i]))
}