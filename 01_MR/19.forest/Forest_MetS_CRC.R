# =============================================================
# Forest_MetS_CRC.R
# 用途：把 MetS->FinnGen CRC 主结果 (17.MR/MetS__FinnGen_CRC.table.MRresult.csv)
#      5 方法画成森林图（参考风格：左文字右图形，OR=1 红虚线）
# 依赖：forestploter / grid / readr
#       install.packages(c("forestploter","grid","readr"))
# 输出：3.MR/19.forest/Forest_MetS_CRC.pdf / .png
# =============================================================

suppressPackageStartupMessages({
  library(grid); library(readr); library(forestploter)
})

# ---- 路径 ----
script_dir <- "D:/BaiduSyncdisk/daixie-colorectal/3.MR/19.forest"
data_path  <- file.path(script_dir, "..", "17.MR",
                        "MetS__FinnGen_CRC.table.MRresult.csv")
if (!file.exists(data_path)) stop("找不到: ", data_path)

# ---- 读 + 筛选 5 主方法 ----
df <- read.csv(data_path, check.names = FALSE)
keep <- c("Inverse variance weighted","MR Egger","Weighted median",
          "Simple mode","Weighted mode")
df <- df[df$method %in% keep, ]
# 方法顺序：IVW 在上
df$method <- factor(df$method, levels = keep)
df <- df[order(df$method), ]
rownames(df) <- NULL
if (nrow(df) != 5) stop("期望 5 行方法，实际 ", nrow(df))

# ---- 整理显示列 ----
df$" "          <- paste(rep(" ", 5), collapse = " ")
df$"OR(95% CI)" <- sprintf("%.2f (%.2f to %.2f)",
                           df$or, df$or_lci95, df$or_uci95)
df$pval <- ifelse(df$pval < 0.001, "<0.001", sprintf("%.3f", df$pval))
df$nsnp <- as.integer(df$nsnp)

# ---- 主题 ----
tm <- forest_theme(base_size = 12,
                   ci_pch = 15, ci_lty = 1, ci_lwd = 1.5, ci_col = "black",
                   ci_fill = "#339900FF", ci_Theight = 0.2,
                   refline_gp = gpar(lty = "dashed", lwd = 1, col = "red"),
                   xaxis_gp  = gpar(fontsize = 12))

# ---- 出图 ----
plot <- forestploter::forest(
  df[, c("exposure","outcome","nsnp","method","pval"," ","OR(95% CI)")],
  est    = df$or,
  lower  = df$or_lci95,
  upper  = df$or_uci95,
  ci_column = 6,
  ref_line  = 1,
  xlim      = c(0.5, 1.5),
  ticks_at  = c(0.5, 1, 1.5),
  theme     = tm
)

# 5 方法 5 色 box（IVW 红 / Egger 青 / Wmedian 绿 / Simple 深蓝 / Wmode 橙）
boxcolor <- c("Inverse variance weighted"="#E64B35",
              "MR Egger"                ="#4DBBD5",
              "Weighted median"         ="#00A087",
              "Simple mode"             ="#3C5488",
              "Weighted mode"           ="#F39B7F")
for (i in 1:nrow(df)) {
  plot <- edit_plot(plot, col = 6, row = i, which = "ci",
                    gp = gpar(fill = boxcolor[as.character(df$method[i])],
                              col  = "black"))
}

# P<0.05 加粗
pos_bold <- which(as.numeric(gsub("<", "", df$pval)) < 0.05)
for (i in pos_bold) {
  plot <- edit_plot(plot, col = 5, row = i, which = "text",
                    gp = gpar(fontface = "bold"))
}

# 表头加粗 + 加顶线
plot <- add_border(plot, part = "header", row = 1, where = "top",
                   gp = gpar(lwd = 2))
plot <- edit_plot(plot, col = 1:ncol(df), row = 1:nrow(df),
                  which = "text", gp = gpar(fontsize = 11))
plot <- edit_plot(plot, col = 1:ncol(df), which = "text",
                  hjust = unit(0.5, "npc"), x = unit(0.5, "npc"),
                  part = "header")

# ---- 输出 ----
out_pdf <- file.path(script_dir, "Forest_MetS_CRC.pdf")
out_png <- file.path(script_dir, "Forest_MetS_CRC.png")
pdf(out_pdf, width = 10, height = 3.5); print(plot); dev.off()
png(out_png, width = 10, height = 3.5, units = "in", res = 300)
print(plot); dev.off()

cat("DONE ->", out_pdf, "and", out_png, "\n\n")
cat(sprintf("%-22s nsnp=%-4d OR=%-7s p=%s\n",
            df$method, df$nsnp,
            sprintf("%.3f (%.3f-%.3f)", df$or, df$or_lci95, df$or_uci95),
            df$pval))