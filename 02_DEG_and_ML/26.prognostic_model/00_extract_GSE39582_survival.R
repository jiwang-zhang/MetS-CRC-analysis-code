# ============================================================
# 附加: 从本地 series matrix 提取 GSE39582 生存数据
# ============================================================
# GSE39582 (Marisa 2013, PLoS Med) 含 OS/RFS 双终点生存数据
# 本脚本用 BASE R 解析本地拷贝的 series matrix (无需联网/GEOquery)
# 输出: GSE39582_survival.csv (供 03_prognostic_GSE39582.R 使用)
#
# 用法: 先确保以下文件存在, 然后 Source 本脚本
#   源文件: 第二部分DEG/7.databases - 副本/GSE39582/GSE39582_series_matrix.txt/GSE39582_series_matrix.txt
# ============================================================

setwd("D:\\BaiduSyncdisk\\daixie-colorectal\\第二部分DEG\\26.prognostic_model")

series_file <- "D:/BaiduSyncdisk/daixie-colorectal/第二部分DEG/7.databases - 副本/GSE39582/GSE39582_series_matrix.txt/GSE39582_series_matrix.txt"
out_file    <- "GSE39582_survival.csv"

cat("=== 解析 series matrix ===\n")
cat("  文件:", series_file, "\n")

samples <- NULL
fields <- list()

# 逐行读取 (内存友好, 支持大文件)
con <- file(series_file, "r", encoding = "UTF-8")
while (length(ln <- readLines(con, n = 1)) > 0) {
  if (startsWith(ln, "!Sample_geo_accession")) {
    parts <- strsplit(ln, "\t")[[1]]
    samples <- gsub('"', '', parts[-1])
  } else if (startsWith(ln, "!Sample_characteristics_ch1")) {
    parts <- strsplit(ln, "\t")[[1]]
    if (length(parts) < 2) next
    first <- gsub('"', '', parts[2])
    # 提取 key: "os.event: 1" -> "os.event" ; "os.delay (months): 50" -> "os.delay (months)"
    key <- sub(":.*$", "", first)
    if (key != "") {
      vals <- gsub('"', '', parts[-1])
      fields[[key]] <- vals
    }
  }
}
close(con)
cat("  样本数:", length(samples), "\n")
cat("  提取到的字段数:", length(fields), "\n")

getf <- function(k) {
  if (k %in% names(fields)) fields[[k]] else rep(NA_character_, length(samples))
}
splitv <- function(v) sub("^[^:]*:\\s*", "", v)   # 取冒号后的部分

os_event  <- getf("os.event")
os_time   <- getf("os.delay(months)")
rfs_event <- getf("rfs.event")
rfs_time  <- getf("rfs.delay")
age       <- getf("age.at.diagnosis(year)")
sex       <- getf("Sex")
stage     <- getf("tnm.stage")
mmr       <- getf("mmr.status")

# 构建数据框
df <- data.frame(
  Sample          = samples,
  os_event        = as.integer(splitv(os_event)),
  os_time_months  = suppressWarnings(as.numeric(splitv(os_time))),
  rfs_event       = as.integer(splitv(rfs_event)),
  rfs_time_months = suppressWarnings(as.numeric(splitv(rfs_time))),
  age_year        = suppressWarnings(as.numeric(splitv(age))),
  sex             = splitv(sex),
  tnm_stage       = splitv(stage),
  mmr_status      = splitv(mmr),
  stringsAsFactors = FALSE
)

write.csv(df, out_file, row.names = FALSE)
cat("  已写入:", out_file, "\n")

# 简单统计
os_valid  <- sum(!is.na(df$os_time_months) & df$os_time_months > 0 & df$os_event %in% c(0, 1), na.rm = TRUE)
rfs_valid <- sum(!is.na(df$rfs_time_months) & df$rfs_time_months > 0 & df$rfs_event %in% c(0, 1), na.rm = TRUE)
cat("  OS  有效样本:", os_valid, " (事件=", sum(df$os_event == 1, na.rm = TRUE), ")\n")
cat("  RFS 有效样本:", rfs_valid, " (事件=", sum(df$rfs_event == 1, na.rm = TRUE), ")\n")
cat("\n提取完成. 接下来运行 03_prognostic_GSE39582.R\n")
