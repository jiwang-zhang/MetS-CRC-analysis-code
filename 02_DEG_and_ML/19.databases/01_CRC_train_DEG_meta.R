# ============================================================
# 01_CRC_train_DEG_meta.R
# MetS x CRC 生信分析 · 第1步（CRC 侧）
# 目标：3 个 GPL570 训练集(GSE21510 / GSE32323 / GSE39582)
#       分别 limma DEG → 基因水平逆方差 meta → 稳健 CRC DEG 清单
#
# 工作目录：请把本脚本放在 19.databases 内，或在 RStudio 里
#          setwd("D:/daixie-colorectal/19.databases")
# 数据格式假设：每个 GSE 已下载 GEO series matrix，文件名
#          GSEXXXX_series_matrix.txt.gz （GEO 官网标准命名）
# 平台注释：GPL570 用 hgu133plus2.db（需安装）
#
# 需安装包（一次性）：
#   install.packages("BiocManager")
#   BiocManager::install(c("limma","hgu133plus2.db","AnnotationDbi"))
# ============================================================

library(limma)
library(hgu133plus2.db)
library(AnnotationDbi)

# ---------- 0. 参数（与两篇模板一致：|log2FC|>0.5 且 adj.P<0.05）----------
DEG_LFC  <- 0.5
DEG_PADJ <- 0.05
WORKDIR  <- "D:/daixie-colorectal/19.databases"   # ← 改成你的实际路径
setwd(WORKDIR)

train_gse <- c("GSE21510", "GSE32323", "GSE39582")

# 每个数据集里区分 Normal / Tumor 的正则（GEO 表型措辞不一，先按常见写法；
# 跑出分组数后如不对，再调这里的词）
NORMAL_RE <- "(?i)normal|control|non-?tumor|adjacent|mucosa"
TUMOR_RE  <- "(?i)tumor|cancer|carcinoma|adenocarcinoma|neoplas|colorectal"

# ---------- 1. 读 GEO series matrix（手动解析表头，透明可控）----------
read_gse_matrix <- function(gse, destdir = ".") {
  f <- file.path(destdir, paste0(gse, "_series_matrix.txt.gz"))
  if (!file.exists(f)) stop("找不到 ", f, " —— 请先确认 series matrix 已下载到本文件夹")
  con <- gzfile(f, "rt"); raw <- readLines(con); close(con)

  begin   <- grep("!series_matrix_table_begin", raw)
  hdr     <- raw[1:begin]                              # 所有 ! 开头的表头行
  samp_ln <- hdr[grepl("^!Sample_", hdr)]              # 每个样本的特征行

  # 把每条 !Sample_xxx 行解析成「字段名 + 按样本排列的取值向量」
  fields <- lapply(samp_ln, function(l) {
    p <- strsplit(l, "\t")[[1]]
    list(name = sub("^!Sample_", "", p[1]), vals = p[-1])
  })

  # 样本顺序 = 表达表里 ID_REF 行之后的列顺序（GSM 编号）
  idrow  <- raw[begin + 1]
  cols   <- strsplit(idrow, "\t")[[1]]
  samp_ids <- cols[-1]
  n <- length(samp_ids)

  anno <- data.frame(sample = samp_ids, stringsAsFactors = FALSE)
  for (fld in fields) {
    anno[[fld$name]] <- if (length(fld$vals) >= n) fld$vals[1:n] else rep(NA_character_, n)
  }

  expr_lines <- raw[(begin + 1):length(raw)]
  expr <- read.delim(text = expr_lines, check.names = FALSE, stringsAsFactors = FALSE)
  rownames(expr) <- expr$ID_REF
  mat  <- as.matrix(expr[, colnames(expr) != "ID_REF", drop = FALSE])
  mode(mat) <- "numeric"
  colnames(mat) <- samp_ids

  list(expr = mat, anno = anno, gse = gse)
}

# ---------- 2. 单数据集 DEG（Normal vs Tumor）----------
deg_one <- function(gse, destdir = ".") {
  g <- read_gse_matrix(gse, destdir)
  anno <- g$anno
  # 把每个样本的所有特征拼成一段文本，用来匹配 Normal/Tumor
  grp_text <- apply(anno, 1, function(r) paste(na.omit(r), collapse = " | "))

  isN <- grepl(NORMAL_RE, grp_text)
  isT <- grepl(TUMOR_RE,  grp_text)
  cat(sprintf("[%s] 样本总数=%d  匹配到 正常=%d  肿瘤=%d\n",
              gse, ncol(g$expr), sum(isN), sum(isT)))
  if (sum(isN) == 0 || sum(isT) == 0)
    stop(sprintf("[%s] 分组失败，请检查 NORMAL_RE / TUMOR_RE 正则或表型字段", gse))

  group <- factor(ifelse(isT, "Tumor", "Normal"), levels = c("Normal", "Tumor"))
  design <- model.matrix(~0 + group); colnames(design) <- c("Normal", "Tumor")
  contr  <- makeContrasts(Tumor - Normal, levels = design)

  fit <- lmFit(g$expr, design) %>% contrasts.fit(contr) %>% eBayes()
  tt  <- topTable(fit, number = Inf, sort.by = "none")
  tt$ID_REF <- rownames(tt)

  # 注释 GPL570 探针 -> 基因符号
  tt$SYMBOL <- mapIds(hgu133plus2.db, keys = tt$ID_REF,
                      column = "SYMBOL", keytype = "PROBEID")
  tt
}

# ---------- 3. 探针 -> 基因（每个基因保留 |logFC| 最大的探针）----------
collapse_gene <- function(tt) {
  tt <- tt[!is.na(tt$SYMBOL) & tt$SYMBOL != "", ]
  tt$abs <- abs(tt$logFC)
  tt <- tt[order(tt$abs, decreasing = TRUE), ]
  tt[!duplicated(tt$SYMBOL), ]
}

# ---------- 4. 基因水平逆方差 meta（透明、不依赖黑盒包）----------
# 思路：每个研究给每个基因的 z = sign(logFC)*qnorm(1-P/2)，SE = |logFC|/|z|
#       跨研究用 1/SE^2 加权合并效应，再算合并 z、合并 P、FDR。
gene_meta <- function(tops, lfc = DEG_LFC, padj = DEG_PADJ, min_studies = 2) {
  genes <- Reduce(union, lapply(tops, function(t) t$SYMBOL))
  out <- data.frame(SYMBOL = genes, stringsAsFactors = FALSE)
  out$meta_logFC <- NA_real_; out$meta_SE <- NA_real_
  out$meta_P <- NA_real_;   out$n_studies <- 0L

  for (i in seq_along(genes)) {
    g <- genes[i]; lfcs <- c(); ses <- c(); ns <- 0L
    for (t in tops) {
      r <- t[t$SYMBOL == g, ]
      if (nrow(r) == 0) next
      z  <- sign(r$logFC) * qnorm(1 - r$P.Value / 2)
      se <- abs(r$logFC) / max(abs(z), 1e-6)
      lfcs <- c(lfcs, r$logFC); ses <- c(ses, se); ns <- ns + 1L
    }
    if (ns < min_studies) next
    w    <- 1 / ses^2
    cLFC <- sum(w * lfcs, na.rm = TRUE) / sum(w, na.rm = TRUE)
    cSE  <- sqrt(1 / sum(w, na.rm = TRUE))
    cZ   <- cLFC / cSE
    cP   <- 2 * pnorm(-abs(cZ))
    out[i, c("meta_logFC", "meta_SE", "meta_P", "n_studies")] <- c(cLFC, cSE, cP, ns)
  }
  out$meta_FDR <- p.adjust(out$meta_P, method = "fdr")
  out[abs(out$meta_logFC) > lfc & out$meta_FDR < padj & out$n_studies >= min_studies, ]
}

# ---------- 5. 主流程 ----------
cat("===== 步骤A：3 个 GPL570 训练集分别 DEG =====\n")
t1 <- collapse_gene(deg_one("GSE21510", WORKDIR))
t2 <- collapse_gene(deg_one("GSE32323", WORKDIR))
t3 <- collapse_gene(deg_one("GSE39582", WORKDIR))
tops <- list(t1, t2, t3)

cat("\n===== 步骤B：基因水平 meta 合并 =====\n")
crc_meta <- gene_meta(tops)
cat(sprintf("CRC 训练集 meta 后显著 DEG 数：%d\n", nrow(crc_meta)))

write.csv(crc_meta, "CRC_train_meta_DEG.csv", row.names = FALSE)
cat("已保存 -> CRC_train_meta_DEG.csv\n")

# ---------- 6. 后续（脚本先留接口，等 MetS 数据到位再接）----------
# 下一步：把 MetS 血液 meta 得到的 DEG 清单 与 crc_meta 取交集 -> co-DEGs
#         co-DEGs 再 ∩ 热点基因集（炎症/胰岛素/Wnt）-> 富集 -> WGCNA -> ML -> ...
cat("\n[提示] 下一步需 MetS 侧 DEG 清单（GSE98895 + GSE200744 血液 meta）。\n")
