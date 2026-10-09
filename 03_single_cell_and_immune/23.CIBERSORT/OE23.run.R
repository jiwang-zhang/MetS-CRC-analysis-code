# ============================================================
# Step 7a: CIBERSORT 免疫细胞去卷积
# 输入: CRC合并表达矩阵 (merged_after_batch_removal.csv)
# 输出: CIBERSORT-Results.txt (22种免疫细胞比例 + P-value)
# 
# 依赖: e1071, parallel, preprocessCore
# ============================================================

inputFile <- "../../第二部分DEG/25.diagnostic_model/merged_after_batch_removal.csv"   # 表达数据文件
setwd("D:/BaiduSyncdisk/daixie-colorectal/4.免疫浸润单细胞/23.CIBERSORT")

source("OE23.CIBERSORT.R")

# --- 准备 CIBERSORT 输入 (制表符分隔, 基因×样本) ---
rt <- read.table(inputFile, header = TRUE, sep = ",", check.names = FALSE)
rt <- as.matrix(rt)
rownames(rt) <- rt[, 1]
exp <- rt[, 2:ncol(rt)]
dn  <- list(rownames(exp), colnames(exp))
exp <- matrix(as.numeric(exp), nrow = nrow(exp), dimnames = dn)
suppressMessages(library(limma))
exp <- avereps(exp)

mixtureFile <- "mixture_input.txt"
write.table(data.frame(Gene = rownames(exp), exp, check.names = FALSE),
            file = mixtureFile, sep = "\t", row.names = FALSE, quote = FALSE)
cat("混合物文件已写:", mixtureFile, "(", nrow(exp), "基因 x", ncol(exp), "样本)\n")

# --- 运行 CIBERSORT ---
cat("运行 CIBERSORT (perm=1000)...\n")
outTab <- CIBERSORT("ref.txt", mixtureFile, perm = 1000)

# 注意：CIBERSORT 的 P-value 是单样本反卷积可信度指标。
# 正常结肠组织与 LM22 训练集差异较大，P 值天然偏高，
# 若按 P<0.05 过滤会系统性删掉正常样本（61 正常仅剩 5 例），
# 造成组间样本量极不平衡、Normal 那根柱不可信。
# 故【取消 P 过滤】，保留全部 748 样本用于免疫浸润构成展示
# （与 ESTIMATE / MCPcounter 所用的全 748 样本对齐）。
# outTab <- outTab[outTab[, "P-value"] < 0.05, ]
outTab <- outTab[, 1:(ncol(outTab) - 3), drop = FALSE]    # 仅保留 22 种免疫细胞
outTab <- rbind(id = colnames(outTab), outTab)
write.table(outTab, file = "CIBERSORT-Results.txt", sep = "\t", quote = FALSE, col.names = FALSE)
cat("CIBERSORT 完成: 保存 CIBERSORT-Results.txt (", nrow(outTab) - 1, "样本)\n")
